from collections import Counter
import numpy as np

def calculate_hmm_parameters_extended(alignments):
    num_sequences = len(alignments)
    sequence_length = len(alignments[0])
    match_threshold = num_sequences / 2

    # Definir estados y transiciones
    states = ['M' if sum(seq[col] != '-' for seq in alignments) > match_threshold else 'I'
              for col in range(sequence_length)]
    transitions = {'MM': 0, 'MI': 0, 'IM': 0, 'II': 0}

    # Contar transiciones
    for i in range(len(states) - 1):
        transition_type = states[i] + states[i + 1]
        transitions[transition_type] += 1

    # Calcular probabilidades de transición
    total_transitions = sum(transitions.values())
    for key in transitions:
        transitions[key] = (transitions[key] + 1) / (total_transitions + len(transitions))  # Laplace smoothing

    # Calcular probabilidades de emisión
    emission_probabilities = {'M': Counter(), 'I': Counter()}
    for i, state in enumerate(states):
        for seq in alignments:
            if seq[i] != '-':
                emission_probabilities[state][seq[i]] += 1

    # Normalizar probabilidades de emisión con Laplace smoothing
    for state in emission_probabilities:
        total = sum(emission_probabilities[state].values())
        for symbol in emission_probabilities[state]:
            emission_probabilities[state][symbol] = (emission_probabilities[state][symbol] + 1) / (total + len(emission_probabilities[state]))

    # Probabilidades de background equiprobables
    symbols = set(''.join(alignments).replace('-', ''))
    background_probabilities = {symbol: 1/len(symbols) for symbol in symbols}

    return states, transitions, emission_probabilities, background_probabilities
alignments = [
    "VGA--HAGEY",
    "V----NVDEV",
    "VEA--DVAGH",
    "VKG------D",
    "VYS--TYETS",
    "FNA--NIPKH",
    "IAGADNGAGY"
]
states, transitions, emission_probabilities, background_probabilities = calculate_hmm_parameters_extended(alignments)
print("States:",states)
print("Transitions",transitions)
print("Emmission Probabilities:",emission_probabilities)
print("Backgroud Probabilities:",background_probabilities)

