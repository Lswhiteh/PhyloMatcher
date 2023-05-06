import regex as re

# Example usage
text = "Pear, apple, xpPle, axple, app e, orange"
sp_key = "apple"
pattern = fr"(?e)({sp_key}){{e<=1}}" # levenshtein + best match
matches = re.findall(pattern, text, re.IGNORECASE)
print(matches)