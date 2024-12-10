original = "48656C6C6F2065766572796F6E652C20496D207372696E69"
decrypted ="dcac31a2e67457d6b5e5bd46b4e23c916dea6efada4f2f673d9a7f82"
decrypted = decrypted.rstrip('0')
# print(decrypted)

if original == decrypted:
    print("Strings match!")
else:
    for i, (original_char, decrypted_char) in enumerate(zip(original, decrypted)):
        if original_char != decrypted_char:
            print("Mismatch at index:", i)
