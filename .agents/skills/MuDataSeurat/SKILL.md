```markdown
# MuDataSeurat Development Patterns

> Auto-generated skill from repository analysis

## Overview
This skill teaches you the core development patterns and conventions used in the MuDataSeurat TypeScript codebase. You'll learn how to structure files, write imports and exports, follow commit message conventions, and write and run tests. The guide also provides suggested commands for common workflows.

## Coding Conventions

### File Naming
- Use **PascalCase** for all file names.
  - Example: `MyComponent.ts`, `DataParser.ts`

### Import Style
- Use **relative imports** for referencing other modules/files.
  - Example:
    ```typescript
    import { MyFunction } from './MyFunction';
    ```

### Export Style
- Use **named exports** rather than default exports.
  - Example:
    ```typescript
    // In MyFunction.ts
    export function MyFunction() { ... }
    ```

### Commit Messages
- Follow the **conventional commit** style.
- Use the `fix` prefix for bug fixes.
- Keep commit message length around 58 characters.
  - Example:
    ```
    fix: correct data parsing error in DataParser
    ```

## Workflows

### Code Contribution
**Trigger:** When adding or updating code in the repository  
**Command:** `/contribute`

1. Create a new file using PascalCase naming.
2. Write code using named exports and relative imports.
3. Write or update corresponding test files (`*.test.ts`).
4. Commit changes using a conventional commit message (e.g., `fix: ...`).
5. Open a pull request for review.

### Testing Code
**Trigger:** When verifying code correctness  
**Command:** `/test`

1. Ensure test files are named with the `.test.ts` suffix.
2. Run the test suite using your preferred TypeScript test runner (framework not specified; common options include Jest or Mocha).
3. Review test results and address any failures.

## Testing Patterns

- Test files are named with the `.test.ts` suffix and are located alongside or near the code they test.
- The specific testing framework is not specified; use a standard TypeScript-compatible test runner (e.g., Jest, Mocha).
- Tests should cover core functionality and edge cases.

## Commands
| Command      | Purpose                                         |
|--------------|-------------------------------------------------|
| /contribute  | Steps to add or update code in the repository   |
| /test        | Steps to run and verify tests                   |
```