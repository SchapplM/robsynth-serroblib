% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRRPRR4
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
%
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRPRR4_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_jacobia_rot_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:05:52
% EndTime: 2019-02-26 20:05:52
% DurationCPUTime: 0.09s
% Computational Cost: add. (256->25), mult. (739->65), div. (57->9), fcn. (1061->11), ass. (0->42)
t55 = sin(pkin(11));
t57 = cos(pkin(11));
t62 = cos(qJ(2));
t58 = cos(pkin(6));
t60 = sin(qJ(2));
t65 = t58 * t60;
t50 = t55 * t62 + t57 * t65;
t59 = sin(qJ(3));
t56 = sin(pkin(6));
t61 = cos(qJ(3));
t67 = t56 * t61;
t41 = t50 * t59 + t57 * t67;
t68 = t56 * t59;
t53 = -t58 * t61 + t60 * t68;
t39 = atan2(-t41, t53);
t35 = sin(t39);
t36 = cos(t39);
t34 = -t35 * t41 + t36 * t53;
t33 = 0.1e1 / t34 ^ 2;
t52 = -t55 * t65 + t57 * t62;
t44 = t52 * t59 - t55 * t67;
t71 = t33 * t44;
t48 = 0.1e1 / t53 ^ 2;
t70 = t41 * t48;
t45 = t52 * t61 + t55 * t68;
t40 = 0.1e1 / t45 ^ 2;
t64 = t58 * t62;
t51 = -t55 * t64 - t57 * t60;
t69 = t51 ^ 2 * t40;
t66 = t56 * t62;
t63 = -t35 * t53 - t36 * t41;
t54 = t58 * t59 + t60 * t67;
t49 = -t55 * t60 + t57 * t64;
t47 = 0.1e1 / t53;
t43 = t50 * t61 - t57 * t68;
t38 = 0.1e1 / (t41 ^ 2 * t48 + 0.1e1);
t37 = 0.1e1 / (0.1e1 + t69);
t32 = 0.1e1 / t34;
t31 = 0.1e1 / (t44 ^ 2 * t33 + 0.1e1);
t30 = (-t47 * t49 + t66 * t70) * t59 * t38;
t29 = (-t43 * t47 + t54 * t70) * t38;
t1 = [0, t30, t29, 0, 0, 0; 0 (t51 * t59 * t32 - ((-t35 * t49 + t36 * t66) * t59 + t63 * t30) * t71) * t31 (t45 * t32 - (t63 * t29 - t35 * t43 + t36 * t54) * t71) * t31, 0, 0, 0; 0 (-t52 / t45 - t61 * t69) * t37, t44 * t51 * t40 * t37, 0, 0, 0;];
Ja_rot  = t1;
