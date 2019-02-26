% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRRR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:57
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPRRR6_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR6_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR6_jacobia_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:56:50
% EndTime: 2019-02-26 19:56:50
% DurationCPUTime: 0.14s
% Computational Cost: add. (334->29), mult. (949->79), div. (65->9), fcn. (1349->13), ass. (0->51)
t61 = sin(pkin(11));
t63 = cos(pkin(11));
t67 = sin(qJ(2));
t64 = cos(pkin(6));
t70 = cos(qJ(2));
t73 = t64 * t70;
t55 = t61 * t67 - t63 * t73;
t69 = cos(qJ(4));
t62 = sin(pkin(6));
t66 = sin(qJ(4));
t78 = t62 * t66;
t50 = t55 * t69 + t63 * t78;
t75 = t62 * t70;
t59 = t64 * t66 + t69 * t75;
t47 = atan2(t50, t59);
t44 = sin(t47);
t45 = cos(t47);
t38 = t44 * t50 + t45 * t59;
t37 = 0.1e1 / t38 ^ 2;
t57 = t61 * t73 + t63 * t67;
t48 = -t57 * t69 + t61 * t78;
t83 = t37 * t48;
t76 = t62 * t69;
t49 = t57 * t66 + t61 * t76;
t68 = cos(qJ(5));
t74 = t64 * t67;
t58 = -t61 * t74 + t63 * t70;
t65 = sin(qJ(5));
t80 = t58 * t65;
t43 = t49 * t68 + t80;
t41 = 0.1e1 / t43 ^ 2;
t79 = t58 * t68;
t42 = t49 * t65 - t79;
t82 = t41 * t42;
t54 = 0.1e1 / t59 ^ 2;
t81 = t50 * t54;
t77 = t62 * t67;
t72 = t42 ^ 2 * t41 + 0.1e1;
t71 = -t44 * t59 + t45 * t50;
t60 = t64 * t69 - t66 * t75;
t56 = t61 * t70 + t63 * t74;
t53 = 0.1e1 / t59;
t51 = -t55 * t66 + t63 * t76;
t46 = 0.1e1 / (t50 ^ 2 * t54 + 0.1e1);
t40 = 0.1e1 / t43;
t39 = 0.1e1 / t72;
t36 = 0.1e1 / t38;
t35 = 0.1e1 / (t48 ^ 2 * t37 + 0.1e1);
t34 = (t53 * t56 + t77 * t81) * t69 * t46;
t33 = (t51 * t53 - t60 * t81) * t46;
t1 = [0, t34, 0, t33, 0, 0; 0 (-t58 * t69 * t36 - ((t44 * t56 - t45 * t77) * t69 + t71 * t34) * t83) * t35, 0 (t49 * t36 - (t71 * t33 + t44 * t51 + t45 * t60) * t83) * t35, 0, 0; 0 ((t57 * t68 + t66 * t80) * t40 - (-t57 * t65 + t66 * t79) * t82) * t39, 0 (-t40 * t65 + t68 * t82) * t48 * t39, t72 * t39, 0;];
Ja_rot  = t1;
