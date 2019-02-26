% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPPR3
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:59
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRPPR3_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:59:23
% EndTime: 2019-02-26 19:59:24
% DurationCPUTime: 0.15s
% Computational Cost: add. (334->30), mult. (949->78), div. (65->9), fcn. (1349->13), ass. (0->50)
t63 = sin(pkin(10));
t65 = cos(pkin(10));
t72 = cos(qJ(2));
t66 = cos(pkin(6));
t69 = sin(qJ(2));
t76 = t66 * t69;
t56 = t63 * t72 + t65 * t76;
t71 = cos(qJ(3));
t64 = sin(pkin(6));
t68 = sin(qJ(3));
t79 = t64 * t68;
t49 = t56 * t71 - t65 * t79;
t78 = t64 * t71;
t60 = t66 * t68 + t69 * t78;
t47 = atan2(-t49, t60);
t44 = sin(t47);
t45 = cos(t47);
t38 = -t44 * t49 + t45 * t60;
t37 = 0.1e1 / t38 ^ 2;
t58 = -t63 * t76 + t65 * t72;
t52 = t58 * t71 + t63 * t79;
t84 = t37 * t52;
t51 = t58 * t68 - t63 * t78;
t70 = cos(qJ(6));
t75 = t66 * t72;
t57 = -t63 * t75 - t65 * t69;
t67 = sin(qJ(6));
t81 = t57 * t67;
t43 = t51 * t70 + t81;
t41 = 0.1e1 / t43 ^ 2;
t80 = t57 * t70;
t42 = t51 * t67 - t80;
t83 = t41 * t42;
t54 = 0.1e1 / t60 ^ 2;
t82 = t49 * t54;
t77 = t64 * t72;
t74 = t42 ^ 2 * t41 + 0.1e1;
t73 = -t44 * t60 - t45 * t49;
t59 = t66 * t71 - t69 * t79;
t55 = -t63 * t69 + t65 * t75;
t53 = 0.1e1 / t60;
t48 = t56 * t68 + t65 * t78;
t46 = 0.1e1 / (t49 ^ 2 * t54 + 0.1e1);
t40 = 0.1e1 / t43;
t39 = 0.1e1 / t74;
t36 = 0.1e1 / t38;
t35 = 0.1e1 / (t52 ^ 2 * t37 + 0.1e1);
t34 = (-t53 * t55 + t77 * t82) * t71 * t46;
t33 = (t48 * t53 + t59 * t82) * t46;
t1 = [0, t34, t33, 0, 0, 0; 0 (t57 * t71 * t36 - ((-t44 * t55 + t45 * t77) * t71 + t73 * t34) * t84) * t35 (-t51 * t36 - (t73 * t33 + t44 * t48 + t45 * t59) * t84) * t35, 0, 0, 0; 0 ((t58 * t70 + t68 * t81) * t40 - (-t58 * t67 + t68 * t80) * t83) * t39 (t40 * t67 - t70 * t83) * t52 * t39, 0, 0, t74 * t39;];
Ja_rot  = t1;
