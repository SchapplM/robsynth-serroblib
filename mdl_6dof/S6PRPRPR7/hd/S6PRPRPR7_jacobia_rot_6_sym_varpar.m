% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRPR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:49
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPRPR7_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:49:47
% EndTime: 2019-02-26 19:49:47
% DurationCPUTime: 0.14s
% Computational Cost: add. (334->29), mult. (949->80), div. (65->9), fcn. (1349->13), ass. (0->50)
t63 = cos(pkin(10));
t66 = sin(qJ(4));
t61 = sin(pkin(10));
t67 = sin(qJ(2));
t64 = cos(pkin(6));
t70 = cos(qJ(2));
t74 = t64 * t70;
t71 = t61 * t67 - t63 * t74;
t62 = sin(pkin(6));
t69 = cos(qJ(4));
t77 = t62 * t69;
t52 = t63 * t77 - t71 * t66;
t76 = t62 * t70;
t60 = t64 * t69 - t66 * t76;
t48 = atan2(t52, t60);
t45 = sin(t48);
t46 = cos(t48);
t39 = t45 * t52 + t46 * t60;
t38 = 0.1e1 / t39 ^ 2;
t57 = t61 * t74 + t63 * t67;
t50 = t57 * t66 + t61 * t77;
t83 = t38 * t50;
t79 = t62 * t66;
t49 = -t57 * t69 + t61 * t79;
t75 = t64 * t67;
t58 = -t61 * t75 + t63 * t70;
t65 = sin(qJ(6));
t68 = cos(qJ(6));
t44 = t49 * t65 + t58 * t68;
t42 = 0.1e1 / t44 ^ 2;
t43 = -t49 * t68 + t58 * t65;
t82 = t42 * t43;
t55 = 0.1e1 / t60 ^ 2;
t81 = t52 * t55;
t80 = t58 * t69;
t78 = t62 * t67;
t73 = t43 ^ 2 * t42 + 0.1e1;
t72 = -t45 * t60 + t46 * t52;
t59 = -t64 * t66 - t69 * t76;
t56 = t61 * t70 + t63 * t75;
t54 = 0.1e1 / t60;
t51 = t63 * t79 + t71 * t69;
t47 = 0.1e1 / (t52 ^ 2 * t55 + 0.1e1);
t41 = 0.1e1 / t44;
t40 = 0.1e1 / t73;
t37 = 0.1e1 / t39;
t36 = 0.1e1 / (t50 ^ 2 * t38 + 0.1e1);
t35 = (-t54 * t56 - t78 * t81) * t66 * t47;
t34 = (-t51 * t54 - t59 * t81) * t47;
t1 = [0, t35, 0, t34, 0, 0; 0 (t58 * t66 * t37 - ((-t45 * t56 + t46 * t78) * t66 + t72 * t35) * t83) * t36, 0 (-t49 * t37 - (t72 * t34 - t45 * t51 + t46 * t59) * t83) * t36, 0, 0; 0 ((-t57 * t65 + t68 * t80) * t41 - (-t57 * t68 - t65 * t80) * t82) * t40, 0 (-t41 * t68 - t65 * t82) * t50 * t40, 0, t73 * t40;];
Ja_rot  = t1;
