% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPPP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha4,d1,d2,d3,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:02
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPPP1_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPP1_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPP1_jacobia_rot_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:02:47
% EndTime: 2019-02-26 22:02:48
% DurationCPUTime: 0.23s
% Computational Cost: add. (423->38), mult. (1241->99), div. (80->9), fcn. (1765->13), ass. (0->55)
t68 = sin(qJ(3));
t67 = cos(pkin(6));
t69 = sin(qJ(2));
t84 = t69 * t67;
t65 = sin(pkin(6));
t72 = cos(qJ(2));
t85 = t65 * t72;
t92 = t68 * t84 + t85;
t71 = cos(qJ(3));
t73 = cos(qJ(1));
t78 = t73 * t71;
t70 = sin(qJ(1));
t80 = t70 * t72;
t58 = t68 * t80 + t78;
t83 = t69 * t70;
t50 = -t58 * t65 - t67 * t83;
t56 = t69 * t68 * t65 - t72 * t67;
t49 = atan2(t50, t56);
t46 = sin(t49);
t47 = cos(t49);
t40 = t46 * t50 + t47 * t56;
t39 = 0.1e1 / t40 ^ 2;
t79 = t73 * t68;
t60 = t70 * t71 - t72 * t79;
t81 = t69 * t73;
t51 = t60 * t65 - t67 * t81;
t91 = t39 * t51;
t61 = t70 * t68 + t72 * t78;
t64 = sin(pkin(10));
t66 = cos(pkin(10));
t74 = t60 * t67 + t65 * t81;
t45 = t61 * t66 + t74 * t64;
t43 = 0.1e1 / t45 ^ 2;
t44 = t61 * t64 - t74 * t66;
t90 = t43 * t44;
t89 = t47 * t50;
t55 = 0.1e1 / t56 ^ 2;
t88 = t50 * t55;
t87 = t51 ^ 2 * t39;
t86 = t61 * t67;
t82 = t69 * t71;
t76 = -t46 * t56 + t89;
t75 = -t58 * t67 + t65 * t83;
t59 = -t71 * t80 + t79;
t57 = t68 * t85 + t84;
t54 = 0.1e1 / t56;
t53 = t56 * t70;
t48 = 0.1e1 / (t50 ^ 2 * t55 + 0.1e1);
t42 = 0.1e1 / t45;
t41 = 0.1e1 / (t44 ^ 2 * t43 + 0.1e1);
t38 = 0.1e1 / t40;
t37 = 0.1e1 / (0.1e1 + t87);
t36 = (t54 * t59 - t82 * t88) * t65 * t48;
t35 = (t53 * t54 - t57 * t88) * t48;
t1 = [t51 * t54 * t48, t35, t36, 0, 0, 0; (t50 * t38 + (t46 + (t54 * t89 - t46) * t48) * t87) * t37 ((t76 * t35 + t46 * t53 + t47 * t57) * t91 - t56 * t38 * t73) * t37 (t61 * t65 * t38 + ((t46 * t59 + t47 * t82) * t65 + t76 * t36) * t91) * t37, 0, 0, 0; ((t59 * t64 + t75 * t66) * t42 - (t59 * t66 - t75 * t64) * t90) * t41 ((-t64 * t82 - t92 * t66) * t42 - (t92 * t64 - t66 * t82) * t90) * t41 * t73 ((t60 * t64 + t66 * t86) * t42 - (t60 * t66 - t64 * t86) * t90) * t41, 0, 0, 0;];
Ja_rot  = t1;
