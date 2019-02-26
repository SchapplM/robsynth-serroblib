% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPPR5
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:05
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPPR5_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR5_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPPR5_jacobia_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:05:45
% EndTime: 2019-02-26 22:05:46
% DurationCPUTime: 0.22s
% Computational Cost: add. (833->38), mult. (1217->89), div. (80->9), fcn. (1746->13), ass. (0->55)
t81 = cos(pkin(6));
t82 = sin(qJ(2));
t85 = cos(qJ(1));
t88 = t85 * t82;
t83 = sin(qJ(1));
t84 = cos(qJ(2));
t89 = t83 * t84;
t70 = t81 * t88 + t89;
t77 = qJ(3) + pkin(11);
t75 = sin(t77);
t76 = cos(t77);
t79 = sin(pkin(6));
t91 = t79 * t85;
t59 = t70 * t75 + t76 * t91;
t94 = t79 * t82;
t67 = t75 * t94 - t81 * t76;
t58 = atan2(-t59, t67);
t53 = sin(t58);
t54 = cos(t58);
t49 = -t53 * t59 + t54 * t67;
t48 = 0.1e1 / t49 ^ 2;
t87 = t85 * t84;
t90 = t83 * t82;
t72 = -t81 * t90 + t87;
t93 = t79 * t83;
t63 = t72 * t75 - t76 * t93;
t101 = t48 * t63;
t64 = t72 * t76 + t75 * t93;
t80 = cos(pkin(12));
t71 = t81 * t89 + t88;
t78 = sin(pkin(12));
t96 = t71 * t78;
t56 = t64 * t80 + t96;
t52 = 0.1e1 / t56 ^ 2;
t95 = t71 * t80;
t55 = t64 * t78 - t95;
t100 = t52 * t55;
t99 = t54 * t59;
t66 = 0.1e1 / t67 ^ 2;
t98 = t59 * t66;
t97 = t63 ^ 2 * t48;
t92 = t79 * t84;
t61 = t70 * t76 - t75 * t91;
t86 = -t53 * t67 - t99;
t69 = t81 * t87 - t90;
t68 = t81 * t75 + t76 * t94;
t65 = 0.1e1 / t67;
t57 = 0.1e1 / (t59 ^ 2 * t66 + 0.1e1);
t51 = 0.1e1 / t56;
t50 = 0.1e1 / (t55 ^ 2 * t52 + 0.1e1);
t47 = 0.1e1 / t49;
t46 = 0.1e1 / (0.1e1 + t97);
t45 = (-t65 * t69 + t92 * t98) * t75 * t57;
t44 = (-t61 * t65 + t68 * t98) * t57;
t1 = [-t63 * t65 * t57, t45, t44, 0, 0, 0; (-t59 * t47 - (-t53 + (t65 * t99 + t53) * t57) * t97) * t46 (-t71 * t75 * t47 - ((-t53 * t69 + t54 * t92) * t75 + t86 * t45) * t101) * t46 (t64 * t47 - (t86 * t44 - t53 * t61 + t54 * t68) * t101) * t46, 0, 0, 0; ((-t61 * t78 - t69 * t80) * t51 - (-t61 * t80 + t69 * t78) * t100) * t50 ((-t72 * t80 - t76 * t96) * t51 - (t72 * t78 - t76 * t95) * t100) * t50 (t80 * t100 - t78 * t51) * t63 * t50, 0, 0, 0;];
Ja_rot  = t1;
