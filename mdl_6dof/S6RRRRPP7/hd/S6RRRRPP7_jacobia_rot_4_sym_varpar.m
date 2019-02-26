% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRRPP7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:29
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRPP7_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP7_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPP7_jacobia_rot_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:29:04
% EndTime: 2019-02-26 22:29:04
% DurationCPUTime: 0.21s
% Computational Cost: add. (459->37), mult. (1300->91), div. (85->9), fcn. (1858->13), ass. (0->54)
t76 = cos(pkin(6));
t79 = sin(qJ(2));
t84 = cos(qJ(1));
t88 = t84 * t79;
t80 = sin(qJ(1));
t83 = cos(qJ(2));
t89 = t80 * t83;
t70 = t76 * t88 + t89;
t78 = sin(qJ(3));
t82 = cos(qJ(3));
t75 = sin(pkin(6));
t91 = t75 * t84;
t59 = t70 * t78 + t82 * t91;
t94 = t75 * t78;
t67 = -t76 * t82 + t79 * t94;
t58 = atan2(-t59, t67);
t55 = sin(t58);
t56 = cos(t58);
t49 = -t55 * t59 + t56 * t67;
t48 = 0.1e1 / t49 ^ 2;
t87 = t84 * t83;
t90 = t80 * t79;
t72 = -t76 * t90 + t87;
t93 = t75 * t82;
t63 = t72 * t78 - t80 * t93;
t100 = t48 * t63;
t64 = t72 * t82 + t80 * t94;
t71 = t76 * t89 + t88;
t77 = sin(qJ(4));
t81 = cos(qJ(4));
t54 = t64 * t81 + t71 * t77;
t52 = 0.1e1 / t54 ^ 2;
t53 = t64 * t77 - t71 * t81;
t99 = t52 * t53;
t98 = t56 * t59;
t66 = 0.1e1 / t67 ^ 2;
t97 = t59 * t66;
t96 = t63 ^ 2 * t48;
t95 = t71 * t82;
t92 = t75 * t83;
t86 = t53 ^ 2 * t52 + 0.1e1;
t61 = t70 * t82 - t78 * t91;
t85 = -t55 * t67 - t98;
t69 = t76 * t87 - t90;
t68 = t76 * t78 + t79 * t93;
t65 = 0.1e1 / t67;
t57 = 0.1e1 / (t59 ^ 2 * t66 + 0.1e1);
t51 = 0.1e1 / t54;
t50 = 0.1e1 / t86;
t47 = 0.1e1 / t49;
t46 = 0.1e1 / (0.1e1 + t96);
t45 = (-t65 * t69 + t92 * t97) * t78 * t57;
t44 = (-t61 * t65 + t68 * t97) * t57;
t1 = [-t63 * t65 * t57, t45, t44, 0, 0, 0; (-t59 * t47 - (-t55 + (t65 * t98 + t55) * t57) * t96) * t46 (-t71 * t78 * t47 - ((-t55 * t69 + t56 * t92) * t78 + t85 * t45) * t100) * t46 (t64 * t47 - (t85 * t44 - t55 * t61 + t56 * t68) * t100) * t46, 0, 0, 0; ((-t61 * t77 - t69 * t81) * t51 - (-t61 * t81 + t69 * t77) * t99) * t50 ((-t72 * t81 - t77 * t95) * t51 - (t72 * t77 - t81 * t95) * t99) * t50 (-t77 * t51 + t81 * t99) * t63 * t50, t86 * t50, 0, 0;];
Ja_rot  = t1;
