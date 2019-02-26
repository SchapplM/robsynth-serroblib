% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:40
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRPR6_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR6_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR6_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:40:37
% EndTime: 2019-02-26 21:40:37
% DurationCPUTime: 0.27s
% Computational Cost: add. (867->39), mult. (2363->96), div. (85->9), fcn. (3329->15), ass. (0->57)
t100 = sin(qJ(4));
t104 = cos(qJ(4));
t106 = cos(qJ(1));
t96 = sin(pkin(6));
t113 = t106 * t96;
t102 = sin(qJ(1));
t101 = sin(qJ(2));
t105 = cos(qJ(2));
t95 = sin(pkin(11));
t97 = cos(pkin(11));
t109 = t101 * t97 + t105 * t95;
t98 = cos(pkin(6));
t90 = t109 * t98;
t91 = t101 * t95 - t105 * t97;
t79 = -t102 * t91 + t106 * t90;
t73 = -t100 * t113 + t79 * t104;
t89 = t109 * t96;
t86 = t98 * t100 + t89 * t104;
t71 = atan2(-t73, t86);
t68 = sin(t71);
t69 = cos(t71);
t62 = -t68 * t73 + t69 * t86;
t61 = 0.1e1 / t62 ^ 2;
t108 = -t102 * t90 - t106 * t91;
t114 = t102 * t96;
t77 = t100 * t114 + t104 * t108;
t120 = t61 * t77;
t103 = cos(qJ(6));
t107 = t91 * t98;
t81 = t102 * t107 - t106 * t109;
t112 = t81 * t103;
t76 = t100 * t108 - t104 * t114;
t99 = sin(qJ(6));
t67 = t76 * t99 - t112;
t65 = 0.1e1 / t67 ^ 2;
t115 = t81 * t99;
t66 = -t76 * t103 - t115;
t119 = t65 * t66;
t118 = t69 * t73;
t84 = 0.1e1 / t86 ^ 2;
t117 = t73 * t84;
t116 = t77 ^ 2 * t61;
t111 = t66 ^ 2 * t65 + 0.1e1;
t110 = -t68 * t86 - t118;
t72 = t79 * t100 + t104 * t113;
t88 = t91 * t96;
t85 = -t89 * t100 + t98 * t104;
t83 = 0.1e1 / t86;
t78 = -t102 * t109 - t106 * t107;
t70 = 0.1e1 / (t73 ^ 2 * t84 + 0.1e1);
t64 = 0.1e1 / t67;
t63 = 0.1e1 / t111;
t60 = 0.1e1 / t62;
t59 = 0.1e1 / (0.1e1 + t116);
t58 = (-t88 * t117 - t78 * t83) * t70 * t104;
t57 = (t85 * t117 + t72 * t83) * t70;
t1 = [-t77 * t83 * t70, t58, 0, t57, 0, 0; (-t73 * t60 - (-t68 + (t83 * t118 + t68) * t70) * t116) * t59 (t81 * t104 * t60 - (t110 * t58 + (-t68 * t78 - t69 * t88) * t104) * t120) * t59, 0 (-t76 * t60 - (t110 * t57 + t68 * t72 + t69 * t85) * t120) * t59, 0, 0; ((t103 * t72 + t78 * t99) * t64 - (t78 * t103 - t72 * t99) * t119) * t63 ((-t100 * t112 + t108 * t99) * t64 - (t100 * t115 + t103 * t108) * t119) * t63, 0 (-t103 * t64 - t99 * t119) * t77 * t63, 0, t111 * t63;];
Ja_rot  = t1;
