% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:59
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRRR10_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR10_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR10_jacobia_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:59:17
% EndTime: 2019-02-26 21:59:17
% DurationCPUTime: 0.25s
% Computational Cost: add. (1003->39), mult. (1383->90), div. (90->9), fcn. (1970->13), ass. (0->58)
t103 = cos(qJ(1));
t98 = sin(pkin(6));
t110 = t103 * t98;
t100 = sin(qJ(2));
t107 = t103 * t100;
t101 = sin(qJ(1));
t102 = cos(qJ(2));
t108 = t101 * t102;
t99 = cos(pkin(6));
t87 = t99 * t107 + t108;
t96 = pkin(12) + qJ(4);
t92 = sin(t96);
t93 = cos(t96);
t76 = t93 * t110 + t87 * t92;
t113 = t100 * t98;
t84 = t92 * t113 - t99 * t93;
t75 = atan2(-t76, t84);
t72 = sin(t75);
t73 = cos(t75);
t66 = -t72 * t76 + t73 * t84;
t65 = 0.1e1 / t66 ^ 2;
t112 = t101 * t98;
t106 = t103 * t102;
t109 = t101 * t100;
t89 = -t99 * t109 + t106;
t80 = -t93 * t112 + t89 * t92;
t120 = t65 * t80;
t88 = t99 * t108 + t107;
t97 = qJ(5) + qJ(6);
t94 = sin(t97);
t115 = t88 * t94;
t81 = t92 * t112 + t89 * t93;
t95 = cos(t97);
t71 = t81 * t95 + t115;
t69 = 0.1e1 / t71 ^ 2;
t114 = t88 * t95;
t70 = t81 * t94 - t114;
t119 = t69 * t70;
t118 = t73 * t76;
t83 = 0.1e1 / t84 ^ 2;
t117 = t76 * t83;
t116 = t80 ^ 2 * t65;
t111 = t102 * t98;
t105 = t70 ^ 2 * t69 + 0.1e1;
t78 = -t92 * t110 + t87 * t93;
t104 = -t72 * t84 - t118;
t86 = t99 * t106 - t109;
t85 = t93 * t113 + t99 * t92;
t82 = 0.1e1 / t84;
t74 = 0.1e1 / (t76 ^ 2 * t83 + 0.1e1);
t68 = 0.1e1 / t71;
t67 = 0.1e1 / t105;
t64 = 0.1e1 / t66;
t63 = 0.1e1 / (0.1e1 + t116);
t62 = (t111 * t117 - t82 * t86) * t92 * t74;
t61 = (t85 * t117 - t78 * t82) * t74;
t60 = t105 * t67;
t1 = [-t80 * t82 * t74, t62, 0, t61, 0, 0; (-t76 * t64 - (-t72 + (t82 * t118 + t72) * t74) * t116) * t63 (-t88 * t92 * t64 - ((t73 * t111 - t72 * t86) * t92 + t104 * t62) * t120) * t63, 0 (t81 * t64 - (t104 * t61 - t72 * t78 + t73 * t85) * t120) * t63, 0, 0; ((-t78 * t94 - t86 * t95) * t68 - (-t78 * t95 + t86 * t94) * t119) * t67 ((-t93 * t115 - t89 * t95) * t68 - (-t93 * t114 + t89 * t94) * t119) * t67, 0 (t95 * t119 - t94 * t68) * t80 * t67, t60, t60;];
Ja_rot  = t1;
