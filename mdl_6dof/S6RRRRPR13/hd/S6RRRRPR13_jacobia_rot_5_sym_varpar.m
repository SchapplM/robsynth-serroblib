% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:37
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRPR13_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR13_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR13_jacobia_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:37:35
% EndTime: 2019-02-26 22:37:36
% DurationCPUTime: 0.32s
% Computational Cost: add. (905->44), mult. (2519->109), div. (107->9), fcn. (3557->13), ass. (0->63)
t103 = sin(qJ(4));
t107 = cos(qJ(4));
t102 = cos(pkin(6));
t109 = cos(qJ(2));
t110 = cos(qJ(1));
t115 = t110 * t109;
t105 = sin(qJ(2));
t106 = sin(qJ(1));
t118 = t106 * t105;
t114 = -t102 * t115 + t118;
t104 = sin(qJ(3));
t108 = cos(qJ(3));
t101 = sin(pkin(6));
t120 = t101 * t110;
t116 = t110 * t105;
t117 = t106 * t109;
t97 = t102 * t116 + t117;
t89 = t104 * t120 - t97 * t108;
t133 = t89 * t103 + t114 * t107;
t112 = t114 * t103;
t132 = t89 * t107 - t112;
t121 = t101 * t108;
t96 = t102 * t104 + t105 * t121;
t85 = t101 * t109 * t107 + t96 * t103;
t73 = atan2(t133, t85);
t69 = sin(t73);
t70 = cos(t73);
t68 = t133 * t69 + t70 * t85;
t67 = 0.1e1 / t68 ^ 2;
t98 = t102 * t117 + t116;
t123 = t98 * t107;
t122 = t101 * t104;
t99 = -t102 * t118 + t115;
t91 = t106 * t122 + t99 * t108;
t79 = t91 * t103 - t123;
t130 = t67 * t79;
t129 = t70 * t133;
t124 = t98 * t103;
t80 = t91 * t107 + t124;
t75 = 0.1e1 / t80 ^ 2;
t90 = -t99 * t104 + t106 * t121;
t128 = t75 * t90;
t83 = 0.1e1 / t85 ^ 2;
t127 = t133 * t83;
t126 = t79 ^ 2 * t67;
t125 = t90 ^ 2 * t75;
t119 = t103 * t109;
t113 = -t69 * t85 + t129;
t111 = t97 * t104 + t108 * t120;
t95 = t102 * t108 - t105 * t122;
t92 = (-t105 * t107 + t108 * t119) * t101;
t86 = -t101 * t119 + t96 * t107;
t82 = 0.1e1 / t85;
t81 = -t97 * t107 - t108 * t112;
t74 = 0.1e1 / t80;
t72 = 0.1e1 / (0.1e1 + t125);
t71 = 0.1e1 / (t133 ^ 2 * t83 + 0.1e1);
t66 = 0.1e1 / t68;
t65 = 0.1e1 / (0.1e1 + t126);
t64 = (t111 * t82 - t95 * t127) * t71 * t103;
t63 = (-t92 * t127 - t81 * t82) * t71;
t62 = (-t86 * t127 + t132 * t82) * t71;
t1 = [-t79 * t82 * t71, t63, t64, t62, 0, 0; (t133 * t66 - (-t69 + (-t82 * t129 + t69) * t71) * t126) * t65 ((-t99 * t107 - t108 * t124) * t66 - (t113 * t63 - t69 * t81 + t70 * t92) * t130) * t65 (t90 * t103 * t66 - (t113 * t64 + (t111 * t69 + t70 * t95) * t103) * t130) * t65 (t80 * t66 - (t113 * t62 + t132 * t69 + t70 * t86) * t130) * t65, 0, 0; (t111 * t74 - t132 * t128) * t72 (t98 * t104 * t74 - (t99 * t103 - t108 * t123) * t128) * t72 (-t107 * t125 - t74 * t91) * t72, t79 * t72 * t128, 0, 0;];
Ja_rot  = t1;
