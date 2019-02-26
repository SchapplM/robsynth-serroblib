% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
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

function Ja_rot = S6RRRRPP7_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP7_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPP7_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:28:55
% EndTime: 2019-02-26 22:28:55
% DurationCPUTime: 0.34s
% Computational Cost: add. (1427->45), mult. (2519->109), div. (107->9), fcn. (3557->13), ass. (0->64)
t109 = qJ(4) + pkin(11);
t107 = sin(t109);
t108 = cos(t109);
t111 = cos(pkin(6));
t116 = cos(qJ(2));
t117 = cos(qJ(1));
t122 = t117 * t116;
t113 = sin(qJ(2));
t114 = sin(qJ(1));
t125 = t114 * t113;
t121 = -t111 * t122 + t125;
t123 = t117 * t113;
t124 = t114 * t116;
t103 = t111 * t123 + t124;
t112 = sin(qJ(3));
t115 = cos(qJ(3));
t110 = sin(pkin(6));
t126 = t110 * t117;
t95 = -t103 * t115 + t112 * t126;
t140 = t95 * t107 + t121 * t108;
t119 = t121 * t107;
t139 = t95 * t108 - t119;
t128 = t110 * t115;
t102 = t111 * t112 + t113 * t128;
t127 = t110 * t116;
t90 = t102 * t107 + t108 * t127;
t78 = atan2(t140, t90);
t75 = sin(t78);
t76 = cos(t78);
t74 = t140 * t75 + t76 * t90;
t73 = 0.1e1 / t74 ^ 2;
t104 = t111 * t124 + t123;
t131 = t104 * t108;
t105 = -t111 * t125 + t122;
t129 = t110 * t112;
t97 = t105 * t115 + t114 * t129;
t85 = t97 * t107 - t131;
t137 = t73 * t85;
t136 = t76 * t140;
t86 = t104 * t107 + t97 * t108;
t81 = 0.1e1 / t86 ^ 2;
t96 = -t105 * t112 + t114 * t128;
t135 = t81 * t96;
t89 = 0.1e1 / t90 ^ 2;
t134 = t140 * t89;
t133 = t85 ^ 2 * t73;
t132 = t96 ^ 2 * t81;
t130 = t107 * t115;
t120 = -t75 * t90 + t136;
t118 = t103 * t112 + t115 * t126;
t101 = t111 * t115 - t113 * t129;
t98 = (-t108 * t113 + t116 * t130) * t110;
t91 = t102 * t108 - t107 * t127;
t88 = 0.1e1 / t90;
t87 = -t103 * t108 - t115 * t119;
t80 = 0.1e1 / t86;
t79 = 0.1e1 / (0.1e1 + t132);
t77 = 0.1e1 / (t140 ^ 2 * t89 + 0.1e1);
t72 = 0.1e1 / t74;
t71 = 0.1e1 / (0.1e1 + t133);
t70 = (-t101 * t134 + t118 * t88) * t77 * t107;
t69 = (-t98 * t134 - t87 * t88) * t77;
t68 = (-t91 * t134 + t139 * t88) * t77;
t1 = [-t85 * t88 * t77, t69, t70, t68, 0, 0; (t140 * t72 - (-t75 + (-t88 * t136 + t75) * t77) * t133) * t71 ((-t104 * t130 - t105 * t108) * t72 - (t120 * t69 - t75 * t87 + t76 * t98) * t137) * t71 (t96 * t107 * t72 - (t120 * t70 + (t101 * t76 + t118 * t75) * t107) * t137) * t71 (t86 * t72 - (t120 * t68 + t139 * t75 + t76 * t91) * t137) * t71, 0, 0; (t118 * t80 - t139 * t135) * t79 (t104 * t112 * t80 - (t105 * t107 - t115 * t131) * t135) * t79 (-t108 * t132 - t80 * t97) * t79, t85 * t79 * t135, 0, 0;];
Ja_rot  = t1;
