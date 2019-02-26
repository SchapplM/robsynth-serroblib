% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:35
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRPR9_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR9_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR9_jacobia_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:35:13
% EndTime: 2019-02-26 22:35:13
% DurationCPUTime: 0.25s
% Computational Cost: add. (1267->39), mult. (1738->90), div. (115->9), fcn. (2485->13), ass. (0->59)
t106 = qJ(3) + qJ(4);
t103 = sin(t106);
t104 = cos(t106);
t107 = sin(pkin(6));
t112 = cos(qJ(1));
t119 = t107 * t112;
t108 = cos(pkin(6));
t109 = sin(qJ(2));
t116 = t112 * t109;
t110 = sin(qJ(1));
t111 = cos(qJ(2));
t117 = t110 * t111;
t96 = t108 * t116 + t117;
t85 = t96 * t103 + t104 * t119;
t122 = t107 * t109;
t93 = t103 * t122 - t108 * t104;
t84 = atan2(-t85, t93);
t81 = sin(t84);
t82 = cos(t84);
t75 = -t81 * t85 + t82 * t93;
t74 = 0.1e1 / t75 ^ 2;
t121 = t107 * t110;
t115 = t112 * t111;
t118 = t110 * t109;
t98 = -t108 * t118 + t115;
t89 = t98 * t103 - t104 * t121;
t129 = t74 * t89;
t105 = pkin(12) + qJ(6);
t102 = cos(t105);
t101 = sin(t105);
t97 = t108 * t117 + t116;
t124 = t97 * t101;
t90 = t103 * t121 + t98 * t104;
t80 = t90 * t102 + t124;
t78 = 0.1e1 / t80 ^ 2;
t123 = t97 * t102;
t79 = t90 * t101 - t123;
t128 = t78 * t79;
t127 = t82 * t85;
t92 = 0.1e1 / t93 ^ 2;
t126 = t85 * t92;
t125 = t89 ^ 2 * t74;
t120 = t107 * t111;
t114 = t79 ^ 2 * t78 + 0.1e1;
t87 = -t103 * t119 + t96 * t104;
t113 = -t81 * t93 - t127;
t95 = t108 * t115 - t118;
t94 = t108 * t103 + t104 * t122;
t91 = 0.1e1 / t93;
t83 = 0.1e1 / (t85 ^ 2 * t92 + 0.1e1);
t77 = 0.1e1 / t80;
t76 = 0.1e1 / t114;
t73 = 0.1e1 / t75;
t72 = 0.1e1 / (0.1e1 + t125);
t71 = (t120 * t126 - t91 * t95) * t83 * t103;
t70 = (t94 * t126 - t87 * t91) * t83;
t69 = (-t101 * t77 + t102 * t128) * t89 * t76;
t68 = (t90 * t73 - (t113 * t70 - t81 * t87 + t82 * t94) * t129) * t72;
t1 = [-t89 * t91 * t83, t71, t70, t70, 0, 0; (-t85 * t73 - (-t81 + (t91 * t127 + t81) * t83) * t125) * t72 (-t97 * t103 * t73 - (t113 * t71 + (t82 * t120 - t81 * t95) * t103) * t129) * t72, t68, t68, 0, 0; ((-t101 * t87 - t95 * t102) * t77 - (t95 * t101 - t102 * t87) * t128) * t76 ((-t98 * t102 - t104 * t124) * t77 - (t98 * t101 - t104 * t123) * t128) * t76, t69, t69, 0, t114 * t76;];
Ja_rot  = t1;
