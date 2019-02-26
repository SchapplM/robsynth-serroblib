% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRRR10
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:20
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRRR10_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR10_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR10_jacobia_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:19:49
% EndTime: 2019-02-26 21:19:50
% DurationCPUTime: 0.20s
% Computational Cost: add. (673->40), mult. (1815->87), div. (60->9), fcn. (2489->15), ass. (0->59)
t105 = sin(pkin(7));
t108 = cos(pkin(7));
t109 = cos(pkin(6));
t107 = cos(pkin(13));
t131 = sin(qJ(1));
t116 = t131 * t107;
t104 = sin(pkin(13));
t112 = cos(qJ(1));
t122 = t112 * t104;
t114 = t109 * t116 + t122;
t106 = sin(pkin(6));
t118 = t106 * t131;
t132 = -t105 * t118 + t114 * t108;
t103 = qJ(4) + qJ(5);
t101 = sin(t103);
t102 = cos(t103);
t110 = sin(qJ(3));
t111 = cos(qJ(3));
t117 = t131 * t104;
t121 = t112 * t107;
t98 = -t109 * t117 + t121;
t87 = -t132 * t110 + t98 * t111;
t93 = t114 * t105 + t108 * t118;
t77 = t93 * t101 + t87 * t102;
t75 = 0.1e1 / t77 ^ 2;
t76 = t87 * t101 - t93 * t102;
t130 = t75 * t76;
t125 = t106 * t112;
t119 = t105 * t125;
t123 = t108 * t111;
t97 = t109 * t122 + t116;
t127 = t97 * t110;
t96 = -t109 * t121 + t117;
t82 = t111 * t119 + t96 * t123 + t127;
t126 = t105 * t109;
t90 = -t111 * t126 + (t104 * t110 - t107 * t123) * t106;
t81 = atan2(-t82, t90);
t79 = cos(t81);
t129 = t79 * t82;
t78 = sin(t81);
t72 = -t78 * t82 + t79 * t90;
t71 = 0.1e1 / t72 ^ 2;
t86 = t98 * t110 + t132 * t111;
t128 = t86 ^ 2 * t71;
t124 = t108 * t110;
t120 = t76 ^ 2 * t75 + 0.1e1;
t85 = t110 * t119 - t97 * t111 + t96 * t124;
t92 = -t96 * t105 + t108 * t125;
t91 = t110 * t126 + (t104 * t111 + t107 * t124) * t106;
t89 = 0.1e1 / t90 ^ 2;
t88 = 0.1e1 / t90;
t80 = 0.1e1 / (t82 ^ 2 * t89 + 0.1e1);
t74 = 0.1e1 / t77;
t73 = 0.1e1 / t120;
t70 = 0.1e1 / t72;
t69 = 0.1e1 / (0.1e1 + t128);
t68 = (t82 * t89 * t91 + t85 * t88) * t80;
t67 = t120 * t73;
t1 = [-t86 * t88 * t80, 0, t68, 0, 0, 0; ((-t127 + (-t108 * t96 - t119) * t111) * t70 - (-t78 + (t88 * t129 + t78) * t80) * t128) * t69, 0 (t87 * t70 - (t78 * t85 + t79 * t91 + (-t78 * t90 - t129) * t68) * t86 * t71) * t69, 0, 0, 0; ((t85 * t101 - t92 * t102) * t74 - (t92 * t101 + t85 * t102) * t130) * t73, 0 (-t101 * t74 + t102 * t130) * t86 * t73, t67, t67, 0;];
Ja_rot  = t1;
