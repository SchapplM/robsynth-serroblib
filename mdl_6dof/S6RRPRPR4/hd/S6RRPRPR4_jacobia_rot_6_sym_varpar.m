% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRPR4_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR4_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR4_jacobia_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:39:32
% EndTime: 2019-02-26 21:39:32
% DurationCPUTime: 0.28s
% Computational Cost: add. (1286->40), mult. (2363->96), div. (85->9), fcn. (3329->15), ass. (0->58)
t107 = qJ(4) + pkin(12);
t105 = sin(t107);
t106 = cos(t107);
t109 = sin(pkin(6));
t117 = cos(qJ(1));
t123 = t109 * t117;
t108 = sin(pkin(11));
t110 = cos(pkin(11));
t113 = sin(qJ(2));
t116 = cos(qJ(2));
t100 = t113 * t108 - t116 * t110;
t114 = sin(qJ(1));
t111 = cos(pkin(6));
t119 = t116 * t108 + t113 * t110;
t99 = t119 * t111;
t88 = -t114 * t100 + t117 * t99;
t81 = t88 * t105 + t106 * t123;
t98 = t119 * t109;
t94 = t98 * t105 - t111 * t106;
t80 = atan2(-t81, t94);
t77 = sin(t80);
t78 = cos(t80);
t71 = -t77 * t81 + t78 * t94;
t70 = 0.1e1 / t71 ^ 2;
t120 = -t117 * t100 - t114 * t99;
t124 = t109 * t114;
t85 = t105 * t120 - t106 * t124;
t131 = t70 * t85;
t115 = cos(qJ(6));
t112 = sin(qJ(6));
t118 = t100 * t111;
t90 = t114 * t118 - t117 * t119;
t126 = t90 * t112;
t86 = t105 * t124 + t106 * t120;
t76 = t86 * t115 - t126;
t74 = 0.1e1 / t76 ^ 2;
t125 = t90 * t115;
t75 = t86 * t112 + t125;
t130 = t74 * t75;
t129 = t78 * t81;
t93 = 0.1e1 / t94 ^ 2;
t128 = t81 * t93;
t127 = t85 ^ 2 * t70;
t122 = t75 ^ 2 * t74 + 0.1e1;
t83 = -t105 * t123 + t88 * t106;
t121 = -t77 * t94 - t129;
t97 = t100 * t109;
t95 = t111 * t105 + t98 * t106;
t92 = 0.1e1 / t94;
t87 = -t114 * t119 - t117 * t118;
t79 = 0.1e1 / (t81 ^ 2 * t93 + 0.1e1);
t73 = 0.1e1 / t76;
t72 = 0.1e1 / t122;
t69 = 0.1e1 / t71;
t68 = 0.1e1 / (0.1e1 + t127);
t67 = (-t97 * t128 - t87 * t92) * t79 * t105;
t66 = (t95 * t128 - t83 * t92) * t79;
t1 = [-t85 * t92 * t79, t67, 0, t66, 0, 0; (-t81 * t69 - (-t77 + (t92 * t129 + t77) * t79) * t127) * t68 (t90 * t105 * t69 - (t121 * t67 + (-t77 * t87 - t78 * t97) * t105) * t131) * t68, 0 (t86 * t69 - (t121 * t66 - t77 * t83 + t78 * t95) * t131) * t68, 0, 0; ((-t112 * t83 - t87 * t115) * t73 - (t87 * t112 - t115 * t83) * t130) * t72 ((t106 * t126 - t115 * t120) * t73 - (t106 * t125 + t112 * t120) * t130) * t72, 0 (-t112 * t73 + t115 * t130) * t85 * t72, 0, t122 * t72;];
Ja_rot  = t1;
