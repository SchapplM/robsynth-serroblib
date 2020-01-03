% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S4PPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
% 
% Output:
% tauB_reg [6x(5*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S4PPRR3_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR3_invdynB_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR3_invdynB_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR3_invdynB_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR3_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR3_invdynB_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:17:27
% EndTime: 2019-12-31 16:17:30
% DurationCPUTime: 1.22s
% Computational Cost: add. (1727->198), mult. (3273->272), div. (0->0), fcn. (2252->6), ass. (0->128)
t169 = cos(qJ(3));
t171 = qJD(3) ^ 2;
t167 = sin(qJ(3));
t179 = qJDD(3) * t167;
t138 = t169 * t171 + t179;
t178 = t169 * qJDD(3);
t139 = t167 * t171 - t178;
t163 = sin(pkin(6));
t164 = cos(pkin(6));
t101 = t164 * t138 + t139 * t163;
t161 = g(3) - qJDD(1);
t121 = pkin(4) * t139 + t161 * t167;
t172 = pkin(4) * t138 + t161 * t169;
t199 = -qJ(1) * t101 + t121 * t163 + t164 * t172;
t175 = -t138 * t163 + t164 * t139;
t197 = -qJ(1) * t175 + t121 * t164 - t163 * t172;
t142 = g(1) * t163 - t164 * g(2);
t136 = -qJDD(2) + t142;
t143 = g(1) * t164 + g(2) * t163;
t104 = t169 * t136 - t143 * t167;
t105 = -t167 * t136 - t169 * t143;
t74 = t104 * t169 - t105 * t167;
t75 = t104 * t167 + t105 * t169;
t55 = t163 * t75 + t164 * t74;
t196 = t163 * t74 - t164 * t75;
t193 = pkin(1) + pkin(2);
t166 = sin(qJ(4));
t95 = -qJDD(3) * pkin(3) - t171 * pkin(5) + t104;
t192 = t166 * t95;
t168 = cos(qJ(4));
t191 = t168 * t95;
t96 = -pkin(3) * t171 + qJDD(3) * pkin(5) + t105;
t81 = t166 * t161 + t168 * t96;
t150 = t166 * t171 * t168;
t144 = qJDD(4) + t150;
t188 = t144 * t166;
t187 = t144 * t168;
t145 = qJDD(4) - t150;
t186 = t145 * t166;
t185 = t145 * t168;
t159 = t166 ^ 2;
t184 = t159 * t171;
t183 = t163 * t161;
t152 = t164 * t161;
t160 = t168 ^ 2;
t182 = t159 + t160;
t181 = qJD(3) * qJD(4);
t180 = qJDD(3) * t166;
t155 = t168 * qJDD(3);
t177 = t166 * t181;
t176 = t168 * t181;
t80 = -t168 * t161 + t166 * t96;
t128 = t163 * t143;
t98 = t136 * t164 - t128;
t108 = t142 * t164 - t128;
t129 = t164 * t143;
t99 = -t136 * t163 - t129;
t109 = -t142 * t163 - t129;
t174 = t167 * t150;
t173 = t169 * t150;
t59 = t166 * t81 - t168 * t80;
t61 = t166 * t80 + t168 * t81;
t170 = qJD(4) ^ 2;
t157 = t160 * t171;
t149 = -t157 - t170;
t148 = t157 - t170;
t147 = -t170 - t184;
t146 = t170 - t184;
t141 = t157 - t184;
t140 = t157 + t184;
t137 = t182 * qJDD(3);
t135 = t155 - 0.2e1 * t177;
t134 = t155 - t177;
t133 = t176 + t180;
t132 = 0.2e1 * t176 + t180;
t131 = t182 * t181;
t119 = qJDD(4) * t167 + t131 * t169;
t118 = t133 * t168 - t159 * t181;
t117 = qJDD(4) * t169 - t131 * t167;
t116 = -t134 * t166 - t160 * t181;
t115 = -t147 * t166 - t185;
t114 = -t146 * t166 + t187;
t113 = t149 * t168 - t188;
t112 = t148 * t168 - t186;
t111 = t147 * t168 - t186;
t110 = t149 * t166 + t187;
t107 = t137 * t169 - t140 * t167;
t106 = t137 * t167 + t140 * t169;
t97 = -t132 * t166 + t135 * t168;
t94 = t118 * t169 - t174;
t93 = t116 * t169 + t174;
t92 = -t118 * t167 - t173;
t91 = -t116 * t167 + t173;
t90 = t114 * t169 + t166 * t179;
t89 = t112 * t169 + t167 * t155;
t88 = -t114 * t167 + t166 * t178;
t87 = -t112 * t167 + t168 * t178;
t85 = t115 * t169 + t132 * t167;
t84 = t113 * t169 - t135 * t167;
t83 = t115 * t167 - t132 * t169;
t82 = t113 * t167 + t135 * t169;
t79 = -t141 * t167 + t169 * t97;
t78 = -t141 * t169 - t167 * t97;
t77 = t106 * t163 + t107 * t164;
t76 = -t106 * t164 + t107 * t163;
t71 = -pkin(5) * t111 + t191;
t70 = -pkin(5) * t110 + t192;
t69 = -pkin(3) * t111 + t81;
t68 = -pkin(3) * t110 + t80;
t67 = pkin(4) * t74 + qJ(2) * t161;
t66 = -pkin(4) * t75 + t193 * t161;
t65 = t163 * t83 + t164 * t85;
t64 = t163 * t82 + t164 * t84;
t63 = t163 * t85 - t164 * t83;
t62 = t163 * t84 - t164 * t82;
t58 = -pkin(4) * t106 - t169 * t59;
t57 = -pkin(4) * t107 + t167 * t59;
t54 = t167 * t95 + t169 * t61;
t53 = t167 * t61 - t169 * t95;
t52 = -pkin(4) * t83 + qJ(2) * t111 - t167 * t69 + t169 * t71;
t51 = -pkin(4) * t82 + qJ(2) * t110 - t167 * t68 + t169 * t70;
t50 = -pkin(4) * t85 + t193 * t111 - t167 * t71 - t169 * t69;
t49 = -pkin(4) * t84 + t193 * t110 - t167 * t70 - t169 * t68;
t48 = t163 * t53 + t164 * t54;
t47 = t163 * t54 - t164 * t53;
t46 = -pkin(4) * t53 + (pkin(3) * t167 - pkin(5) * t169 + qJ(2)) * t59;
t45 = -pkin(4) * t54 + (pkin(3) * t169 + pkin(5) * t167 + t193) * t59;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t109, 0, 0, 0, 0, 0, 0, 0, 0, 0, t99, 0, 0, 0, 0, 0, 0, -t101, t175, 0, -t196, 0, 0, 0, 0, 0, 0, t64, t65, t77, t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t108, 0, 0, 0, 0, 0, 0, 0, 0, 0, t98, 0, 0, 0, 0, 0, 0, t175, t101, 0, t55, 0, 0, 0, 0, 0, 0, t62, t63, t76, t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t161, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t161, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t161, 0, 0, 0, 0, 0, 0, -t110, -t111, 0, -t59; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t183, -t152, -t108, -qJ(1) * t108, 0, 0, 0, 0, 0, 0, -t183, -t98, t152, -qJ(1) * t98 + (-pkin(1) * t163 + qJ(2) * t164) * t161, 0, 0, -t175, 0, -t101, 0, t197, t199, t55, -qJ(1) * t55 - t163 * t66 + t164 * t67, -t163 * t92 + t164 * t94, -t163 * t78 + t164 * t79, -t163 * t88 + t164 * t90, -t163 * t91 + t164 * t93, -t163 * t87 + t164 * t89, -t117 * t163 + t119 * t164, -qJ(1) * t62 - t163 * t49 + t164 * t51, -qJ(1) * t63 - t163 * t50 + t164 * t52, -qJ(1) * t76 - t163 * t57 + t164 * t58, -qJ(1) * t47 - t163 * t45 + t164 * t46; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t152, -t183, t109, qJ(1) * t109, 0, 0, 0, 0, 0, 0, t152, t99, t183, qJ(1) * t99 + (pkin(1) * t164 + qJ(2) * t163) * t161, 0, 0, -t101, 0, t175, 0, t199, -t197, t196, -qJ(1) * t196 + t163 * t67 + t164 * t66, t163 * t94 + t164 * t92, t163 * t79 + t164 * t78, t163 * t90 + t164 * t88, t163 * t93 + t164 * t91, t163 * t89 + t164 * t87, t117 * t164 + t119 * t163, qJ(1) * t64 + t163 * t51 + t164 * t49, qJ(1) * t65 + t163 * t52 + t164 * t50, qJ(1) * t77 + t163 * t58 + t164 * t57, qJ(1) * t48 + t163 * t46 + t164 * t45; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t142, t143, 0, 0, 0, 0, 0, 0, 0, 0, t136, 0, -t143, pkin(1) * t136 - qJ(2) * t143, 0, 0, 0, 0, 0, -qJDD(3), -qJ(2) * t138 + t193 * t139 + t104, qJ(2) * t139 + t193 * t138 + t105, 0, qJ(2) * t75 + t193 * t74, (-t133 - t176) * t166, -t132 * t168 - t135 * t166, -t146 * t168 - t188, (-t134 + t177) * t168, -t148 * t166 - t185, 0, -pkin(3) * t135 - pkin(5) * t113 + qJ(2) * t84 - t193 * t82 + t191, pkin(3) * t132 - pkin(5) * t115 + qJ(2) * t85 - t193 * t83 - t192, -pkin(3) * t140 - pkin(5) * t137 + qJ(2) * t107 - t193 * t106 - t61, pkin(3) * t95 - pkin(5) * t61 + qJ(2) * t54 - t193 * t53;];
tauB_reg = t1;
