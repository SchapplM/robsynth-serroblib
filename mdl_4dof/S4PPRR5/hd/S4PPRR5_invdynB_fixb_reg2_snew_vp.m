% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S4PPRR5
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
% Datum: 2019-12-31 16:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S4PPRR5_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR5_invdynB_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR5_invdynB_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR5_invdynB_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR5_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR5_invdynB_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:19:53
% EndTime: 2019-12-31 16:19:55
% DurationCPUTime: 1.02s
% Computational Cost: add. (1680->192), mult. (3072->269), div. (0->0), fcn. (2061->6), ass. (0->129)
t195 = pkin(2) + qJ(1);
t161 = sin(qJ(3));
t163 = cos(qJ(3));
t174 = t163 * qJDD(3);
t194 = qJD(3) ^ 2;
t133 = -t161 * t194 + t174;
t157 = sin(pkin(6));
t158 = cos(pkin(6));
t137 = t158 * g(1) + t157 * g(2);
t193 = pkin(1) + pkin(4);
t85 = t193 * t133 + t161 * t137;
t160 = sin(qJ(4));
t136 = t157 * g(1) - t158 * g(2);
t130 = -qJDD(2) + t136;
t155 = g(3) - qJDD(1);
t168 = -t163 * t130 + t161 * t155;
t96 = -qJDD(3) * pkin(3) - t194 * pkin(5) - t168;
t192 = t160 * t96;
t162 = cos(qJ(4));
t191 = t162 * t96;
t180 = t161 * t130 + t163 * t155;
t97 = -t194 * pkin(3) + qJDD(3) * pkin(5) - t180;
t80 = -t160 * t137 + t162 * t97;
t175 = t161 * qJDD(3);
t132 = t163 * t194 + t175;
t190 = t157 * t132;
t189 = t157 * t133;
t188 = t157 * t155;
t187 = t158 * t132;
t186 = t158 * t155;
t144 = t160 * t194 * t162;
t138 = qJDD(4) + t144;
t185 = t160 * t138;
t139 = qJDD(4) - t144;
t184 = t160 * t139;
t182 = t162 * t138;
t181 = t162 * t139;
t153 = t160 ^ 2;
t154 = t162 ^ 2;
t179 = t153 + t154;
t178 = t153 * t194;
t177 = qJD(3) * qJD(4);
t176 = t160 * qJDD(3);
t149 = t162 * qJDD(3);
t173 = t160 * t177;
t172 = t162 * t177;
t131 = t179 * qJDD(3);
t151 = t154 * t194;
t134 = t151 + t178;
t79 = t162 * t137 + t160 * t97;
t61 = t160 * t79 + t162 * t80;
t98 = t161 * t131 + t163 * t134;
t99 = t163 * t131 - t161 * t134;
t171 = pkin(3) * t134 + pkin(5) * t131 - qJ(2) * t99 + t195 * t98 + t61;
t170 = qJ(2) * t133 - t195 * t132 + t180;
t169 = qJ(2) * t132 + t195 * t133 + t168;
t120 = t157 * t137;
t94 = t158 * t130 - t120;
t100 = t158 * t136 - t120;
t121 = t158 * t137;
t95 = -t157 * t130 - t121;
t101 = -t157 * t136 - t121;
t167 = t161 * t144;
t166 = t163 * t144;
t59 = t160 * t80 - t162 * t79;
t165 = t161 * t180 - t163 * t168;
t74 = -t161 * t168 - t163 * t180;
t86 = t193 * t132 - t163 * t137;
t164 = qJD(4) ^ 2;
t143 = -t151 - t164;
t142 = t151 - t164;
t141 = -t164 - t178;
t140 = t164 - t178;
t135 = t151 - t178;
t129 = qJ(2) * t137;
t128 = t149 - 0.2e1 * t173;
t127 = t149 - t173;
t126 = t172 + t176;
t125 = 0.2e1 * t172 + t176;
t124 = t179 * t177;
t119 = t158 * t133;
t115 = t162 * t126 - t153 * t177;
t113 = t163 * qJDD(4) - t161 * t124;
t112 = -t160 * t127 - t154 * t177;
t111 = -t160 * t141 - t181;
t110 = -t160 * t140 + t182;
t109 = t162 * t143 - t185;
t108 = t162 * t142 - t184;
t107 = t162 * t141 - t184;
t106 = t162 * t140 + t185;
t105 = t160 * t143 + t182;
t104 = t160 * t142 + t181;
t103 = (t126 + t172) * t160;
t102 = (t127 - t173) * t162;
t93 = -t160 * t125 + t162 * t128;
t92 = t162 * t125 + t160 * t128;
t90 = -t161 * t115 - t166;
t89 = -t161 * t112 + t166;
t88 = -t161 * t110 + t160 * t174;
t87 = -t161 * t108 + t162 * t174;
t84 = t163 * t111 + t161 * t125;
t83 = t163 * t109 - t161 * t128;
t82 = t161 * t111 - t163 * t125;
t81 = t161 * t109 + t163 * t128;
t76 = -t163 * t135 - t161 * t93;
t71 = -pkin(5) * t107 + t191;
t70 = -pkin(5) * t105 + t192;
t69 = -t157 * t165 - t121;
t68 = t158 * t165 - t120;
t67 = -pkin(3) * t107 + t80;
t66 = -pkin(3) * t105 + t79;
t65 = t158 * t107 + t157 * t82;
t64 = t158 * t105 + t157 * t81;
t63 = t157 * t107 - t158 * t82;
t62 = t157 * t105 - t158 * t81;
t58 = -pkin(2) * t165 - qJ(2) * t74;
t57 = -pkin(2) * t137 - t193 * t74;
t56 = t161 * t96 + t163 * t61;
t55 = t161 * t61 - t163 * t96;
t54 = t161 * t59 - t193 * t99;
t53 = pkin(2) * t82 - pkin(3) * t125 + pkin(5) * t111 - qJ(2) * t84 + t192;
t52 = pkin(2) * t81 + pkin(3) * t128 + pkin(5) * t109 - qJ(2) * t83 - t191;
t50 = t157 * t55 + t158 * t59;
t49 = t157 * t59 - t158 * t55;
t48 = pkin(2) * t107 - t161 * t71 - t163 * t67 - t193 * t84;
t47 = pkin(2) * t105 - t161 * t70 - t163 * t66 - t193 * t83;
t46 = pkin(2) * t55 - pkin(3) * t96 + pkin(5) * t61 - qJ(2) * t56;
t45 = -t193 * t56 + (pkin(3) * t163 + pkin(5) * t161 + pkin(2)) * t59;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95, 0, 0, 0, 0, 0, 0, t189, -t190, 0, t69, 0, 0, 0, 0, 0, 0, t64, t65, t157 * t98, t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t100, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, 0, 0, 0, 0, 0, 0, -t119, t187, 0, t68, 0, 0, 0, 0, 0, 0, t62, t63, -t158 * t98, t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t155, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t155, 0, 0, 0, 0, 0, 0, -t132, -t133, 0, t74, 0, 0, 0, 0, 0, 0, t83, t84, t99, t56; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t188, -t186, -t100, -qJ(1) * t100, 0, 0, 0, 0, 0, 0, -t94, t188, t186, -qJ(1) * t94 + (-pkin(1) * t157 + qJ(2) * t158) * t155, 0, 0, t190, 0, t189, t158 * qJDD(3), -t157 * t86 + t169 * t158, -t157 * t85 + t170 * t158, t157 * t74, -qJ(1) * t68 - t157 * t57 + t158 * t58, t158 * t103 - t157 * t90, -t157 * t76 + t158 * t92, t158 * t106 - t157 * t88, t158 * t102 - t157 * t89, t158 * t104 - t157 * t87, -t157 * t113, -qJ(1) * t62 - t157 * t47 + t158 * t52, -qJ(1) * t63 - t157 * t48 + t158 * t53, -t157 * t54 + t171 * t158, -qJ(1) * t49 - t157 * t45 + t158 * t46; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t186, -t188, t101, qJ(1) * t101, 0, 0, 0, 0, 0, 0, t95, -t186, t188, qJ(1) * t95 + (pkin(1) * t158 + qJ(2) * t157) * t155, 0, 0, -t187, 0, -t119, t157 * qJDD(3), t169 * t157 + t158 * t86, t170 * t157 + t158 * t85, -t158 * t74, qJ(1) * t69 + t157 * t58 + t158 * t57, t157 * t103 + t158 * t90, t157 * t92 + t158 * t76, t157 * t106 + t158 * t88, t157 * t102 + t158 * t89, t157 * t104 + t158 * t87, t158 * t113, qJ(1) * t64 + t157 * t52 + t158 * t47, qJ(1) * t65 + t157 * t53 + t158 * t48, t171 * t157 + t158 * t54, qJ(1) * t50 + t157 * t46 + t158 * t45; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t136, t137, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t130, -t137, pkin(1) * t130 - t129, 0, 0, t133, 0, -t132, 0, -t85, t86, t165, t165 * t193 - t129, t163 * t115 - t167, -t161 * t135 + t163 * t93, t163 * t110 + t160 * t175, t163 * t112 + t167, t163 * t108 + t161 * t149, t161 * qJDD(4) + t163 * t124, qJ(2) * t105 - t161 * t66 + t163 * t70 - t193 * t81, qJ(2) * t107 - t161 * t67 + t163 * t71 - t193 * t82, -t163 * t59 - t193 * t98, -t193 * t55 + (pkin(3) * t161 - pkin(5) * t163 + qJ(2)) * t59;];
tauB_reg = t1;
