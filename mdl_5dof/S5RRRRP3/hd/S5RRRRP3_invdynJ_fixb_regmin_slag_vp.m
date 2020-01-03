% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% tau_reg [5x20]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRRP3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP3_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP3_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP3_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:49:28
% EndTime: 2019-12-31 21:49:31
% DurationCPUTime: 1.25s
% Computational Cost: add. (2156->256), mult. (3073->285), div. (0->0), fcn. (1662->12), ass. (0->164)
t98 = cos(qJ(2));
t162 = qJD(1) * t98;
t137 = qJD(2) * t162;
t94 = sin(qJ(2));
t151 = qJDD(1) * t94;
t106 = (t137 + t151) * pkin(1);
t167 = pkin(1) * qJD(1);
t144 = t94 * t167;
t192 = t98 * pkin(1);
t79 = qJDD(1) * t192;
t87 = qJDD(1) + qJDD(2);
t36 = t87 * pkin(2) - qJD(2) * t144 + t79;
t141 = pkin(1) * t162;
t88 = qJD(1) + qJD(2);
t44 = t88 * pkin(2) + t141;
t128 = qJD(3) * t144;
t93 = sin(qJ(3));
t56 = t93 * t128;
t97 = cos(qJ(3));
t201 = -(qJD(3) * t44 + t106) * t97 - t93 * t36 + t56;
t92 = sin(qJ(4));
t96 = cos(qJ(4));
t122 = t96 * pkin(4) + t92 * qJ(5);
t45 = -pkin(3) - t122;
t148 = qJDD(4) * qJ(5);
t30 = t97 * t144 + t93 * t44;
t83 = qJD(3) + t88;
t25 = t83 * pkin(8) + t30;
t180 = t92 * t25;
t82 = qJDD(3) + t87;
t9 = t82 * pkin(8) - t201;
t7 = t96 * t9;
t4 = t148 + t7 + (qJD(5) - t180) * qJD(4);
t155 = qJDD(4) * pkin(4);
t156 = qJD(4) * t96;
t198 = t25 * t156 - t155;
t6 = t92 * t9;
t5 = qJDD(5) + t6 + t198;
t199 = t4 * t96 + t5 * t92;
t91 = qJ(1) + qJ(2);
t86 = qJ(3) + t91;
t74 = sin(t86);
t75 = cos(t86);
t172 = g(1) * t75 + g(2) * t74;
t66 = g(1) * t74;
t84 = sin(t91);
t71 = g(1) * t84;
t197 = g(2) * t75;
t196 = t82 * pkin(3);
t195 = t83 * pkin(3);
t193 = t97 * pkin(2);
t100 = qJD(4) ^ 2;
t191 = pkin(8) * t100;
t159 = qJD(3) * t97;
t160 = qJD(3) * t93;
t179 = t93 * t94;
t78 = pkin(2) + t192;
t17 = t78 * t159 + (-t94 * t160 + (t97 * t98 - t179) * qJD(2)) * pkin(1);
t190 = t17 * t83;
t178 = t94 * t97;
t119 = t93 * t98 + t178;
t18 = t78 * t160 + (t119 * qJD(2) + t94 * t159) * pkin(1);
t189 = t18 * t83;
t61 = t93 * t144;
t29 = t97 * t44 - t61;
t188 = t29 * t83;
t187 = t30 * t83;
t186 = t45 * t82;
t185 = t45 * t83;
t184 = t74 * t92;
t183 = t75 * t92;
t182 = t83 * t92;
t181 = t83 * t96;
t177 = t96 * t25;
t176 = t96 * t82;
t142 = pkin(2) * t160;
t152 = qJD(4) * qJ(5);
t153 = t92 * qJD(5);
t157 = qJD(4) * t92;
t41 = pkin(4) * t157 - t96 * t152 - t153;
t31 = t41 + t142;
t39 = t119 * t167;
t175 = t31 - t39;
t174 = t41 - t30;
t173 = g(1) * t184 - g(2) * t183;
t171 = pkin(1) * t178 + t93 * t78;
t85 = cos(t91);
t170 = g(1) * t85 + g(2) * t84;
t89 = t92 ^ 2;
t90 = t96 ^ 2;
t169 = t89 - t90;
t168 = t89 + t90;
t38 = pkin(8) + t171;
t166 = t100 * t38;
t76 = t93 * pkin(2) + pkin(8);
t165 = t100 * t76;
t163 = pkin(8) * qJDD(4);
t158 = qJD(4) * t83;
t16 = t152 + t177;
t154 = t16 * qJD(4);
t150 = qJDD(4) * t38;
t149 = qJDD(4) * t76;
t68 = pkin(1) * t179;
t81 = t83 ^ 2;
t147 = t92 * t81 * t96;
t53 = t96 * t66;
t146 = t29 * t157 + t30 * t181 + t53;
t40 = t97 * t141 - t61;
t145 = t40 * t157 + t39 * t181 + t53;
t143 = pkin(2) * t159;
t140 = t83 * t160;
t126 = t93 * pkin(1) * t137 + qJDD(1) * t68 + t44 * t160 + (t128 - t36) * t97;
t10 = t126 - t196;
t139 = -t10 - t197;
t138 = t168 * t82;
t134 = t97 * t78 - t68;
t27 = -t134 + t45;
t136 = t27 * t83 - t17;
t24 = -t29 - t195;
t133 = t10 * t92 + t24 * t156 - t173;
t132 = t74 * pkin(8) - t45 * t75;
t131 = -qJD(4) * pkin(4) + qJD(5);
t130 = qJD(1) * (-qJD(2) + t88);
t129 = qJD(2) * (-qJD(1) - t88);
t127 = -g(2) * t85 + t71 + t79;
t125 = pkin(2) * t85 + t132;
t123 = t191 - t196;
t121 = pkin(4) * t92 - qJ(5) * t96;
t15 = t131 + t180;
t120 = t15 * t92 + t16 * t96;
t118 = g(1) * t183 + g(2) * t184 - g(3) * t96 - t6;
t77 = -pkin(3) - t193;
t117 = t77 * t82 + t165;
t1 = (t121 * qJD(4) - t153) * t83 + t186 + t126;
t116 = -t1 - t186 - t191;
t43 = t45 - t193;
t115 = -t43 * t82 - t1 - t165;
t114 = t43 * t83 - t143;
t113 = t77 * t83 - t143;
t112 = -qJDD(5) + t118;
t37 = -pkin(3) - t134;
t110 = t37 * t82 + t166 + t189;
t109 = -t126 + t66 - t197;
t108 = t45 * t66;
t107 = t15 * t156 - t92 * t154 - t172 + t199;
t12 = t18 + t41;
t105 = -t12 * t83 - t27 * t82 - t1 - t166;
t104 = -t150 + (t37 * t83 - t17) * qJD(4);
t63 = t75 * pkin(8);
t103 = -g(1) * t63 - t108;
t102 = (t15 * t96 - t16 * t92) * qJD(4) + t199;
t101 = t172 + t201;
t99 = cos(qJ(1));
t95 = sin(qJ(1));
t59 = t92 * t82;
t47 = qJDD(4) * t96 - t100 * t92;
t46 = qJDD(4) * t92 + t100 * t96;
t35 = t121 * t83;
t34 = 0.2e1 * t156 * t182 + t89 * t82;
t22 = -0.2e1 * t169 * t158 + 0.2e1 * t92 * t176;
t20 = t24 * t157;
t13 = -t29 + t185;
t11 = t13 * t157;
t2 = [qJDD(1), g(1) * t95 - g(2) * t99, g(1) * t99 + g(2) * t95, t87, (t94 * t129 + t87 * t98) * pkin(1) + t127, ((-qJDD(1) - t87) * t94 + t98 * t129) * pkin(1) + t170, t82, t134 * t82 + t109 - t189, -t171 * t82 + t101 - t190, t34, t22, t46, t47, 0, t20 + t53 + t104 * t92 + (-t110 + t139) * t96, t104 * t96 + t110 * t92 + t133, t11 + t53 + (t136 * qJD(4) - t150) * t92 + (t105 - t197) * t96, t38 * t138 + t168 * t190 + t107, (t150 + (-t13 - t136) * qJD(4)) * t96 + t105 * t92 + t173, t1 * t27 + t13 * t12 - g(1) * (-t95 * pkin(1) - pkin(2) * t84 + t63) - g(2) * (t99 * pkin(1) + t125) - t108 + t120 * t17 + t102 * t38; 0, 0, 0, t87, t94 * pkin(1) * t130 + t127, (t98 * t130 - t151) * pkin(1) + t170, t82, t39 * t83 + (t82 * t97 - t140) * pkin(2) + t109, t40 * t83 + t56 + (-pkin(2) * t82 - t36) * t93 + ((-pkin(2) * t83 - t44) * qJD(3) - t106) * t97 + t172, t34, t22, t46, t47, 0, t20 + (t113 * qJD(4) - t149) * t92 + (-pkin(2) * t140 - t117 + t139) * t96 + t145, (-t149 + (t113 + t40) * qJD(4)) * t96 + ((-t39 + t142) * t83 + t117) * t92 + t133, t11 + (qJD(4) * t114 - t149) * t92 + (-t31 * t83 + t115 - t197) * t96 + t145, t76 * t138 + t107 + (t143 - t40) * t83 * t168, (t149 + (-t114 - t13 - t40) * qJD(4)) * t96 + (-t175 * t83 + t115) * t92 + t173, t1 * t43 - g(2) * t125 - t120 * t40 + t175 * t13 + (t120 * t159 + t71) * pkin(2) + t102 * t76 + t103; 0, 0, 0, 0, 0, 0, t82, t109 + t187, t101 + t188, t34, t22, t46, t47, 0, t20 + (-pkin(3) * t158 - t163) * t92 + (-t123 + t139) * t96 + t146, (-t163 + (t29 - t195) * qJD(4)) * t96 + (t123 - t187) * t92 + t133, t11 + (t158 * t45 - t163) * t92 + (-t41 * t83 + t116 - t197) * t96 + t146, pkin(8) * t138 - t168 * t188 + t107, (t163 + (-t13 - t29 - t185) * qJD(4)) * t96 + (-t174 * t83 + t116) * t92 + t173, t102 * pkin(8) - g(2) * t132 + t1 * t45 - t120 * t29 + t174 * t13 + t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t147, t169 * t81, t59, t176, qJDD(4), -t24 * t182 + t118, g(3) * t92 - t7 + (-t24 * t83 + t172) * t96, 0.2e1 * t155 + (-t13 * t92 + t35 * t96) * t83 + t112, -t121 * t82 + ((t16 - t152) * t92 + (t131 - t15) * t96) * t83, 0.2e1 * t148 + 0.2e1 * qJD(4) * qJD(5) + t7 + (t35 * t83 - g(3)) * t92 + (t13 * t83 - t172) * t96, t4 * qJ(5) - t5 * pkin(4) - t13 * t35 - t15 * t177 - g(3) * t122 + (qJD(5) + t180) * t16 + t172 * t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(4) - t147, t59, -t89 * t81 - t100, t13 * t182 - t112 - t154 + t198;];
tau_reg = t2;
