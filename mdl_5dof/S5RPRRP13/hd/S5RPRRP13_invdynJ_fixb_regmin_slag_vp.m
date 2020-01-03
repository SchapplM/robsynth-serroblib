% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRRP13
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% tau_reg [5x24]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRP13_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP13_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP13_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP13_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP13_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP13_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:59:40
% EndTime: 2019-12-31 18:59:45
% DurationCPUTime: 1.85s
% Computational Cost: add. (1945->311), mult. (3662->395), div. (0->0), fcn. (2150->6), ass. (0->148)
t81 = sin(qJ(1));
t184 = g(1) * t81;
t84 = cos(qJ(1));
t192 = -g(2) * t84 + t184;
t80 = sin(qJ(3));
t64 = t80 * qJD(1) + qJD(4);
t183 = g(3) * t80;
t83 = cos(qJ(3));
t194 = -t192 * t83 + t183;
t120 = qJD(4) * t80 + qJD(1);
t82 = cos(qJ(4));
t138 = t82 * qJD(3);
t126 = t83 * t138;
t79 = sin(qJ(4));
t140 = t79 * qJD(3);
t149 = qJD(1) * t83;
t50 = t82 * t149 + t140;
t147 = qJD(3) * t50;
t134 = t83 * qJDD(1);
t143 = qJD(4) * t79;
t98 = t80 * t138 + t83 * t143;
t13 = qJD(1) * t98 - qJD(4) * t138 - t79 * qJDD(3) - t82 * t134;
t163 = t83 * t13;
t133 = qJD(1) * qJD(3);
t122 = t83 * t133;
t135 = t80 * qJDD(1);
t46 = qJDD(4) + t122 + t135;
t166 = t82 * t46;
t193 = (t120 * t79 - t126) * t64 + (t147 - t166) * t80 + t163;
t169 = t80 * t82;
t117 = t80 * pkin(3) - t83 * pkin(7);
t55 = qJ(2) + t117;
t85 = -pkin(1) - pkin(6);
t156 = t85 * t169 + t79 * t55;
t87 = qJD(1) ^ 2;
t191 = -t87 * qJ(2) - t192;
t144 = qJD(4) * t50;
t171 = t79 * t80;
t14 = -t82 * qJDD(3) - t133 * t171 + t79 * t134 + t144;
t62 = t85 * qJD(1) + qJD(2);
t162 = t83 * t62;
t42 = -qJD(3) * pkin(3) - t162;
t48 = t79 * t149 - t138;
t10 = t48 * pkin(4) - t50 * qJ(5) + t42;
t185 = pkin(7) * t46;
t190 = t10 * t64 - t185;
t152 = pkin(7) * qJD(4);
t127 = t64 * t152;
t189 = t127 - t194;
t186 = t50 ^ 2;
t182 = g(3) * t83;
t181 = t46 * pkin(4);
t35 = t55 * qJD(1);
t54 = t80 * t62;
t41 = qJD(3) * pkin(7) + t54;
t12 = t79 * t35 + t82 * t41;
t9 = t64 * qJ(5) + t12;
t178 = t9 * t64;
t177 = t12 * t64;
t176 = t13 * t79;
t175 = t48 * t64;
t174 = t50 * t48;
t173 = t50 * t64;
t172 = t50 * t82;
t170 = t79 * t85;
t168 = t81 * t79;
t167 = t81 * t82;
t118 = pkin(3) * t83 + pkin(7) * t80;
t47 = qJD(3) * t118 + qJD(2);
t165 = t82 * t47;
t164 = t82 * t55;
t161 = t84 * t79;
t160 = t84 * t82;
t109 = pkin(4) * t79 - qJ(5) * t82;
t159 = -t79 * qJD(5) + t64 * t109 - t54;
t53 = t118 * qJD(1);
t158 = t82 * t162 + t79 * t53;
t157 = g(2) * t83 * t160 + g(3) * t169;
t155 = t84 * pkin(1) + t81 * qJ(2);
t78 = t83 ^ 2;
t154 = t80 ^ 2 - t78;
t86 = qJD(3) ^ 2;
t153 = -t86 - t87;
t151 = t46 * qJ(5);
t150 = pkin(1) * qJDD(1);
t148 = qJD(3) * t48;
t146 = qJD(3) * t80;
t145 = qJD(3) * t83;
t142 = qJD(4) * t82;
t141 = t42 * qJD(4);
t11 = t82 * t35 - t79 * t41;
t137 = qJD(5) - t11;
t136 = qJDD(3) * t80;
t132 = qJDD(1) * qJ(2);
t131 = t83 * t184;
t19 = qJD(1) * t47 + qJDD(1) * t55;
t60 = t85 * qJDD(1) + qJDD(2);
t25 = qJDD(3) * pkin(7) + t62 * t145 + t80 * t60;
t130 = -t35 * t142 - t79 * t19 - t82 * t25;
t129 = t85 * t126 + t55 * t142 + t79 * t47;
t125 = 0.2e1 * qJD(1) * qJD(2);
t124 = -pkin(4) + t170;
t121 = -t41 * t142 - t35 * t143 + t82 * t19 - t79 * t25;
t119 = qJDD(2) - t150;
t36 = t80 * t168 - t160;
t38 = t80 * t161 + t167;
t116 = g(1) * t38 + g(2) * t36;
t37 = t80 * t167 + t161;
t39 = t80 * t160 - t168;
t115 = -g(1) * t39 - g(2) * t37;
t114 = g(1) * t84 + g(2) * t81;
t8 = -t64 * pkin(4) + t137;
t112 = t79 * t9 - t8 * t82;
t111 = t79 * t8 + t82 * t9;
t110 = t82 * pkin(4) + t79 * qJ(5);
t106 = pkin(3) + t110;
t103 = t109 - t85;
t102 = t64 * t142 + t79 * t46;
t101 = -t64 * t143 + t166;
t24 = -qJDD(3) * pkin(3) + t62 * t146 - t83 * t60;
t100 = -t41 * t143 - t130;
t99 = 0.2e1 * qJ(2) * t133 + qJDD(3) * t85;
t97 = -t60 - t191;
t96 = t42 * t64 - t185;
t95 = -t192 * t80 - t182;
t94 = g(1) * t36 - g(2) * t38 + t79 * t182 + t121;
t93 = -t114 + t125 + 0.2e1 * t132;
t1 = t64 * qJD(5) + t100 + t151;
t2 = qJDD(5) - t121 - t181;
t92 = -qJD(4) * t112 + t1 * t82 + t2 * t79;
t91 = t10 * t50 + qJDD(5) - t94;
t90 = -t85 * t86 + t93;
t89 = -g(1) * t37 + g(2) * t39 - t82 * t182 + t100;
t88 = -t46 * t171 - t83 * t14 + t48 * t146 + (-t120 * t82 - t83 * t140) * t64;
t74 = t84 * qJ(2);
t71 = qJDD(3) * t83;
t28 = t103 * t83;
t27 = t124 * t80 - t164;
t26 = t80 * qJ(5) + t156;
t22 = t50 * pkin(4) + t48 * qJ(5);
t18 = -t82 * t53 + (-pkin(4) * qJD(1) + t62 * t79) * t83;
t17 = qJ(5) * t149 + t158;
t7 = (qJD(4) * t110 - qJD(5) * t82) * t83 - t103 * t146;
t6 = -t13 + t175;
t5 = t156 * qJD(4) + t124 * t145 - t165;
t4 = qJ(5) * t145 + (-t85 * t143 + qJD(5)) * t80 + t129;
t3 = t14 * pkin(4) + t13 * qJ(5) - t50 * qJD(5) + t24;
t15 = [qJDD(1), t192, t114, qJDD(2) - t192 - 0.2e1 * t150, t93, -t119 * pkin(1) - g(1) * (-t81 * pkin(1) + t74) - g(2) * t155 + (t125 + t132) * qJ(2), t78 * qJDD(1) - 0.2e1 * t122 * t80, 0.2e1 * t154 * t133 - 0.2e1 * t80 * t134, -t86 * t80 + t71, -t86 * t83 - t136, 0, t80 * t90 + t83 * t99, -t80 * t99 + t83 * t90, -t82 * t163 - t50 * t98, (t48 * t82 + t50 * t79) * t146 + (t176 - t14 * t82 + (t48 * t79 - t172) * qJD(4)) * t83, (-t64 * t138 - t13) * t80 + (t101 + t147) * t83, (t64 * t140 - t14) * t80 + (-t102 - t148) * t83, t64 * t145 + t46 * t80, (-t55 * t143 + t165) * t64 + t46 * t164 + (-t42 * t140 + (-t102 + t148) * t85 + t121) * t80 + (t82 * t141 - t85 * t14 + t24 * t79 + (-t64 * t170 + t11) * qJD(3)) * t83 + t115, -t129 * t64 - t156 * t46 + ((t64 * t85 + t41) * t143 + (-t42 * t82 + t50 * t85) * qJD(3) + t130) * t80 + (-t12 * qJD(3) + t85 * t13 - t79 * t141 + t24 * t82) * t83 + t116, t28 * t14 - t27 * t46 + t7 * t48 - t5 * t64 + (-t10 * t140 - t2) * t80 + (-qJD(3) * t8 + t10 * t142 + t3 * t79) * t83 + t115, -t27 * t13 - t26 * t14 - t4 * t48 + t5 * t50 + t112 * t146 + (-qJD(4) * t111 - t1 * t79 + t2 * t82 + t114) * t83, t28 * t13 + t26 * t46 + t4 * t64 - t7 * t50 + (t10 * t138 + t1) * t80 + (qJD(3) * t9 + t10 * t143 - t3 * t82) * t83 - t116, t1 * t26 + t9 * t4 + t3 * t28 + t10 * t7 + t2 * t27 + t8 * t5 - g(1) * (t39 * pkin(4) + t38 * qJ(5) + t117 * t84 + t74) - g(2) * (t37 * pkin(4) + t84 * pkin(6) + t36 * qJ(5) + t155) + (-g(1) * t85 - g(2) * t117) * t81; 0, 0, 0, qJDD(1), -t87, t119 + t191, 0, 0, 0, 0, 0, t153 * t80 + t71, t153 * t83 - t136, 0, 0, 0, 0, 0, t88, t193, t88, (t120 * t50 - t14 * t80 - t48 * t145) * t82 + (t120 * t48 - t13 * t80 + t50 * t145) * t79, -t193, -t112 * qJD(1) + (qJD(3) * t111 - t3) * t83 + (qJD(3) * t10 + t92) * t80 - t192; 0, 0, 0, 0, 0, 0, t83 * t87 * t80, -t154 * t87, t134, -t135, qJDD(3), -t83 * t97 + t183, t80 * t97 + t182, t64 * t172 - t176, (-t13 - t175) * t82 + (-t14 - t173) * t79, (t64 * t169 - t50 * t83) * qJD(1) + t102, (-t64 * t171 + t48 * t83) * qJD(1) + t101, -t64 * t149, -t11 * t149 - t48 * t54 - pkin(3) * t14 + (-t131 - t24 + (-t53 - t152) * t64) * t82 + (t64 * t162 + t96) * t79 + t157, pkin(3) * t13 + t158 * t64 + t12 * t149 - t50 * t54 + t96 * t82 + (t189 + t24) * t79, t8 * t149 - t106 * t14 + t18 * t64 + t159 * t48 + (-t127 - t3 - t131) * t82 + t190 * t79 + t157, t17 * t48 - t18 * t50 + (t1 + t64 * t8 + (-t14 + t144) * pkin(7)) * t82 + (t2 - t178 + (qJD(4) * t48 - t13) * pkin(7)) * t79 + t95, -t9 * t149 - t106 * t13 - t17 * t64 - t159 * t50 - t190 * t82 + (-t189 - t3) * t79, -t9 * t17 - t8 * t18 + t159 * t10 + (t92 + t95) * pkin(7) + (-t3 + t194) * t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t174, -t48 ^ 2 + t186, t6, t173 - t14, t46, -t42 * t50 + t177 + t94, t11 * t64 + t42 * t48 - t89, -t22 * t48 + t177 + 0.2e1 * t181 - t91, pkin(4) * t13 - t14 * qJ(5) + (-t12 + t9) * t50 + (t8 - t137) * t48, 0.2e1 * t151 - t10 * t48 + t22 * t50 + (0.2e1 * qJD(5) - t11) * t64 + t89, t1 * qJ(5) - t2 * pkin(4) - t10 * t22 - t8 * t12 - g(1) * (-pkin(4) * t36 + qJ(5) * t37) - g(2) * (pkin(4) * t38 - qJ(5) * t39) + t137 * t9 + t109 * t182; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46 + t174, t6, -t64 ^ 2 - t186, -t178 + t91 - t181;];
tau_reg = t15;
