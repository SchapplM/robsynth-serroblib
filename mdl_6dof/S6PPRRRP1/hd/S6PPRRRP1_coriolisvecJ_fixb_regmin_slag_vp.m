% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PPRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% 
% Output:
% tauc_reg [6x21]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PPRRRP1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRP1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:54:24
% EndTime: 2019-03-08 18:54:31
% DurationCPUTime: 1.94s
% Computational Cost: add. (2422->266), mult. (6486->405), div. (0->0), fcn. (5476->12), ass. (0->151)
t100 = sin(pkin(12));
t102 = sin(pkin(6));
t108 = sin(qJ(3));
t111 = cos(qJ(3));
t103 = cos(pkin(12));
t104 = cos(pkin(7));
t167 = t103 * t104;
t116 = (t100 * t111 + t108 * t167) * t102;
t101 = sin(pkin(7));
t170 = t101 * t108;
t105 = cos(pkin(6));
t93 = qJD(1) * t105 + qJD(2);
t49 = qJD(1) * t116 + t93 * t170;
t107 = sin(qJ(4));
t110 = cos(qJ(4));
t129 = pkin(4) * t107 - pkin(10) * t110;
t86 = t129 * qJD(4);
t204 = -t49 + t86;
t169 = t101 * t111;
t197 = (-t100 * t108 + t111 * t167) * t102;
t203 = t105 * t169 + t197;
t202 = qJD(1) * t197 + t93 * t169;
t109 = cos(qJ(5));
t157 = qJD(4) * t107;
t106 = sin(qJ(5));
t166 = t106 * t110;
t194 = pkin(9) * t106;
t163 = qJD(1) * t102;
t143 = t103 * t163;
t132 = t104 * t143;
t144 = t100 * t163;
t48 = -t108 * t144 + (t101 * t93 + t132) * t111;
t201 = t204 * t109 + t157 * t194 + t48 * t166;
t154 = qJD(5) * t109;
t164 = t109 * t110;
t88 = -pkin(4) * t110 - pkin(10) * t107 - pkin(3);
t200 = t204 * t106 + t88 * t154 - t48 * t164;
t47 = qJD(3) * pkin(9) + t49;
t64 = -t101 * t143 + t104 * t93;
t199 = -t107 * t47 + t110 * t64;
t161 = qJD(3) * t108;
t142 = t101 * t161;
t159 = qJD(3) * t111;
t45 = t132 * t161 + t93 * t142 + t144 * t159;
t198 = qJD(3) * t49 - t45;
t138 = t107 * t154;
t150 = qJD(4) * qJD(5);
t156 = qJD(4) * t110;
t63 = (t106 * t156 + t138) * qJD(3) + t106 * t150;
t158 = qJD(4) * t106;
t162 = qJD(3) * t107;
t82 = t109 * t162 + t158;
t196 = t82 ^ 2;
t31 = t107 * t64 + t110 * t47;
t28 = qJD(4) * pkin(10) + t31;
t39 = t88 * qJD(3) - t48;
t9 = -t106 * t28 + t109 * t39;
t6 = -qJ(6) * t82 + t9;
t160 = qJD(3) * t110;
t94 = -qJD(5) + t160;
t5 = -pkin(5) * t94 + t6;
t195 = t5 - t6;
t153 = t109 * qJD(4);
t80 = t106 * t162 - t153;
t193 = t80 * t94;
t192 = t82 * t94;
t191 = -qJ(6) - pkin(10);
t124 = pkin(5) * t107 - qJ(6) * t164;
t152 = t109 * qJD(6);
t172 = qJ(6) * t107;
t95 = pkin(9) * t164;
t190 = -t107 * t152 + t124 * qJD(4) + (-t95 + (-t88 + t172) * t106) * qJD(5) + t201;
t85 = t129 * qJD(3);
t189 = t106 * t85 + t109 * t199;
t165 = t107 * t109;
t188 = (-pkin(9) * qJD(4) - qJ(6) * qJD(5)) * t165 + (-qJD(6) * t107 + (-pkin(9) * qJD(5) - qJ(6) * qJD(4)) * t110) * t106 + t200;
t133 = qJD(5) * t191;
t187 = t152 - t189 + (qJ(6) * t160 + t133) * t106;
t134 = -t106 * t199 + t109 * t85;
t186 = -t124 * qJD(3) - t106 * qJD(6) + t109 * t133 - t134;
t183 = t106 * t88 + t95;
t98 = t107 ^ 2;
t182 = -t110 ^ 2 + t98;
t181 = qJD(3) * pkin(3);
t27 = -qJD(4) * pkin(4) - t199;
t180 = t106 * t27;
t179 = t106 * t39;
t178 = t106 * t94;
t177 = t109 * t94;
t44 = t202 * qJD(3);
t13 = t107 * t44 + t47 * t156 + t64 * t157;
t175 = t13 * t106;
t174 = t13 * t109;
t155 = qJD(5) * t106;
t139 = t107 * t155;
t140 = t110 * t153;
t62 = -t109 * t150 + (t139 - t140) * qJD(3);
t173 = t62 * t106;
t113 = qJD(3) ^ 2;
t168 = t101 * t113;
t151 = qJD(3) * qJD(4);
t148 = t94 * t155;
t147 = t108 * t168;
t145 = pkin(5) * t106 + pkin(9);
t141 = t101 * t159;
t136 = t107 * t151;
t12 = qJD(4) * t199 + t110 * t44;
t35 = qJD(3) * t86 + t45;
t135 = t106 * t12 - t109 * t35;
t131 = t110 * t141;
t130 = t107 * t141;
t8 = pkin(5) * t63 + t13;
t10 = t109 * t28 + t179;
t54 = t105 * t170 + t116;
t67 = -t101 * t102 * t103 + t104 * t105;
t38 = t107 * t67 + t110 * t54;
t19 = -t106 * t203 + t109 * t38;
t18 = -t106 * t38 - t109 * t203;
t37 = t107 * t54 - t110 * t67;
t127 = qJD(3) * t98 - t110 * t94;
t112 = qJD(4) ^ 2;
t126 = pkin(9) * t112 - t198;
t46 = -t48 - t181;
t125 = qJD(4) * (t46 + t48 - t181);
t69 = t104 * t107 + t110 * t170;
t57 = -t106 * t69 - t109 * t169;
t123 = t106 * t169 - t109 * t69;
t68 = -t104 * t110 + t107 * t170;
t121 = t106 * t35 + t109 * t12 + t39 * t154 - t28 * t155;
t114 = -t10 * qJD(5) - t135;
t90 = t191 * t109;
t89 = t191 * t106;
t79 = t109 * t88;
t77 = t80 ^ 2;
t59 = -t106 * t172 + t183;
t56 = t69 * qJD(4) + t130;
t55 = -t68 * qJD(4) + t131;
t52 = -qJ(6) * t165 + t79 + (-pkin(5) - t194) * t110;
t51 = t54 * qJD(3);
t50 = t203 * qJD(3);
t26 = t123 * qJD(5) - t106 * t55 + t109 * t142;
t25 = t57 * qJD(5) + t106 * t142 + t109 * t55;
t20 = pkin(5) * t80 + qJD(6) + t27;
t17 = -t37 * qJD(4) + t50 * t110;
t16 = t38 * qJD(4) + t50 * t107;
t7 = -qJ(6) * t80 + t10;
t4 = t18 * qJD(5) + t51 * t106 + t17 * t109;
t3 = -t19 * qJD(5) - t17 * t106 + t51 * t109;
t2 = -qJ(6) * t63 - qJD(6) * t80 + t121;
t1 = pkin(5) * t136 + qJ(6) * t62 - qJD(6) * t82 + t114;
t11 = [0, 0, 0, -t51 * qJD(3), -t50 * qJD(3), 0, 0, 0, 0, 0, -t16 * qJD(4) + (-t110 * t51 - t157 * t203) * qJD(3), -t17 * qJD(4) + (t107 * t51 - t156 * t203) * qJD(3), 0, 0, 0, 0, 0, t18 * t136 + t16 * t80 - t3 * t94 + t37 * t63, -t136 * t19 + t16 * t82 - t37 * t62 + t4 * t94, t18 * t62 - t19 * t63 - t3 * t82 - t4 * t80, t1 * t18 + t16 * t20 + t19 * t2 + t3 * t5 + t37 * t8 + t4 * t7; 0, 0, 0, -t147, -t111 * t168, 0, 0, 0, 0, 0, -t110 * t147 + (-t56 - t130) * qJD(4), t107 * t147 + (-t55 - t131) * qJD(4), 0, 0, 0, 0, 0, t136 * t57 - t26 * t94 + t56 * t80 + t63 * t68, t123 * t136 + t25 * t94 + t56 * t82 - t62 * t68, t123 * t63 - t25 * t80 - t26 * t82 + t57 * t62, t1 * t57 - t123 * t2 + t20 * t56 + t25 * t7 + t26 * t5 + t68 * t8; 0, 0, 0, t198 (-t202 + t48) * qJD(3), 0.2e1 * t110 * t136, -0.2e1 * t182 * t151, t112 * t110, -t112 * t107, 0, t107 * t125 - t126 * t110, t126 * t107 + t110 * t125, t82 * t140 + (-t109 * t62 - t82 * t155) * t107 (-t106 * t82 - t109 * t80) * t156 + (t173 - t109 * t63 + (t106 * t80 - t109 * t82) * qJD(5)) * t107, t94 * t139 + t110 * t62 + (t107 * t82 + t127 * t109) * qJD(4), t94 * t138 + t110 * t63 + (-t127 * t106 - t107 * t80) * qJD(4) (-t94 - t160) * t157 (t155 * t88 - t201) * t94 + ((pkin(9) * t80 + t180) * qJD(4) + (t179 + (pkin(9) * t94 + t28) * t109) * qJD(5) + t135) * t110 + (t27 * t154 + pkin(9) * t63 + t175 - t48 * t80 + ((-pkin(9) * t166 + t79) * qJD(3) + t9) * qJD(4)) * t107, t200 * t94 + (t27 * t153 + (qJD(4) * t82 - t148) * pkin(9) + t121) * t110 + (-t27 * t155 - pkin(9) * t62 + t174 - t48 * t82 + (-pkin(9) * t177 - t183 * qJD(3) - t10) * qJD(4)) * t107, t52 * t62 - t59 * t63 - t190 * t82 - t188 * t80 + (-t106 * t7 - t109 * t5) * t156 + (-t1 * t109 - t106 * t2 + (t106 * t5 - t109 * t7) * qJD(5)) * t107, t1 * t52 + t2 * t59 + t188 * t7 + t190 * t5 + t20 * t145 * t156 + (t8 * t145 + (pkin(5) * t154 - t48) * t20) * t107; 0, 0, 0, 0, 0, -t107 * t113 * t110, t182 * t113, 0, 0, 0, qJD(4) * t31 - t46 * t162 - t13 (-qJD(3) * t46 - t44) * t110, -t82 * t177 - t173 (-t62 + t193) * t109 + (-t63 + t192) * t106, -t94 * t154 + (t94 * t164 + (-t82 + t158) * t107) * qJD(3), t148 + (-t94 * t166 + (t80 + t153) * t107) * qJD(3), t94 * t162, -pkin(4) * t63 - t174 + t134 * t94 - t31 * t80 + (pkin(10) * t177 + t180) * qJD(5) + (-t9 * t107 + (-pkin(10) * t157 - t110 * t27) * t106) * qJD(3), pkin(4) * t62 + t175 - t189 * t94 - t31 * t82 + (-pkin(10) * t178 + t109 * t27) * qJD(5) + (-t27 * t164 + (-pkin(10) * t153 + t10) * t107) * qJD(3), t62 * t89 + t63 * t90 - t186 * t82 - t187 * t80 + (t5 * t94 + t2) * t109 + (t7 * t94 - t1) * t106, -t2 * t90 + t1 * t89 + t8 * (-pkin(5) * t109 - pkin(4)) + t187 * t7 + t186 * t5 + (-pkin(5) * t178 - t31) * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82 * t80, -t77 + t196, -t62 - t193, -t192 - t63, t136, -t10 * t94 - t27 * t82 + t114, t27 * t80 - t9 * t94 - t121, pkin(5) * t62 - t195 * t80, t195 * t7 + (-t20 * t82 + t1) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t77 - t196, t5 * t82 + t7 * t80 + t8;];
tauc_reg  = t11;
