% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
% 
% Output:
% tauc_reg [6x25]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPPRPR4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:47:08
% EndTime: 2019-03-09 01:47:12
% DurationCPUTime: 1.18s
% Computational Cost: add. (1943->214), mult. (4062->315), div. (0->0), fcn. (2705->8), ass. (0->134)
t93 = sin(qJ(4));
t151 = qJD(1) * t93;
t152 = cos(pkin(10));
t95 = cos(qJ(4));
t130 = t152 * t95;
t69 = qJD(1) * t130;
t88 = sin(pkin(10));
t55 = t88 * t151 - t69;
t143 = qJD(6) - t55;
t94 = cos(qJ(6));
t144 = t94 * qJD(4);
t131 = t152 * t93;
t61 = t88 * t95 + t131;
t57 = t61 * qJD(1);
t92 = sin(qJ(6));
t35 = -t92 * t57 - t144;
t180 = t143 * t35;
t128 = t143 * t94;
t56 = t61 * qJD(4);
t50 = qJD(1) * t56;
t163 = t92 * t50;
t179 = t128 * t143 - t163;
t89 = sin(pkin(9));
t97 = qJD(4) ^ 2;
t98 = qJD(1) ^ 2;
t177 = (t97 + t98) * t89;
t90 = cos(pkin(9));
t150 = qJD(2) * t90;
t123 = -qJD(5) + t150;
t114 = t123 * t93;
t140 = qJ(5) * qJD(4);
t145 = t93 * qJD(3);
t142 = qJ(2) * qJD(1);
t96 = -pkin(1) - pkin(2);
t67 = t96 * qJD(1) + qJD(2);
t53 = t90 * t142 + t89 * t67;
t46 = -qJD(1) * pkin(7) + t53;
t176 = (-t95 * t46 - t145) * qJD(4) + (t95 * t140 - t114) * qJD(1);
t175 = t95 * qJD(3) - t93 * t46;
t174 = -qJD(6) + t143;
t149 = qJD(4) * t93;
t58 = qJD(4) * t130 - t88 * t149;
t113 = t123 * t95;
t16 = t175 * qJD(4) + (t93 * t140 + t113) * qJD(1);
t3 = -t152 * t176 + t88 * t16;
t74 = t88 * pkin(4) + pkin(8);
t173 = (-pkin(4) * t151 - t57 * pkin(5) - t55 * pkin(8) + qJD(6) * t74) * t143 + t3;
t148 = qJD(6) * t61;
t172 = -(t94 * t148 + t92 * t58) * t143 + t61 * t163;
t164 = t88 * t93;
t106 = t130 - t164;
t156 = t90 * qJ(2) + t89 * t96;
t63 = -pkin(7) + t156;
t153 = qJ(5) - t63;
t125 = qJD(4) * t153;
t101 = t95 * t125 - t114;
t31 = t93 * t125 + t113;
t11 = t88 * t101 + t152 * t31;
t47 = t153 * t95;
t19 = -t152 * t47 + t153 * t164;
t119 = t19 * t50 - t3 * t61;
t52 = -t89 * t142 + t90 * t67;
t45 = qJD(1) * pkin(3) - t52;
t85 = t95 * pkin(4);
t34 = qJD(1) * t85 + qJD(5) + t45;
t13 = -t55 * pkin(5) + t57 * pkin(8) + t34;
t132 = -t89 * qJ(2) + t90 * t96;
t62 = pkin(3) - t132;
t115 = t62 + t85;
t23 = pkin(5) * t106 + t61 * pkin(8) + t115;
t4 = t152 * t16 + t176 * t88;
t141 = qJ(5) * qJD(1);
t33 = t145 + (t46 - t141) * t95;
t165 = t88 * t33;
t32 = t93 * t141 + t175;
t29 = qJD(4) * pkin(4) + t32;
t7 = t152 * t29 - t165;
t5 = -qJD(4) * pkin(5) - t7;
t171 = -(qJD(6) * t23 + t11) * t143 - (qJD(6) * t13 + t4) * t106 - t5 * t58 + t119;
t139 = qJD(1) * qJD(2);
t70 = t89 * t139;
t170 = 0.2e1 * t70;
t147 = qJD(6) * t92;
t138 = qJD(1) * qJD(4);
t133 = t93 * t138;
t51 = qJD(4) * t69 - t88 * t133;
t20 = qJD(6) * t144 + t57 * t147 - t94 * t51;
t169 = t20 * t92;
t168 = t23 * t50;
t37 = t92 * qJD(4) - t94 * t57;
t167 = t37 * t57;
t166 = t57 * t35;
t162 = t92 * t51;
t39 = t94 * t50;
t160 = t97 * t93;
t159 = t97 * t95;
t27 = t152 * t33;
t8 = t88 * t29 + t27;
t158 = -t90 * t57 + t58 * t89;
t157 = -t106 * t90 * qJD(1) - t89 * t56;
t155 = t93 ^ 2 - t95 ^ 2;
t146 = t89 * qJD(1);
t135 = 0.2e1 * t139;
t134 = 0.2e1 * t138;
t124 = t95 * t134;
t122 = -pkin(4) * t149 + t89 * qJD(2);
t121 = -qJD(6) * t90 + t157;
t6 = qJD(4) * pkin(8) + t8;
t1 = t94 * t13 - t92 * t6;
t2 = t92 * t13 + t94 * t6;
t118 = t106 * t20 - t56 * t37;
t21 = t37 * qJD(6) - t162;
t117 = -t106 * t21 + t56 * t35;
t116 = t52 * t89 - t53 * t90;
t112 = (t45 - t150) * qJD(1);
t49 = t106 * t89;
t111 = qJD(6) * t49 + t146;
t110 = -t39 + (t55 * t92 - t147) * t143;
t59 = -pkin(4) * t133 + t70;
t108 = t61 * t147 - t94 * t58;
t107 = -t63 * t97 + t170;
t103 = qJD(4) * (-qJD(1) * t62 - t150 - t45);
t12 = t152 * t32 - t165;
t102 = t74 * t50 + (t12 + t5) * t143;
t100 = t108 * t143 + t61 * t39;
t75 = -t152 * pkin(4) - pkin(5);
t48 = t61 * t89;
t22 = -t56 * pkin(5) + t58 * pkin(8) + t122;
t18 = -t153 * t131 - t88 * t47;
t17 = -t50 * pkin(5) + t51 * pkin(8) + t59;
t15 = t94 * t17;
t10 = t88 * t32 + t27;
t9 = -t152 * t101 + t88 * t31;
t14 = [0, 0, 0, 0, t135, qJ(2) * t135, t170, t90 * t135 ((-t89 * t132 + t90 * t156) * qJD(1) - t116) * qJD(2), t93 * t124, -t155 * t134, -t159, t160, 0, t93 * t103 + t107 * t95, t95 * t103 - t107 * t93, -t106 * t4 + t11 * t55 - t18 * t51 + t8 * t56 - t9 * t57 + t7 * t58 + t119, t8 * t11 + t59 * t115 + t34 * t122 + t3 * t18 + t4 * t19 - t7 * t9, -t20 * t94 * t61 + t108 * t37 (t35 * t94 + t37 * t92) * t58 + (t169 + t21 * t94 + (-t35 * t92 + t37 * t94) * qJD(6)) * t61, t100 + t118, t117 - t172, -t106 * t50 - t143 * t56, -t1 * t56 + t15 * t106 + t18 * t21 + t9 * t35 + (t22 * t143 - t168 + (-t106 * t6 - t143 * t19 - t5 * t61) * qJD(6)) * t94 + t171 * t92, t18 * t20 + t2 * t56 + t9 * t37 + (-(-qJD(6) * t19 + t22) * t143 + t168 - (-qJD(6) * t6 + t17) * t106 + t5 * t148) * t92 + t171 * t94; 0, 0, 0, 0, -t98, -t98 * qJ(2), -t89 * t98, -t90 * t98, t116 * qJD(1), 0, 0, 0, 0, 0, 0.2e1 * t90 * t133 - t95 * t177, t90 * t124 + t93 * t177, t157 * t55 - t158 * t57 - t48 * t51 + t49 * t50, -t34 * t146 + t157 * t8 - t158 * t7 + t3 * t48 + t4 * t49 - t59 * t90, 0, 0, 0, 0, 0 -(-t92 * t49 - t94 * t90) * t50 + t48 * t21 - (t111 * t94 + t121 * t92) * t143 + t158 * t35 (t94 * t49 - t92 * t90) * t50 + t48 * t20 - (-t111 * t92 + t121 * t94) * t143 + t158 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t160, -t159, t106 * t51 + t61 * t50 + t58 * t55 - t56 * t57, -t106 * t3 + t4 * t61 - t7 * t56 + t8 * t58, 0, 0, 0, 0, 0, t117 + t172, t100 - t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t93 * t98 * t95, t155 * t98, 0, 0, 0, t93 * t112, t95 * t112 (t10 - t8) * t57 + (-t12 + t7) * t55 + (t152 * t51 + t50 * t88) * pkin(4), t7 * t10 - t8 * t12 + (t34 * t151 - t152 * t3 + t4 * t88) * pkin(4), t37 * t128 + t169 (t20 - t180) * t94 + (-t143 * t37 - t21) * t92, t167 + t179, t110 - t166, t143 * t57, t1 * t57 - t10 * t35 + t102 * t92 - t173 * t94 + t75 * t21, -t10 * t37 + t102 * t94 + t173 * t92 - t2 * t57 + t75 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55 ^ 2 - t57 ^ 2, -t8 * t55 - t7 * t57 + t59, 0, 0, 0, 0, 0, t110 + t166, t167 - t179; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37 * t35, -t35 ^ 2 + t37 ^ 2, t20 + t180, t174 * t37 + t162, -t50, t174 * t2 - t5 * t37 - t92 * t4 + t15, t174 * t1 - t92 * t17 + t5 * t35 - t94 * t4;];
tauc_reg  = t14;
