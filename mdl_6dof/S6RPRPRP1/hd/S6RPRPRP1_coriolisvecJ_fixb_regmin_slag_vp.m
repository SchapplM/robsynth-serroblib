% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% 
% Output:
% tauc_reg [6x22]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRPRP1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:03:03
% EndTime: 2019-03-09 03:03:07
% DurationCPUTime: 1.45s
% Computational Cost: add. (3003->257), mult. (7132->352), div. (0->0), fcn. (4902->8), ass. (0->143)
t100 = sin(pkin(9)) * pkin(1) + pkin(7);
t161 = qJ(4) + t100;
t112 = sin(qJ(5));
t108 = sin(pkin(10));
t110 = cos(pkin(10));
t115 = cos(qJ(3));
t158 = qJD(1) * t115;
t113 = sin(qJ(3));
t159 = qJD(1) * t113;
t82 = -t108 * t159 + t110 * t158;
t80 = qJD(5) - t82;
t142 = t112 * t80;
t114 = cos(qJ(5));
t92 = t108 * t115 + t110 * t113;
t84 = t92 * qJD(1);
t68 = qJD(3) * t112 + t114 * t84;
t190 = t68 * t142;
t141 = t161 * qJD(1);
t70 = t115 * qJD(2) - t141 * t113;
t152 = qJD(1) * qJD(4);
t71 = t113 * qJD(2) + t141 * t115;
t189 = -t71 * qJD(3) - t113 * t152;
t169 = t110 * t71;
t170 = qJD(3) * pkin(3);
t63 = t70 + t170;
t25 = t108 * t63 + t169;
t22 = qJD(3) * pkin(8) + t25;
t102 = -cos(pkin(9)) * pkin(1) - pkin(2);
t127 = -pkin(3) * t115 + t102;
t122 = t127 * qJD(1);
t81 = qJD(4) + t122;
t34 = -pkin(4) * t82 - pkin(8) * t84 + t81;
t12 = t112 * t34 + t114 * t22;
t57 = qJD(3) * t70 + t115 * t152;
t18 = t108 * t189 + t110 * t57;
t153 = qJD(1) * qJD(3);
t146 = t115 * t153;
t147 = t113 * t153;
t125 = -t108 * t147 + t110 * t146;
t83 = t92 * qJD(3);
t78 = qJD(1) * t83;
t97 = pkin(3) * t147;
t39 = t78 * pkin(4) - t125 * pkin(8) + t97;
t36 = t114 * t39;
t119 = -t12 * qJD(5) - t112 * t18 + t36;
t154 = t114 * qJD(3);
t157 = qJD(5) * t112;
t37 = -qJD(5) * t154 - t114 * t125 + t84 * t157;
t1 = pkin(5) * t78 + qJ(6) * t37 - qJD(6) * t68 + t119;
t66 = t112 * t84 - t154;
t7 = -qJ(6) * t66 + t12;
t188 = t80 * t7 + t1;
t177 = t92 * t78;
t91 = t108 * t113 - t110 * t115;
t86 = t91 * qJD(3);
t131 = -t80 * t86 + t177;
t149 = t92 * t157;
t187 = -t131 * t114 + t80 * t149;
t186 = t68 ^ 2;
t11 = -t112 * t22 + t114 * t34;
t6 = -qJ(6) * t68 + t11;
t5 = pkin(5) * t80 + t6;
t185 = t5 - t6;
t99 = pkin(3) * t108 + pkin(8);
t167 = qJ(6) + t99;
t144 = qJD(5) * t167;
t165 = qJ(6) * t114;
t60 = t108 * t71;
t27 = t110 * t70 - t60;
t48 = pkin(3) * t159 + pkin(4) * t84 - pkin(8) * t82;
t42 = t114 * t48;
t184 = -pkin(5) * t84 - t114 * t144 + t82 * t165 - t42 + (-qJD(6) + t27) * t112;
t17 = t108 * t57 - t110 * t189;
t183 = t17 * t92;
t24 = t110 * t63 - t60;
t21 = -qJD(3) * pkin(4) - t24;
t182 = t21 * t86;
t181 = t66 * t82;
t180 = t66 * t84;
t179 = t68 * t84;
t178 = t68 * t86;
t121 = t112 * t125;
t38 = qJD(5) * t68 + t121;
t176 = (-t38 * t92 + t66 * t86) * t114;
t175 = t112 * t48 + t114 * t27;
t156 = qJD(5) * t114;
t174 = -t112 * t38 - t66 * t156;
t173 = -t37 * t91 + t68 * t83;
t88 = t161 * t113;
t90 = t161 * t115;
t55 = -t108 * t88 + t110 * t90;
t51 = t114 * t55;
t52 = pkin(4) * t91 - pkin(8) * t92 + t127;
t172 = t112 * t52 + t51;
t166 = qJ(6) * t112;
t171 = qJD(6) * t114 - t112 * t144 + t82 * t166 - t175;
t168 = t112 * t78;
t164 = qJD(5) * t66;
t116 = qJD(3) ^ 2;
t163 = t116 * t113;
t162 = t116 * t115;
t160 = t113 ^ 2 - t115 ^ 2;
t94 = qJD(1) * t102;
t140 = qJD(3) * t161;
t72 = qJD(4) * t115 - t113 * t140;
t73 = -qJD(4) * t113 - t115 * t140;
t33 = t108 * t73 + t110 * t72;
t150 = t113 * t170;
t49 = pkin(4) * t83 + pkin(8) * t86 + t150;
t151 = t112 * t49 + t114 * t33 + t52 * t156;
t148 = t92 * t156;
t101 = -pkin(3) * t110 - pkin(4);
t26 = t108 * t70 + t169;
t32 = t108 * t72 - t110 * t73;
t54 = t108 * t90 + t110 * t88;
t143 = t114 * t80;
t138 = t68 * t148;
t137 = qJD(5) * t80 * t99 + t17;
t124 = t112 * t39 + t114 * t18 + t34 * t156 - t22 * t157;
t2 = -qJ(6) * t38 - qJD(6) * t66 + t124;
t136 = -t80 * t5 + t2;
t135 = -t112 * t7 - t114 * t5;
t134 = t112 * t5 - t114 * t7;
t133 = -t55 * t78 + t183;
t132 = -t38 * t91 - t66 * t83;
t130 = qJ(6) * t86 - qJD(6) * t92;
t128 = 0.2e1 * qJD(3) * t94;
t126 = t114 * t78 + t82 * t142 - t80 * t157;
t9 = pkin(5) * t38 + t17;
t123 = t80 * t21 - t99 * t78;
t118 = -t131 * t112 - t80 * t148;
t117 = qJD(1) ^ 2;
t89 = t167 * t114;
t87 = t167 * t112;
t65 = t66 ^ 2;
t47 = t114 * t52;
t43 = t114 * t49;
t15 = pkin(5) * t66 + qJD(6) + t21;
t14 = -t92 * t166 + t172;
t13 = pkin(5) * t91 - t112 * t55 - t92 * t165 + t47;
t4 = -qJ(6) * t148 + (-qJD(5) * t55 + t130) * t112 + t151;
t3 = pkin(5) * t83 - t112 * t33 + t43 + t130 * t114 + (-t51 + (qJ(6) * t92 - t52) * t112) * qJD(5);
t8 = [0, 0, 0, 0, 0.2e1 * t113 * t146, -0.2e1 * t160 * t153, t162, -t163, 0, -t100 * t162 + t113 * t128, t100 * t163 + t115 * t128, t54 * t125 - t18 * t91 + t24 * t86 - t25 * t83 + t32 * t84 + t33 * t82 + t133, t17 * t54 + t18 * t55 - t24 * t32 + t25 * t33 + (t81 + t122) * t150, -t68 * t149 + (-t37 * t92 - t178) * t114, -t138 + (t178 + (t37 + t164) * t92) * t112 + t176, t173 - t187, t118 + t132, t78 * t91 + t80 * t83 (-t55 * t156 + t43) * t80 + t47 * t78 + (-t22 * t156 + t36) * t91 + t11 * t83 + t32 * t66 + t54 * t38 + t21 * t148 + ((-qJD(5) * t52 - t33) * t80 + (-qJD(5) * t34 - t18) * t91 - t182 + t133) * t112 -(-t55 * t157 + t151) * t80 - t172 * t78 - t124 * t91 - t12 * t83 + t32 * t68 - t54 * t37 - t21 * t149 + (-t182 + t183) * t114, t13 * t37 - t14 * t38 - t3 * t68 - t4 * t66 - t135 * t86 + (t134 * qJD(5) - t1 * t114 - t112 * t2) * t92, t2 * t14 + t7 * t4 + t1 * t13 + t5 * t3 + t9 * (pkin(5) * t112 * t92 + t54) + t15 * ((-t112 * t86 + t148) * pkin(5) + t32); 0, 0, 0, 0, 0, 0, 0, 0, 0, -t163, -t162, t91 * t125 - t86 * t82 + t83 * t84 - t177, t17 * t91 + t18 * t92 - t24 * t83 - t25 * t86, 0, 0, 0, 0, 0, t118 - t132, t173 + t187, t138 + (-t178 + (-t37 + t164) * t92) * t112 + t176, t15 * t83 + t9 * t91 + t134 * t86 + (t135 * qJD(5) - t1 * t112 + t114 * t2) * t92; 0, 0, 0, 0, -t113 * t117 * t115, t160 * t117, 0, 0, 0, -t94 * t159, -t94 * t158 (t25 - t26) * t84 + (-t27 + t24) * t82 + (-t108 * t78 - t110 * t125) * pkin(3), t24 * t26 - t25 * t27 + (t108 * t18 - t110 * t17 - t81 * t159) * pkin(3), -t112 * t37 + t68 * t143 (-t37 + t181) * t114 - t190 + t174, t80 * t143 + t168 - t179, t126 + t180, -t80 * t84, t101 * t38 - t11 * t84 - t26 * t66 - t42 * t80 - t137 * t114 + (t27 * t80 + t123) * t112, -t101 * t37 + t137 * t112 + t123 * t114 + t12 * t84 + t175 * t80 - t26 * t68, -t112 * t188 + t136 * t114 - t171 * t66 - t184 * t68 - t37 * t87 - t38 * t89, t2 * t89 - t1 * t87 + t9 * (-pkin(5) * t114 + t101) + t171 * t7 + t184 * t5 + (pkin(5) * t142 - t26) * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t82 ^ 2 - t84 ^ 2, t24 * t84 - t25 * t82 + t97, 0, 0, 0, 0, 0, t126 - t180, -t80 ^ 2 * t114 - t168 - t179 (t37 + t181) * t114 + t190 + t174, t136 * t112 + t114 * t188 - t15 * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68 * t66, -t65 + t186, t66 * t80 - t37, -t121 + (-qJD(5) + t80) * t68, t78, t12 * t80 - t21 * t68 + t119, t11 * t80 + t21 * t66 - t124, pkin(5) * t37 - t185 * t66, t185 * t7 + (-t15 * t68 + t1) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65 - t186, t5 * t68 + t66 * t7 + t9;];
tauc_reg  = t8;
