% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% 
% Output:
% tauc_reg [6x26]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPPRRP7_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP7_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP7_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP7_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:13:54
% EndTime: 2019-03-09 02:13:58
% DurationCPUTime: 1.55s
% Computational Cost: add. (2703->250), mult. (5999->346), div. (0->0), fcn. (4203->6), ass. (0->130)
t101 = sin(pkin(9));
t102 = cos(pkin(9));
t105 = sin(qJ(4));
t107 = cos(qJ(4));
t77 = t101 * t107 + t102 * t105;
t112 = qJD(1) * t77;
t106 = cos(qJ(5));
t142 = qJD(5) * t106;
t104 = sin(qJ(5));
t146 = qJD(4) * t104;
t147 = qJD(1) * t101;
t133 = t105 * t147;
t150 = t107 * t102;
t135 = qJD(1) * t150;
t180 = t133 - t135;
t24 = -t180 * t142 + (qJD(5) - t112) * t146;
t177 = qJD(5) + t112;
t128 = t106 * t177;
t144 = qJD(4) * t107;
t134 = t102 * t144;
t84 = qJD(4) * t133;
t63 = qJD(1) * t134 - t84;
t156 = t104 * t63;
t179 = -t128 * t177 - t156;
t103 = -pkin(1) - qJ(3);
t174 = t103 * qJD(1);
t85 = qJD(2) + t174;
t132 = -pkin(7) * qJD(1) + t85;
t64 = t132 * t101;
t65 = t132 * t102;
t175 = -t105 * t64 + t107 * t65;
t164 = -pkin(7) + t103;
t78 = t164 * t101;
t79 = t164 * t102;
t48 = t105 * t78 - t107 * t79;
t158 = t101 ^ 2 + t102 ^ 2;
t173 = t158 * qJD(3);
t35 = t105 * t65 + t107 * t64;
t31 = qJD(4) * pkin(8) + t35;
t100 = qJD(1) * qJ(2);
t94 = qJD(3) + t100;
t80 = pkin(3) * t147 + t94;
t32 = pkin(4) * t112 + pkin(8) * t180 + t80;
t13 = t104 * t32 + t106 * t31;
t114 = t77 * qJD(3);
t17 = -qJD(1) * t114 + qJD(4) * t175;
t113 = t77 * qJD(4);
t110 = qJD(1) * t113;
t99 = qJD(1) * qJD(2);
t33 = t63 * pkin(4) + pkin(8) * t110 + t99;
t28 = t106 * t33;
t109 = -qJD(5) * t13 - t104 * t17 + t28;
t141 = t106 * qJD(4);
t143 = qJD(5) * t104;
t23 = -qJD(5) * t141 + t106 * t110 - t143 * t180;
t54 = -t106 * t180 + t146;
t1 = pkin(5) * t63 + qJ(6) * t23 - qJD(6) * t54 + t109;
t52 = -t104 * t180 - t141;
t7 = -qJ(6) * t52 + t13;
t172 = t177 * t7 + t1;
t171 = t54 ^ 2;
t139 = 0.2e1 * t99;
t12 = -t104 * t31 + t106 * t32;
t6 = -qJ(6) * t54 + t12;
t5 = pkin(5) * t177 + t6;
t170 = t5 - t6;
t163 = -qJ(6) - pkin(8);
t130 = qJD(5) * t163;
t151 = qJ(6) * t106;
t46 = -pkin(4) * t180 + pkin(8) * t112;
t42 = t106 * t46;
t169 = pkin(5) * t180 + t106 * t130 - t112 * t151 - t42 + (-qJD(6) + t175) * t104;
t168 = t52 * t112;
t167 = t52 * t180;
t166 = t54 * t180;
t76 = t101 * t105 - t150;
t165 = t76 * t63;
t162 = -t104 * t24 - t142 * t52;
t161 = t104 * t46 + t106 * t175;
t90 = pkin(3) * t101 + qJ(2);
t45 = pkin(4) * t77 + pkin(8) * t76 + t90;
t49 = t105 * t79 + t107 * t78;
t47 = t106 * t49;
t160 = t104 * t45 + t47;
t152 = qJ(6) * t104;
t159 = t106 * qJD(6) + t104 * t130 - t112 * t152 - t161;
t157 = t104 * t54;
t59 = t106 * t63;
t153 = t23 * t104;
t145 = qJD(4) * t105;
t74 = -t101 * t145 + t134;
t149 = t74 * qJD(4);
t25 = -qJD(4) * t48 - t114;
t73 = -t101 * t144 - t102 * t145;
t43 = pkin(4) * t74 - pkin(8) * t73 + qJD(2);
t140 = t104 * t43 + t106 * t25 + t142 * t45;
t138 = t76 * t143;
t136 = t76 * t142;
t131 = qJD(1) * t158;
t129 = t104 * t177;
t127 = qJD(5) * t77 + qJD(1);
t111 = qJD(3) * t180 - t64 * t144 - t65 * t145;
t126 = pkin(8) * qJD(5) * t177 - t111;
t115 = t104 * t33 + t106 * t17 + t142 * t32 - t143 * t31;
t2 = -qJ(6) * t24 - qJD(6) * t52 + t115;
t125 = -t177 * t5 + t2;
t30 = -qJD(4) * pkin(4) - t175;
t124 = t111 * t76 + t30 * t73;
t123 = t23 * t76 + t54 * t73;
t122 = -t23 * t77 + t54 * t74;
t121 = -t24 * t77 - t52 * t74;
t120 = -t177 * t73 + t165;
t118 = -qJ(6) * t73 + qJD(6) * t76;
t117 = t59 + (-t104 * t112 - t143) * t177;
t116 = -pkin(8) * t63 + t177 * t30;
t9 = pkin(5) * t24 - t111;
t26 = -qJD(3) * t76 + qJD(4) * t49;
t108 = qJD(1) ^ 2;
t82 = t163 * t106;
t81 = t163 * t104;
t66 = t73 * qJD(4);
t51 = t52 ^ 2;
t40 = t106 * t45;
t38 = t106 * t43;
t15 = pkin(5) * t52 + qJD(6) + t30;
t14 = t152 * t76 + t160;
t11 = pkin(5) * t77 - t104 * t49 + t151 * t76 + t40;
t4 = qJ(6) * t136 + (-qJD(5) * t49 + t118) * t104 + t140;
t3 = pkin(5) * t74 - t104 * t25 + t38 + t118 * t106 + (-t47 + (-qJ(6) * t76 - t45) * t104) * qJD(5);
t8 = [0, 0, 0, 0, t139, qJ(2) * t139, t101 * t139, t102 * t139, 0.2e1 * qJD(3) * t131 (t94 + t100) * qJD(2) + (-t85 - t174) * t173, t110 * t76 - t180 * t73, t110 * t77 - t112 * t73 + t180 * t74 + t165, t66, -t149, 0, 0.2e1 * qJD(2) * t112 - t26 * qJD(4) + t90 * t63 + t80 * t74, -qJD(2) * t180 - t25 * qJD(4) + t80 * t73 + (-qJD(2) * t76 - t113 * t90) * qJD(1), t106 * t123 + t138 * t54 (-t106 * t52 - t157) * t73 + (-t153 + t106 * t24 + (-t104 * t52 + t106 * t54) * qJD(5)) * t76, -t106 * t120 + t138 * t177 + t122, t104 * t120 + t136 * t177 + t121, t177 * t74 + t63 * t77 (-t49 * t142 + t38) * t177 + t40 * t63 + (-t31 * t142 + t28) * t77 + t12 * t74 + t26 * t52 + t48 * t24 - t30 * t136 + ((-qJD(5) * t45 - t25) * t177 - t49 * t63 + (-qJD(5) * t32 - t17) * t77 + t124) * t104 -(-t143 * t49 + t140) * t177 - t160 * t63 - t115 * t77 - t13 * t74 + t26 * t54 - t48 * t23 + t30 * t138 + t124 * t106, t11 * t23 - t14 * t24 - t3 * t54 - t4 * t52 + (-t104 * t7 - t106 * t5) * t73 + (t1 * t106 + t104 * t2 + (-t104 * t5 + t106 * t7) * qJD(5)) * t76, t2 * t14 + t7 * t4 + t1 * t11 + t5 * t3 + t9 * (-pkin(5) * t104 * t76 + t48) + t15 * ((t104 * t73 - t136) * pkin(5) + t26); 0, 0, 0, 0, -t108, -t108 * qJ(2), -t108 * t101, -t108 * t102, 0 (-t94 - t173) * qJD(1), 0, 0, 0, 0, 0, -qJD(1) * t112 + t66, qJD(1) * t180 - t149, 0, 0, 0, 0, 0, -t77 * t156 + t24 * t76 - t52 * t73 + (-t104 * t74 - t106 * t127) * t177, -t77 * t59 + (t104 * t127 - t106 * t74) * t177 - t123 (t127 * t54 + t121) * t106 + (t127 * t52 + t122) * t104, -t15 * t73 + t76 * t9 + (-t127 * t5 + t2 * t77 + t7 * t74) * t106 + (-t1 * t77 - t127 * t7 - t5 * t74) * t104; 0, 0, 0, 0, 0, 0, 0, 0, -t158 * t108, t131 * t85 + t99, 0, 0, 0, 0, 0, -t84 + (-t180 + t135) * qJD(4), -0.2e1 * t112 * qJD(4), 0, 0, 0, 0, 0, t117 + t167, t166 + t179 (t23 - t168) * t106 + t54 * t129 + t162, t104 * t125 + t106 * t172 + t15 * t180; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t180 * t112, -t112 ^ 2 + t180 ^ 2, 0, t84 + (-t180 - t135) * qJD(4), 0, qJD(4) * t35 + t180 * t80 + t111 (qJD(3) + t80) * t112, t128 * t54 - t153 (-t23 - t168) * t106 - t177 * t157 + t162, t166 - t179, t117 - t167, t177 * t180, -pkin(4) * t24 + t12 * t180 - t35 * t52 - t42 * t177 - t126 * t106 + (t175 * t177 + t116) * t104, pkin(4) * t23 + t104 * t126 + t106 * t116 - t13 * t180 + t161 * t177 - t35 * t54, -t104 * t172 + t106 * t125 - t159 * t52 - t169 * t54 + t23 * t81 + t24 * t82, -t2 * t82 + t1 * t81 + t9 * (-pkin(5) * t106 - pkin(4)) + t159 * t7 + t169 * t5 + (pkin(5) * t129 - t35) * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54 * t52, -t51 + t171, t177 * t52 - t23, t177 * t54 - t24, t63, t13 * t177 - t30 * t54 + t109, t12 * t177 + t30 * t52 - t115, pkin(5) * t23 - t170 * t52, t170 * t7 + (-t15 * t54 + t1) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51 - t171, t5 * t54 + t52 * t7 + t9;];
tauc_reg  = t8;
