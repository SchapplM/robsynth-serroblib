% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% tauc_reg [6x21]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRPRRP1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:58:29
% EndTime: 2019-03-08 19:58:34
% DurationCPUTime: 1.68s
% Computational Cost: add. (1952->254), mult. (4929->366), div. (0->0), fcn. (3686->10), ass. (0->146)
t106 = cos(qJ(2));
t98 = sin(pkin(6));
t159 = qJD(1) * t98;
t139 = t106 * t159;
t103 = sin(qJ(2));
t140 = t103 * t159;
t99 = cos(pkin(11));
t84 = t99 * t140;
t97 = sin(pkin(11));
t53 = t97 * t139 + t84;
t102 = sin(qJ(4));
t105 = cos(qJ(4));
t123 = pkin(4) * t102 - pkin(9) * t105;
t80 = t123 * qJD(4);
t193 = -t53 + t80;
t58 = (t103 * t97 - t106 * t99) * t98;
t104 = cos(qJ(5));
t149 = qJD(4) * t102;
t101 = sin(qJ(5));
t157 = t101 * t105;
t91 = pkin(2) * t97 + pkin(8);
t169 = t101 * t91;
t83 = t97 * t140;
t56 = t99 * t139 - t83;
t192 = t104 * t193 + t149 * t169 + t56 * t157;
t146 = qJD(5) * t104;
t155 = t104 * t105;
t116 = -pkin(4) * t105 - pkin(9) * t102 - pkin(3);
t186 = pkin(2) * t99;
t71 = t116 - t186;
t191 = t101 * t193 + t71 * t146 - t56 * t155;
t82 = qJD(2) * pkin(2) + t139;
t48 = t97 * t82 + t84;
t45 = qJD(2) * pkin(8) + t48;
t100 = cos(pkin(6));
t89 = qJD(1) * t100 + qJD(3);
t190 = -t102 * t45 + t105 * t89;
t151 = qJD(2) * t105;
t90 = -qJD(5) + t151;
t95 = t102 ^ 2;
t118 = qJD(2) * t95 - t105 * t90;
t147 = qJD(5) * t101;
t136 = t102 * t147;
t145 = t104 * qJD(4);
t189 = -t118 * t145 - t90 * t136;
t135 = t102 * t146;
t142 = qJD(4) * qJD(5);
t148 = qJD(4) * t105;
t52 = (t101 * t148 + t135) * qJD(2) + t101 * t142;
t150 = qJD(4) * t101;
t152 = qJD(2) * t102;
t76 = t104 * t152 + t150;
t188 = t76 ^ 2;
t30 = t102 * t89 + t105 * t45;
t27 = qJD(4) * pkin(9) + t30;
t47 = t99 * t82 - t83;
t35 = t116 * qJD(2) - t47;
t9 = -t101 * t27 + t104 * t35;
t6 = -qJ(6) * t76 + t9;
t5 = -pkin(5) * t90 + t6;
t187 = t5 - t6;
t74 = t101 * t152 - t145;
t185 = t74 * t90;
t184 = t76 * t90;
t183 = -qJ(6) - pkin(9);
t115 = pkin(5) * t102 - qJ(6) * t155;
t144 = t104 * qJD(6);
t160 = qJ(6) * t102;
t78 = t91 * t155;
t182 = -t102 * t144 + t115 * qJD(4) + (-t78 + (-t71 + t160) * t101) * qJD(5) + t192;
t156 = t102 * t104;
t181 = (-qJ(6) * qJD(5) - qJD(4) * t91) * t156 + (-qJD(6) * t102 + (-qJ(6) * qJD(4) - qJD(5) * t91) * t105) * t101 + t191;
t79 = t123 * qJD(2);
t180 = t101 * t79 + t104 * t190;
t137 = t105 * t145;
t179 = -t74 * t137 - t52 * t156;
t127 = qJD(5) * t183;
t177 = t144 - t180 + (qJ(6) * t151 + t127) * t101;
t129 = -t101 * t190 + t104 * t79;
t176 = -t115 * qJD(2) - t101 * qJD(6) + t104 * t127 - t129;
t174 = t101 * t71 + t78;
t173 = -t105 ^ 2 + t95;
t26 = -qJD(4) * pkin(4) - t190;
t172 = t101 * t26;
t171 = t101 * t35;
t170 = t101 * t90;
t168 = t102 * t74;
t167 = t104 * t26;
t166 = t104 * t90;
t165 = t105 * t52;
t108 = qJD(2) ^ 2;
t163 = t108 * t98;
t55 = qJD(2) * t58;
t50 = qJD(1) * t55;
t16 = -t102 * t50 + t45 * t148 + t89 * t149;
t162 = t16 * t101;
t161 = t16 * t104;
t158 = qJD(5) * t74;
t107 = qJD(4) ^ 2;
t154 = t107 * t102;
t153 = t107 * t105;
t143 = qJD(2) * qJD(4);
t15 = qJD(4) * t190 - t105 * t50;
t59 = (t103 * t99 + t106 * t97) * t98;
t36 = (qJD(1) * t59 + t80) * qJD(2);
t141 = t101 * t36 + t104 * t15 + t35 * t146;
t138 = t76 * t148;
t134 = pkin(5) * t101 + t91;
t133 = t90 * t91 + t27;
t131 = t102 * t143;
t130 = t101 * t15 - t104 * t36;
t51 = -t104 * t142 + (t136 - t137) * qJD(2);
t128 = t51 * t105 + t76 * t149;
t125 = t90 * t135;
t124 = t76 * t135;
t10 = t104 * t27 + t171;
t7 = -qJ(6) * t74 + t10;
t122 = -t101 * t7 - t104 * t5;
t121 = t101 * t5 - t104 * t7;
t8 = pkin(5) * t52 + t16;
t43 = t100 * t102 + t105 * t59;
t22 = t101 * t58 + t104 * t43;
t21 = -t101 * t43 + t104 * t58;
t117 = t100 * t105 - t102 * t59;
t54 = qJD(2) * t59;
t49 = qJD(1) * t54;
t114 = qJD(2) * t53 - t107 * t91 - t49;
t44 = -qJD(2) * pkin(3) - t47;
t113 = qJD(4) * (qJD(2) * (-pkin(3) - t186) + t44 + t56);
t112 = -t27 * t147 + t141;
t111 = t118 * t101;
t109 = -t10 * qJD(5) - t130;
t87 = t183 * t104;
t86 = t183 * t101;
t73 = t74 ^ 2;
t62 = t104 * t71;
t37 = -t101 * t160 + t174;
t34 = -qJ(6) * t156 + t62 + (-pkin(5) - t169) * t105;
t20 = t117 * qJD(4) - t55 * t105;
t19 = t43 * qJD(4) - t55 * t102;
t18 = pkin(5) * t74 + qJD(6) + t26;
t4 = t21 * qJD(5) + t54 * t101 + t20 * t104;
t3 = -t22 * qJD(5) - t20 * t101 + t54 * t104;
t2 = -qJ(6) * t52 - qJD(6) * t74 + t112;
t1 = pkin(5) * t131 + t51 * qJ(6) - t76 * qJD(6) + t109;
t11 = [0, 0, -t103 * t163, -t106 * t163, -t47 * t54 - t48 * t55 + t49 * t58 - t50 * t59, 0, 0, 0, 0, 0, -t19 * qJD(4) + (-t105 * t54 + t58 * t149) * qJD(2), -t20 * qJD(4) + (t102 * t54 + t58 * t148) * qJD(2), 0, 0, 0, 0, 0, -t117 * t52 + t21 * t131 + t19 * t74 - t3 * t90, t117 * t51 - t22 * t131 + t19 * t76 + t4 * t90, t21 * t51 - t22 * t52 - t3 * t76 - t4 * t74, t1 * t21 - t117 * t8 + t18 * t19 + t2 * t22 + t3 * t5 + t4 * t7; 0, 0, 0, 0, t47 * t53 - t48 * t56 + (-t49 * t99 - t50 * t97) * pkin(2), 0.2e1 * t105 * t131, -0.2e1 * t173 * t143, t153, -t154, 0, t102 * t113 + t114 * t105, -t114 * t102 + t105 * t113, t76 * t137 + (-t51 * t104 - t76 * t147) * t102, -t124 + (-t138 + (t51 + t158) * t102) * t101 + t179, t128 - t189, t125 + t165 + (-t111 - t168) * qJD(4) (-t90 - t151) * t149 (t71 * t147 - t192) * t90 + ((t74 * t91 + t172) * qJD(4) + (t133 * t104 + t171) * qJD(5) + t130) * t105 + (t26 * t146 + t162 + t91 * t52 - t56 * t74 + ((-t91 * t157 + t62) * qJD(2) + t9) * qJD(4)) * t102, t191 * t90 + (-t133 * t147 + (t76 * t91 + t167) * qJD(4) + t141) * t105 + (-t26 * t147 + t161 - t91 * t51 - t56 * t76 + (-t174 * qJD(2) - t91 * t166 - t10) * qJD(4)) * t102, t34 * t51 - t37 * t52 - t182 * t76 - t181 * t74 + t122 * t148 + (t121 * qJD(5) - t1 * t104 - t101 * t2) * t102, t1 * t34 + t2 * t37 + t181 * t7 + t182 * t5 + t18 * t134 * t148 + (t8 * t134 + (pkin(5) * t146 - t56) * t18) * t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t154, -t153, 0, 0, 0, 0, 0, t125 - t165 + (-t111 + t168) * qJD(4), t128 + t189, t124 + (t138 + (-t51 + t158) * t102) * t101 + t179 (-qJD(4) * t121 - t8) * t105 + (qJD(4) * t18 + t122 * qJD(5) - t1 * t101 + t2 * t104) * t102; 0, 0, 0, 0, 0, -t102 * t108 * t105, t173 * t108, 0, 0, 0, qJD(4) * t30 - t44 * t152 - t16 (-qJD(2) * t44 + t50) * t105, -t51 * t101 - t76 * t166 (-t51 + t185) * t104 + (-t52 + t184) * t101, -t90 * t146 + (t90 * t155 + (-t76 + t150) * t102) * qJD(2), t90 * t147 + (-t90 * t157 + (t74 + t145) * t102) * qJD(2), t90 * t152, -pkin(4) * t52 - t161 + t129 * t90 - t30 * t74 + (pkin(9) * t166 + t172) * qJD(5) + (-t9 * t102 + (-pkin(9) * t149 - t105 * t26) * t101) * qJD(2), pkin(4) * t51 + t162 - t180 * t90 - t30 * t76 + (-pkin(9) * t170 + t167) * qJD(5) + (-t26 * t155 + (-pkin(9) * t145 + t10) * t102) * qJD(2), t86 * t51 + t87 * t52 - t176 * t76 - t177 * t74 + (t90 * t5 + t2) * t104 + (t90 * t7 - t1) * t101, -t2 * t87 + t1 * t86 + t8 * (-pkin(5) * t104 - pkin(4)) + t177 * t7 + t176 * t5 + (-pkin(5) * t170 - t30) * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76 * t74, -t73 + t188, -t51 - t185, -t184 - t52, t131, -t10 * t90 - t26 * t76 + t109, t26 * t74 - t9 * t90 - t112, pkin(5) * t51 - t187 * t74, t187 * t7 + (-t18 * t76 + t1) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t73 - t188, t5 * t76 + t7 * t74 + t8;];
tauc_reg  = t11;
