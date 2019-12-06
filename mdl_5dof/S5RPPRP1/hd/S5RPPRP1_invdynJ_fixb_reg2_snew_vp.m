% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPPRP1
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPPRP1_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP1_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP1_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP1_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:36:25
% EndTime: 2019-12-05 17:36:30
% DurationCPUTime: 1.16s
% Computational Cost: add. (2257->181), mult. (4967->241), div. (0->0), fcn. (2902->8), ass. (0->128)
t106 = sin(pkin(8));
t100 = t106 ^ 2;
t108 = cos(pkin(8));
t101 = t108 ^ 2;
t114 = qJD(1) ^ 2;
t80 = (t100 + t101) * t114;
t169 = 2 * qJD(3);
t110 = sin(qJ(4));
t112 = cos(qJ(4));
t155 = t100 * t114;
t141 = t112 * t155;
t81 = t110 * t141;
t147 = t108 * qJDD(1);
t89 = -qJDD(4) + t147;
t61 = -t81 - t89;
t178 = pkin(4) * t61;
t152 = qJD(1) * t112;
t172 = (qJD(4) * t152 + qJDD(1) * t110) * t106;
t136 = t106 * t152;
t151 = t108 * qJD(1);
t90 = -qJD(4) + t151;
t71 = t90 * t136;
t51 = t172 - t71;
t177 = t106 * t51;
t107 = sin(pkin(7));
t109 = cos(pkin(7));
t111 = sin(qJ(1));
t113 = cos(qJ(1));
t159 = g(2) * t113 + g(3) * t111;
t75 = qJDD(1) * pkin(1) + t159;
t125 = -g(2) * t111 + g(3) * t113;
t76 = -pkin(1) * t114 - t125;
t160 = t107 * t75 + t109 * t76;
t50 = -pkin(2) * t114 + qJDD(1) * qJ(3) + t160;
t176 = qJD(1) * t169 + t50;
t154 = -g(1) + qJDD(2);
t38 = t106 * t154 + t108 * t176;
t161 = t108 * pkin(3);
t126 = -pkin(6) * t106 - t161;
t74 = t126 * qJD(1);
t30 = t151 * t74 + t38;
t133 = -t107 * t76 + t109 * t75;
t49 = -qJDD(1) * pkin(2) - qJ(3) * t114 + qJDD(3) - t133;
t41 = qJDD(1) * t126 + t49;
t12 = t110 * t41 + t112 * t30;
t135 = qJD(1) * qJD(5) * t106;
t59 = -pkin(4) * t90 - qJ(5) * t136;
t122 = -qJ(5) * t172 - 0.2e1 * t110 * t135 + t90 * t59 + t12;
t103 = t110 ^ 2;
t140 = t103 * t155;
t104 = t112 ^ 2;
t139 = t104 * t155;
t88 = t90 ^ 2;
t55 = -t139 - t88;
t175 = -t122 + (t55 + t140) * pkin(4);
t134 = pkin(1) * t107 + qJ(3);
t138 = pkin(1) * t109 + pkin(2);
t174 = -t138 * qJDD(1) + t134 * t80 + t49;
t95 = t108 * t154;
t131 = pkin(4) * t172 - qJ(5) * t140 + qJDD(5) - t95;
t132 = -t112 * t59 - t74;
t14 = (t50 + (t169 - t132) * qJD(1)) * t106 + t131;
t11 = t110 * t30 - t112 * t41;
t153 = qJD(1) * t110;
t137 = t106 * t153;
t130 = t90 * t137;
t148 = t106 * qJDD(1);
t85 = qJD(4) * t137;
t63 = t112 * t148 - t85;
t120 = -t11 + (t130 - t63) * qJ(5) + t178;
t128 = t112 * t135;
t84 = -0.2e1 * t128;
t7 = t120 + t84;
t168 = pkin(4) * t7;
t60 = -t81 + t89;
t158 = t110 * t60;
t43 = t112 * t55 + t158;
t167 = pkin(3) * t43;
t157 = t112 * t61;
t64 = -t88 - t140;
t47 = t110 * t64 + t157;
t166 = pkin(3) * t47;
t142 = t90 * t153;
t149 = qJDD(1) * t112;
t54 = -t85 + (-t142 + t149) * t106;
t165 = pkin(4) * t54;
t52 = t71 + t172;
t32 = -t110 * t52 - t112 * t54;
t164 = pkin(6) * t32;
t163 = pkin(6) * t43;
t162 = pkin(6) * t47;
t156 = qJ(5) * t112;
t33 = t110 * t54 - t112 * t52;
t65 = (t103 + t104) * t155;
t19 = -t106 * t65 + t108 * t33;
t146 = qJ(3) * t19 + pkin(1) * (t107 * t19 - t109 * t32) - pkin(2) * t32;
t44 = -t110 * t55 + t112 * t60;
t53 = -t85 + (t142 + t149) * t106;
t22 = t106 * t53 + t108 * t44;
t145 = pkin(1) * (t107 * t22 - t109 * t43) + qJ(3) * t22 - pkin(2) * t43;
t48 = -t110 * t61 + t112 * t64;
t27 = t108 * t48 + t177;
t144 = pkin(1) * (t107 * t27 - t109 * t47) + qJ(3) * t27 - pkin(2) * t47;
t37 = t106 * t176 - t95;
t15 = t106 * t37 + t108 * t38;
t5 = -t11 * t112 + t110 * t12;
t118 = t120 + t178;
t97 = t101 * qJDD(1);
t96 = t100 * qJDD(1);
t83 = 0.2e1 * t128;
t79 = t108 * t89;
t77 = t97 + t96;
t66 = (-t103 + t104) * t155;
t36 = (t106 * (t63 + t130) - t108 * t110 * t155) * t112;
t35 = (t108 * t141 + t177) * t110;
t29 = -t95 + (t50 + (t169 + t74) * qJD(1)) * t106;
t26 = t106 * t48 - t108 * t51;
t25 = t106 * (t112 * (-t88 + t140) + t158) + t108 * t52;
t24 = t106 * (t157 - t110 * (t88 - t139)) - t108 * t54;
t21 = t106 * t44 - t108 * t53;
t18 = t106 * (-t110 * t53 - t112 * t51) - t108 * t66;
t17 = t106 * t33 + t108 * t65;
t8 = -pkin(4) * t140 + t122;
t6 = t11 * t110 + t112 * t12;
t3 = -t110 * t7 + t112 * t8;
t2 = t110 * t8 + t112 * t7;
t1 = t106 * t14 + t108 * t3;
t4 = [0, 0, 0, 0, 0, qJDD(1), t159, t125, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (qJDD(1) * t109 - t107 * t114) + t133, pkin(1) * (-qJDD(1) * t107 - t109 * t114) - t160, 0, pkin(1) * (t107 * t160 + t109 * t133), t96, 0.2e1 * t106 * t147, 0, t97, 0, 0, -t174 * t108, t174 * t106, pkin(2) * t80 + qJ(3) * t77 + pkin(1) * (t107 * t77 + t109 * t80) + t15, -pkin(2) * t49 + qJ(3) * t15 + pkin(1) * (t107 * t15 - t109 * t49), t36, t18, t24, t35, t25, t79, t106 * (t110 * t29 - t162) + t108 * (t11 - t166) + t144, t106 * (t112 * t29 - t163) + t108 * (t12 - t167) + t145, t106 * (-t5 - t164) - t32 * t161 + t146, t134 * (t106 * t29 + t108 * t6) + (t126 - t138) * t5, t36, t18, t24, t35, t25, t79, t106 * (-t61 * t156 - t162 + (pkin(4) * t51 - qJ(5) * t64 + t131 + (-qJD(1) * t132 + t176) * t106) * t110) + t108 * (-t118 + t83 - t166) + t144, t106 * (-t110 * (-pkin(4) * t53 + qJ(5) * t60) - t163 + (-qJ(5) * t55 + t14) * t112) + t108 * (-t167 - t175) + t145, t106 * (t112 * (qJ(5) * t54 - t120 + t83) - t110 * (-qJ(5) * t52 + (t65 - t140) * pkin(4) + t122) - t164) + t108 * (-pkin(3) * t32 + t165) + t146, t106 * (-t7 * t156 - t110 * (-pkin(4) * t14 + qJ(5) * t8) - pkin(6) * t2) + t108 * (-pkin(3) * t2 - t168) - pkin(2) * t2 + qJ(3) * t1 + pkin(1) * (t1 * t107 - t109 * t2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t154, 0, 0, 0, 0, 0, 0, 0, 0, 0, t106 * t38 - t108 * t37, 0, 0, 0, 0, 0, 0, t26, t21, t17, t106 * t6 - t108 * t29, 0, 0, 0, 0, 0, 0, t26, t21, t17, t106 * t3 - t108 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t147, t148, -t80, t49, 0, 0, 0, 0, 0, 0, t47, t43, t32, t5, 0, 0, 0, 0, 0, 0, t47, t43, t32, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, t66, t54, -t81, -t52, -t89, -t11, -t12, 0, 0, t81, t66, t54, -t81, -t52, -t89, t118 + t84, t175, -t165, t168; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, t53, -t65, t14;];
tauJ_reg = t4;
