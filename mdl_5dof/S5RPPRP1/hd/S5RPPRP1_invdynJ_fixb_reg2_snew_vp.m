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
% Datum: 2022-01-23 09:13
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-23 09:12:46
% EndTime: 2022-01-23 09:12:48
% DurationCPUTime: 1.11s
% Computational Cost: add. (2257->181), mult. (4967->240), div. (0->0), fcn. (2902->8), ass. (0->128)
t107 = cos(pkin(8));
t100 = t107 ^ 2;
t113 = qJD(1) ^ 2;
t105 = sin(pkin(8));
t99 = t105 ^ 2;
t80 = (t100 + t99) * t113;
t168 = 2 * qJD(3);
t109 = sin(qJ(4));
t111 = cos(qJ(4));
t156 = t99 * t113;
t81 = t109 * t111 * t156;
t147 = t107 * qJDD(1);
t89 = -qJDD(4) + t147;
t61 = -t81 - t89;
t177 = pkin(4) * t61;
t152 = qJD(1) * t111;
t170 = (qJD(4) * t152 + qJDD(1) * t109) * t105;
t136 = t105 * t152;
t151 = t107 * qJD(1);
t90 = -qJD(4) + t151;
t71 = t90 * t136;
t51 = t170 - t71;
t176 = t105 * t51;
t106 = sin(pkin(7));
t108 = cos(pkin(7));
t110 = sin(qJ(1));
t112 = cos(qJ(1));
t135 = g(1) * t110 - g(2) * t112;
t75 = qJDD(1) * pkin(1) + t135;
t124 = g(1) * t112 + g(2) * t110;
t76 = -pkin(1) * t113 - t124;
t159 = t106 * t75 + t108 * t76;
t50 = -pkin(2) * t113 + qJDD(1) * qJ(3) + t159;
t175 = qJD(1) * t168 + t50;
t154 = -g(3) + qJDD(2);
t38 = t105 * t154 + t107 * t175;
t160 = t107 * pkin(3);
t125 = -pkin(6) * t105 - t160;
t74 = t125 * qJD(1);
t30 = t151 * t74 + t38;
t132 = -t106 * t76 + t108 * t75;
t49 = -qJDD(1) * pkin(2) - qJ(3) * t113 + qJDD(3) - t132;
t41 = qJDD(1) * t125 + t49;
t12 = t109 * t41 + t111 * t30;
t134 = qJD(1) * qJD(5) * t105;
t59 = -pkin(4) * t90 - qJ(5) * t136;
t121 = -qJ(5) * t170 - 0.2e1 * t109 * t134 + t90 * t59 + t12;
t102 = t109 ^ 2;
t143 = t102 * t156;
t103 = t111 ^ 2;
t142 = t103 * t156;
t88 = t90 ^ 2;
t55 = -t142 - t88;
t174 = -t121 + (t55 + t143) * pkin(4);
t133 = pkin(1) * t106 + qJ(3);
t138 = pkin(1) * t108 + pkin(2);
t173 = -t138 * qJDD(1) + t133 * t80 + t49;
t95 = t107 * t154;
t130 = pkin(4) * t170 - qJ(5) * t143 + qJDD(5) - t95;
t131 = -t111 * t59 - t74;
t14 = (t50 + (t168 - t131) * qJD(1)) * t105 + t130;
t11 = t109 * t30 - t111 * t41;
t153 = qJD(1) * t109;
t137 = t105 * t153;
t129 = t90 * t137;
t148 = t105 * qJDD(1);
t85 = qJD(4) * t137;
t63 = t111 * t148 - t85;
t119 = -t11 + (t129 - t63) * qJ(5) + t177;
t127 = t111 * t134;
t84 = -0.2e1 * t127;
t7 = t119 + t84;
t167 = pkin(4) * t7;
t60 = -t81 + t89;
t158 = t109 * t60;
t43 = t111 * t55 + t158;
t166 = pkin(3) * t43;
t157 = t111 * t61;
t64 = -t88 - t143;
t47 = t109 * t64 + t157;
t165 = pkin(3) * t47;
t139 = t90 * t153;
t149 = qJDD(1) * t111;
t54 = -t85 + (-t139 + t149) * t105;
t164 = pkin(4) * t54;
t52 = t71 + t170;
t32 = -t109 * t52 - t111 * t54;
t163 = pkin(6) * t32;
t162 = pkin(6) * t43;
t161 = pkin(6) * t47;
t155 = qJ(5) * t111;
t33 = t109 * t54 - t111 * t52;
t65 = (t102 + t103) * t156;
t19 = -t105 * t65 + t107 * t33;
t146 = qJ(3) * t19 + pkin(1) * (t106 * t19 - t108 * t32) - pkin(2) * t32;
t44 = -t109 * t55 + t111 * t60;
t53 = -t85 + (t139 + t149) * t105;
t22 = t105 * t53 + t107 * t44;
t145 = pkin(1) * (t106 * t22 - t108 * t43) + qJ(3) * t22 - pkin(2) * t43;
t48 = -t109 * t61 + t111 * t64;
t27 = t107 * t48 + t176;
t144 = pkin(1) * (t106 * t27 - t108 * t47) + qJ(3) * t27 - pkin(2) * t47;
t141 = t107 * t156;
t37 = t105 * t175 - t95;
t15 = t105 * t37 + t107 * t38;
t5 = t109 * t12 - t11 * t111;
t117 = t119 + t177;
t97 = t100 * qJDD(1);
t96 = t99 * qJDD(1);
t83 = 0.2e1 * t127;
t79 = t107 * t89;
t77 = t97 + t96;
t66 = (-t102 + t103) * t156;
t36 = (t105 * (t63 + t129) - t109 * t141) * t111;
t35 = (t111 * t141 + t176) * t109;
t29 = -t95 + (t50 + (t168 + t74) * qJD(1)) * t105;
t26 = t105 * t48 - t107 * t51;
t25 = t105 * (t111 * (-t88 + t143) + t158) + t107 * t52;
t24 = t105 * (t157 - t109 * (t88 - t142)) - t107 * t54;
t21 = t105 * t44 - t107 * t53;
t18 = t105 * (-t109 * t53 - t111 * t51) - t107 * t66;
t17 = t105 * t33 + t107 * t65;
t8 = -pkin(4) * t143 + t121;
t6 = t109 * t11 + t111 * t12;
t3 = -t109 * t7 + t111 * t8;
t2 = t109 * t8 + t111 * t7;
t1 = t105 * t14 + t107 * t3;
t4 = [0, 0, 0, 0, 0, qJDD(1), t135, t124, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (qJDD(1) * t108 - t106 * t113) + t132, pkin(1) * (-qJDD(1) * t106 - t108 * t113) - t159, 0, pkin(1) * (t106 * t159 + t108 * t132), t96, 0.2e1 * t105 * t147, 0, t97, 0, 0, -t173 * t107, t173 * t105, pkin(2) * t80 + qJ(3) * t77 + pkin(1) * (t106 * t77 + t108 * t80) + t15, -pkin(2) * t49 + qJ(3) * t15 + pkin(1) * (t106 * t15 - t108 * t49), t36, t18, t24, t35, t25, t79, t105 * (t109 * t29 - t161) + t107 * (t11 - t165) + t144, t105 * (t111 * t29 - t162) + t107 * (t12 - t166) + t145, t105 * (-t5 - t163) - t32 * t160 + t146, t133 * (t105 * t29 + t107 * t6) + (t125 - t138) * t5, t36, t18, t24, t35, t25, t79, t105 * (-t61 * t155 - t161 + (pkin(4) * t51 - qJ(5) * t64 + t130 + (-qJD(1) * t131 + t175) * t105) * t109) + t107 * (-t117 + t83 - t165) + t144, t105 * (-t109 * (-pkin(4) * t53 + qJ(5) * t60) - t162 + (-qJ(5) * t55 + t14) * t111) + t107 * (-t166 - t174) + t145, t105 * (t111 * (qJ(5) * t54 - t119 + t83) - t109 * (-qJ(5) * t52 + (t65 - t143) * pkin(4) + t121) - t163) + t107 * (-pkin(3) * t32 + t164) + t146, t105 * (-t7 * t155 - t109 * (-pkin(4) * t14 + qJ(5) * t8) - pkin(6) * t2) + t107 * (-pkin(3) * t2 - t167) - pkin(2) * t2 + qJ(3) * t1 + pkin(1) * (t1 * t106 - t108 * t2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t154, 0, 0, 0, 0, 0, 0, 0, 0, 0, t105 * t38 - t107 * t37, 0, 0, 0, 0, 0, 0, t26, t21, t17, t105 * t6 - t107 * t29, 0, 0, 0, 0, 0, 0, t26, t21, t17, t105 * t3 - t107 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t147, t148, -t80, t49, 0, 0, 0, 0, 0, 0, t47, t43, t32, t5, 0, 0, 0, 0, 0, 0, t47, t43, t32, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, t66, t54, -t81, -t52, -t89, -t11, -t12, 0, 0, t81, t66, t54, -t81, -t52, -t89, t117 + t84, t174, -t164, t167; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, t53, -t65, t14;];
tauJ_reg = t4;
