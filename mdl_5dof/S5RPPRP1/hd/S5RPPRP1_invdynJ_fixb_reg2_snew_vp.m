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
% Datum: 2020-01-03 11:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:25:49
% EndTime: 2020-01-03 11:25:56
% DurationCPUTime: 1.12s
% Computational Cost: add. (2257->181), mult. (4967->240), div. (0->0), fcn. (2902->8), ass. (0->128)
t112 = qJD(1) ^ 2;
t104 = sin(pkin(8));
t98 = t104 ^ 2;
t106 = cos(pkin(8));
t99 = t106 ^ 2;
t80 = (t98 + t99) * t112;
t167 = 2 * qJD(3);
t108 = sin(qJ(4));
t110 = cos(qJ(4));
t155 = t98 * t112;
t81 = t108 * t110 * t155;
t146 = t106 * qJDD(1);
t89 = -qJDD(4) + t146;
t61 = -t81 - t89;
t176 = pkin(4) * t61;
t151 = qJD(1) * t110;
t169 = (qJD(4) * t151 + qJDD(1) * t108) * t104;
t135 = t104 * t151;
t150 = t106 * qJD(1);
t90 = -qJD(4) + t150;
t71 = t90 * t135;
t51 = t169 - t71;
t175 = t104 * t51;
t105 = sin(pkin(7));
t107 = cos(pkin(7));
t109 = sin(qJ(1));
t111 = cos(qJ(1));
t124 = -t111 * g(2) - t109 * g(3);
t75 = qJDD(1) * pkin(1) + t124;
t123 = t109 * g(2) - t111 * g(3);
t76 = -t112 * pkin(1) - t123;
t158 = t105 * t75 + t107 * t76;
t50 = -t112 * pkin(2) + qJDD(1) * qJ(3) + t158;
t174 = qJD(1) * t167 + t50;
t153 = -g(1) + qJDD(2);
t38 = t104 * t153 + t174 * t106;
t159 = t106 * pkin(3);
t125 = -t104 * pkin(6) - t159;
t74 = t125 * qJD(1);
t30 = t74 * t150 + t38;
t132 = -t105 * t76 + t107 * t75;
t49 = -qJDD(1) * pkin(2) - t112 * qJ(3) + qJDD(3) - t132;
t41 = t125 * qJDD(1) + t49;
t12 = t108 * t41 + t110 * t30;
t134 = qJD(1) * qJD(5) * t104;
t59 = -t90 * pkin(4) - qJ(5) * t135;
t120 = -qJ(5) * t169 - 0.2e1 * t108 * t134 + t90 * t59 + t12;
t101 = t108 ^ 2;
t142 = t101 * t155;
t102 = t110 ^ 2;
t141 = t102 * t155;
t88 = t90 ^ 2;
t55 = -t141 - t88;
t173 = -t120 + (t55 + t142) * pkin(4);
t133 = pkin(1) * t105 + qJ(3);
t137 = pkin(1) * t107 + pkin(2);
t172 = -t137 * qJDD(1) + t133 * t80 + t49;
t95 = t106 * t153;
t130 = pkin(4) * t169 - qJ(5) * t142 + qJDD(5) - t95;
t131 = -t110 * t59 - t74;
t14 = (t50 + (t167 - t131) * qJD(1)) * t104 + t130;
t11 = t108 * t30 - t110 * t41;
t152 = qJD(1) * t108;
t136 = t104 * t152;
t129 = t90 * t136;
t147 = t104 * qJDD(1);
t85 = qJD(4) * t136;
t63 = t110 * t147 - t85;
t118 = -t11 + (t129 - t63) * qJ(5) + t176;
t127 = t110 * t134;
t84 = -0.2e1 * t127;
t7 = t118 + t84;
t166 = pkin(4) * t7;
t60 = -t81 + t89;
t157 = t108 * t60;
t43 = t110 * t55 + t157;
t165 = pkin(3) * t43;
t156 = t110 * t61;
t64 = -t88 - t142;
t47 = t108 * t64 + t156;
t164 = pkin(3) * t47;
t138 = t90 * t152;
t148 = qJDD(1) * t110;
t54 = -t85 + (-t138 + t148) * t104;
t163 = pkin(4) * t54;
t52 = t71 + t169;
t32 = -t108 * t52 - t110 * t54;
t162 = pkin(6) * t32;
t161 = pkin(6) * t43;
t160 = pkin(6) * t47;
t154 = qJ(5) * t110;
t33 = t108 * t54 - t110 * t52;
t65 = (t101 + t102) * t155;
t19 = -t104 * t65 + t106 * t33;
t145 = qJ(3) * t19 + pkin(1) * (t105 * t19 - t107 * t32) - pkin(2) * t32;
t44 = -t108 * t55 + t110 * t60;
t53 = -t85 + (t138 + t148) * t104;
t22 = t104 * t53 + t106 * t44;
t144 = pkin(1) * (t105 * t22 - t107 * t43) + qJ(3) * t22 - pkin(2) * t43;
t48 = -t108 * t61 + t110 * t64;
t27 = t106 * t48 + t175;
t143 = pkin(1) * (t105 * t27 - t107 * t47) + qJ(3) * t27 - pkin(2) * t47;
t140 = t106 * t155;
t37 = t174 * t104 - t95;
t15 = t104 * t37 + t106 * t38;
t5 = t108 * t12 - t110 * t11;
t116 = t118 + t176;
t97 = t99 * qJDD(1);
t96 = t98 * qJDD(1);
t83 = 0.2e1 * t127;
t79 = t106 * t89;
t77 = t97 + t96;
t66 = (-t101 + t102) * t155;
t36 = (t104 * (t63 + t129) - t108 * t140) * t110;
t35 = (t110 * t140 + t175) * t108;
t29 = -t95 + (t50 + (t167 + t74) * qJD(1)) * t104;
t26 = t104 * t48 - t106 * t51;
t25 = t104 * (t110 * (-t88 + t142) + t157) + t106 * t52;
t24 = t104 * (t156 - t108 * (t88 - t141)) - t106 * t54;
t21 = t104 * t44 - t106 * t53;
t18 = t104 * (-t108 * t53 - t110 * t51) - t106 * t66;
t17 = t104 * t33 + t106 * t65;
t8 = -pkin(4) * t142 + t120;
t6 = t108 * t11 + t110 * t12;
t3 = -t108 * t7 + t110 * t8;
t2 = t108 * t8 + t110 * t7;
t1 = t104 * t14 + t106 * t3;
t4 = [0, 0, 0, 0, 0, qJDD(1), t124, t123, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (t107 * qJDD(1) - t105 * t112) + t132, pkin(1) * (-t105 * qJDD(1) - t107 * t112) - t158, 0, pkin(1) * (t105 * t158 + t107 * t132), t96, 0.2e1 * t104 * t146, 0, t97, 0, 0, -t172 * t106, t172 * t104, pkin(2) * t80 + qJ(3) * t77 + pkin(1) * (t105 * t77 + t107 * t80) + t15, -pkin(2) * t49 + qJ(3) * t15 + pkin(1) * (t105 * t15 - t107 * t49), t36, t18, t24, t35, t25, t79, t104 * (t108 * t29 - t160) + t106 * (t11 - t164) + t143, t104 * (t110 * t29 - t161) + t106 * (t12 - t165) + t144, t104 * (-t5 - t162) - t32 * t159 + t145, t133 * (t104 * t29 + t106 * t6) + (t125 - t137) * t5, t36, t18, t24, t35, t25, t79, t104 * (-t61 * t154 - t160 + (pkin(4) * t51 - qJ(5) * t64 + t130 + (-t131 * qJD(1) + t174) * t104) * t108) + t106 * (-t116 + t83 - t164) + t143, t104 * (-t108 * (-pkin(4) * t53 + qJ(5) * t60) - t161 + (-qJ(5) * t55 + t14) * t110) + t106 * (-t165 - t173) + t144, t104 * (t110 * (qJ(5) * t54 - t118 + t83) - t108 * (-qJ(5) * t52 + (t65 - t142) * pkin(4) + t120) - t162) + t106 * (-pkin(3) * t32 + t163) + t145, t104 * (-t7 * t154 - t108 * (-pkin(4) * t14 + qJ(5) * t8) - pkin(6) * t2) + t106 * (-pkin(3) * t2 - t166) - pkin(2) * t2 + qJ(3) * t1 + pkin(1) * (t105 * t1 - t107 * t2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t153, 0, 0, 0, 0, 0, 0, 0, 0, 0, t104 * t38 - t106 * t37, 0, 0, 0, 0, 0, 0, t26, t21, t17, t104 * t6 - t106 * t29, 0, 0, 0, 0, 0, 0, t26, t21, t17, t104 * t3 - t106 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t146, t147, -t80, t49, 0, 0, 0, 0, 0, 0, t47, t43, t32, t5, 0, 0, 0, 0, 0, 0, t47, t43, t32, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, t66, t54, -t81, -t52, -t89, -t11, -t12, 0, 0, t81, t66, t54, -t81, -t52, -t89, t116 + t84, t173, -t163, t166; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, t53, -t65, t14;];
tauJ_reg = t4;
