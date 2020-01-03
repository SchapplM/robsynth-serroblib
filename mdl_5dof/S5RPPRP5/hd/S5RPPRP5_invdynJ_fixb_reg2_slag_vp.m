% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPPRP5
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPRP5_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP5_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP5_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP5_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP5_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP5_invdynJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:53:48
% EndTime: 2019-12-31 17:53:50
% DurationCPUTime: 1.30s
% Computational Cost: add. (1276->255), mult. (2983->299), div. (0->0), fcn. (2133->6), ass. (0->135)
t111 = sin(qJ(4));
t109 = cos(pkin(7));
t164 = qJD(1) * t109;
t153 = t111 * t164;
t108 = sin(pkin(7));
t113 = cos(qJ(4));
t57 = t108 * t111 + t109 * t113;
t129 = -qJD(4) * t153 + qJDD(1) * t57;
t165 = qJD(1) * t108;
t152 = t113 * t165;
t51 = t152 - t153;
t186 = t129 + (-t51 + t152) * qJD(4);
t60 = -t109 * pkin(2) - t108 * qJ(3) - pkin(1);
t53 = t109 * pkin(3) - t60;
t49 = t57 * qJD(1);
t112 = sin(qJ(1));
t114 = cos(qJ(1));
t166 = g(1) * t114 + g(2) * t112;
t173 = t108 * t113;
t185 = t109 * t111 - t173;
t157 = t109 * qJDD(1);
t159 = qJD(1) * qJD(3);
t81 = t108 * t159;
t93 = t108 * qJDD(1);
t183 = -pkin(2) * t157 - qJ(3) * t93 - t81;
t182 = t49 ^ 2;
t181 = t51 ^ 2;
t180 = g(1) * t112;
t103 = g(2) * t114;
t179 = t51 * t49;
t178 = -pkin(6) + qJ(2);
t117 = qJ(2) ^ 2;
t106 = t109 ^ 2;
t160 = qJD(1) * qJD(2);
t155 = 0.2e1 * t160;
t75 = t106 * t155;
t94 = t106 * qJDD(1);
t177 = qJ(2) * t75 + t117 * t94;
t63 = t178 * t109;
t59 = qJD(1) * t63;
t176 = t111 * t59;
t175 = pkin(1) * t114 + qJ(2) * t112;
t107 = qJDD(1) * pkin(1);
t174 = qJDD(4) * pkin(4);
t172 = t108 * t114;
t170 = t109 * t112;
t169 = t109 * t114;
t67 = qJ(2) * t165 + qJD(3);
t55 = -pkin(6) * t165 + t67;
t24 = t111 * t55 + t113 * t59;
t168 = t24 * qJD(4);
t23 = t113 * t55 - t176;
t167 = qJD(5) - t23;
t163 = qJD(4) * t111;
t162 = qJD(4) * t113;
t161 = t108 * qJD(3);
t158 = qJDD(4) * qJ(5);
t91 = qJDD(2) - t107;
t54 = qJ(2) * t93 + t108 * t160 + qJDD(3);
t36 = -pkin(6) * t93 + t54;
t38 = (qJDD(1) * t178 + t160) * t109;
t156 = t111 * t36 + t113 * t38 + t162 * t55;
t154 = -t181 + t182;
t99 = t114 * qJ(2);
t151 = -t112 * pkin(1) + t99;
t150 = -t114 * pkin(6) + t99;
t149 = t108 * t157;
t148 = -t103 + t180;
t147 = t23 + t176;
t146 = t111 * t38 - t113 * t36 + t162 * t59 + t163 * t55;
t105 = t108 ^ 2;
t116 = qJD(1) ^ 2;
t61 = (-t105 - t106) * t116;
t145 = pkin(2) * t169 + qJ(3) * t172 + t175;
t144 = 0.2e1 * qJ(2) * t94 - t166 + t75;
t46 = -qJD(1) * pkin(1) - pkin(2) * t164 - qJ(3) * t165 + qJD(2);
t143 = pkin(3) * t169 + t145;
t142 = -pkin(4) * t57 - qJ(5) * t185;
t22 = qJD(4) * t152 + t129;
t48 = t108 * t162 - t109 * t163;
t139 = t22 * t57 + t48 * t49;
t33 = pkin(3) * t164 - t46;
t138 = t111 * t157 - t113 * t93;
t32 = t91 + t183;
t62 = t178 * t108;
t28 = t111 * t63 - t113 * t62;
t29 = t111 * t62 + t113 * t63;
t137 = qJD(4) * t48 + qJDD(4) * t57;
t136 = -t148 + t91;
t87 = pkin(3) * t157;
t27 = t87 - t32;
t135 = t107 - t91 - t103;
t134 = -qJDD(1) * t60 - t103 - t32;
t47 = t57 * qJD(4);
t133 = (qJ(2) * qJDD(1) + t160) * t105;
t21 = qJD(1) * t47 + t138;
t132 = -t111 * t22 + t113 * t21 - t162 * t49 + t163 * t51;
t14 = qJD(2) * t57 - qJD(4) * t28;
t15 = qJD(2) * t185 + t29 * qJD(4);
t131 = -t14 * t49 + t15 * t51 - t21 * t28 - t22 * t29 + t166;
t130 = t136 + t183;
t115 = qJD(4) ^ 2;
t127 = qJDD(4) * t111 + t113 * t115 + t165 * t51;
t42 = t57 * t112;
t44 = t57 * t114;
t126 = g(1) * t44 + g(2) * t42 - g(3) * t185 - t156;
t125 = t185 * t22 + t21 * t57 + t47 * t49 - t48 * t51;
t41 = t111 * t170 - t112 * t173;
t43 = t111 * t169 - t113 * t172;
t124 = g(1) * t41 - g(2) * t43 + qJD(4) * t14 + qJDD(4) * t29;
t123 = g(1) * t42 - g(2) * t44 - qJD(4) * t15 - qJDD(4) * t28;
t122 = g(1) * t43 + g(2) * t41 + g(3) * t57 - t146;
t121 = -t22 * pkin(4) - t21 * qJ(5) - t27;
t120 = (g(2) * pkin(6) + g(1) * t53) * t112;
t10 = pkin(4) * t49 - qJ(5) * t51 + t33;
t119 = t10 * t51 + qJDD(5) - t122;
t118 = (-t51 - t152) * qJD(4) - t129;
t100 = g(3) * t109;
t92 = t105 * qJDD(1);
t78 = g(1) * t170;
t26 = qJDD(4) * t113 - t111 * t115 - t165 * t49;
t25 = -qJD(4) * t47 - qJDD(4) * t185;
t20 = pkin(4) * t51 + qJ(5) * t49;
t18 = qJD(4) * qJ(5) + t24;
t17 = -qJD(4) * pkin(4) + t167;
t16 = -t142 + t53;
t12 = 0.2e1 * qJD(4) * t49 + t138;
t11 = pkin(4) * t48 + qJ(5) * t47 + qJD(5) * t185 + t161;
t7 = t181 + t182;
t6 = t185 * t21 - t47 * t51;
t4 = -t163 * t59 + t156;
t3 = qJDD(5) + t146 - t174;
t2 = t158 + (qJD(5) - t176) * qJD(4) + t156;
t1 = -qJD(5) * t51 - t121;
t5 = [0, 0, 0, 0, 0, qJDD(1), t148, t166, 0, 0, t92, 0.2e1 * t149, 0, t94, 0, 0, t109 * t135 + t78, (-t135 - t180) * t108, 0.2e1 * t133 + t144, -t91 * pkin(1) - g(1) * t151 - g(2) * t175 + (qJ(2) * t155 + qJDD(1) * t117) * t105 + t177, t92, 0, -0.2e1 * t149, 0, 0, t94, t78 + (t134 + t81) * t109, t54 * t108 + t133 + t144, t105 * t159 + (t134 + t180) * t108, t32 * t60 - g(1) * (-pkin(2) * t170 + t151) - g(2) * t145 + (qJ(2) * t54 + qJ(3) * t180 + qJD(2) * t67 - qJD(3) * t46) * t108 + t177, t6, t125, t25, t139, -t137, 0, t161 * t49 + t22 * t53 + t27 * t57 + t33 * t48 + t123, t161 * t51 - t185 * t27 - t21 * t53 - t33 * t47 - t124, -t146 * t185 + t23 * t47 - t24 * t48 - t4 * t57 + t131, -g(1) * t150 - g(2) * t143 + t24 * t14 + t146 * t28 - t23 * t15 + t161 * t33 + t27 * t53 + t4 * t29 + t120, t6, t25, -t125, 0, t137, t139, t1 * t57 + t10 * t48 + t11 * t49 + t16 * t22 + t123, -t17 * t47 - t18 * t48 - t185 * t3 - t2 * t57 + t131, t1 * t185 + t10 * t47 - t11 * t51 + t16 * t21 + t124, t2 * t29 + t18 * t14 + t1 * t16 + t10 * t11 + t3 * t28 + t17 * t15 - g(1) * (-t42 * pkin(4) - t41 * qJ(5) + t150) - g(2) * (pkin(4) * t44 + qJ(5) * t43 + t143) + t120; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t157, t93, t61, qJ(2) * t61 + t136, 0, 0, 0, 0, 0, 0, -t157, t61, -t93, -qJ(2) * t106 * t116 - t165 * t67 + t130, 0, 0, 0, 0, 0, 0, t118, t12, t7, -t23 * t51 - t24 * t49 + t130 - t87, 0, 0, 0, 0, 0, 0, t118, t7, -t12, -t18 * t49 + (qJD(5) + t17) * t51 + t121 - t148; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t108 * t116 * t109, t93, -t105 * t116, t100 + (qJD(1) * t46 - t166) * t108 + t54, 0, 0, 0, 0, 0, 0, t26, -t127, t132, t4 * t111 - t146 * t113 + t100 + (-t111 * t23 + t113 * t24) * qJD(4) + (-qJD(1) * t33 - t166) * t108, 0, 0, 0, 0, 0, 0, t26, t132, t127, t2 * t111 - t3 * t113 + t100 + (t111 * t17 + t113 * t18) * qJD(4) + (-qJD(1) * t10 - t166) * t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t179, -t154, -t138, -t179, -t186, qJDD(4), -t33 * t51 + t122 + t168, qJD(4) * t147 + t33 * t49 + t126, 0, 0, t179, -t138, t154, qJDD(4), t186, -t179, -t20 * t49 - t119 + t168 + 0.2e1 * t174, pkin(4) * t21 - t22 * qJ(5) + (t18 - t24) * t51 + (t17 - t167) * t49, 0.2e1 * t158 - t10 * t49 + t20 * t51 + (0.2e1 * qJD(5) - t147) * qJD(4) - t126, t2 * qJ(5) - t3 * pkin(4) - t10 * t20 - t17 * t24 - g(1) * (-pkin(4) * t43 + qJ(5) * t44) - g(2) * (-pkin(4) * t41 + qJ(5) * t42) - g(3) * t142 + t167 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(4) + t179, -t138, -t115 - t181, -qJD(4) * t18 + t119 - t174;];
tau_reg = t5;
