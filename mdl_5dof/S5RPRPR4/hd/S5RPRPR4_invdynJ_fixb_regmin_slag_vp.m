% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRPR4
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% tau_reg [5x22]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:23
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRPR4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR4_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:23:10
% EndTime: 2022-01-23 09:23:13
% DurationCPUTime: 1.28s
% Computational Cost: add. (1636->229), mult. (3529->314), div. (0->0), fcn. (2526->16), ass. (0->140)
t125 = sin(qJ(5));
t128 = cos(qJ(5));
t120 = sin(pkin(9));
t122 = cos(pkin(9));
t129 = cos(qJ(3));
t169 = t122 * t129;
t159 = qJD(1) * t169;
t126 = sin(qJ(3));
t165 = qJD(1) * t126;
t74 = t120 * t165 - t159;
t83 = t120 * t129 + t122 * t126;
t77 = t83 * qJD(1);
t144 = t125 * t74 - t128 * t77;
t161 = t129 * qJDD(1);
t162 = t126 * qJDD(1);
t149 = t120 * t162 - t122 * t161;
t76 = t83 * qJD(3);
t42 = qJD(1) * t76 + t149;
t163 = qJD(1) * qJD(3);
t158 = t126 * t163;
t138 = qJDD(1) * t83 - t120 * t158;
t157 = t129 * t163;
t43 = t122 * t157 + t138;
t135 = qJD(5) * t144 - t125 * t43 - t128 * t42;
t115 = qJD(3) + qJD(5);
t171 = t144 * t115;
t189 = t135 - t171;
t188 = 0.2e1 * qJD(3);
t66 = t128 * t74;
t35 = -t125 * t77 - t66;
t170 = t35 * t115;
t164 = qJD(5) * t125;
t7 = -qJD(5) * t66 - t125 * t42 + t128 * t43 - t164 * t77;
t187 = t7 - t170;
t186 = t144 * t35;
t121 = sin(pkin(8));
t103 = pkin(1) * t121 + pkin(6);
t168 = qJ(4) + t103;
t117 = qJ(1) + pkin(8);
t108 = sin(t117);
t110 = cos(t117);
t152 = g(1) * t110 + g(2) * t108;
t185 = t144 ^ 2 - t35 ^ 2;
t116 = qJ(3) + pkin(9);
t112 = qJ(5) + t116;
t101 = sin(t112);
t102 = cos(t112);
t181 = t74 * pkin(7);
t154 = t168 * qJD(1);
t60 = t126 * qJD(2) + t129 * t154;
t172 = t122 * t60;
t173 = qJD(3) * pkin(3);
t59 = qJD(2) * t129 - t126 * t154;
t53 = t59 + t173;
t24 = t120 * t53 + t172;
t11 = t24 - t181;
t106 = pkin(3) * t129 + pkin(2);
t123 = cos(pkin(8));
t174 = t123 * pkin(1);
t88 = -t106 - t174;
t72 = qJD(1) * t88 + qJD(4);
t44 = t74 * pkin(4) + t72;
t184 = g(3) * t101 + t102 * t152 + t11 * t164 - t44 * t35;
t111 = t129 * qJDD(2);
t91 = t103 * qJDD(1);
t136 = qJ(4) * qJDD(1) + qJD(1) * qJD(4) + qJD(2) * qJD(3) + t91;
t143 = t154 * qJD(3);
t21 = qJDD(3) * pkin(3) - t126 * t136 - t129 * t143 + t111;
t25 = (qJDD(2) - t143) * t126 + t136 * t129;
t4 = -t120 * t25 + t122 * t21;
t2 = qJDD(3) * pkin(4) - pkin(7) * t43 + t4;
t5 = t120 * t21 + t122 * t25;
t3 = -pkin(7) * t42 + t5;
t183 = -g(3) * t102 + t101 * t152 - t125 * t3 + t128 * t2 + t44 * t144;
t182 = qJD(5) - t115;
t180 = t77 * pkin(7);
t179 = pkin(3) * t120;
t178 = pkin(3) * t126;
t175 = g(3) * t129;
t49 = t120 * t60;
t27 = t122 * t59 - t49;
t153 = qJD(3) * t168;
t63 = t129 * qJD(4) - t126 * t153;
t64 = -t126 * qJD(4) - t129 * t153;
t29 = t120 * t64 + t122 * t63;
t80 = t168 * t126;
t81 = t168 * t129;
t41 = -t120 * t80 + t122 * t81;
t167 = qJDD(2) - g(3);
t118 = t126 ^ 2;
t166 = -t129 ^ 2 + t118;
t105 = -pkin(2) - t174;
t94 = qJD(1) * t105;
t160 = t126 * t173;
t23 = t122 * t53 - t49;
t26 = -t120 * t59 - t172;
t28 = -t120 * t63 + t122 * t64;
t40 = -t120 * t81 - t122 * t80;
t151 = g(1) * t108 - g(2) * t110;
t127 = sin(qJ(1));
t130 = cos(qJ(1));
t150 = g(1) * t127 - g(2) * t130;
t10 = qJD(3) * pkin(4) - t180 + t23;
t148 = -t125 * t10 - t128 * t11;
t114 = qJDD(3) + qJDD(5);
t82 = t120 * t126 - t169;
t45 = t125 * t83 + t128 * t82;
t79 = t82 * qJD(3);
t12 = -qJD(5) * t45 - t125 * t76 - t128 * t79;
t46 = -t125 * t82 + t128 * t83;
t147 = t114 * t46 + t115 * t12;
t30 = -pkin(7) * t83 + t40;
t31 = -pkin(7) * t82 + t41;
t146 = -t125 * t31 + t128 * t30;
t145 = t125 * t30 + t128 * t31;
t104 = pkin(3) * t122 + pkin(4);
t142 = t104 * t125 + t128 * t179;
t141 = t104 * t128 - t125 * t179;
t139 = -qJD(1) * t94 + t152 - t91;
t137 = -qJDD(3) * t103 + t188 * t94;
t58 = pkin(3) * t158 + qJDD(1) * t88 + qJDD(4);
t131 = qJD(3) ^ 2;
t134 = -0.2e1 * qJDD(1) * t105 - t103 * t131 + t151;
t132 = qJD(1) ^ 2;
t124 = -qJ(4) - pkin(6);
t109 = cos(t116);
t107 = sin(t116);
t90 = qJDD(3) * t129 - t126 * t131;
t89 = qJDD(3) * t126 + t129 * t131;
t62 = pkin(4) * t76 + t160;
t61 = pkin(3) * t165 + pkin(4) * t77;
t57 = pkin(4) * t82 + t88;
t22 = t42 * pkin(4) + t58;
t17 = -pkin(7) * t76 + t29;
t16 = pkin(7) * t79 + t28;
t15 = t27 - t180;
t14 = t26 + t181;
t13 = qJD(5) * t46 - t125 * t79 + t128 * t76;
t6 = -t114 * t45 - t115 * t13;
t1 = [qJDD(1), t150, g(1) * t130 + g(2) * t127, (t150 + (t121 ^ 2 + t123 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), qJDD(1) * t118 + 0.2e1 * t126 * t157, 0.2e1 * t126 * t161 - 0.2e1 * t163 * t166, t89, t90, 0, t126 * t137 + t129 * t134, -t126 * t134 + t129 * t137, t40 * qJDD(3) + t88 * t42 + t58 * t82 + t72 * t76 + t151 * t109 + (t178 * t74 + t28) * qJD(3), -t41 * qJDD(3) + t88 * t43 + t58 * t83 - t72 * t79 - t151 * t107 + (t178 * t77 - t29) * qJD(3), t23 * t79 - t24 * t76 - t28 * t77 - t29 * t74 - t4 * t83 - t40 * t43 - t41 * t42 - t5 * t82 - t152, t5 * t41 + t24 * t29 + t4 * t40 + t23 * t28 + t58 * t88 + t72 * t160 - g(1) * (-pkin(1) * t127 - t106 * t108 - t110 * t124) - g(2) * (pkin(1) * t130 + t106 * t110 - t108 * t124), -t12 * t144 + t46 * t7, t12 * t35 + t13 * t144 + t135 * t46 - t45 * t7, t147, t6, 0, -t62 * t35 - t57 * t135 + t22 * t45 + t44 * t13 + (-qJD(5) * t145 - t125 * t17 + t128 * t16) * t115 + t146 * t114 + t151 * t102, -t62 * t144 + t57 * t7 + t22 * t46 + t44 * t12 - (qJD(5) * t146 + t125 * t16 + t128 * t17) * t115 - t145 * t114 - t151 * t101; 0, 0, 0, t167, 0, 0, 0, 0, 0, t90, -t89, -qJD(3) * t76 - qJDD(3) * t82, qJD(3) * t79 - qJDD(3) * t83, -t42 * t83 + t43 * t82 + t74 * t79 + t76 * t77, -t23 * t76 - t24 * t79 - t4 * t82 + t5 * t83 - g(3), 0, 0, 0, 0, 0, t6, -t147; 0, 0, 0, 0, -t126 * t132 * t129, t166 * t132, t162, t161, qJDD(3), t126 * t139 + t111 - t175, -t126 * t167 + t129 * t139, -g(3) * t109 - t26 * qJD(3) - t72 * t77 + t152 * t107 + (qJDD(3) * t122 - t165 * t74) * pkin(3) + t4, g(3) * t107 + t27 * qJD(3) + t72 * t74 + t152 * t109 + (-qJDD(3) * t120 - t165 * t77) * pkin(3) - t5, (t24 + t26) * t77 + (-t23 + t27) * t74 + (-t120 * t42 - t122 * t43) * pkin(3), -t23 * t26 - t24 * t27 + (-t175 + t120 * t5 + t122 * t4 + (-qJD(1) * t72 + t152) * t126) * pkin(3), t186, t185, t187, t189, t114, t141 * t114 + t61 * t35 - (-t125 * t15 + t128 * t14) * t115 + (-t115 * t142 + t148) * qJD(5) + t183, -t142 * t114 - t128 * t3 - t125 * t2 + t61 * t144 + (t125 * t14 + t128 * t15) * t115 + (-t128 * t10 - t115 * t141) * qJD(5) + t184; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t188 * t77 + t149, (-t74 + t159) * qJD(3) + t138, -t74 ^ 2 - t77 ^ 2, t23 * t77 + t24 * t74 - t151 + t58, 0, 0, 0, 0, 0, -t135 - t171, t7 + t170; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t186, t185, t187, t189, t114, t148 * t182 + t183, (-t11 * t115 - t2) * t125 + (-t10 * t182 - t3) * t128 + t184;];
tau_reg = t1;
