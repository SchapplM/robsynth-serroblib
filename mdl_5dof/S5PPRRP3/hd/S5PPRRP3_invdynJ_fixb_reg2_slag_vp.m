% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PPRRP3
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
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:11
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PPRRP3_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP3_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP3_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRP3_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP3_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP3_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:11:22
% EndTime: 2019-12-05 15:11:24
% DurationCPUTime: 1.22s
% Computational Cost: add. (946->210), mult. (2251->274), div. (0->0), fcn. (1736->8), ass. (0->137)
t82 = cos(pkin(8));
t152 = qJD(1) * t82;
t80 = sin(pkin(8));
t153 = qJD(1) * t80;
t85 = sin(qJ(3));
t87 = cos(qJ(3));
t47 = t85 * qJD(2) + t87 * t153;
t41 = qJD(3) * pkin(6) + t47;
t84 = sin(qJ(4));
t86 = cos(qJ(4));
t26 = -t84 * t152 + t86 * t41;
t181 = t26 * qJD(4);
t104 = t86 * t152 + t84 * t41;
t180 = qJD(5) + t104;
t12 = -qJD(4) * pkin(4) + t180;
t169 = t80 * t85;
t129 = g(3) * t169;
t179 = t47 * qJD(3) + t129;
t178 = qJD(2) * qJD(3) + qJDD(1) * t80;
t15 = qJD(4) * qJ(5) + t26;
t89 = qJD(3) ^ 2;
t176 = qJDD(3) * t85 + t89 * t87;
t175 = t86 * pkin(4) + t84 * qJ(5);
t144 = qJDD(4) * pkin(4);
t174 = qJDD(5) - t144;
t167 = t80 * t87;
t128 = g(3) * t167;
t79 = t86 ^ 2;
t135 = t79 * qJDD(3);
t78 = t84 ^ 2;
t136 = t78 * qJDD(3);
t158 = t78 + t79;
t64 = t85 * t153;
t46 = t87 * qJD(2) - t64;
t173 = -t158 * t46 * qJD(3) - t128 + (t135 + t136) * pkin(6);
t171 = t82 ^ 2 * qJDD(1) - g(3);
t170 = t80 * t84;
t168 = t80 * t86;
t81 = sin(pkin(7));
t166 = t81 * t85;
t165 = t81 * t87;
t83 = cos(pkin(7));
t164 = t83 * t85;
t163 = t83 * t87;
t109 = pkin(4) * t84 - qJ(5) * t86;
t35 = t109 * qJD(4) - t84 * qJD(5);
t160 = t35 - t47;
t159 = -t78 + t79;
t88 = qJD(4) ^ 2;
t157 = t88 + t89;
t156 = qJD(3) * pkin(3);
t154 = pkin(6) * qJDD(4);
t151 = qJD(3) * t35;
t50 = -pkin(3) - t175;
t150 = qJD(3) * t50;
t149 = qJD(3) * t80;
t148 = qJD(3) * t84;
t147 = qJD(3) * t86;
t146 = qJD(3) * t87;
t145 = qJDD(3) * pkin(3);
t140 = qJDD(1) * t82;
t139 = qJDD(3) * t50;
t137 = qJDD(4) * t84;
t74 = t84 * qJDD(3);
t134 = t86 * qJDD(3);
t133 = t87 * qJDD(3);
t131 = qJD(3) * qJD(4);
t130 = qJDD(4) * qJ(5);
t127 = t86 * t167;
t126 = t84 * t46 * qJD(4) + t179 * t86;
t125 = -t85 * qJDD(2) - t178 * t87;
t124 = t85 * t149;
t123 = -g(1) * t81 + g(2) * t83;
t122 = qJD(1) * t149;
t121 = t84 * t131;
t120 = t86 * t131;
t24 = -t85 * t122 - t125;
t14 = qJDD(3) * pkin(6) + t24;
t119 = t84 * t14 + t86 * t140 + t181;
t40 = -t46 - t156;
t118 = t40 - t156;
t117 = (-qJDD(2) + t122) * t87 + t178 * t85;
t27 = -t46 + t150;
t116 = t27 + t150;
t114 = t84 * t120;
t36 = -t82 * t166 - t163;
t38 = -t82 * t164 + t165;
t112 = -g(1) * t38 - g(2) * t36;
t37 = t82 * t165 - t164;
t39 = t82 * t163 + t166;
t111 = -g(1) * t39 - g(2) * t37;
t110 = t86 * t14 - t84 * t140;
t108 = t12 * t84 + t15 * t86;
t107 = -t104 * t84 - t26 * t86;
t106 = -t89 * t85 + t133;
t42 = t84 * t167 + t82 * t86;
t103 = t176 * t80;
t102 = -pkin(6) * t88 + t112;
t17 = -t81 * t168 + t37 * t84;
t19 = -t83 * t168 + t39 * t84;
t100 = g(1) * t19 + g(2) * t17 + g(3) * t42 - t119;
t13 = t117 - t145;
t99 = t102 - t13 + t145;
t6 = t117 + t139 + t151;
t98 = t102 - t6 - t139;
t18 = t81 * t170 + t37 * t86;
t20 = t83 * t170 + t39 * t86;
t43 = -t82 * t84 + t127;
t97 = g(1) * t20 + g(2) * t18 + g(3) * t43 - t110;
t2 = t130 + (qJD(5) - t104) * qJD(4) + t110;
t3 = t119 + t174;
t96 = t2 * t86 + t3 * t84 + (t12 * t86 - t15 * t84) * qJD(4);
t4 = -t104 * qJD(4) + t110;
t95 = t4 * t86 + (t104 * t86 - t26 * t84) * qJD(4) + t119 * t84;
t94 = t100 + t181;
t21 = -qJD(4) * t127 + (qJD(4) * t82 + t124) * t84;
t93 = t21 * qJD(4) - t42 * qJDD(4) - t86 * t103 + t121 * t169;
t92 = t111 + t96;
t91 = t111 + t95;
t22 = -t42 * qJD(4) - t86 * t124;
t90 = t42 * t74 + t22 * t147 + t43 * t134 + (-t21 * t84 + (t42 * t86 - t43 * t84) * qJD(4)) * qJD(3);
t66 = pkin(6) * t167;
t65 = t84 * t89 * t86;
t54 = t159 * t89;
t52 = qJDD(4) * t86 - t88 * t84;
t51 = t88 * t86 + t137;
t48 = t109 * qJD(3);
t45 = -0.2e1 * t114 + t135;
t44 = 0.2e1 * t114 + t136;
t34 = t38 * pkin(3);
t33 = t36 * pkin(3);
t30 = t159 * t131 + t84 * t134;
t16 = t176 * t158;
t8 = (-0.2e1 * t121 + t134) * t87 + (-t157 * t86 - t137) * t85;
t7 = (-qJDD(4) * t85 - 0.2e1 * t87 * t131) * t86 + (t157 * t85 - t133) * t84;
t1 = t22 * qJD(4) + t43 * qJDD(4) + (-t85 * t120 - t176 * t84) * t80;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1) - g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t80 ^ 2 * qJDD(1) + t171, 0, 0, 0, 0, 0, 0, -t103, -t106 * t80, 0, (t117 * t85 + t24 * t87 + (-t46 * t87 - t47 * t85) * qJD(3)) * t80 + t171, 0, 0, 0, 0, 0, 0, t93, -t1, t90, -t104 * t21 + t26 * t22 + t4 * t43 + t119 * t42 - g(3) + (t13 * t85 + t40 * t146) * t80, 0, 0, 0, 0, 0, 0, t93, t90, t1, -t12 * t21 + t15 * t22 + t2 * t43 + t3 * t42 - g(3) + (t146 * t27 + t6 * t85) * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2) + t123, 0, 0, 0, 0, 0, 0, t106, -t176, 0, -t117 * t87 + t24 * t85 + (-t46 * t85 + t47 * t87) * qJD(3) + t123, 0, 0, 0, 0, 0, 0, t8, t7, t16, (-t107 * qJD(3) - t13) * t87 + (qJD(3) * t40 + t95) * t85 + t123, 0, 0, 0, 0, 0, 0, t8, t16, -t7, (qJD(3) * t108 - t6) * t87 + (qJD(3) * t27 + t96) * t85 + t123; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), t112 - t117 + t179, t128 + (t46 + t64) * qJD(3) - t111 + t125, 0, 0, t44, 0.2e1 * t30, t51, t45, t52, 0, (t118 * qJD(4) - t154) * t84 + t99 * t86 + t126, (-t154 + (t118 + t46) * qJD(4)) * t86 + (-t99 - t179) * t84, t91 + t173, -t13 * pkin(3) - t40 * t47 - g(1) * t34 - g(2) * t33 - g(3) * (-pkin(3) * t169 + t66) + t107 * t46 + t91 * pkin(6), t44, t51, -0.2e1 * t30, 0, -t52, t45, (t116 * qJD(4) - t154) * t84 + (t98 - t151) * t86 + t126, t92 + t173, (t154 + (-t116 - t46) * qJD(4)) * t86 + (-qJD(3) * t160 + t129 + t98) * t84, -g(1) * (t175 * t38 + t34) - g(2) * (t175 * t36 + t33) - g(3) * t66 - t108 * t46 + t160 * t27 + t92 * pkin(6) + (t6 - t129) * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65, -t54, t74, t65, t134, qJDD(4), -t40 * t148 + t94, -t40 * t147 + t97, 0, 0, -t65, t74, t54, qJDD(4), -t134, t65, 0.2e1 * t144 - qJDD(5) + (-t27 * t84 + t48 * t86) * qJD(3) + t94, -t109 * qJDD(3), 0.2e1 * t130 + (t27 * t86 + t48 * t84) * qJD(3) + 0.2e1 * qJD(4) * qJD(5) - t97, t2 * qJ(5) - t3 * pkin(4) - t27 * t48 - t12 * t26 - g(1) * (-t19 * pkin(4) + t20 * qJ(5)) - g(2) * (-t17 * pkin(4) + t18 * qJ(5)) - g(3) * (-t42 * pkin(4) + t43 * qJ(5)) + t180 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(4) - t65, t74, -t78 * t89 - t88, -t15 * qJD(4) + t148 * t27 - t100 + t174;];
tau_reg = t5;
