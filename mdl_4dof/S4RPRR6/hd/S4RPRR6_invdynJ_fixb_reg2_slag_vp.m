% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4RPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RPRR6_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR6_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR6_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR6_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR6_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR6_invdynJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:52:44
% EndTime: 2019-12-31 16:52:46
% DurationCPUTime: 1.21s
% Computational Cost: add. (1998->218), mult. (5037->287), div. (0->0), fcn. (3811->12), ass. (0->123)
t114 = sin(qJ(4));
t117 = cos(qJ(4));
t112 = cos(pkin(7));
t118 = cos(qJ(3));
t154 = t118 * t112;
t142 = qJD(1) * t154;
t111 = sin(pkin(7));
t115 = sin(qJ(3));
t155 = t115 * t111;
t143 = qJD(1) * t155;
t65 = -t142 + t143;
t76 = t111 * t118 + t112 * t115;
t67 = t76 * qJD(1);
t131 = t114 * t65 - t117 * t67;
t146 = t112 * qJDD(1);
t147 = t111 * qJDD(1);
t144 = qJD(3) * t142 + t115 * t146 + t118 * t147;
t39 = qJD(3) * t143 - t144;
t132 = t115 * t147 - t118 * t146;
t70 = t76 * qJD(3);
t40 = qJD(1) * t70 + t132;
t125 = qJD(4) * t131 + t114 * t39 - t117 * t40;
t110 = qJD(3) + qJD(4);
t158 = t131 * t110;
t180 = t125 - t158;
t149 = qJD(4) * t117;
t150 = qJD(4) * t114;
t129 = -t114 * t40 - t117 * t39 - t149 * t65 - t150 * t67;
t34 = -t114 * t67 - t117 * t65;
t157 = t34 * t110;
t179 = t129 - t157;
t165 = t131 ^ 2;
t166 = t34 ^ 2;
t178 = t165 - t166;
t164 = t34 * t131;
t156 = qJDD(1) * pkin(1);
t116 = sin(qJ(1));
t119 = cos(qJ(1));
t175 = g(1) * t116 - g(2) * t119;
t130 = -qJDD(2) + t156 + t175;
t134 = g(1) * t119 + g(2) * t116;
t162 = pkin(5) + qJ(2);
t87 = t162 * t112;
t78 = qJD(1) * t87;
t160 = t115 * t78;
t86 = t162 * t111;
t77 = qJD(1) * t86;
t43 = -t118 * t77 - t160;
t24 = -pkin(6) * t67 + t43;
t23 = qJD(3) * pkin(3) + t24;
t44 = -t115 * t77 + t118 * t78;
t25 = -pkin(6) * t65 + t44;
t148 = qJD(1) * qJD(2);
t173 = qJDD(1) * t162 + t148;
t53 = t173 * t111;
t54 = t173 * t112;
t137 = -t115 * t54 - t118 * t53;
t17 = -t44 * qJD(3) + t137;
t6 = qJDD(3) * pkin(3) + t39 * pkin(6) + t17;
t151 = qJD(3) * t118;
t145 = t115 * t53 - t118 * t54 + t151 * t77;
t152 = qJD(3) * t115;
t16 = -t152 * t78 - t145;
t7 = -pkin(6) * t40 + t16;
t1 = (qJD(4) * t23 + t7) * t117 + t114 * t6 - t25 * t150;
t99 = pkin(2) * t112 + pkin(1);
t81 = -qJD(1) * t99 + qJD(2);
t45 = t65 * pkin(3) + t81;
t109 = pkin(7) + qJ(3);
t103 = qJ(4) + t109;
t97 = sin(t103);
t98 = cos(t103);
t177 = g(3) * t97 + t134 * t98 - t45 * t34 - t1;
t159 = t117 * t25;
t11 = t114 * t23 + t159;
t2 = -qJD(4) * t11 - t114 * t7 + t117 * t6;
t176 = -g(3) * t98 + t45 * t131 + t134 * t97 + t2;
t174 = qJ(2) * qJDD(1);
t172 = t67 ^ 2;
t171 = pkin(3) * t70;
t170 = t40 * pkin(3);
t163 = t67 * t65;
t47 = -t115 * t86 + t118 * t87;
t161 = t114 * t25;
t107 = t111 ^ 2;
t108 = t112 ^ 2;
t153 = t107 + t108;
t46 = -t115 * t87 - t118 * t86;
t136 = t153 * qJD(1) ^ 2;
t135 = 0.2e1 * t153;
t29 = -pkin(6) * t76 + t46;
t75 = -t154 + t155;
t30 = -pkin(6) * t75 + t47;
t14 = -t114 * t30 + t117 * t29;
t15 = t114 * t29 + t117 * t30;
t42 = -t114 * t75 + t117 * t76;
t80 = -qJDD(1) * t99 + qJDD(2);
t128 = t130 + t156;
t101 = sin(t109);
t102 = cos(t109);
t127 = -g(3) * t102 + t101 * t134;
t27 = -t86 * t151 + qJD(2) * t154 + (-qJD(2) * t111 - qJD(3) * t87) * t115;
t126 = -t175 + t80;
t123 = t135 * t148 - t134;
t28 = -qJD(2) * t76 - qJD(3) * t47;
t106 = -pkin(6) - t162;
t105 = qJDD(3) + qJDD(4);
t79 = pkin(3) * t102 + t99;
t69 = t111 * t152 - t112 * t151;
t62 = t65 ^ 2;
t51 = pkin(3) * t75 - t99;
t41 = t114 * t76 + t117 * t75;
t26 = t80 + t170;
t21 = t69 * pkin(6) + t28;
t20 = -t70 * pkin(6) + t27;
t19 = qJD(4) * t42 - t114 * t69 + t117 * t70;
t18 = t114 * t70 + t117 * t69 + t149 * t75 + t150 * t76;
t13 = t117 * t24 - t161;
t12 = -t114 * t24 - t159;
t10 = t117 * t23 - t161;
t4 = -qJD(4) * t15 - t114 * t20 + t117 * t21;
t3 = qJD(4) * t14 + t114 * t21 + t117 * t20;
t5 = [0, 0, 0, 0, 0, qJDD(1), t175, t134, 0, 0, t107 * qJDD(1), 0.2e1 * t111 * t146, 0, t108 * qJDD(1), 0, 0, t128 * t112, -t128 * t111, t135 * t174 + t123, t130 * pkin(1) + (t153 * t174 + t123) * qJ(2), -t39 * t76 - t67 * t69, t39 * t75 - t40 * t76 + t65 * t69 - t67 * t70, -qJD(3) * t69 + qJDD(3) * t76, t40 * t75 + t65 * t70, -qJD(3) * t70 - qJDD(3) * t75, 0, qJD(3) * t28 + qJDD(3) * t46 + t102 * t175 - t40 * t99 + t70 * t81 + t75 * t80, -t27 * qJD(3) - t47 * qJDD(3) - t101 * t175 + t99 * t39 - t81 * t69 + t80 * t76, -t16 * t75 - t17 * t76 - t27 * t65 - t28 * t67 + t39 * t46 - t40 * t47 + t43 * t69 - t44 * t70 - t134, t16 * t47 + t44 * t27 + t17 * t46 + t43 * t28 - t80 * t99 - g(1) * (-t116 * t99 + t119 * t162) - g(2) * (t116 * t162 + t119 * t99), t129 * t42 + t131 * t18, t125 * t42 - t129 * t41 + t131 * t19 - t18 * t34, t105 * t42 - t110 * t18, -t125 * t41 - t19 * t34, -t105 * t41 - t110 * t19, 0, t14 * t105 + t4 * t110 - t125 * t51 - t171 * t34 + t175 * t98 + t45 * t19 + t26 * t41, -t15 * t105 - t3 * t110 + t129 * t51 - t131 * t171 - t175 * t97 - t45 * t18 + t26 * t42, -t1 * t41 + t10 * t18 - t11 * t19 + t125 * t15 - t129 * t14 + t131 * t4 - t2 * t42 + t3 * t34 - t134, t1 * t15 + t11 * t3 + t2 * t14 + t10 * t4 + t26 * t51 + t45 * t171 - g(1) * (-t106 * t119 - t116 * t79) - g(2) * (-t106 * t116 + t119 * t79); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t146, t147, -t136, -qJ(2) * t136 - t130, 0, 0, 0, 0, 0, 0, 0.2e1 * t67 * qJD(3) + t132, (-t65 - t143) * qJD(3) + t144, -t62 - t172, t43 * t67 + t44 * t65 + t126, 0, 0, 0, 0, 0, 0, -t125 - t158, t129 + t157, -t165 - t166, -t10 * t131 - t11 * t34 + t126 + t170; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t163, -t62 + t172, (t65 - t143) * qJD(3) + t144, -t163, -t132, qJDD(3), -t81 * t67 + t127 + t137, g(3) * t101 + t81 * t65 + t134 * t102 + (t43 + t160) * qJD(3) + t145, 0, 0, t164, t178, t179, -t164, t180, t105, -t12 * t110 + (t105 * t117 - t110 * t150 + t34 * t67) * pkin(3) + t176, t13 * t110 + (-t105 * t114 - t110 * t149 + t131 * t67) * pkin(3) + t177, t10 * t34 - t11 * t131 - t12 * t131 - t13 * t34 + (t114 * t125 - t117 * t129 + (-t114 * t131 + t117 * t34) * qJD(4)) * pkin(3), -t10 * t12 - t11 * t13 + (t1 * t114 + t117 * t2 - t45 * t67 + (-t10 * t114 + t11 * t117) * qJD(4) + t127) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t164, t178, t179, -t164, t180, t105, t11 * t110 + t176, t10 * t110 + t177, 0, 0;];
tau_reg = t5;
