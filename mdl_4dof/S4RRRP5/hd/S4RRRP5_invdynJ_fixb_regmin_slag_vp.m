% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4RRRP5
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% tau_reg [4x21]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RRRP5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP5_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP5_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP5_invdynJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:17:07
% EndTime: 2019-12-31 17:17:09
% DurationCPUTime: 0.78s
% Computational Cost: add. (1164->196), mult. (2603->257), div. (0->0), fcn. (1674->8), ass. (0->118)
t81 = qJD(2) + qJD(3);
t149 = pkin(6) + pkin(5);
t84 = qJ(2) + qJ(3);
t78 = sin(t84);
t79 = cos(t84);
t136 = t79 * pkin(3) + t78 * qJ(4);
t80 = qJDD(2) + qJDD(3);
t75 = t80 * qJ(4);
t76 = t81 * qJD(4);
t153 = t75 + t76;
t147 = cos(qJ(3));
t85 = sin(qJ(3));
t86 = sin(qJ(2));
t88 = cos(qJ(2));
t43 = t147 * t86 + t85 * t88;
t39 = t43 * qJD(1);
t77 = t80 * pkin(3);
t152 = qJDD(4) - t77;
t87 = sin(qJ(1));
t89 = cos(qJ(1));
t110 = g(1) * t87 - g(2) * t89;
t49 = t149 * t86;
t50 = t149 * t88;
t31 = t147 * t50 - t49 * t85;
t105 = -t147 * t49 - t50 * t85;
t122 = qJD(2) * t149;
t45 = t86 * t122;
t47 = t88 * t122;
t7 = qJD(3) * t105 - t147 * t45 - t85 * t47;
t151 = t110 * t78 + t31 * t80 + t7 * t81;
t150 = t39 ^ 2;
t148 = pkin(2) * t88;
t46 = qJD(1) * t50;
t139 = t85 * t46;
t133 = qJD(2) * pkin(2);
t44 = qJD(1) * t49;
t41 = -t44 + t133;
t18 = t147 * t41 - t139;
t146 = t18 * t81;
t124 = t147 * t46;
t19 = t41 * t85 + t124;
t145 = t19 * t81;
t123 = t147 * t88;
t114 = qJD(1) * t123;
t132 = qJD(1) * t86;
t125 = t85 * t132;
t37 = -t114 + t125;
t144 = t39 * t37;
t143 = t78 * t87;
t142 = t78 * t89;
t141 = t79 * t87;
t140 = t79 * t89;
t138 = t85 * t86;
t119 = qJD(3) * t147;
t23 = -t147 * t44 - t139;
t137 = pkin(2) * t119 + qJD(4) - t23;
t82 = t86 ^ 2;
t135 = -t88 ^ 2 + t82;
t134 = qJ(4) * t79;
t131 = qJD(3) * t85;
t130 = qJD(4) - t18;
t129 = t86 * qJDD(1);
t128 = t88 * qJDD(1);
t127 = qJD(1) * qJD(2);
t126 = t86 * t133;
t74 = pkin(1) + t148;
t121 = t86 * t127;
t120 = t88 * t127;
t118 = qJDD(1) * t147;
t28 = qJDD(2) * pkin(2) + t149 * (-t120 - t129);
t29 = t149 * (-t121 + t128);
t117 = t119 * t41 - t131 * t46 + t147 * t29 + t28 * t85;
t116 = t119 * t46 + t131 * t41 - t147 * t28 + t29 * t85;
t115 = -t114 * t81 - t86 * t118 - t85 * t128;
t22 = -t44 * t85 + t124;
t113 = pkin(2) * t131 - t22;
t112 = -pkin(2) * t86 - pkin(3) * t78;
t111 = g(1) * t89 + g(2) * t87;
t109 = -t118 * t88 + t129 * t85;
t48 = t74 * qJD(1);
t107 = t81 * t138;
t16 = pkin(3) * t39 + qJ(4) * t37;
t106 = t74 + t136;
t104 = g(1) * t140 + g(2) * t141 + g(3) * t78 - t117;
t34 = pkin(2) * t121 - qJDD(1) * t74;
t103 = -0.2e1 * pkin(1) * t127 - pkin(5) * qJDD(2);
t100 = g(1) * t142 + g(2) * t143 - g(3) * t79 - t116;
t8 = qJD(3) * t31 + t147 * t47 - t85 * t45;
t99 = g(1) * t141 - g(2) * t140 + t105 * t80 - t8 * t81;
t12 = pkin(3) * t37 - qJ(4) * t39 - t48;
t98 = -t12 * t37 - t104;
t97 = -t48 * t37 + t104;
t91 = qJD(2) ^ 2;
t96 = 0.2e1 * qJDD(1) * pkin(1) - pkin(5) * t91 + t110;
t92 = qJD(1) ^ 2;
t95 = pkin(1) * t92 - pkin(5) * qJDD(1) + t111;
t94 = t48 * t39 + t100;
t93 = t12 * t39 - t100 + t152;
t25 = t81 * t43;
t73 = -pkin(2) * t147 - pkin(3);
t69 = pkin(2) * t85 + qJ(4);
t52 = t89 * t134;
t51 = t87 * t134;
t42 = -t123 + t138;
t24 = -qJD(2) * t123 - t119 * t88 + t107;
t17 = pkin(3) * t42 - qJ(4) * t43 - t74;
t15 = qJ(4) * t81 + t19;
t14 = pkin(2) * t132 + t16;
t13 = -pkin(3) * t81 + t130;
t11 = -t37 ^ 2 + t150;
t10 = qJD(1) * t25 + t109;
t9 = qJD(1) * t107 + t115;
t5 = -t115 + (-t125 + t37) * t81;
t4 = pkin(3) * t25 + qJ(4) * t24 - qJD(4) * t43 + t126;
t3 = t116 + t152;
t2 = t117 + t153;
t1 = pkin(3) * t10 + qJ(4) * t9 - qJD(4) * t39 + t34;
t6 = [qJDD(1), t110, t111, qJDD(1) * t82 + 0.2e1 * t120 * t86, -0.2e1 * t127 * t135 + 0.2e1 * t128 * t86, qJDD(2) * t86 + t88 * t91, qJDD(2) * t88 - t86 * t91, 0, t103 * t86 + t88 * t96, t103 * t88 - t86 * t96, -t24 * t39 - t43 * t9, -t10 * t43 + t24 * t37 - t25 * t39 + t42 * t9, -t24 * t81 + t43 * t80, -t25 * t81 - t42 * t80, 0, -t10 * t74 + t126 * t37 - t25 * t48 + t34 * t42 + t99, t126 * t39 + t24 * t48 + t34 * t43 + t74 * t9 - t151, t1 * t42 + t10 * t17 + t12 * t25 + t37 * t4 + t99, -t10 * t31 + t105 * t9 - t13 * t24 - t15 * t25 - t2 * t42 + t3 * t43 - t37 * t7 + t39 * t8 - t111, -t1 * t43 + t12 * t24 + t17 * t9 - t39 * t4 + t151, t1 * t17 + t12 * t4 + t13 * t8 + t15 * t7 + t2 * t31 - t3 * t105 + (-g(1) * t149 - g(2) * t106) * t89 + (g(1) * t106 - g(2) * t149) * t87; 0, 0, 0, -t86 * t92 * t88, t135 * t92, t129, t128, qJDD(2), -g(3) * t88 + t86 * t95, g(3) * t86 + t88 * t95, t144, t11, t5, -t109, t80, t22 * t81 + (-t131 * t81 - t132 * t37 + t147 * t80) * pkin(2) + t94, t23 * t81 + (-t119 * t81 - t132 * t39 - t80 * t85) * pkin(2) + t97, -t113 * t81 - t14 * t37 - t73 * t80 - t93, -t10 * t69 - t73 * t9 + (t113 + t15) * t39 + (t13 - t137) * t37, t137 * t81 + t14 * t39 + t69 * t80 + t153 + t98, t2 * t69 + t3 * t73 - t12 * t14 - g(1) * (t112 * t89 + t52) - g(2) * (t112 * t87 + t51) - g(3) * (t136 + t148) + t137 * t15 + t113 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t144, t11, t5, -t109, t80, t94 + t145, t97 + t146, -t16 * t37 + t145 + t77 - t93, pkin(3) * t9 - qJ(4) * t10 + (t15 - t19) * t39 + (t13 - t130) * t37, t16 * t39 - t146 + 0.2e1 * t75 + 0.2e1 * t76 + t98, t2 * qJ(4) - t3 * pkin(3) - t12 * t16 - t13 * t19 - g(1) * (-pkin(3) * t142 + t52) - g(2) * (-pkin(3) * t143 + t51) - g(3) * t136 + t130 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t80 + t144, t5, -t81 ^ 2 - t150, -t15 * t81 + t93;];
tau_reg = t6;
