% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRRP1
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% tau_reg [5x22]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:00
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRP1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP1_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:00:00
% EndTime: 2019-12-05 18:00:04
% DurationCPUTime: 0.96s
% Computational Cost: add. (1327->197), mult. (2536->250), div. (0->0), fcn. (1587->8), ass. (0->118)
t104 = qJD(1) ^ 2;
t101 = cos(qJ(1));
t98 = sin(qJ(1));
t149 = g(1) * t98 - g(2) * t101;
t111 = -t104 * qJ(2) - t149;
t95 = qJ(3) + qJ(4);
t79 = sin(t95);
t80 = cos(t95);
t165 = g(3) * t79 - t149 * t80;
t100 = cos(qJ(3));
t89 = qJD(3) + qJD(4);
t99 = cos(qJ(4));
t123 = t89 * t99;
t134 = qJD(1) * qJD(3);
t97 = sin(qJ(3));
t130 = t97 * t134;
t145 = qJD(1) * t97;
t96 = sin(qJ(4));
t132 = t96 * t145;
t135 = t97 * qJDD(1);
t13 = -qJD(4) * t132 + (qJD(1) * t123 + qJDD(1) * t96) * t100 - t96 * t130 + t99 * t135;
t125 = t100 * t134;
t164 = t125 + t135;
t102 = -pkin(1) - pkin(6);
t154 = pkin(7) - t102;
t55 = t154 * t97;
t56 = t154 * t100;
t151 = -t99 * t55 - t96 * t56;
t63 = t102 * qJD(1) + qJD(2);
t31 = -pkin(7) * t145 + t97 * t63;
t26 = t96 * t31;
t139 = qJD(1) * t100;
t32 = -pkin(7) * t139 + t100 * t63;
t29 = qJD(3) * pkin(3) + t32;
t128 = t99 * t29 - t26;
t41 = t99 * t139 - t132;
t33 = t41 * qJ(5);
t162 = t33 - t128;
t49 = t96 * t100 + t99 * t97;
t39 = t49 * qJD(1);
t120 = g(1) * t101 + g(2) * t98;
t91 = qJD(1) * qJD(2);
t131 = 0.2e1 * t91;
t90 = qJDD(1) * qJ(2);
t161 = -t120 + 0.2e1 * t90 + t131;
t143 = qJD(4) * t96;
t133 = t100 * qJDD(1);
t144 = qJD(3) * t97;
t62 = t102 * qJDD(1) + qJDD(2);
t51 = t100 * t62;
t18 = -t63 * t144 + qJDD(3) * pkin(3) + t51 + (t130 - t133) * pkin(7);
t138 = qJD(3) * t100;
t19 = -t164 * pkin(7) + t63 * t138 + t97 * t62;
t160 = (qJD(4) * t29 + t19) * t99 - t31 * t143 + t96 * t18;
t159 = t41 ^ 2;
t7 = t89 * pkin(4) - t162;
t157 = t7 + t162;
t57 = pkin(3) * t145 + qJD(1) * qJ(2);
t23 = t39 * pkin(4) + qJD(5) + t57;
t156 = t23 * t41;
t155 = t41 * t39;
t27 = t99 * t31;
t21 = t89 * t49;
t50 = t99 * t100 - t96 * t97;
t87 = qJDD(3) + qJDD(4);
t153 = -t21 * t89 + t50 * t87;
t152 = t99 * t32 - t26;
t150 = t101 * pkin(1) + t98 * qJ(2);
t94 = t100 ^ 2;
t148 = t97 ^ 2 - t94;
t147 = t39 * qJ(5);
t83 = t97 * pkin(3);
t70 = qJ(2) + t83;
t146 = pkin(1) * qJDD(1);
t142 = qJD(4) * t99;
t64 = pkin(3) * t138 + qJD(2);
t103 = qJD(3) ^ 2;
t140 = -t103 - t104;
t136 = qJDD(3) * t97;
t127 = -t96 * t32 - t27;
t126 = t96 * t55 - t99 * t56;
t30 = t164 * pkin(3) + t90 + t91;
t122 = qJDD(2) - t146;
t119 = -t99 * t133 + t96 * t135;
t12 = t21 * qJD(1) + t119;
t118 = -t50 * t12 - t41 * t21;
t22 = t100 * t123 - t97 * t143 - t96 * t144;
t117 = -t22 * t89 - t49 * t87;
t116 = -t96 * t29 - t27;
t47 = t154 * t144;
t48 = qJD(3) * t56;
t114 = -t56 * t142 + t55 * t143 + t96 * t47 - t99 * t48;
t113 = 0.2e1 * qJ(2) * t134 + qJDD(3) * t102;
t112 = t13 * pkin(4) + qJDD(5) + t30;
t110 = t116 * qJD(4) + t99 * t18 - t96 * t19;
t109 = -t151 * qJD(4) + t99 * t47 + t96 * t48;
t1 = t87 * pkin(4) + t12 * qJ(5) - t41 * qJD(5) + t110;
t2 = -t13 * qJ(5) - t39 * qJD(5) + t160;
t9 = -t116 - t147;
t108 = t1 * t50 + t2 * t49 - t7 * t21 + t9 * t22 - t149;
t107 = -t102 * t103 + t161;
t106 = g(3) * t80 + t149 * t79 + t57 * t39 - t160;
t105 = -t57 * t41 + t110 + t165;
t88 = -qJ(5) - pkin(7) - pkin(6);
t82 = t101 * qJ(2);
t78 = qJDD(3) * t100;
t73 = t99 * pkin(3) + pkin(4);
t53 = pkin(4) * t79 + t83;
t38 = t39 ^ 2;
t17 = -t49 * qJ(5) + t151;
t16 = -t50 * qJ(5) + t126;
t14 = -t38 + t159;
t11 = -t33 + t152;
t10 = t127 + t147;
t6 = t41 * t89 - t13;
t4 = t21 * qJ(5) - t50 * qJD(5) + t109;
t3 = -t22 * qJ(5) - t49 * qJD(5) + t114;
t5 = [qJDD(1), t149, t120, qJDD(2) - 0.2e1 * t146 - t149, t161, -t122 * pkin(1) - g(1) * (-t98 * pkin(1) + t82) - g(2) * t150 + (t131 + t90) * qJ(2), t94 * qJDD(1) - 0.2e1 * t97 * t125, -0.2e1 * t97 * t133 + 0.2e1 * t148 * t134, -t103 * t97 + t78, -t103 * t100 - t136, 0, t113 * t100 + t107 * t97, t107 * t100 - t113 * t97, t118, t12 * t49 - t50 * t13 + t21 * t39 - t41 * t22, t153, t117, 0, t109 * t89 - t120 * t79 + t126 * t87 + t70 * t13 + t57 * t22 + t30 * t49 + t64 * t39, -t114 * t89 - t70 * t12 - t120 * t80 - t151 * t87 - t57 * t21 + t30 * t50 + t64 * t41, t16 * t12 - t17 * t13 - t3 * t39 - t4 * t41 - t108, t2 * t17 + t9 * t3 + t1 * t16 + t7 * t4 + t112 * (t49 * pkin(4) + t70) + t23 * (t22 * pkin(4) + t64) - g(1) * (t101 * t53 + t82 + (-pkin(1) + t88) * t98) - g(2) * (-t101 * t88 + t98 * t53 + t150); 0, 0, 0, qJDD(1), -t104, t122 + t111, 0, 0, 0, 0, 0, t140 * t97 + t78, t140 * t100 - t136, 0, 0, 0, 0, 0, -qJD(1) * t39 + t153, -qJD(1) * t41 + t117, -t49 * t13 - t22 * t39 - t118, -t23 * qJD(1) + t108; 0, 0, 0, 0, 0, 0, t100 * t104 * t97, -t148 * t104, t133, -t135, qJDD(3), g(3) * t97 + t111 * t100 + t51, g(3) * t100 + (-t111 - t62) * t97, t155, t14, -t119, t6, t87, -t127 * t89 + (-t39 * t139 - t89 * t143 + t99 * t87) * pkin(3) + t105, t152 * t89 + (-t41 * t139 - t89 * t142 - t96 * t87) * pkin(3) + t106, t73 * t12 + (t10 + t9) * t41 + (t11 - t7) * t39 + (-t13 * t96 + (-t39 * t99 + t41 * t96) * qJD(4)) * pkin(3), -pkin(4) * t156 + g(3) * t53 + t1 * t73 - t7 * t10 - t9 * t11 - t149 * (t100 * pkin(3) + pkin(4) * t80) + (-t23 * t139 + t2 * t96 + (-t7 * t96 + t9 * t99) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t155, t14, -t119, t6, t87, -t116 * t89 + t105, t128 * t89 + t106, pkin(4) * t12 - t157 * t39, t157 * t9 + (t1 - t156 + t165) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38 - t159, t9 * t39 + t7 * t41 + t112 - t120;];
tau_reg = t5;
