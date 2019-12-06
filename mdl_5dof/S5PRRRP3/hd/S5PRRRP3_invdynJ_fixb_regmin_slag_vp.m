% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PRRRP3
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% tau_reg [5x20]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:44
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRRP3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP3_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP3_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP3_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:44:16
% EndTime: 2019-12-05 16:44:19
% DurationCPUTime: 0.90s
% Computational Cost: add. (1151->178), mult. (2512->235), div. (0->0), fcn. (1702->8), ass. (0->115)
t147 = cos(qJ(4));
t92 = sin(qJ(3));
t132 = t92 * qJD(1);
t154 = pkin(6) + pkin(7);
t93 = cos(qJ(3));
t62 = t154 * t93;
t34 = qJD(2) * t62 + t132;
t91 = sin(qJ(4));
t28 = t91 * t34;
t137 = qJD(3) * pkin(3);
t123 = qJD(2) * t154;
t33 = t93 * qJD(1) - t92 * t123;
t31 = t33 + t137;
t120 = t147 * t31 - t28;
t46 = t147 * t92 + t91 * t93;
t39 = t46 * qJD(2);
t135 = t39 * qJ(5);
t158 = t135 - t120;
t86 = pkin(8) + qJ(2);
t77 = sin(t86);
t78 = cos(t86);
t112 = g(1) * t78 + g(2) * t77;
t90 = qJ(3) + qJ(4);
t81 = sin(t90);
t82 = cos(t90);
t157 = -g(3) * t82 + t112 * t81;
t87 = qJD(3) + qJD(4);
t83 = t93 * pkin(3);
t149 = pkin(2) + t83;
t61 = t154 * t92;
t140 = t147 * t62 - t91 * t61;
t155 = t39 ^ 2;
t7 = t87 * pkin(4) - t158;
t153 = t7 + t158;
t116 = qJDD(2) * t147;
t130 = t92 * qJDD(2);
t110 = -t93 * t116 + t91 * t130;
t24 = t87 * t46;
t15 = qJD(2) * t24 + t110;
t142 = t91 * t92;
t109 = t87 * t142;
t117 = t147 * qJD(4);
t124 = t147 * t93;
t23 = -qJD(3) * t124 - t93 * t117 + t109;
t113 = qJD(2) * t124;
t134 = qJD(2) * t92;
t126 = t91 * t134;
t37 = -t113 + t126;
t148 = -t46 * t15 + t23 * t37;
t60 = t149 * qJD(2);
t25 = t37 * pkin(4) + qJD(5) - t60;
t146 = t25 * t39;
t145 = t39 * t37;
t144 = t77 * t82;
t143 = t78 * t82;
t141 = t147 * t33 - t28;
t139 = pkin(4) * t82 + t83;
t88 = t92 ^ 2;
t138 = -t93 ^ 2 + t88;
t136 = t37 * qJ(5);
t133 = qJD(4) * t91;
t131 = qJDD(1) - g(3);
t129 = t93 * qJDD(2);
t128 = qJD(2) * qJD(3);
t127 = t92 * t137;
t30 = t147 * t34;
t122 = qJD(3) * t154;
t121 = t92 * t128;
t119 = -t91 * t33 - t30;
t118 = -t147 * t61 - t91 * t62;
t115 = -t87 * t113 - t92 * t116 - t91 * t129;
t111 = g(1) * t77 - g(2) * t78;
t14 = qJD(2) * t109 + t115;
t45 = -t124 + t142;
t108 = -t45 * t14 + t39 * t24;
t84 = qJDD(3) + qJDD(4);
t107 = t23 * t87 - t46 * t84;
t106 = -t91 * t31 - t30;
t51 = t92 * t122;
t52 = t93 * t122;
t105 = -t61 * t117 - t62 * t133 - t147 * t51 - t91 * t52;
t35 = pkin(3) * t121 - qJDD(2) * t149;
t104 = -0.2e1 * pkin(2) * t128 - pkin(6) * qJDD(3);
t94 = qJD(3) ^ 2;
t103 = 0.2e1 * qJDD(2) * pkin(2) - pkin(6) * t94 + t111;
t95 = qJD(2) ^ 2;
t102 = t95 * pkin(2) - qJDD(2) * pkin(6) + t112;
t101 = t15 * pkin(4) + qJDD(5) + t35;
t79 = t93 * qJDD(1);
t19 = qJDD(3) * pkin(3) + t79 - qJDD(2) * t61 + (-t93 * t123 - t132) * qJD(3);
t22 = t33 * qJD(3) + t92 * qJDD(1) + qJDD(2) * t62;
t100 = t106 * qJD(4) + t147 * t19 - t91 * t22;
t99 = -t140 * qJD(4) - t147 * t52 + t91 * t51;
t98 = t31 * t117 - t34 * t133 + t147 * t22 + t91 * t19;
t97 = g(1) * t143 + g(2) * t144 + g(3) * t81 - t60 * t37 - t98;
t96 = t60 * t39 + t100 + t157;
t85 = -qJ(5) - t154;
t75 = t147 * pkin(3) + pkin(4);
t59 = qJDD(3) * t93 - t94 * t92;
t58 = qJDD(3) * t92 + t94 * t93;
t50 = pkin(2) + t139;
t36 = t37 ^ 2;
t21 = -t45 * qJ(5) + t140;
t20 = -t46 * qJ(5) + t118;
t16 = -t36 + t155;
t13 = -t24 * t87 - t45 * t84;
t12 = -t135 + t141;
t11 = t119 + t136;
t10 = -t106 - t136;
t5 = -t115 + (-t126 + t37) * t87;
t4 = t23 * qJ(5) - t46 * qJD(5) + t99;
t3 = -t24 * qJ(5) - t45 * qJD(5) + t105;
t2 = -t15 * qJ(5) - t37 * qJD(5) + t98;
t1 = t84 * pkin(4) + t14 * qJ(5) - t39 * qJD(5) + t100;
t6 = [t131, 0, 0, 0, 0, 0, 0, 0, 0, t59, -t58, 0, 0, 0, 0, 0, t13, t107, t108 + t148, -t1 * t45 - t10 * t23 + t2 * t46 - t7 * t24 - g(3); 0, qJDD(2), t111, t112, t88 * qJDD(2) + 0.2e1 * t93 * t121, -0.2e1 * t138 * t128 + 0.2e1 * t92 * t129, t58, t59, 0, t103 * t93 + t104 * t92, -t103 * t92 + t104 * t93, -t14 * t46 - t39 * t23, -t108 + t148, -t107, t13, 0, g(1) * t144 - g(2) * t143 + t118 * t84 + t37 * t127 - t149 * t15 - t60 * t24 + t35 * t45 + t99 * t87, -t105 * t87 - t111 * t81 + t39 * t127 + t14 * t149 - t140 * t84 + t60 * t23 + t35 * t46, -t1 * t46 - t10 * t24 + t20 * t14 - t21 * t15 - t2 * t45 + t7 * t23 - t3 * t37 - t4 * t39 - t112, t2 * t21 + t10 * t3 + t1 * t20 + t7 * t4 + t101 * (t45 * pkin(4) - t149) + t25 * (t24 * pkin(4) + t127) - g(1) * (-t77 * t50 - t78 * t85) - g(2) * (t78 * t50 - t77 * t85); 0, 0, 0, 0, -t92 * t95 * t93, t138 * t95, t130, t129, qJDD(3), -g(3) * t93 + t102 * t92 + t79, t102 * t93 - t131 * t92, t145, t16, t5, -t110, t84, -t119 * t87 + (-t87 * t133 - t37 * t134 + t147 * t84) * pkin(3) + t96, t141 * t87 + (-t87 * t117 - t39 * t134 - t91 * t84) * pkin(3) + t97, t75 * t14 + (t10 + t11) * t39 + (t12 - t7) * t37 + (-t15 * t91 + (-t147 * t37 + t39 * t91) * qJD(4)) * pkin(3), t1 * t75 - t10 * t12 - t7 * t11 - pkin(4) * t146 - g(3) * t139 - t112 * (-t92 * pkin(3) - pkin(4) * t81) + (-t25 * t134 + t2 * t91 + (t147 * t10 - t7 * t91) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t145, t16, t5, -t110, t84, -t106 * t87 + t96, t120 * t87 + t97, pkin(4) * t14 - t153 * t37, t153 * t10 + (t1 - t146 + t157) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36 - t155, t10 * t37 + t7 * t39 + t101 - t111;];
tau_reg = t6;
