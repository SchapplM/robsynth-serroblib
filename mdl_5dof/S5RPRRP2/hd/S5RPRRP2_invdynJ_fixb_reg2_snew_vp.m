% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPRRP2
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPRRP2_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP2_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP2_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP2_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:45:29
% EndTime: 2020-01-03 11:45:32
% DurationCPUTime: 0.76s
% Computational Cost: add. (2246->143), mult. (3374->184), div. (0->0), fcn. (1805->8), ass. (0->97)
t92 = qJD(1) + qJD(3);
t123 = (qJD(5) * t92);
t142 = 2 * t123;
t102 = cos(qJ(4));
t90 = t92 ^ 2;
t99 = sin(qJ(4));
t79 = t102 * t90 * t99;
t73 = qJDD(4) + t79;
t141 = pkin(4) * t73;
t100 = sin(qJ(3));
t103 = cos(qJ(3));
t101 = sin(qJ(1));
t104 = cos(qJ(1));
t115 = -t104 * g(2) - t101 * g(3);
t70 = qJDD(1) * pkin(1) + t115;
t106 = qJD(1) ^ 2;
t114 = t101 * g(2) - t104 * g(3);
t71 = -t106 * pkin(1) - t114;
t97 = sin(pkin(8));
t98 = cos(pkin(8));
t118 = t98 * t70 - t97 * t71;
t110 = qJDD(1) * pkin(2) + t118;
t128 = t97 * t70 + t98 * t71;
t40 = -t106 * pkin(2) + t128;
t31 = t100 * t110 + t103 * t40;
t91 = qJDD(1) + qJDD(3);
t28 = -t90 * pkin(3) + t91 * pkin(7) + t31;
t95 = -g(1) + qJDD(2);
t22 = t102 * t28 + t99 * t95;
t124 = qJD(4) * t92;
t121 = t99 * t124;
t85 = t102 * t91;
t59 = t85 - t121;
t126 = qJ(5) * t99;
t72 = qJD(4) * pkin(4) - t92 * t126;
t107 = t59 * qJ(5) - qJD(4) * t72 + t102 * t142 + t22;
t119 = t102 * t124;
t132 = t99 * t91;
t58 = t119 + t132;
t87 = t102 * t95;
t108 = -t87 + (-t119 + t58) * qJ(5) - t141;
t94 = t102 ^ 2;
t135 = t94 * t90;
t93 = t99 ^ 2;
t127 = t93 + t94;
t68 = t127 * t90;
t140 = t99 * ((qJ(5) * t91 + t142 + t28) * t99 + t108) + t102 * (qJ(5) * t85 + (t68 - t135) * pkin(4) + t107);
t30 = -t100 * t40 + t103 * t110;
t27 = -t91 * pkin(3) - t90 * pkin(7) - t30;
t16 = t99 * t92 * t72 - t59 * pkin(4) - qJ(5) * t135 + qJDD(5) + t27;
t57 = 0.2e1 * t119 + t132;
t74 = qJDD(4) - t79;
t105 = qJD(4) ^ 2;
t136 = t93 * t90;
t76 = -t105 - t136;
t139 = t99 * (-qJ(5) * t76 + t16) + t102 * (-pkin(4) * t57 - qJ(5) * t74);
t133 = t99 * t73;
t77 = -t105 - t135;
t49 = t102 * t77 - t133;
t60 = t85 - 0.2e1 * t121;
t34 = t100 * t49 + t103 * t60;
t138 = pkin(1) * (t97 * (-t100 * t60 + t103 * t49) + t98 * t34) + pkin(2) * t34;
t134 = t99 * t28;
t21 = -t87 + t134;
t7 = t102 * t22 + t99 * t21;
t137 = -pkin(3) * t27 + pkin(7) * t7;
t131 = pkin(3) * t60 + pkin(7) * t49;
t125 = t102 * t74;
t51 = -t99 * t76 - t125;
t130 = -pkin(3) * t57 + pkin(7) * t51;
t65 = t127 * t91;
t129 = pkin(3) * t68 + pkin(7) * t65;
t122 = -t102 * t27 + t131;
t35 = t100 * t51 - t103 * t57;
t117 = pkin(1) * (t97 * (t100 * t57 + t103 * t51) + t98 * t35) + pkin(2) * t35 + t130;
t39 = t100 * t65 + t103 * t68;
t116 = pkin(1) * (t97 * (-t100 * t68 + t103 * t65) + t98 * t39) + pkin(2) * t39 + t129;
t66 = -t100 * t91 - t103 * t90;
t113 = t100 * t90 - t103 * t91;
t13 = -0.2e1 * t99 * t123 - t108 - t134;
t14 = -pkin(4) * t135 + t107;
t3 = t102 * t14 - t99 * t13;
t111 = pkin(7) * t3 - t13 * t126 - pkin(3) * t16 + t102 * (-pkin(4) * t16 + qJ(5) * t14);
t109 = -t73 * t126 + t131 + t102 * (pkin(4) * t60 + qJ(5) * t77 - t16);
t69 = (t93 - t94) * t90;
t50 = t125 + t99 * (-t105 + t135);
t48 = t102 * (t105 - t136) + t133;
t47 = t102 * t76 - t99 * t74;
t46 = t102 * t73 + t99 * t77;
t43 = (t59 - t121) * t102;
t42 = (t58 + t119) * t99;
t36 = t102 * t57 + t99 * t60;
t23 = t99 * t27;
t10 = t100 * t31 + t103 * t30;
t4 = t100 * t7 - t103 * t27;
t1 = t100 * t3 - t103 * t16;
t2 = [0, 0, 0, 0, 0, qJDD(1), t115, t114, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (t98 * qJDD(1) - t97 * t106) + t118, pkin(1) * (-t97 * qJDD(1) - t98 * t106) - t128, 0, pkin(1) * (t98 * t118 + t97 * t128), 0, 0, 0, 0, 0, t91, pkin(1) * (-t113 * t98 + t97 * t66) - pkin(2) * t113 + t30, pkin(1) * (t97 * t113 + t98 * t66) + pkin(2) * t66 - t31, 0, pkin(1) * (t97 * (-t100 * t30 + t103 * t31) + t98 * t10) + pkin(2) * t10, t42, t36, t48, t43, t50, 0, t122 + t138, t23 + t117, t116 + t7, pkin(1) * (t97 * (t100 * t27 + t103 * t7) + t98 * t4) + pkin(2) * t4 + t137, t42, t36, t48, t43, t50, 0, t109 + t138, t117 + t139, t116 + t140, pkin(1) * (t97 * (t100 * t16 + t103 * t3) + t98 * t1) + pkin(2) * t1 + t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95, 0, 0, 0, 0, 0, 0, t46, t47, 0, -t102 * t21 + t99 * t22, 0, 0, 0, 0, 0, 0, t46, t47, 0, t102 * t13 + t99 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t91, t30, -t31, 0, 0, t42, t36, t48, t43, t50, 0, t122, t23 + t130, t129 + t7, t137, t42, t36, t48, t43, t50, 0, t109, t130 + t139, t129 + t140, t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79, t69, t132, t79, t85, qJDD(4), -t21, -t22, 0, 0, -t79, t69, t132, t79, t85, qJDD(4), t13 + t141, (t76 + t135) * pkin(4) - t107, -pkin(4) * t132, pkin(4) * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60, t57, -t68, t16;];
tauJ_reg = t2;
