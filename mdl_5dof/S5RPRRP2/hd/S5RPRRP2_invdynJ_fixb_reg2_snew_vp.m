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
% Datum: 2022-01-23 09:28
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-23 09:28:12
% EndTime: 2022-01-23 09:28:14
% DurationCPUTime: 0.67s
% Computational Cost: add. (2246->143), mult. (3374->184), div. (0->0), fcn. (1805->8), ass. (0->97)
t93 = qJD(1) + qJD(3);
t124 = (qJD(5) * t93);
t143 = 2 * t124;
t100 = sin(qJ(4));
t103 = cos(qJ(4));
t91 = t93 ^ 2;
t79 = t103 * t91 * t100;
t73 = qJDD(4) + t79;
t142 = pkin(4) * t73;
t101 = sin(qJ(3));
t104 = cos(qJ(3));
t102 = sin(qJ(1));
t105 = cos(qJ(1));
t119 = t102 * g(1) - t105 * g(2);
t70 = qJDD(1) * pkin(1) + t119;
t107 = qJD(1) ^ 2;
t115 = t105 * g(1) + t102 * g(2);
t71 = -t107 * pkin(1) - t115;
t98 = sin(pkin(8));
t99 = cos(pkin(8));
t118 = t99 * t70 - t98 * t71;
t112 = qJDD(1) * pkin(2) + t118;
t132 = t98 * t70 + t99 * t71;
t40 = -t107 * pkin(2) + t132;
t31 = t101 * t112 + t104 * t40;
t92 = qJDD(1) + qJDD(3);
t28 = -t91 * pkin(3) + t92 * pkin(7) + t31;
t96 = -g(3) + qJDD(2);
t22 = t100 * t96 + t103 * t28;
t125 = qJD(4) * t93;
t121 = t100 * t125;
t85 = t103 * t92;
t59 = t85 - t121;
t126 = qJ(5) * t100;
t72 = qJD(4) * pkin(4) - t93 * t126;
t108 = t59 * qJ(5) - qJD(4) * t72 + t103 * t143 + t22;
t120 = t103 * t125;
t128 = t100 * t92;
t58 = t120 + t128;
t87 = t103 * t96;
t109 = -t87 + (-t120 + t58) * qJ(5) - t142;
t95 = t103 ^ 2;
t136 = t95 * t91;
t94 = t100 ^ 2;
t131 = t94 + t95;
t68 = t131 * t91;
t141 = t100 * ((qJ(5) * t92 + t143 + t28) * t100 + t109) + t103 * (qJ(5) * t85 + (t68 - t136) * pkin(4) + t108);
t30 = -t101 * t40 + t104 * t112;
t27 = -t92 * pkin(3) - t91 * pkin(7) - t30;
t16 = t100 * t93 * t72 - t59 * pkin(4) - qJ(5) * t136 + qJDD(5) + t27;
t57 = 0.2e1 * t120 + t128;
t74 = qJDD(4) - t79;
t106 = qJD(4) ^ 2;
t137 = t94 * t91;
t76 = -t106 - t137;
t140 = t100 * (-qJ(5) * t76 + t16) + t103 * (-pkin(4) * t57 - qJ(5) * t74);
t129 = t100 * t73;
t77 = -t106 - t136;
t50 = t103 * t77 - t129;
t60 = t85 - 0.2e1 * t121;
t34 = t101 * t50 + t104 * t60;
t139 = pkin(1) * (t98 * (-t101 * t60 + t104 * t50) + t99 * t34) + pkin(2) * t34;
t130 = t100 * t28;
t21 = -t87 + t130;
t7 = t100 * t21 + t103 * t22;
t138 = -pkin(3) * t27 + pkin(7) * t7;
t135 = pkin(3) * t60 + pkin(7) * t50;
t127 = t103 * t74;
t51 = -t100 * t76 - t127;
t134 = -pkin(3) * t57 + pkin(7) * t51;
t65 = t131 * t92;
t133 = pkin(3) * t68 + pkin(7) * t65;
t123 = -t103 * t27 + t135;
t35 = t101 * t51 - t104 * t57;
t117 = pkin(1) * (t98 * (t101 * t57 + t104 * t51) + t99 * t35) + pkin(2) * t35 + t134;
t39 = t101 * t65 + t104 * t68;
t116 = pkin(1) * (t98 * (-t101 * t68 + t104 * t65) + t99 * t39) + pkin(2) * t39 + t133;
t66 = -t101 * t92 - t104 * t91;
t114 = t101 * t91 - t104 * t92;
t13 = -0.2e1 * t100 * t124 - t109 - t130;
t14 = -pkin(4) * t136 + t108;
t3 = -t100 * t13 + t103 * t14;
t111 = pkin(7) * t3 - t13 * t126 - pkin(3) * t16 + t103 * (-pkin(4) * t16 + qJ(5) * t14);
t110 = -t73 * t126 + t135 + t103 * (pkin(4) * t60 + qJ(5) * t77 - t16);
t69 = (t94 - t95) * t91;
t49 = -t100 * t74 + t103 * t76;
t48 = t129 + t103 * (t106 - t137);
t47 = t100 * t77 + t103 * t73;
t46 = t100 * (-t106 + t136) + t127;
t43 = (t58 + t120) * t100;
t42 = (t59 - t121) * t103;
t36 = t100 * t60 + t103 * t57;
t23 = t100 * t27;
t10 = t101 * t31 + t104 * t30;
t4 = t101 * t7 - t104 * t27;
t1 = t101 * t3 - t104 * t16;
t2 = [0, 0, 0, 0, 0, qJDD(1), t119, t115, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (t99 * qJDD(1) - t98 * t107) + t118, pkin(1) * (-t98 * qJDD(1) - t99 * t107) - t132, 0, pkin(1) * (t99 * t118 + t98 * t132), 0, 0, 0, 0, 0, t92, pkin(1) * (-t114 * t99 + t98 * t66) - pkin(2) * t114 + t30, pkin(1) * (t98 * t114 + t99 * t66) + pkin(2) * t66 - t31, 0, pkin(1) * (t98 * (-t101 * t30 + t104 * t31) + t99 * t10) + pkin(2) * t10, t43, t36, t48, t42, t46, 0, t123 + t139, t23 + t117, t116 + t7, pkin(1) * (t98 * (t101 * t27 + t104 * t7) + t99 * t4) + pkin(2) * t4 + t138, t43, t36, t48, t42, t46, 0, t110 + t139, t117 + t140, t116 + t141, pkin(1) * (t98 * (t101 * t16 + t104 * t3) + t99 * t1) + pkin(2) * t1 + t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96, 0, 0, 0, 0, 0, 0, t47, t49, 0, t100 * t22 - t103 * t21, 0, 0, 0, 0, 0, 0, t47, t49, 0, t100 * t14 + t103 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t92, t30, -t31, 0, 0, t43, t36, t48, t42, t46, 0, t123, t23 + t134, t133 + t7, t138, t43, t36, t48, t42, t46, 0, t110, t134 + t140, t133 + t141, t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79, t69, t128, t79, t85, qJDD(4), -t21, -t22, 0, 0, -t79, t69, t128, t79, t85, qJDD(4), t13 + t142, (t76 + t136) * pkin(4) - t108, -pkin(4) * t128, pkin(4) * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60, t57, -t68, t16;];
tauJ_reg = t2;
