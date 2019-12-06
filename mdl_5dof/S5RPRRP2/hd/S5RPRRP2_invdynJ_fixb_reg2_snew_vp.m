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
% Datum: 2019-12-05 18:02
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:01:44
% EndTime: 2019-12-05 18:01:48
% DurationCPUTime: 0.66s
% Computational Cost: add. (2246->143), mult. (3374->184), div. (0->0), fcn. (1805->8), ass. (0->97)
t94 = qJD(1) + qJD(3);
t124 = (qJD(5) * t94);
t144 = 2 * t124;
t101 = sin(qJ(4));
t104 = cos(qJ(4));
t92 = t94 ^ 2;
t79 = t104 * t92 * t101;
t73 = qJDD(4) + t79;
t143 = pkin(4) * t73;
t102 = sin(qJ(3));
t105 = cos(qJ(3));
t100 = cos(pkin(8));
t103 = sin(qJ(1));
t106 = cos(qJ(1));
t132 = t106 * g(2) + t103 * g(3);
t70 = qJDD(1) * pkin(1) + t132;
t108 = qJD(1) ^ 2;
t116 = -t103 * g(2) + t106 * g(3);
t71 = -t108 * pkin(1) - t116;
t99 = sin(pkin(8));
t119 = t100 * t70 - t99 * t71;
t113 = qJDD(1) * pkin(2) + t119;
t133 = t100 * t71 + t99 * t70;
t40 = -t108 * pkin(2) + t133;
t31 = t102 * t113 + t105 * t40;
t93 = qJDD(1) + qJDD(3);
t28 = -t92 * pkin(3) + t93 * pkin(7) + t31;
t97 = -g(1) + qJDD(2);
t22 = t101 * t97 + t104 * t28;
t125 = qJD(4) * t94;
t121 = t101 * t125;
t85 = t104 * t93;
t59 = t85 - t121;
t126 = qJ(5) * t101;
t72 = qJD(4) * pkin(4) - t94 * t126;
t109 = t59 * qJ(5) - qJD(4) * t72 + t104 * t144 + t22;
t120 = t104 * t125;
t128 = t101 * t93;
t58 = t120 + t128;
t87 = t104 * t97;
t110 = -t87 + (-t120 + t58) * qJ(5) - t143;
t96 = t104 ^ 2;
t137 = t96 * t92;
t95 = t101 ^ 2;
t131 = t95 + t96;
t68 = t131 * t92;
t142 = t101 * ((qJ(5) * t93 + t144 + t28) * t101 + t110) + t104 * (qJ(5) * t85 + (t68 - t137) * pkin(4) + t109);
t30 = -t102 * t40 + t105 * t113;
t27 = -t93 * pkin(3) - t92 * pkin(7) - t30;
t16 = t101 * t94 * t72 - t59 * pkin(4) - qJ(5) * t137 + qJDD(5) + t27;
t57 = 0.2e1 * t120 + t128;
t74 = qJDD(4) - t79;
t107 = qJD(4) ^ 2;
t138 = t95 * t92;
t76 = -t107 - t138;
t141 = t101 * (-qJ(5) * t76 + t16) + t104 * (-pkin(4) * t57 - qJ(5) * t74);
t129 = t101 * t73;
t77 = -t107 - t137;
t50 = t104 * t77 - t129;
t60 = t85 - 0.2e1 * t121;
t34 = t102 * t50 + t105 * t60;
t140 = pkin(1) * (t99 * (-t102 * t60 + t105 * t50) + t100 * t34) + pkin(2) * t34;
t130 = t101 * t28;
t21 = -t87 + t130;
t7 = t101 * t21 + t104 * t22;
t139 = -pkin(3) * t27 + pkin(7) * t7;
t136 = pkin(3) * t60 + pkin(7) * t50;
t127 = t104 * t74;
t51 = -t101 * t76 - t127;
t135 = -pkin(3) * t57 + pkin(7) * t51;
t65 = t131 * t93;
t134 = pkin(3) * t68 + pkin(7) * t65;
t123 = -t104 * t27 + t136;
t35 = t102 * t51 - t105 * t57;
t118 = pkin(1) * (t99 * (t102 * t57 + t105 * t51) + t100 * t35) + pkin(2) * t35 + t135;
t39 = t102 * t65 + t105 * t68;
t117 = pkin(1) * (t99 * (-t102 * t68 + t105 * t65) + t100 * t39) + pkin(2) * t39 + t134;
t66 = -t102 * t93 - t105 * t92;
t115 = t102 * t92 - t105 * t93;
t13 = -0.2e1 * t101 * t124 - t110 - t130;
t14 = -pkin(4) * t137 + t109;
t3 = -t101 * t13 + t104 * t14;
t112 = pkin(7) * t3 - t13 * t126 - pkin(3) * t16 + t104 * (-pkin(4) * t16 + qJ(5) * t14);
t111 = -t73 * t126 + t136 + t104 * (pkin(4) * t60 + qJ(5) * t77 - t16);
t69 = (t95 - t96) * t92;
t49 = -t101 * t74 + t104 * t76;
t48 = t129 + t104 * (t107 - t138);
t47 = t101 * t77 + t104 * t73;
t46 = t101 * (-t107 + t137) + t127;
t43 = (t58 + t120) * t101;
t42 = (t59 - t121) * t104;
t36 = t101 * t60 + t104 * t57;
t23 = t101 * t27;
t10 = t102 * t31 + t105 * t30;
t4 = t102 * t7 - t105 * t27;
t1 = t102 * t3 - t105 * t16;
t2 = [0, 0, 0, 0, 0, qJDD(1), t132, t116, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (t100 * qJDD(1) - t99 * t108) + t119, pkin(1) * (-t99 * qJDD(1) - t100 * t108) - t133, 0, pkin(1) * (t100 * t119 + t99 * t133), 0, 0, 0, 0, 0, t93, pkin(1) * (-t100 * t115 + t99 * t66) - pkin(2) * t115 + t30, pkin(1) * (t100 * t66 + t99 * t115) + pkin(2) * t66 - t31, 0, pkin(1) * (t99 * (-t102 * t30 + t105 * t31) + t100 * t10) + pkin(2) * t10, t43, t36, t48, t42, t46, 0, t123 + t140, t23 + t118, t117 + t7, pkin(1) * (t99 * (t102 * t27 + t105 * t7) + t100 * t4) + pkin(2) * t4 + t139, t43, t36, t48, t42, t46, 0, t111 + t140, t118 + t141, t117 + t142, pkin(1) * (t99 * (t102 * t16 + t105 * t3) + t100 * t1) + pkin(2) * t1 + t112; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97, 0, 0, 0, 0, 0, 0, t47, t49, 0, t101 * t22 - t104 * t21, 0, 0, 0, 0, 0, 0, t47, t49, 0, t101 * t14 + t104 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t93, t30, -t31, 0, 0, t43, t36, t48, t42, t46, 0, t123, t23 + t135, t134 + t7, t139, t43, t36, t48, t42, t46, 0, t111, t135 + t141, t134 + t142, t112; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79, t69, t128, t79, t85, qJDD(4), -t21, -t22, 0, 0, -t79, t69, t128, t79, t85, qJDD(4), t13 + t143, (t76 + t137) * pkin(4) - t109, -pkin(4) * t128, pkin(4) * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60, t57, -t68, t16;];
tauJ_reg = t2;
