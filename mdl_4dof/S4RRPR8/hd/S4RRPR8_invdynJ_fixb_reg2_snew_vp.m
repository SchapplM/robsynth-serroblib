% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S4RRPR8
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
%   pkin=[a2,a3,a4,d1,d2,d4]';
% 
% Output:
% tauJ_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S4RRPR8_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR8_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR8_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR8_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR8_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR8_invdynJ_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:08:23
% EndTime: 2019-12-31 17:08:26
% DurationCPUTime: 0.89s
% Computational Cost: add. (1655->188), mult. (3603->219), div. (0->0), fcn. (2010->6), ass. (0->122)
t92 = sin(qJ(4));
t93 = sin(qJ(2));
t95 = cos(qJ(4));
t96 = cos(qJ(2));
t49 = (-t93 * t92 - t96 * t95) * qJD(1);
t123 = qJD(1) * t96;
t124 = qJD(1) * t93;
t51 = -t92 * t123 + t95 * t124;
t142 = t51 * t49;
t84 = qJDD(2) - qJDD(4);
t156 = -t84 + t142;
t160 = t156 * t92;
t159 = t156 * t95;
t98 = qJD(1) ^ 2;
t130 = t96 * t98;
t75 = t93 * t130;
t69 = qJDD(2) - t75;
t131 = t96 * t69;
t89 = t93 ^ 2;
t139 = t89 * t98;
t97 = qJD(2) ^ 2;
t71 = t97 + t139;
t158 = pkin(5) * (-t93 * t71 + t131);
t85 = qJD(2) - qJD(4);
t143 = t49 * t85;
t118 = qJD(1) * qJD(2);
t114 = t96 * t118;
t120 = t93 * qJDD(1);
t62 = t114 + t120;
t115 = t93 * t118;
t119 = t96 * qJDD(1);
t63 = -t115 + t119;
t25 = t49 * qJD(4) + t95 * t62 - t92 * t63;
t17 = t25 - t143;
t157 = t63 - t115;
t144 = cos(qJ(1));
t94 = sin(qJ(1));
t104 = t144 * g(1) + t94 * g(2);
t122 = qJDD(1) * pkin(5);
t53 = -t98 * pkin(1) - t104 + t122;
t125 = t93 * qJ(3);
t147 = t96 * pkin(2);
t108 = -t125 - t147;
t59 = t108 * qJD(1);
t110 = qJD(1) * t59 + t53;
t155 = t110 * t93;
t90 = t96 ^ 2;
t138 = t90 * t98;
t154 = t131 + t93 * (-t97 + t138);
t112 = t92 * t62 + t95 * t63;
t14 = (qJD(4) + t85) * t51 + t112;
t47 = t49 ^ 2;
t48 = t51 ^ 2;
t83 = t85 ^ 2;
t152 = 2 * qJD(3);
t151 = pkin(2) + pkin(3);
t116 = t94 * g(1) - t144 * g(2);
t52 = qJDD(1) * pkin(1) + t98 * pkin(5) + t116;
t100 = t157 * pkin(2) + t52;
t126 = qJ(3) * t96;
t70 = -qJD(2) * pkin(3) - pkin(6) * t124;
t7 = t62 * qJ(3) + t63 * pkin(3) - pkin(6) * t138 + (qJD(2) * t126 + (t152 + t70) * t93) * qJD(1) + t100;
t150 = t92 * t7;
t149 = t93 * g(3);
t148 = t95 * t7;
t146 = t96 * g(3);
t145 = t97 * pkin(2);
t141 = t85 * t92;
t140 = t85 * t95;
t28 = t142 + t84;
t137 = t92 * t28;
t64 = -0.2e1 * t115 + t119;
t136 = t93 * t64;
t68 = qJDD(2) + t75;
t135 = t93 * t68;
t132 = t95 * t28;
t73 = -t97 - t138;
t129 = pkin(5) * (t96 * t73 - t135) + pkin(1) * t64;
t127 = t89 + t90;
t66 = t127 * t98;
t128 = pkin(1) * t66 + t127 * t122;
t117 = -qJDD(2) * pkin(2) - t97 * qJ(3) + qJDD(3);
t105 = t117 + t146;
t10 = -qJDD(2) * pkin(3) + (-t62 + t114) * pkin(6) + (-pkin(3) * t130 + t110) * t93 + t105;
t45 = t96 * t53;
t111 = qJDD(2) * qJ(3) + qJD(2) * t152 + t59 * t123 + t45;
t103 = t111 - t149;
t22 = t103 - t145;
t9 = -pkin(3) * t138 - t63 * pkin(6) + qJD(2) * t70 + t22;
t3 = -t95 * t10 + t92 * t9;
t39 = t93 * t53 + t146;
t40 = t45 - t149;
t113 = t93 * t39 + t96 * t40;
t4 = t92 * t10 + t95 * t9;
t1 = -t95 * t3 + t92 * t4;
t2 = t92 * t3 + t95 * t4;
t107 = t62 + t114;
t61 = 0.2e1 * t114 + t120;
t106 = t96 * t61 + t136;
t102 = t96 * t151 + pkin(1) + t125;
t101 = t117 + t155;
t99 = t124 * t152 + t100;
t67 = (t89 - t90) * t98;
t42 = -t48 + t83;
t41 = t47 - t83;
t38 = -t48 - t83;
t37 = t135 + t96 * (t97 - t139);
t36 = t107 * t93;
t35 = t157 * t96;
t31 = t48 - t47;
t27 = -t83 - t47;
t24 = -t51 * qJD(4) - t112;
t23 = t101 + t146;
t20 = -t92 * t38 + t132;
t19 = t95 * t38 + t137;
t18 = t25 + t143;
t13 = (qJD(4) - t85) * t51 + t112;
t12 = t95 * t27 - t160;
t11 = t92 * t27 + t159;
t6 = -t14 * t95 + t92 * t18;
t5 = -t14 * t92 - t95 * t18;
t8 = [0, 0, 0, 0, 0, qJDD(1), t116, t104, 0, 0, t36, t106, t37, t35, t154, 0, t96 * t52 + t129, -pkin(1) * t61 - t93 * t52 - t158, t113 + t128, pkin(1) * t52 + pkin(5) * t113, t36, t37, -t106, 0, -t154, t35, t96 * (pkin(2) * t64 + t99) + (t96 * t107 + t136) * qJ(3) + t129, t96 * (pkin(2) * t66 + t111 - t145) + (qJ(3) * t66 + t101) * t93 + t128, t93 * t99 + t158 + (pkin(1) + t147) * t61 + (t107 + t61) * t125, pkin(5) * (t96 * t22 + t93 * t23) + (pkin(1) - t108) * (t107 * qJ(3) + t99), t93 * (t51 * t141 + t95 * t25) + t96 * (t51 * t140 - t92 * t25), t93 * (-t95 * t13 - t92 * t17) + t96 * (t92 * t13 - t95 * t17), t93 * (-t92 * t42 + t159) + t96 * (-t95 * t42 - t160), t93 * (t49 * t140 - t92 * t24) + t96 * (-t49 * t141 - t95 * t24), t93 * (t95 * t41 + t137) + t96 * (-t92 * t41 + t132), (t93 * (-t49 * t95 - t51 * t92) + t96 * (t49 * t92 - t51 * t95)) * t85, t93 * (-pkin(6) * t11 + t150) + t96 * (-pkin(6) * t12 + t148) + pkin(5) * (t93 * t11 + t96 * t12) + t102 * t13, t93 * (-pkin(6) * t19 + t148) + t96 * (-pkin(6) * t20 - t150) + pkin(5) * (t93 * t19 + t96 * t20) + t102 * t17, t93 * (-pkin(6) * t5 - t1) + t96 * (-pkin(6) * t6 - t2) + pkin(5) * (t93 * t5 + t96 * t6) + t102 * (-t47 - t48), t102 * t7 + (pkin(5) - pkin(6)) * (t93 * t1 + t96 * t2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t75, t67, t120, t75, t119, qJDD(2), -t39, -t40, 0, 0, -t75, t120, -t67, qJDD(2), -t119, t75, pkin(2) * t68 + qJ(3) * t73 - t105 - t155, (-pkin(2) * t93 + t126) * qJDD(1), qJ(3) * t69 + (t71 - t97) * pkin(2) + t103, -pkin(2) * t23 + qJ(3) * t22, t142, -t31, -t18, -t142, t14, t84, qJ(3) * t12 - t151 * t11 + t3, qJ(3) * t20 - t151 * t19 + t4, qJ(3) * t6 - t151 * t5, qJ(3) * t2 - t151 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68, t120, -t71, t23, 0, 0, 0, 0, 0, 0, t11, t19, t5, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t142, t31, t18, t142, -t14, -t84, -t3, -t4, 0, 0;];
tauJ_reg = t8;
