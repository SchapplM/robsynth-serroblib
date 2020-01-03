% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S4RRPR10
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
% Datum: 2019-12-31 17:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S4RRPR10_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR10_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR10_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR10_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR10_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR10_invdynJ_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:11:55
% EndTime: 2019-12-31 17:11:58
% DurationCPUTime: 0.84s
% Computational Cost: add. (1721->186), mult. (3674->218), div. (0->0), fcn. (2002->6), ass. (0->126)
t91 = sin(qJ(2));
t97 = qJD(1) ^ 2;
t135 = t91 * t97;
t94 = cos(qJ(2));
t121 = t94 * t135;
t65 = qJDD(2) - t121;
t132 = t94 * t65;
t86 = t91 ^ 2;
t81 = t86 * t97;
t96 = qJD(2) ^ 2;
t68 = -t81 - t96;
t161 = pkin(5) * (t91 * t68 + t132);
t129 = qJD(1) * t94;
t90 = sin(qJ(4));
t93 = cos(qJ(4));
t52 = t90 * qJD(2) + t93 * t129;
t54 = t93 * qJD(2) - t90 * t129;
t34 = t54 * t52;
t124 = qJD(1) * qJD(2);
t78 = t94 * t124;
t80 = t91 * qJDD(1);
t58 = t80 + t78;
t48 = qJDD(4) + t58;
t158 = -t34 + t48;
t160 = t158 * t90;
t159 = t158 * t93;
t150 = pkin(3) + pkin(5);
t104 = -pkin(3) * t124 - pkin(6) * t135 + g(3);
t127 = t91 * qJD(1);
t92 = sin(qJ(1));
t95 = cos(qJ(1));
t116 = t95 * g(1) + t92 * g(2);
t128 = qJDD(1) * pkin(5);
t44 = -t97 * pkin(1) - t116 + t128;
t140 = t91 * t44;
t130 = t91 * qJ(3);
t146 = t94 * pkin(2);
t113 = -t130 - t146;
t55 = t113 * qJD(1);
t107 = -qJDD(2) * pkin(2) - t96 * qJ(3) + t55 * t127 + qJDD(3) + t140;
t99 = t58 * pkin(3) - qJDD(2) * pkin(6) + t107;
t157 = t104 * t94 + t99;
t112 = t58 + t78;
t156 = t112 * qJ(3);
t87 = t94 ^ 2;
t142 = t87 * t97;
t155 = t132 + t91 * (-t96 + t142);
t64 = qJDD(2) + t121;
t138 = t91 * t64;
t125 = t94 * qJDD(1);
t77 = t91 * t124;
t60 = -0.2e1 * t77 + t125;
t70 = -t96 - t142;
t153 = pkin(5) * (-t94 * t70 + t138) - pkin(1) * t60;
t46 = t52 ^ 2;
t47 = t54 ^ 2;
t74 = qJD(4) + t127;
t71 = t74 ^ 2;
t152 = 2 * qJD(3);
t151 = -pkin(2) - pkin(6);
t83 = t91 * g(3);
t59 = -t77 + t125;
t117 = t90 * qJDD(2) + t93 * t59;
t103 = (-qJD(4) + t74) * t54 - t117;
t109 = t93 * qJDD(2) - t90 * t59;
t28 = -t52 * qJD(4) + t109;
t41 = t74 * t52;
t20 = t28 + t41;
t5 = t103 * t90 - t93 * t20;
t148 = t91 * t5;
t120 = t92 * g(1) - t95 * g(2);
t108 = -qJDD(1) * pkin(1) - t120;
t66 = pkin(3) * t127 - qJD(2) * pkin(6);
t73 = pkin(2) * t77;
t119 = qJD(3) * t127;
t76 = -0.2e1 * t119;
t7 = -t66 * t127 + t73 + t76 + (-pkin(3) * t87 - pkin(5)) * t97 + t151 * t59 - t156 + t108;
t147 = t93 * t7;
t145 = t94 * g(3);
t144 = t74 * t90;
t143 = t74 * t93;
t25 = t34 + t48;
t141 = t90 * t25;
t139 = t91 * t60;
t105 = (qJD(1) * t55 + t44) * t94 - t83 - t96 * pkin(2);
t123 = qJDD(2) * qJ(3);
t11 = t123 + t59 * pkin(3) - pkin(6) * t142 + (t152 + t66) * qJD(2) + t105;
t134 = t93 * t11;
t133 = t93 * t25;
t62 = t81 + t142;
t131 = pkin(1) * t62 + (t86 + t87) * t128;
t126 = qJD(4) + t74;
t122 = t91 * t34;
t3 = -t93 * t157 + t90 * t7;
t37 = t140 + t145;
t38 = t94 * t44 - t83;
t118 = t91 * t37 + t94 * t38;
t4 = t157 * t90 + t147;
t1 = -t93 * t3 + t90 * t4;
t115 = t90 * t3 + t93 * t4;
t111 = t94 * (-t81 + t96) + t138;
t106 = t94 * t151 - pkin(1) - t130;
t43 = t97 * pkin(5) - t108;
t22 = t107 + t145;
t101 = qJD(2) * t152 + t105;
t100 = t59 * pkin(2) + t43 - t73;
t21 = t101 + t123;
t98 = t100 + 0.2e1 * t119;
t63 = t81 - t142;
t57 = t80 + 0.2e1 * t78;
t40 = -t47 + t71;
t39 = t46 - t71;
t36 = t112 * t91;
t35 = (t59 - t77) * t94;
t32 = t47 - t46;
t31 = -t47 - t71;
t30 = t94 * t57 + t139;
t29 = -t71 - t46;
t27 = -t54 * qJD(4) - t117;
t23 = -t46 - t47;
t19 = t28 - t41;
t18 = -t126 * t52 + t109;
t15 = t126 * t54 + t117;
t12 = t93 * t31 - t141;
t9 = t90 * t29 + t159;
t2 = [0, 0, 0, 0, 0, qJDD(1), t120, t116, 0, 0, t36, t30, t111, t35, t155, 0, t94 * t43 - t153, -pkin(1) * t57 - t91 * t43 - t161, t118 + t131, pkin(1) * t43 + pkin(5) * t118, 0, -t111, -t155, t36, t30, t35, t91 * (qJ(3) * t62 + t107) + (pkin(2) * t62 + t21 + t83) * t94 + t131, t94 * (-pkin(2) * t60 - t100 + t76) + (-t94 * t112 - t139) * qJ(3) + t153, t91 * t98 + t161 + (pkin(1) + t146) * t57 + (t112 + t57) * t130, pkin(5) * (t94 * t21 + t91 * t22) + (pkin(1) - t113) * (t98 + t156), t122 + t94 * (-t54 * t143 - t90 * t28), t91 * t32 + t94 * (t90 * t15 - t93 * t19), t91 * t20 + t94 * (-t93 * t40 - t160), -t122 + t94 * (-t52 * t144 - t93 * t27), t91 * t103 + t94 * (-t90 * t39 - t133), t91 * t48 + t94 * (t52 * t90 + t54 * t93) * t74, t91 * (pkin(3) * t9 - t3) + t94 * (pkin(3) * t15 + t134) + pkin(5) * (t94 * t15 + t91 * t9) + t106 * (t93 * t29 - t160), (t150 * t12 - t90 * t99 - t147) * t91 + ((-t91 * t104 - t11) * t90 + t150 * t18) * t94 + t106 * (-t90 * t31 - t133), pkin(3) * t148 + t94 * (pkin(3) * t23 - t115) + pkin(5) * (t94 * t23 + t148) + t106 * (t103 * t93 + t90 * t20), t106 * t115 + t150 * (t91 * t1 + t94 * t11); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t121, t63, t80, t121, t125, qJDD(2), -t37, -t38, 0, 0, qJDD(2), -t80, -t125, -t121, t63, t121, (-pkin(2) * t91 + qJ(3) * t94) * qJDD(1), -pkin(2) * t64 - qJ(3) * t70 + t22, -pkin(2) * t68 + (qJDD(2) + t65) * qJ(3) + t101, -pkin(2) * t22 + qJ(3) * t21, -t54 * t144 + t93 * t28, -t93 * t15 - t90 * t19, -t90 * t40 + t159, t52 * t143 - t90 * t27, t93 * t39 - t141, (-t52 * t93 + t54 * t90) * t74, qJ(3) * t15 + t90 * t11 + t151 * t9, qJ(3) * t18 + t151 * t12 + t134, qJ(3) * t23 + t151 * t5 - t1, qJ(3) * t11 + t151 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, t64, t68, t22, 0, 0, 0, 0, 0, 0, t9, t12, t5, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, t32, t20, -t34, t103, t48, -t3, -t4, 0, 0;];
tauJ_reg = t2;
