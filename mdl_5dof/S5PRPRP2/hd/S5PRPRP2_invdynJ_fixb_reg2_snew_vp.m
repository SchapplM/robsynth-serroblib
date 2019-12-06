% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PRPRP2
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
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PRPRP2_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP2_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP2_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP2_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP2_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP2_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:31:05
% EndTime: 2019-12-05 15:31:09
% DurationCPUTime: 0.92s
% Computational Cost: add. (1635->158), mult. (3746->208), div. (0->0), fcn. (2357->8), ass. (0->122)
t90 = sin(pkin(8));
t84 = t90 ^ 2;
t92 = cos(pkin(8));
t85 = t92 ^ 2;
t160 = t84 + t85;
t152 = 2 * qJD(3);
t98 = qJD(2) ^ 2;
t143 = t84 * t98;
t94 = sin(qJ(4));
t96 = cos(qJ(4));
t65 = t96 * t94 * t143;
t130 = t92 * qJDD(2);
t75 = -qJDD(4) + t130;
t52 = -t65 - t75;
t159 = pkin(4) * t52;
t135 = qJD(2) * t96;
t154 = (qJD(4) * t135 + qJDD(2) * t94) * t90;
t121 = t90 * t135;
t134 = t92 * qJD(2);
t76 = -qJD(4) + t134;
t59 = t76 * t121;
t42 = t154 - t59;
t158 = t90 * t42;
t91 = sin(pkin(7));
t93 = cos(pkin(7));
t66 = t91 * g(1) - t93 * g(2);
t67 = -t93 * g(1) - t91 * g(2);
t95 = sin(qJ(2));
t97 = cos(qJ(2));
t108 = -t95 * t66 - t97 * t67;
t41 = -t98 * pkin(2) + qJDD(2) * qJ(3) - t108;
t157 = qJD(2) * t152 + t41;
t119 = qJD(2) * qJD(5) * t90;
t50 = -t76 * pkin(4) - qJ(5) * t121;
t133 = -g(3) + qJDD(1);
t29 = t90 * t133 + t157 * t92;
t144 = t92 * pkin(3);
t109 = -t90 * pkin(6) - t144;
t61 = t109 * qJD(2);
t24 = t61 * t134 + t29;
t116 = t97 * t66 - t95 * t67;
t86 = t98 * qJ(3);
t89 = qJDD(2) * pkin(2);
t40 = qJDD(3) - t116 - t86 - t89;
t35 = t109 * qJDD(2) + t40;
t9 = t96 * t24 + t94 * t35;
t106 = -qJ(5) * t154 - 0.2e1 * t94 * t119 + t76 * t50 + t9;
t87 = t94 ^ 2;
t127 = t87 * t143;
t88 = t96 ^ 2;
t126 = t88 * t143;
t74 = t76 ^ 2;
t46 = -t126 - t74;
t156 = -t106 + (t46 + t127) * pkin(4);
t81 = t92 * t133;
t113 = pkin(4) * t154 - qJ(5) * t127 + qJDD(5) - t81;
t117 = -t50 * t96 - t61;
t7 = (t41 + (t152 - t117) * qJD(2)) * t90 + t113;
t153 = t160 * t86 + t40 - t89;
t136 = qJD(2) * t94;
t122 = t90 * t136;
t114 = t76 * t122;
t131 = t90 * qJDD(2);
t71 = qJD(4) * t122;
t54 = t96 * t131 - t71;
t8 = t94 * t24 - t96 * t35;
t104 = -t8 + (t114 - t54) * qJ(5) + t159;
t111 = t96 * t119;
t70 = -0.2e1 * t111;
t5 = t104 + t70;
t151 = pkin(4) * t5;
t51 = -t65 + t75;
t142 = t94 * t51;
t33 = t96 * t46 + t142;
t150 = pkin(3) * t33;
t141 = t96 * t52;
t55 = -t74 - t127;
t37 = t94 * t55 + t141;
t149 = pkin(3) * t37;
t123 = t76 * t136;
t132 = qJDD(2) * t96;
t45 = -t71 + (-t123 + t132) * t90;
t148 = pkin(4) * t45;
t43 = t59 + t154;
t20 = -t94 * t43 - t96 * t45;
t147 = pkin(6) * t20;
t146 = pkin(6) * t33;
t145 = pkin(6) * t37;
t21 = -t96 * t43 + t94 * t45;
t56 = (t87 + t88) * t143;
t140 = qJ(3) * (t92 * t21 - t90 * t56) - pkin(2) * t20;
t34 = -t94 * t46 + t96 * t51;
t44 = -t71 + (t123 + t132) * t90;
t139 = qJ(3) * (t92 * t34 + t90 * t44) - pkin(2) * t33;
t38 = -t94 * t52 + t96 * t55;
t138 = qJ(3) * (t92 * t38 + t158) - pkin(2) * t37;
t137 = t96 * qJ(5);
t128 = t92 * t143;
t28 = t157 * t90 - t81;
t118 = t90 * t28 + t92 * t29;
t3 = -t96 * t8 + t94 * t9;
t102 = t104 + t159;
t83 = t85 * qJDD(2);
t82 = t84 * qJDD(2);
t69 = 0.2e1 * t111;
t64 = t160 * t98;
t63 = t92 * t75;
t57 = (-t87 + t88) * t143;
t26 = (t90 * (t54 + t114) - t94 * t128) * t96;
t25 = (t96 * t128 + t158) * t94;
t23 = -t81 + (t41 + (t152 + t61) * qJD(2)) * t90;
t18 = t90 * t38 - t92 * t42;
t17 = t90 * (t96 * (-t74 + t127) + t142) + t92 * t43;
t16 = t90 * (t141 - t94 * (t74 - t126)) - t92 * t45;
t14 = t90 * t34 - t92 * t44;
t12 = t90 * (-t96 * t42 - t94 * t44) - t92 * t57;
t11 = t90 * t21 + t92 * t56;
t6 = -pkin(4) * t127 + t106;
t4 = t94 * t8 + t96 * t9;
t2 = -t94 * t5 + t96 * t6;
t1 = t96 * t5 + t94 * t6;
t10 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t133, 0, 0, 0, 0, 0, 0, 0, 0, 0, t133, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t92 * t28 + t90 * t29, 0, 0, 0, 0, 0, 0, t18, t14, t11, -t92 * t23 + t90 * t4, 0, 0, 0, 0, 0, 0, t18, t14, t11, t90 * t2 - t92 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t116, t108, 0, 0, t82, 0.2e1 * t90 * t130, 0, t83, 0, 0, -t153 * t92, t153 * t90, pkin(2) * t64 + qJ(3) * (t83 + t82) + t118, -pkin(2) * t40 + qJ(3) * t118, t26, t12, t16, t25, t17, t63, t90 * (t94 * t23 - t145) + t92 * (t8 - t149) + t138, t90 * (t96 * t23 - t146) + t92 * (t9 - t150) + t139, t90 * (-t3 - t147) - t20 * t144 + t140, qJ(3) * (t90 * t23 + t92 * t4) + (-pkin(2) + t109) * t3, t26, t12, t16, t25, t17, t63, t90 * (-t52 * t137 - t145 + (pkin(4) * t42 - qJ(5) * t55 + t113 + (-t117 * qJD(2) + t157) * t90) * t94) + t92 * (-t102 + t69 - t149) + t138, t90 * (-t94 * (-pkin(4) * t44 + qJ(5) * t51) - t146 + (-qJ(5) * t46 + t7) * t96) + t92 * (-t150 - t156) + t139, t90 * (t96 * (qJ(5) * t45 - t104 + t69) - t94 * (-qJ(5) * t43 + (t56 - t127) * pkin(4) + t106) - t147) + t92 * (-pkin(3) * t20 + t148) + t140, t90 * (-t5 * t137 - t94 * (-pkin(4) * t7 + qJ(5) * t6) - pkin(6) * t1) + t92 * (-pkin(3) * t1 - t151) - pkin(2) * t1 + qJ(3) * (t92 * t2 + t90 * t7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t130, t131, -t64, t40, 0, 0, 0, 0, 0, 0, t37, t33, t20, t3, 0, 0, 0, 0, 0, 0, t37, t33, t20, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, t57, t45, -t65, -t43, -t75, -t8, -t9, 0, 0, t65, t57, t45, -t65, -t43, -t75, t102 + t70, t156, -t148, t151; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, t44, -t56, t7;];
tauJ_reg = t10;
