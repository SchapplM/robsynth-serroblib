% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 13:32
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPPPRR1_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR1_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR1_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPPRR1_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR1_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR1_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 13:32:24
% EndTime: 2019-05-05 13:32:27
% DurationCPUTime: 0.92s
% Computational Cost: add. (2578->182), mult. (4701->219), div. (0->0), fcn. (2556->8), ass. (0->133)
t98 = cos(qJ(5));
t137 = qJD(1) * t98;
t94 = sin(qJ(6));
t97 = cos(qJ(6));
t57 = -t97 * qJD(5) + t94 * t137;
t59 = t94 * qJD(5) + t97 * t137;
t44 = t59 * t57;
t132 = qJD(1) * qJD(5);
t125 = t98 * t132;
t95 = sin(qJ(5));
t134 = t95 * qJDD(1);
t65 = -t125 - t134;
t55 = qJDD(6) - t65;
t154 = -t44 + t55;
t156 = t154 * t94;
t155 = t154 * t97;
t96 = sin(qJ(1));
t99 = cos(qJ(1));
t116 = t96 * g(1) - t99 * g(2);
t107 = qJDD(1) * pkin(1) + t116;
t101 = qJD(1) ^ 2;
t117 = t99 * g(1) + t96 * g(2);
t62 = -t101 * pkin(1) - t117;
t91 = sin(pkin(9));
t92 = cos(pkin(9));
t123 = t92 * t107 - t91 * t62;
t114 = qJDD(3) - t123;
t61 = pkin(1) * (-t92 * qJDD(1) + t91 * t101);
t90 = qJDD(1) * pkin(2);
t153 = t114 - 0.2e1 * t90 + t61;
t142 = -pkin(2) - qJ(4);
t152 = t142 * t101 + qJDD(4);
t126 = t95 * t132;
t133 = t98 * qJDD(1);
t66 = -t126 + t133;
t122 = -t97 * qJDD(5) + t94 * t66;
t77 = t95 * qJD(1) + qJD(6);
t23 = (qJD(6) - t77) * t59 + t122;
t53 = t57 ^ 2;
t54 = t59 ^ 2;
t76 = t77 ^ 2;
t151 = 2 * qJD(3);
t60 = pkin(1) * (t91 * qJDD(1) + t92 * t101);
t150 = t101 * pkin(7);
t149 = t77 * t94;
t148 = t77 * t97;
t37 = t44 + t55;
t147 = t94 * t37;
t128 = t95 * t101 * t98;
t72 = qJDD(5) + t128;
t146 = t95 * t72;
t145 = t97 * t37;
t100 = qJD(5) ^ 2;
t141 = t91 * t107 + t92 * t62;
t103 = (-pkin(7) + qJ(3)) * qJDD(1) + t141 + t152;
t127 = qJD(1) * t151;
t102 = t103 + t127;
t88 = -g(3) + qJDD(2);
t21 = -t98 * t102 + t95 * t88;
t118 = t95 * pkin(5) - t98 * pkin(8);
t63 = t118 * qJD(1);
t15 = -qJDD(5) * pkin(5) - t100 * pkin(8) + t63 * t137 + t21;
t144 = t98 * t15;
t143 = t98 * t88;
t86 = t95 ^ 2;
t87 = t98 ^ 2;
t140 = t86 + t87;
t139 = t86 * t101;
t138 = t87 * t101;
t135 = qJD(6) + t77;
t131 = qJD(4) * qJD(1);
t130 = qJDD(1) * qJ(3);
t129 = t95 * t44;
t34 = -t101 * qJ(3) + t114 - t90;
t84 = qJDD(1) * qJ(4);
t106 = -t34 + t84;
t112 = -t66 + t126;
t113 = -t65 + t125;
t83 = 0.2e1 * t131;
t14 = t113 * pkin(5) + t112 * pkin(8) + t106 - t150 + t83;
t16 = -t100 * pkin(5) + qJDD(5) * pkin(8) + t143 + ((t151 - t63) * qJD(1) + t103) * t95;
t5 = -t97 * t14 + t94 * t16;
t6 = t94 * t14 + t97 * t16;
t3 = t94 * t5 + t97 * t6;
t124 = pkin(1) * t91 + qJ(3);
t120 = pkin(1) * t92 - t142;
t119 = -pkin(7) + t124;
t115 = t97 * t5 - t94 * t6;
t22 = t95 * t102 + t143;
t10 = -t98 * t21 + t95 * t22;
t111 = t127 + t141;
t110 = -t94 * qJDD(5) - t97 * t66;
t108 = t111 + t60 + 0.2e1 * t130;
t105 = t111 + t130;
t40 = -t57 * qJD(6) - t110;
t104 = t118 + t120;
t31 = -t106 - 0.2e1 * t131;
t75 = -t100 - t138;
t74 = -t100 - t139;
t73 = qJDD(5) - t128;
t71 = t140 * t101;
t70 = t140 * qJDD(1);
t67 = -0.2e1 * t126 + t133;
t64 = 0.2e1 * t125 + t134;
t56 = t98 * t73;
t49 = t77 * t57;
t48 = -t54 + t76;
t47 = t53 - t76;
t46 = t98 * t75 - t146;
t45 = t95 * t74 + t56;
t43 = t54 - t53;
t42 = -t54 - t76;
t41 = -t76 - t53;
t39 = -t59 * qJD(6) - t122;
t35 = t53 + t54;
t33 = -t101 * pkin(2) + t105;
t32 = t105 + t152;
t30 = t31 + t150;
t28 = t135 * t57 + t110;
t27 = t40 + t49;
t26 = t40 - t49;
t24 = -t135 * t59 - t122;
t20 = -t94 * t42 - t145;
t19 = t97 * t42 - t147;
t18 = t97 * t41 - t156;
t17 = t94 * t41 + t155;
t12 = -t23 * t97 + t94 * t27;
t11 = -t23 * t94 - t97 * t27;
t9 = t95 * t20 + t98 * t28;
t8 = t95 * t18 + t98 * t24;
t7 = t95 * t12 + t98 * t35;
t1 = t95 * t3 - t144;
t2 = [0, 0, 0, 0, 0, qJDD(1), t116, t117, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t123 - t61, -t141 - t60, 0, pkin(1) * (t92 * t123 + t91 * t141), qJDD(1), 0, 0, 0, 0, 0, 0, t153, t108, pkin(1) * (t91 * t33 - t92 * t34) - pkin(2) * t34 + qJ(3) * t33, qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(4) + t108, t83 + 0.2e1 * t84 - t153, -t120 * t31 + t124 * t32, -t112 * t98, -t98 * t64 - t95 * t67, t56 - t95 * (t100 - t138), t113 * t95, t98 * (-t100 + t139) - t146, 0, t119 * t45 + t120 * t64 - t95 * t30, t119 * t46 + t120 * t67 - t98 * t30, -t119 * t70 - t120 * t71 - t10, t119 * t10 - t120 * t30, t98 * (-t59 * t149 + t97 * t40) + t129, t98 * (t97 * t24 - t94 * t26) + t95 * t43, t98 * (-t94 * t48 + t155) + t95 * t27, t98 * (t57 * t148 - t94 * t39) - t129, t98 * (t97 * t47 - t147) - t95 * t23, t95 * t55 + t98 * (-t57 * t97 + t59 * t94) * t77, t104 * t17 + t119 * t8 + t94 * t144 - t95 * t5, t104 * t19 + t119 * t9 + t97 * t144 - t95 * t6, t104 * t11 + t98 * t115 + t119 * t7, t119 * t1 - t104 * t115; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88, 0, 0, 0, 0, 0, 0, -t95 * t73 + t98 * t74, -t98 * t72 - t95 * t75, 0, t95 * t21 + t98 * t22, 0, 0, 0, 0, 0, 0, t98 * t18 - t95 * t24, t98 * t20 - t95 * t28, t98 * t12 - t95 * t35, t95 * t15 + t98 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t101, t34, 0, 0, 0, 0, 0, 0, 0, -t101, -qJDD(1), t31, 0, 0, 0, 0, 0, 0, -t64, -t67, t71, t30, 0, 0, 0, 0, 0, 0, -t17, -t19, -t11, t115; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t101, t32, 0, 0, 0, 0, 0, 0, t45, t46, -t70, t10, 0, 0, 0, 0, 0, 0, t8, t9, t7, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t128, (-t86 + t87) * t101, t133, -t128, -t134, qJDD(5), -t21, -t22, 0, 0, t59 * t148 + t94 * t40, t94 * t24 + t97 * t26, t97 * t48 + t156, t57 * t149 + t97 * t39, t94 * t47 + t145, (-t57 * t94 - t59 * t97) * t77, pkin(5) * t24 + pkin(8) * t18 - t97 * t15, pkin(5) * t28 + pkin(8) * t20 + t94 * t15, pkin(5) * t35 + pkin(8) * t12 + t3, -pkin(5) * t15 + pkin(8) * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, t43, t27, -t44, -t23, t55, -t5, -t6, 0, 0;];
tauJ_reg  = t2;
