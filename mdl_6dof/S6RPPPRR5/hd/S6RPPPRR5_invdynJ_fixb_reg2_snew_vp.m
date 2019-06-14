% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPPPRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta4]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 13:51
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPPPRR5_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR5_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR5_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPPRR5_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR5_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR5_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 13:50:39
% EndTime: 2019-05-05 13:50:42
% DurationCPUTime: 1.19s
% Computational Cost: add. (4418->212), mult. (7655->268), div. (0->0), fcn. (3566->8), ass. (0->144)
t120 = sin(qJ(1));
t123 = cos(qJ(1));
t142 = g(1) * t120 - t123 * g(2);
t135 = -qJDD(2) + t142;
t164 = 2 * qJD(1);
t127 = (qJD(3) * t164) + t135;
t162 = (pkin(1) + qJ(3));
t138 = t162 * qJDD(1);
t167 = t127 + t138;
t125 = qJD(1) ^ 2;
t118 = sin(qJ(6));
t121 = cos(qJ(6));
t119 = sin(qJ(5));
t147 = qJD(1) * t119;
t76 = -t121 * qJD(5) + t118 * t147;
t78 = qJD(5) * t118 + t121 * t147;
t163 = t78 * t76;
t122 = cos(qJ(5));
t101 = t122 * qJDD(1);
t146 = qJD(1) * qJD(5);
t98 = t119 * t146;
t82 = t101 - t98;
t75 = -qJDD(6) + t82;
t130 = -t75 - t163;
t166 = t118 * t130;
t165 = t121 * t130;
t143 = t122 * t146;
t145 = t119 * qJDD(1);
t81 = t143 + t145;
t140 = -t121 * qJDD(5) + t118 * t81;
t94 = qJD(1) * t122 - qJD(6);
t40 = (qJD(6) + t94) * t78 + t140;
t73 = t76 ^ 2;
t74 = t78 ^ 2;
t92 = t94 ^ 2;
t161 = (qJ(2) + pkin(3));
t114 = sin(pkin(9));
t115 = cos(pkin(9));
t109 = qJDD(1) * qJ(2);
t136 = g(1) * t123 + g(2) * t120;
t129 = (qJD(2) * t164) - t136;
t128 = qJDD(3) + t129;
t68 = -(t162 * t125) + t109 + t128;
t126 = qJDD(1) * pkin(3) + t68;
t64 = -(t161 * t125) - t167;
t39 = t114 * t126 + t115 * t64;
t124 = qJD(5) ^ 2;
t137 = -pkin(5) * t122 - pkin(8) * t119;
t37 = -(pkin(4) * t125) + qJDD(1) * pkin(7) + t39;
t141 = t125 * t137 + t37;
t150 = -g(3) + qJDD(4);
t97 = t122 * t150;
t27 = -qJDD(5) * pkin(5) - t124 * pkin(8) + t141 * t119 - t97;
t160 = t118 * t27;
t51 = t75 - t163;
t159 = t118 * t51;
t158 = t118 * t94;
t93 = t122 * t125 * t119;
t88 = qJDD(5) + t93;
t157 = t119 * t88;
t156 = t121 * t27;
t155 = t121 * t51;
t154 = t121 * t94;
t89 = qJDD(5) - t93;
t153 = t122 * t89;
t152 = qJ(2) * t125;
t151 = qJDD(1) * pkin(1);
t149 = qJD(6) - t94;
t144 = t122 * t163;
t38 = -t114 * t64 + t115 * t126;
t133 = t81 + t143;
t134 = -t82 + t98;
t36 = -qJDD(1) * pkin(4) - t125 * pkin(7) - t38;
t26 = t134 * pkin(5) - t133 * pkin(8) + t36;
t139 = t119 * t150;
t28 = -t124 * pkin(5) + qJDD(5) * pkin(8) + t141 * t122 + t139;
t15 = t118 * t28 - t121 * t26;
t16 = t118 * t26 + t121 * t28;
t5 = t118 * t15 + t121 * t16;
t32 = t119 * t37 - t97;
t33 = t122 * t37 + t139;
t18 = t119 * t32 + t122 * t33;
t4 = t118 * t16 - t121 * t15;
t132 = -qJDD(5) * t118 - t121 * t81;
t131 = -pkin(4) + t137;
t55 = -qJD(6) * t76 - t132;
t111 = t122 ^ 2;
t110 = t119 ^ 2;
t105 = 0.2e1 * t109;
t104 = t111 * t125;
t102 = t110 * t125;
t91 = -t104 - t124;
t90 = -t102 - t124;
t87 = t102 + t104;
t86 = (t110 + t111) * qJDD(1);
t85 = qJDD(1) * t114 + t115 * t125;
t84 = qJDD(1) * t115 - t114 * t125;
t83 = t101 - 0.2e1 * t98;
t80 = 0.2e1 * t143 + t145;
t72 = t135 + t151 + t152;
t71 = t76 * t94;
t70 = -t74 + t92;
t69 = t73 - t92;
t67 = t152 + t167;
t66 = -t119 * t90 - t153;
t65 = t122 * t91 - t157;
t63 = t74 - t73;
t62 = -t74 - t92;
t61 = -t114 * t87 + t115 * t86;
t60 = t114 * t86 + t115 * t87;
t56 = -t92 - t73;
t54 = -qJD(6) * t78 - t140;
t50 = t73 + t74;
t49 = t114 * t80 + t115 * t66;
t48 = -t114 * t83 + t115 * t65;
t47 = t114 * t66 - t115 * t80;
t46 = t114 * t65 + t115 * t83;
t45 = t149 * t76 + t132;
t44 = t55 - t71;
t43 = t55 + t71;
t41 = -t149 * t78 - t140;
t35 = -t118 * t62 + t155;
t34 = t121 * t62 + t159;
t31 = t121 * t56 - t166;
t30 = t118 * t56 + t165;
t24 = t118 * t44 - t121 * t40;
t23 = -t118 * t40 - t121 * t44;
t22 = -t114 * t38 + t115 * t39;
t21 = t114 * t39 + t115 * t38;
t20 = -t119 * t45 + t122 * t35;
t19 = -t119 * t41 + t122 * t31;
t17 = -t119 * t50 + t122 * t24;
t13 = t114 * t34 + t115 * t20;
t12 = t114 * t20 - t115 * t34;
t11 = t114 * t30 + t115 * t19;
t10 = t114 * t19 - t115 * t30;
t9 = t114 * t36 + t115 * t18;
t8 = t114 * t18 - t115 * t36;
t7 = t114 * t23 + t115 * t17;
t6 = t114 * t17 - t115 * t23;
t3 = t119 * t27 + t122 * t5;
t2 = t114 * t4 + t115 * t3;
t1 = t114 * t3 - t115 * t4;
t14 = [0, 0, 0, 0, 0, qJDD(1), t142, t136, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, -t135 - 0.2e1 * t151, t105 + t129, pkin(1) * t72 + qJ(2) * (-pkin(1) * t125 + t109 + t129), 0, 0, 0, qJDD(1), 0, 0, t105 + t128, 0, t127 + 0.2e1 * t138, qJ(2) * t68 + t162 * t67, 0, 0, 0, 0, 0, qJDD(1), t161 * t84 + t162 * t85 + t38, -t161 * t85 + t162 * t84 - t39, 0, t161 * t21 - t162 * t22, t133 * t119, t119 * t83 + t122 * t80, t157 + t122 * (-t102 + t124), -t134 * t122, t119 * (t104 - t124) + t153, 0, pkin(4) * t83 + pkin(7) * t65 - t122 * t36 + t161 * t46 - t162 * t48, -pkin(4) * t80 + pkin(7) * t66 + t119 * t36 + t161 * t47 - t162 * t49, pkin(4) * t87 + pkin(7) * t86 + t161 * t60 - t162 * t61 + t18, -pkin(4) * t36 + pkin(7) * t18 + t161 * t8 - t162 * t9, t119 * (t121 * t55 + t78 * t158) - t144, t119 * (-t118 * t43 + t121 * t41) - t122 * t63, t119 * (-t118 * t70 + t165) - t122 * t44, t119 * (-t118 * t54 - t76 * t154) + t144, t119 * (t121 * t69 + t159) + t122 * t40, t122 * t75 + t119 * (-t118 * t78 + t121 * t76) * t94, t119 * (-pkin(8) * t30 + t160) + t122 * (-pkin(5) * t30 + t15) - pkin(4) * t30 + pkin(7) * t19 - t162 * t11 + t161 * t10, t119 * (-pkin(8) * t34 + t156) + t122 * (-pkin(5) * t34 + t16) - pkin(4) * t34 + pkin(7) * t20 - t162 * t13 + t161 * t12, pkin(7) * t17 - t119 * t4 + t131 * t23 + t161 * t6 - t162 * t7, pkin(7) * t3 + t161 * t1 + t131 * t4 - t162 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t125, -t72, 0, 0, 0, 0, 0, 0, -t125, 0, -qJDD(1), -t67, 0, 0, 0, 0, 0, 0, -t85, -t84, 0, t22, 0, 0, 0, 0, 0, 0, t48, t49, t61, t9, 0, 0, 0, 0, 0, 0, t11, t13, t7, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t125, t68, 0, 0, 0, 0, 0, 0, t84, -t85, 0, t21, 0, 0, 0, 0, 0, 0, t46, t47, t60, t8, 0, 0, 0, 0, 0, 0, t10, t12, t6, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t150, 0, 0, 0, 0, 0, 0, t119 * t91 + t122 * t88, -t119 * t89 + t122 * t90, 0, t119 * t33 - t122 * t32, 0, 0, 0, 0, 0, 0, t119 * t31 + t122 * t41, t119 * t35 + t122 * t45, t119 * t24 + t122 * t50, t119 * t5 - t122 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t93, t102 - t104, t145, t93, t101, qJDD(5), -t32, -t33, 0, 0, t118 * t55 - t78 * t154, t118 * t41 + t121 * t43, t121 * t70 + t166, t121 * t54 - t76 * t158, t118 * t69 - t155, (t118 * t76 + t121 * t78) * t94, pkin(5) * t41 + pkin(8) * t31 - t156, pkin(5) * t45 + pkin(8) * t35 + t160, pkin(5) * t50 + pkin(8) * t24 + t5, -pkin(5) * t27 + pkin(8) * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t163, t63, t44, -t163, -t40, -t75, -t15, -t16, 0, 0;];
tauJ_reg  = t14;
