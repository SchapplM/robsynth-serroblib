% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PRPRR4
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PRPRR4_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR4_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR4_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR4_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR4_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR4_invdynJ_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:51:36
% EndTime: 2019-12-05 15:51:40
% DurationCPUTime: 1.17s
% Computational Cost: add. (3322->204), mult. (6218->301), div. (0->0), fcn. (4581->12), ass. (0->136)
t107 = -g(3) + qJDD(1);
t112 = sin(pkin(5));
t115 = cos(pkin(5));
t111 = sin(pkin(9));
t114 = cos(pkin(9));
t87 = t111 * g(1) - t114 * g(2);
t155 = t115 * t87;
t159 = t107 * t112 + t155;
t125 = qJD(2) ^ 2;
t118 = sin(qJ(5));
t121 = cos(qJ(5));
t119 = sin(qJ(4));
t143 = qJD(2) * t119;
t75 = -t121 * qJD(4) + t118 * t143;
t77 = t118 * qJD(4) + t121 * t143;
t156 = t77 * t75;
t122 = cos(qJ(4));
t101 = t122 * qJDD(2);
t142 = qJD(2) * qJD(4);
t99 = t119 * t142;
t81 = t101 - t99;
t74 = -qJDD(5) + t81;
t129 = -t74 - t156;
t158 = t118 * t129;
t157 = t121 * t129;
t139 = t122 * t142;
t141 = t119 * qJDD(2);
t80 = t139 + t141;
t136 = -t121 * qJDD(4) + t118 * t80;
t96 = t122 * qJD(2) - qJD(5);
t39 = (qJD(5) + t96) * t77 + t136;
t72 = t75 ^ 2;
t73 = t77 ^ 2;
t94 = t96 ^ 2;
t110 = sin(pkin(10));
t113 = cos(pkin(10));
t120 = sin(qJ(2));
t123 = cos(qJ(2));
t88 = -t114 * g(1) - t111 * g(2);
t55 = -t120 * t88 + t159 * t123;
t127 = qJDD(2) * pkin(2) + t55;
t56 = t159 * t120 + t123 * t88;
t50 = -t125 * pkin(2) + t56;
t34 = t110 * t127 + t113 * t50;
t124 = qJD(4) ^ 2;
t134 = -t122 * pkin(4) - t119 * pkin(8);
t32 = -t125 * pkin(3) + qJDD(2) * pkin(7) + t34;
t137 = t125 * t134 + t32;
t135 = t115 * t107 + qJDD(3);
t128 = -t112 * t87 + t135;
t63 = t122 * t128;
t18 = -qJDD(4) * pkin(4) - t124 * pkin(8) + t137 * t119 - t63;
t154 = t118 * t18;
t52 = t74 - t156;
t153 = t118 * t52;
t152 = t118 * t96;
t95 = t119 * t125 * t122;
t89 = qJDD(4) + t95;
t151 = t119 * t89;
t150 = t121 * t18;
t149 = t121 * t52;
t148 = t121 * t96;
t90 = qJDD(4) - t95;
t147 = t122 * t90;
t145 = qJD(5) - t96;
t140 = t122 * t156;
t126 = t119 * t128;
t19 = -t124 * pkin(4) + qJDD(4) * pkin(8) + t137 * t122 + t126;
t132 = t80 + t139;
t133 = -t81 + t99;
t138 = t110 * t50 - t113 * t127;
t31 = -qJDD(2) * pkin(3) - t125 * pkin(7) + t138;
t23 = t133 * pkin(4) - t132 * pkin(8) + t31;
t8 = t118 * t19 - t121 * t23;
t9 = t118 * t23 + t121 * t19;
t5 = t118 * t8 + t121 * t9;
t27 = t119 * t32 - t63;
t28 = t122 * t32 + t126;
t12 = t119 * t27 + t122 * t28;
t4 = t118 * t9 - t121 * t8;
t131 = -t118 * qJDD(4) - t121 * t80;
t130 = -pkin(3) + t134;
t58 = -t75 * qJD(5) - t131;
t106 = t122 ^ 2;
t105 = t119 ^ 2;
t104 = t106 * t125;
t102 = t105 * t125;
t93 = -t104 - t124;
t92 = -t102 - t124;
t86 = t102 + t104;
t85 = (t105 + t106) * qJDD(2);
t84 = -t110 * qJDD(2) - t113 * t125;
t83 = t113 * qJDD(2) - t110 * t125;
t82 = t101 - 0.2e1 * t99;
t79 = 0.2e1 * t139 + t141;
t70 = t75 * t96;
t69 = -t73 + t94;
t68 = t72 - t94;
t67 = -t119 * t92 - t147;
t66 = t122 * t93 - t151;
t65 = -t119 * t90 + t122 * t92;
t64 = t119 * t93 + t122 * t89;
t62 = t73 - t72;
t61 = -t73 - t94;
t60 = t110 * t85 + t113 * t86;
t59 = -t94 - t72;
t57 = -t77 * qJD(5) - t136;
t51 = t72 + t73;
t49 = t110 * t67 - t113 * t79;
t48 = t110 * t66 + t113 * t82;
t44 = t145 * t75 + t131;
t43 = t58 - t70;
t42 = t58 + t70;
t40 = -t145 * t77 - t136;
t38 = -t118 * t61 + t149;
t37 = t121 * t61 + t153;
t36 = t121 * t59 - t158;
t35 = t118 * t59 + t157;
t30 = t118 * t43 - t121 * t39;
t29 = -t118 * t39 - t121 * t43;
t25 = -t119 * t44 + t122 * t38;
t24 = t119 * t38 + t122 * t44;
t21 = -t119 * t40 + t122 * t36;
t20 = t119 * t36 + t122 * t40;
t17 = -t119 * t51 + t122 * t30;
t16 = t119 * t30 + t122 * t51;
t15 = t110 * t34 - t113 * t138;
t14 = t110 * t25 - t113 * t37;
t13 = t110 * t21 - t113 * t35;
t11 = t119 * t28 - t122 * t27;
t10 = t110 * t17 - t113 * t29;
t6 = t110 * t12 - t113 * t31;
t3 = t119 * t18 + t122 * t5;
t2 = t119 * t5 - t122 * t18;
t1 = t110 * t3 - t113 * t4;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t107, 0, 0, 0, 0, 0, 0, (qJDD(2) * t123 - t120 * t125) * t112, (-qJDD(2) * t120 - t123 * t125) * t112, 0, t115 ^ 2 * t107 + (t120 * t56 + t123 * t55 - t155) * t112, 0, 0, 0, 0, 0, 0, (t120 * t84 + t123 * t83) * t112, (-t120 * t83 + t123 * t84) * t112, 0, t115 * t135 + (t120 * (t110 * t138 + t113 * t34) + t123 * t15 - t155) * t112, 0, 0, 0, 0, 0, 0, t115 * t64 + (t120 * (-t110 * t82 + t113 * t66) + t123 * t48) * t112, t115 * t65 + (t120 * (t110 * t79 + t113 * t67) + t123 * t49) * t112, (t120 * (-t110 * t86 + t113 * t85) + t123 * t60) * t112, t115 * t11 + (t120 * (t110 * t31 + t113 * t12) + t123 * t6) * t112, 0, 0, 0, 0, 0, 0, t115 * t20 + (t120 * (t110 * t35 + t113 * t21) + t123 * t13) * t112, t115 * t24 + (t120 * (t110 * t37 + t113 * t25) + t123 * t14) * t112, t115 * t16 + (t120 * (t110 * t29 + t113 * t17) + t123 * t10) * t112, t115 * t2 + (t120 * (t110 * t4 + t113 * t3) + t123 * t1) * t112; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t55, -t56, 0, 0, 0, 0, 0, 0, 0, qJDD(2), pkin(2) * t83 - t138, pkin(2) * t84 - t34, 0, pkin(2) * t15, t132 * t119, t119 * t82 + t122 * t79, t151 + t122 * (-t102 + t124), -t133 * t122, t119 * (t104 - t124) + t147, 0, pkin(2) * t48 + pkin(3) * t82 + pkin(7) * t66 - t122 * t31, pkin(2) * t49 - pkin(3) * t79 + pkin(7) * t67 + t119 * t31, pkin(2) * t60 + pkin(3) * t86 + pkin(7) * t85 + t12, pkin(2) * t6 - pkin(3) * t31 + pkin(7) * t12, t119 * (t121 * t58 + t77 * t152) - t140, t119 * (-t118 * t42 + t121 * t40) - t122 * t62, t119 * (-t118 * t69 + t157) - t122 * t43, t119 * (-t118 * t57 - t75 * t148) + t140, t119 * (t121 * t68 + t153) + t122 * t39, t122 * t74 + t119 * (-t118 * t77 + t121 * t75) * t96, t119 * (-pkin(8) * t35 + t154) + t122 * (-pkin(4) * t35 + t8) - pkin(3) * t35 + pkin(7) * t21 + pkin(2) * t13, t119 * (-pkin(8) * t37 + t150) + t122 * (-pkin(4) * t37 + t9) - pkin(3) * t37 + pkin(7) * t25 + pkin(2) * t14, pkin(2) * t10 + pkin(7) * t17 - t119 * t4 + t130 * t29, pkin(2) * t1 + pkin(7) * t3 + t130 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t128, 0, 0, 0, 0, 0, 0, t64, t65, 0, t11, 0, 0, 0, 0, 0, 0, t20, t24, t16, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t95, t102 - t104, t141, t95, t101, qJDD(4), -t27, -t28, 0, 0, t118 * t58 - t77 * t148, t118 * t40 + t121 * t42, t121 * t69 + t158, t121 * t57 - t75 * t152, t118 * t68 - t149, (t118 * t75 + t121 * t77) * t96, pkin(4) * t40 + pkin(8) * t36 - t150, pkin(4) * t44 + pkin(8) * t38 + t154, pkin(4) * t51 + pkin(8) * t30 + t5, -pkin(4) * t18 + pkin(8) * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t156, t62, t43, -t156, -t39, -t74, -t8, -t9, 0, 0;];
tauJ_reg = t7;
