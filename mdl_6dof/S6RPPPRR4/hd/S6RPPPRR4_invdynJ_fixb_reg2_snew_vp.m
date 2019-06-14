% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPPPRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 13:46
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPPPRR4_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR4_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR4_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPPRR4_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR4_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR4_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 13:46:23
% EndTime: 2019-05-05 13:46:26
% DurationCPUTime: 1.14s
% Computational Cost: add. (3916->194), mult. (6850->253), div. (0->0), fcn. (3346->8), ass. (0->138)
t103 = sin(qJ(6));
t106 = cos(qJ(6));
t107 = cos(qJ(5));
t144 = qJD(1) * t107;
t68 = qJD(5) * t106 + t103 * t144;
t69 = -t103 * qJD(5) + t106 * t144;
t159 = t68 * t69;
t143 = qJD(1) * qJD(5);
t134 = t107 * t143;
t104 = sin(qJ(5));
t140 = t104 * qJDD(1);
t72 = t134 + t140;
t67 = -qJDD(6) + t72;
t112 = -t67 - t159;
t165 = t103 * t112;
t164 = t106 * t112;
t100 = sin(pkin(9));
t101 = cos(pkin(9));
t110 = qJD(1) ^ 2;
t105 = sin(qJ(1));
t108 = cos(qJ(1));
t136 = t105 * g(1) - g(2) * t108;
t128 = qJDD(2) - t136;
t115 = -qJ(2) * t110 + t128;
t161 = pkin(1) + pkin(2);
t111 = -t161 * qJDD(1) + t115;
t126 = g(1) * t108 + g(2) * t105;
t121 = 0.2e1 * qJD(2) * qJD(1) - t126;
t94 = qJDD(1) * qJ(2);
t116 = t121 + t94;
t59 = -t161 * t110 + t116;
t133 = t100 * t59 - t101 * t111;
t75 = -t100 * qJDD(1) + t101 * t110;
t76 = qJDD(1) * t101 + t100 * t110;
t163 = -qJ(2) * t75 + t161 * t76 + t133;
t41 = t100 * t111 + t101 * t59;
t162 = qJ(2) * t76 + t161 * t75 + t41;
t135 = t104 * t143;
t139 = t107 * qJDD(1);
t74 = t135 - t139;
t132 = -t106 * qJDD(5) + t103 * t74;
t84 = qJD(1) * t104 - qJD(6);
t35 = (qJD(6) + t84) * t69 - t132;
t65 = t68 ^ 2;
t66 = t69 ^ 2;
t83 = t84 ^ 2;
t160 = pkin(3) + pkin(7);
t96 = t104 ^ 2;
t97 = t107 ^ 2;
t158 = t96 + t97;
t45 = t67 - t159;
t157 = t103 * t45;
t156 = t103 * t84;
t137 = t104 * t110 * t107;
t79 = qJDD(5) + t137;
t155 = t104 * t79;
t154 = t106 * t45;
t153 = t106 * t84;
t109 = qJD(5) ^ 2;
t127 = -pkin(5) * t104 + pkin(8) * t107;
t117 = t110 * t127;
t119 = -t110 * qJ(4) + qJDD(4) + t133;
t113 = t160 * qJDD(1) + t119;
t98 = g(3) + qJDD(3);
t28 = t104 * t98 - t107 * t113;
t22 = -qJDD(5) * pkin(5) - t109 * pkin(8) - t107 * t117 + t28;
t152 = t107 * t22;
t80 = qJDD(5) - t137;
t151 = t107 * t80;
t150 = t110 * t96;
t149 = t110 * t97;
t148 = qJDD(1) * pkin(1);
t147 = qJDD(1) * pkin(3);
t146 = qJD(6) - t84;
t142 = qJD(4) * qJD(1);
t141 = qJDD(1) * qJ(4);
t138 = t104 * t159;
t124 = t72 + t134;
t125 = -t74 - t135;
t120 = -0.2e1 * t142 - t141 + t41;
t31 = -t160 * t110 + t120;
t21 = -pkin(5) * t124 + pkin(8) * t125 + t31;
t29 = t104 * t113 + t107 * t98;
t23 = -t109 * pkin(5) + qJDD(5) * pkin(8) + t104 * t117 + t29;
t10 = t103 * t23 - t106 * t21;
t11 = t103 * t21 + t106 * t23;
t4 = t10 * t103 + t106 * t11;
t131 = qJ(2) * t101 - qJ(4);
t130 = qJ(2) * t100 + t160;
t123 = t10 * t106 - t103 * t11;
t13 = t104 * t29 - t107 * t28;
t122 = -t103 * qJDD(5) - t106 * t74;
t49 = qJD(6) * t68 - t122;
t114 = t127 + t131;
t82 = -t109 - t149;
t81 = -t109 - t150;
t78 = t158 * t110;
t77 = t158 * qJDD(1);
t73 = 0.2e1 * t135 - t139;
t71 = 0.2e1 * t134 + t140;
t63 = -t115 + t148;
t62 = t68 * t84;
t61 = -t66 + t83;
t60 = t65 - t83;
t56 = t107 * t82 - t155;
t55 = t104 * t81 + t151;
t53 = t66 - t65;
t52 = -t66 - t83;
t51 = -t100 * t78 - t101 * t77;
t50 = -t83 - t65;
t48 = qJD(6) * t69 - t132;
t44 = t65 + t66;
t43 = t100 * t73 - t101 * t56;
t42 = -t100 * t71 - t101 * t55;
t39 = -t146 * t68 + t122;
t38 = t49 + t62;
t37 = t49 - t62;
t36 = t146 * t69 - t132;
t33 = t119 + t147;
t32 = -pkin(3) * t110 + t120;
t27 = -t103 * t52 + t154;
t26 = t106 * t52 + t157;
t25 = t106 * t50 - t165;
t24 = t103 * t50 + t164;
t19 = t100 * t41 - t101 * t133;
t18 = t103 * t38 + t106 * t35;
t17 = t103 * t35 - t106 * t38;
t16 = t100 * t32 - t101 * t33;
t15 = t104 * t27 + t107 * t39;
t14 = t104 * t25 + t107 * t36;
t12 = t104 * t18 + t107 * t44;
t8 = t100 * t26 - t101 * t15;
t7 = t100 * t31 - t101 * t13;
t6 = t100 * t24 - t101 * t14;
t5 = t100 * t17 - t101 * t12;
t2 = t104 * t4 - t152;
t1 = -t100 * t123 - t101 * t2;
t3 = [0, 0, 0, 0, 0, qJDD(1), t136, t126, 0, 0, 0, 0, 0, qJDD(1), 0, 0, -t128 + 0.2e1 * t148, 0, t121 + 0.2e1 * t94, qJ(2) * (-pkin(1) * t110 + t116) + pkin(1) * t63, 0, 0, 0, 0, 0, qJDD(1), t163, t162, 0, qJ(2) * (t100 * t133 + t101 * t41) - t161 * t19, qJDD(1), 0, 0, 0, 0, 0, 0, -qJDD(4) - 0.2e1 * t147 - t163, 0.2e1 * t141 + 0.2e1 * t142 - t162, qJ(2) * (t100 * t33 + t101 * t32) + pkin(3) * t33 - qJ(4) * t32 - t161 * t16, t125 * t107, t104 * t73 - t107 * t71, -t151 + t104 * (t109 - t149), t124 * t104, -t107 * (-t109 + t150) + t155, 0, -t104 * t31 + t130 * t55 - t131 * t71 - t161 * t42, -t107 * t31 + t130 * t56 + t131 * t73 - t161 * t43, t130 * t77 - t131 * t78 - t161 * t51 + t13, t130 * t13 + t131 * t31 - t161 * t7, -t107 * (t106 * t49 - t156 * t69) - t138, -t107 * (-t103 * t37 + t106 * t36) - t104 * t53, -t107 * (-t103 * t61 + t164) - t104 * t38, -t107 * (-t103 * t48 + t153 * t68) + t138, -t107 * (t106 * t60 + t157) - t104 * t35, t104 * t67 - t107 * (t103 * t69 - t106 * t68) * t84, t104 * t10 - t103 * t152 + t114 * t24 + t130 * t14 - t161 * t6, t104 * t11 - t106 * t152 + t114 * t26 + t130 * t15 - t161 * t8, -t107 * t123 + t114 * t17 + t130 * t12 - t161 * t5, -t161 * t1 - t114 * t123 + t130 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), 0, -t110, -t63, 0, 0, 0, 0, 0, 0, -t76, -t75, 0, t19, 0, 0, 0, 0, 0, 0, 0, t76, t75, t16, 0, 0, 0, 0, 0, 0, t42, t43, t51, t7, 0, 0, 0, 0, 0, 0, t6, t8, t5, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t98, 0, 0, 0, 0, 0, 0, 0, 0, 0, t98, 0, 0, 0, 0, 0, 0, -t104 * t80 + t107 * t81, -t104 * t82 - t107 * t79, 0, t104 * t28 + t107 * t29, 0, 0, 0, 0, 0, 0, -t104 * t36 + t107 * t25, -t104 * t39 + t107 * t27, -t104 * t44 + t107 * t18, t104 * t22 + t107 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), -t110, t33, 0, 0, 0, 0, 0, 0, t55, t56, t77, t13, 0, 0, 0, 0, 0, 0, t14, t15, t12, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t137, (-t96 + t97) * t110, -t139, -t137, t140, qJDD(5), -t28, -t29, 0, 0, t103 * t49 + t153 * t69, t103 * t36 + t106 * t37, t106 * t61 + t165, t106 * t48 + t156 * t68, t103 * t60 - t154, (-t103 * t68 - t106 * t69) * t84, pkin(5) * t36 + pkin(8) * t25 - t106 * t22, pkin(5) * t39 + pkin(8) * t27 + t103 * t22, pkin(5) * t44 + pkin(8) * t18 + t4, -pkin(5) * t22 + pkin(8) * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t159, t53, t38, -t159, t35, -t67, -t10, -t11, 0, 0;];
tauJ_reg  = t3;
