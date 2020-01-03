% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S4RPRR7
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% tauJ_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S4RPRR7_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR7_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR7_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR7_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR7_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR7_invdynJ_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:54:05
% EndTime: 2019-12-31 16:54:09
% DurationCPUTime: 1.65s
% Computational Cost: add. (3641->239), mult. (8793->338), div. (0->0), fcn. (6236->8), ass. (0->148)
t119 = sin(qJ(1));
t173 = cos(qJ(1));
t129 = t173 * g(1) + t119 * g(2);
t185 = (2 * qJD(2) * qJD(1)) - t129;
t117 = sin(qJ(4));
t114 = sin(pkin(7));
t121 = cos(qJ(3));
t115 = cos(pkin(7));
t118 = sin(qJ(3));
t154 = t115 * t118;
t131 = t114 * t121 + t154;
t101 = t131 * qJD(1);
t120 = cos(qJ(4));
t85 = -t120 * qJD(3) + t101 * t117;
t87 = qJD(3) * t117 + t101 * t120;
t63 = t87 * t85;
t147 = t115 * qJDD(1);
t148 = t114 * qJDD(1);
t132 = t118 * t148 - t121 * t147;
t150 = t101 * qJD(3);
t77 = -t132 - t150;
t69 = qJDD(4) - t77;
t178 = -t63 + t69;
t184 = t117 * t178;
t155 = t114 * t118;
t99 = (-t115 * t121 + t155) * qJD(1);
t81 = t101 * t99;
t177 = qJDD(3) - t81;
t183 = t118 * t177;
t182 = t120 * t178;
t181 = t121 * t177;
t123 = qJD(1) ^ 2;
t180 = -(t123 * pkin(1)) + qJDD(1) * qJ(2) + t185;
t139 = g(1) * t119 - t173 * g(2);
t133 = -qJDD(2) + t139;
t141 = pkin(2) * t115 + pkin(1);
t111 = t114 ^ 2;
t112 = t115 ^ 2;
t151 = t111 + t112;
t72 = t141 * qJDD(1) + (t151 * pkin(5) + qJ(2)) * t123 + t133;
t179 = pkin(5) + qJ(2);
t156 = t112 * t123;
t157 = qJDD(1) * pkin(1);
t159 = qJ(2) * t123;
t95 = t133 + t157 + t159;
t176 = qJ(2) * t156 + t111 * t159 - t157 - t95;
t158 = qJD(3) * t99;
t98 = t131 * qJDD(1);
t79 = t98 - t158;
t136 = -t120 * qJDD(3) + t117 * t79;
t94 = qJD(4) + t99;
t36 = (qJD(4) - t94) * t87 + t136;
t83 = t85 ^ 2;
t84 = t87 ^ 2;
t93 = t94 ^ 2;
t96 = t99 ^ 2;
t97 = t101 ^ 2;
t172 = t115 * g(3);
t127 = (-t179 * qJDD(1) + t141 * t123 - t185) * t114;
t126 = t127 - t172;
t135 = -g(3) * t114 + t180 * t115;
t67 = -pkin(2) * t156 + pkin(5) * t147 + t135;
t43 = t118 * t67 - t121 * t126;
t162 = t121 * t67;
t44 = -g(3) * t154 + t118 * t127 + t162;
t22 = t118 * t44 - t121 * t43;
t171 = t114 * t22;
t122 = qJD(3) ^ 2;
t70 = pkin(3) * t99 - pkin(6) * t101;
t28 = -qJDD(3) * pkin(3) - t122 * pkin(6) + t101 * t70 + t43;
t170 = t117 * t28;
t47 = t63 + t69;
t169 = t117 * t47;
t168 = t117 * t94;
t167 = t118 * t72;
t74 = qJDD(3) + t81;
t166 = t118 * t74;
t165 = t120 * t28;
t164 = t120 * t47;
t163 = t120 * t94;
t161 = t121 * t72;
t160 = t121 * t74;
t152 = qJD(4) + t94;
t146 = t118 * t63;
t145 = t121 * t63;
t140 = -pkin(3) * t121 - pkin(2);
t138 = t114 * (t180 * t114 + t172) + t115 * t135;
t29 = -t122 * pkin(3) + qJDD(3) * pkin(6) + t118 * t126 - t99 * t70 + t162;
t31 = (-t79 + t158) * pkin(6) + (-t77 + t150) * pkin(3) - t72;
t13 = t117 * t29 - t120 * t31;
t14 = t117 * t31 + t120 * t29;
t5 = t117 * t13 + t120 * t14;
t23 = t118 * t43 + t121 * t44;
t4 = t117 * t14 - t120 * t13;
t130 = -qJDD(3) * t117 - t120 * t79;
t55 = -qJD(4) * t85 - t130;
t107 = t112 * qJDD(1);
t106 = t111 * qJDD(1);
t102 = t151 * t123;
t90 = -t97 - t122;
t89 = -t97 + t122;
t88 = t96 - t122;
t78 = t98 - 0.2e1 * t158;
t76 = t132 + 0.2e1 * t150;
t71 = -t122 - t96;
t68 = t94 * t85;
t66 = -t84 + t93;
t65 = t83 - t93;
t62 = -t96 - t97;
t61 = t84 - t83;
t59 = -t84 - t93;
t58 = -t118 * t90 - t160;
t57 = t121 * t90 - t166;
t56 = -t93 - t83;
t54 = -qJD(4) * t87 - t136;
t53 = t83 + t84;
t52 = t118 * t98 - t121 * t132;
t51 = -t118 * t132 - t121 * t98;
t50 = t121 * t71 - t183;
t49 = t118 * t71 + t181;
t45 = (t117 * t87 - t120 * t85) * t94;
t41 = t152 * t85 + t130;
t40 = t55 + t68;
t39 = t55 - t68;
t37 = -t152 * t87 - t136;
t35 = t120 * t55 - t87 * t168;
t34 = -t117 * t54 + t85 * t163;
t33 = t120 * t65 - t169;
t32 = -t117 * t66 + t182;
t27 = -t117 * t59 - t164;
t26 = t120 * t59 - t169;
t25 = t120 * t56 - t184;
t24 = t117 * t56 + t182;
t21 = t117 * t40 - t120 * t36;
t20 = -t117 * t39 + t120 * t37;
t19 = -t117 * t36 - t120 * t40;
t18 = -t118 * t41 + t121 * t27;
t17 = t118 * t27 + t121 * t41;
t16 = -t118 * t37 + t121 * t25;
t15 = t118 * t25 + t121 * t37;
t11 = -pkin(6) * t26 + t165;
t10 = -t118 * t53 + t121 * t21;
t9 = t118 * t21 + t121 * t53;
t8 = -pkin(6) * t24 + t170;
t7 = -pkin(3) * t26 + t14;
t6 = -pkin(3) * t24 + t13;
t1 = -pkin(6) * t19 - t4;
t2 = [0, 0, 0, 0, 0, qJDD(1), t139, t129, 0, 0, t106, 0.2e1 * t114 * t147, 0, t107, 0, 0, -t176 * t115, t176 * t114, pkin(1) * t102 + qJ(2) * (t107 + t106) + t138, pkin(1) * t95 + qJ(2) * t138, t114 * (-t118 * t150 + t121 * t79) + t115 * (t118 * t79 + t121 * t150), t114 * (-t118 * t78 - t121 * t76) + t115 * (-t118 * t76 + t121 * t78), t114 * (-t118 * t89 + t181) + t115 * (t121 * t89 + t183), t114 * (-t118 * t77 + t121 * t158) + t115 * (t118 * t158 + t121 * t77), t114 * (t121 * t88 - t166) + t115 * (t118 * t88 + t160), (t114 * (t101 * t118 - t121 * t99) + t115 * (-t101 * t121 - t118 * t99)) * qJD(3), t114 * (-pkin(5) * t49 - t167) + t115 * (-pkin(2) * t76 + pkin(5) * t50 + t161) - pkin(1) * t76 + qJ(2) * (-t114 * t49 + t115 * t50), t114 * (-pkin(5) * t57 - t161) + t115 * (-pkin(2) * t78 + pkin(5) * t58 - t167) - pkin(1) * t78 + qJ(2) * (-t114 * t57 + t115 * t58), t114 * (-pkin(5) * t51 - t22) + t115 * (-pkin(2) * t62 + pkin(5) * t52 + t23) - pkin(1) * t62 + qJ(2) * (-t114 * t51 + t115 * t52), -pkin(5) * t171 + t115 * (pkin(2) * t72 + pkin(5) * t23) + pkin(1) * t72 + qJ(2) * (t115 * t23 - t171), t114 * (t121 * t35 + t146) + t115 * (t118 * t35 - t145), t114 * (t118 * t61 + t121 * t20) + t115 * (t118 * t20 - t121 * t61), t114 * (t118 * t40 + t121 * t32) + t115 * (t118 * t32 - t121 * t40), t114 * (t121 * t34 - t146) + t115 * (t118 * t34 + t145), t114 * (-t118 * t36 + t121 * t33) + t115 * (t118 * t33 + t121 * t36), t114 * (t118 * t69 + t121 * t45) + t115 * (t118 * t45 - t121 * t69), t114 * (-pkin(5) * t15 - t118 * t6 + t121 * t8) + t115 * (-pkin(2) * t24 + pkin(5) * t16 + t118 * t8 + t121 * t6) - pkin(1) * t24 + qJ(2) * (-t114 * t15 + t115 * t16), t114 * (-pkin(5) * t17 + t11 * t121 - t118 * t7) + t115 * (-pkin(2) * t26 + pkin(5) * t18 + t11 * t118 + t121 * t7) - pkin(1) * t26 + qJ(2) * (-t114 * t17 + t115 * t18), t114 * (-pkin(5) * t9 + t1 * t121) + t115 * (pkin(5) * t10 + t1 * t118) + qJ(2) * (t10 * t115 - t114 * t9) + (pkin(3) * t155 + t115 * t140 - pkin(1)) * t19, (t114 * (pkin(3) * t118 - pkin(6) * t121) + t115 * (-pkin(6) * t118 + t140) - pkin(1)) * t4 + t179 * (-t114 * (t118 * t5 - t121 * t28) + t115 * (t118 * t28 + t121 * t5)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t147, t148, -t102, -t95, 0, 0, 0, 0, 0, 0, t76, t78, t62, -t72, 0, 0, 0, 0, 0, 0, t24, t26, t19, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, t97 - t96, t98, -t81, -t132, qJDD(3), -t43, -t44, 0, 0, t117 * t55 + t87 * t163, t117 * t37 + t120 * t39, t120 * t66 + t184, t120 * t54 + t85 * t168, t117 * t65 + t164, (-t117 * t85 - t120 * t87) * t94, pkin(3) * t37 + pkin(6) * t25 - t165, pkin(3) * t41 + pkin(6) * t27 + t170, pkin(3) * t53 + pkin(6) * t21 + t5, -pkin(3) * t28 + pkin(6) * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, t61, t40, -t63, -t36, t69, -t13, -t14, 0, 0;];
tauJ_reg = t2;
