% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPPPR6
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPPPR6_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR6_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR6_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR6_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR6_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR6_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:47:51
% EndTime: 2019-12-31 17:47:56
% DurationCPUTime: 1.73s
% Computational Cost: add. (3426->230), mult. (9163->325), div. (0->0), fcn. (5555->8), ass. (0->154)
t123 = qJD(1) ^ 2;
t118 = cos(pkin(7));
t116 = sin(pkin(7));
t169 = t116 * qJ(3);
t131 = -t169 + t118 * (-pkin(2) - qJ(4)) - pkin(1);
t120 = sin(qJ(1));
t122 = cos(qJ(1));
t151 = g(1) * t120 - t122 * g(2);
t141 = qJDD(2) - t151;
t163 = t118 * qJD(1);
t112 = t116 ^ 2;
t114 = t118 ^ 2;
t165 = t112 + t114;
t164 = qJD(1) * t116;
t154 = qJD(3) * t164;
t99 = -0.2e1 * t154;
t200 = t99 + (-t165 * pkin(3) - qJ(2)) * t123 + t131 * qJDD(1) + t141 - 0.2e1 * qJD(4) * t163;
t192 = 2 * qJD(2);
t142 = g(1) * t122 + g(2) * t120;
t78 = -pkin(1) * t123 + qJDD(1) * qJ(2) - t142;
t146 = qJD(1) * t192 + t78;
t119 = sin(qJ(5));
t115 = sin(pkin(8));
t121 = cos(qJ(5));
t70 = (t118 * t115 * t119 + t116 * t121) * qJD(1);
t71 = t115 * t121 * t163 - t119 * t164;
t55 = t70 * t71;
t117 = cos(pkin(8));
t161 = t118 * qJDD(1);
t152 = t117 * t161;
t91 = qJDD(5) + t152;
t195 = -t55 + t91;
t198 = t119 * t195;
t197 = t121 * t195;
t139 = -t118 * pkin(2) - t169;
t84 = t139 * qJD(1);
t137 = t78 + (t192 + t84) * qJD(1);
t190 = g(3) * t118;
t148 = qJDD(3) + t190;
t166 = t118 * t123;
t126 = (pkin(3) * qJDD(1) - qJ(4) * t166 + t137) * t116 + t148;
t20 = t200 * t115 - t117 * t126;
t196 = pkin(3) + qJ(2);
t143 = pkin(4) * t117 + pkin(6) * t115;
t136 = pkin(3) + t143;
t170 = t114 * t123;
t88 = t165 * t123;
t21 = t115 * t126 + t200 * t117;
t135 = pkin(1) - t139;
t128 = t135 * qJDD(1);
t193 = qJ(2) * t88;
t67 = t70 ^ 2;
t68 = t71 ^ 2;
t92 = t117 * t163 + qJD(5);
t90 = t92 ^ 2;
t189 = t70 * t92;
t106 = t112 * qJDD(1);
t108 = t114 * qJDD(1);
t187 = qJ(2) * (t108 + t106) + pkin(1) * t88;
t107 = t116 * qJDD(1);
t153 = t115 * t161;
t186 = t121 * t107 + t119 * t153;
t185 = t116 * t193;
t184 = t118 * t193;
t167 = t117 * t123;
t87 = t114 * t115 * t167;
t75 = t107 + t87;
t183 = t115 * t75;
t157 = t116 * t167;
t162 = qJDD(1) * t115;
t64 = (-t157 + t162) * t118;
t156 = t116 * t166;
t86 = t115 * t156;
t66 = -t86 - t152;
t45 = t115 * t66 + t117 * t64;
t182 = t116 * t45;
t76 = t107 - t87;
t181 = t117 * t76;
t129 = t143 * t170;
t172 = t112 * t123;
t14 = -pkin(4) * t107 - pkin(6) * t172 - t115 * t129 + t20;
t180 = t119 * t14;
t43 = t55 + t91;
t179 = t119 * t43;
t178 = t119 * t92;
t177 = t121 * t14;
t176 = t121 * t43;
t175 = t121 * t92;
t174 = qJDD(1) * pkin(1);
t111 = t115 ^ 2;
t173 = t111 * t114;
t113 = t117 ^ 2;
t171 = t113 * t114;
t168 = t116 * t123;
t160 = t117 * t55;
t158 = -t116 * g(3) + t146 * t118;
t15 = -pkin(4) * t172 + pkin(6) * t107 - t117 * t129 + t21;
t147 = t84 * t163 + t158;
t133 = -qJ(4) * t170 + qJDD(4) + t147;
t25 = ((-pkin(4) * t115 + pkin(6) * t117) * t168 + t136 * qJDD(1)) * t118 + t133;
t6 = t119 * t15 - t121 * t25;
t7 = t119 * t25 + t121 * t15;
t4 = t119 * t6 + t121 * t7;
t150 = t116 * (t146 * t116 + t190) + t118 * t158;
t130 = qJ(2) * t123 - t141;
t74 = t130 + t174;
t149 = t74 + t174;
t145 = t115 * t157;
t3 = t119 * t7 - t121 * t6;
t9 = t115 * t20 + t117 * t21;
t138 = t116 * t87;
t134 = (qJD(5) - t92) * t71 + t186;
t132 = -t70 * qJD(5) - t119 * t107 + t121 * t153;
t31 = t132 - t189;
t127 = t130 + 0.2e1 * t154;
t95 = 0.2e1 * t116 * t161;
t81 = (-t112 - t171) * t123;
t79 = (-t112 - t173) * t123;
t77 = (t111 + t113) * t170;
t65 = -t86 + t152;
t63 = (t157 + t162) * t118;
t59 = -t68 + t90;
t58 = t67 - t90;
t56 = t128 + t127;
t54 = -t115 * t79 - t181;
t53 = t117 * t81 - t183;
t52 = -t115 * t76 + t117 * t79;
t51 = t115 * t81 + t117 * t75;
t50 = t68 - t67;
t48 = qJD(5) * t71 + t186;
t47 = -t68 - t90;
t46 = -t115 * t64 + t117 * t66;
t41 = -t90 - t67;
t40 = t137 * t116 + t148;
t37 = t67 + t68;
t36 = pkin(3) * t161 + t133;
t32 = t132 + t189;
t30 = (qJD(5) + t92) * t71 + t186;
t27 = -t119 * t47 - t176;
t26 = t121 * t47 - t179;
t24 = t121 * t41 - t198;
t23 = t119 * t41 + t197;
t19 = -t119 * t32 + t121 * t134;
t18 = t119 * t134 + t121 * t32;
t17 = -t115 * t31 + t117 * t27;
t16 = t115 * t27 + t117 * t31;
t13 = -t115 * t30 + t117 * t24;
t12 = t115 * t24 + t117 * t30;
t11 = -t115 * t37 + t117 * t19;
t10 = t115 * t19 + t117 * t37;
t8 = t115 * t21 - t117 * t20;
t2 = t115 * t14 + t117 * t4;
t1 = t115 * t4 - t117 * t14;
t5 = [0, 0, 0, 0, 0, qJDD(1), t151, t142, 0, 0, t106, t95, 0, t108, 0, 0, t149 * t118 - t184, -t116 * t149 + t185, t150 + t187, pkin(1) * t74 + qJ(2) * t150, 0, 0, 0, t106, t95, t108, t118 * (pkin(2) * t88 + t147) + (qJ(3) * t88 + (qJD(1) * t84 + t146) * t116 + t148) * t116 + t187, t184 + (-t130 + t99 - 0.2e1 * t128) * t118, -t185 + (t127 + 0.2e1 * t128) * t116, qJ(2) * (t116 * t40 + t118 * t147) + t135 * t56, -t138 + (qJDD(1) * t111 + t145) * t114, t118 * (t115 * t65 + t117 * t63) + (t111 - t113) * t114 * t168, -t116 * t64 + t118 * (-t183 - (t112 - t173) * t167), t138 + (qJDD(1) * t113 - t145) * t114, t116 * t66 + t118 * (-t181 - t115 * (-t112 + t171) * t123), t106, t116 * (pkin(3) * t51 - t20) + t118 * (pkin(3) * t65 + t117 * t36) + qJ(2) * (t116 * t51 + t118 * t65) + t131 * t53, t116 * (pkin(3) * t52 - t21) + t118 * (-pkin(3) * t63 - t115 * t36) + qJ(2) * (t116 * t52 - t118 * t63) + t131 * t54, pkin(3) * t182 + t118 * (-pkin(3) * t77 - t9) + qJ(2) * (-t118 * t77 + t182) + t131 * t46, t131 * t9 + t196 * (t116 * t8 + t118 * t36), t116 * (-t119 * t132 - t175 * t71) + t118 * (-t115 * (-t121 * t132 + t71 * t178) + t160), t116 * (t119 * t30 - t121 * t31) + t118 * (-t115 * (t119 * t31 + t121 * t30) + t117 * t50), t116 * (t121 * t59 + t198) + t118 * (-t115 * (-t119 * t59 + t197) - t117 * t32), t116 * (t121 * t48 - t70 * t178) + t118 * (-t115 * (-t119 * t48 - t175 * t70) - t160), t116 * (t119 * t58 + t176) + t118 * (-t115 * (t121 * t58 - t179) + t117 * t134), t116 * (t119 * t70 + t121 * t71) * t92 + t118 * (t117 * t91 - t115 * (-t119 * t71 + t121 * t70) * t92), t116 * (pkin(3) * t12 + pkin(4) * t30 + pkin(6) * t24 - t177) + t118 * (-t115 * (-pkin(6) * t23 + t180) - t117 * (-pkin(4) * t23 + t6) + pkin(3) * t23) + qJ(2) * (t116 * t12 + t118 * t23) + t131 * t13, t116 * (pkin(3) * t16 + pkin(4) * t31 + pkin(6) * t27 + t180) + t118 * (-t115 * (-pkin(6) * t26 + t177) - t117 * (-pkin(4) * t26 + t7) + pkin(3) * t26) + qJ(2) * (t116 * t16 + t118 * t26) + t131 * t17, t116 * (pkin(3) * t10 + pkin(4) * t37 + pkin(6) * t19 + t4) + t118 * (t115 * t3 + t136 * t18) + qJ(2) * (t10 * t116 + t118 * t18) + t131 * t11, (-pkin(4) * t14 + pkin(6) * t4 + t196 * t1) * t116 + (qJ(2) + t136) * t118 * t3 + t131 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t161, t107, -t88, -t74, 0, 0, 0, 0, 0, 0, -t88, t161, -t107, -t56, 0, 0, 0, 0, 0, 0, t53, t54, t46, t9, 0, 0, 0, 0, 0, 0, t13, t17, t11, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t107, t156, -t172, t40, 0, 0, 0, 0, 0, 0, t51, t52, t45, t8, 0, 0, 0, 0, 0, 0, t12, t16, t10, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, -t63, -t77, t36, 0, 0, 0, 0, 0, 0, t23, t26, t18, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, t50, -t32, -t55, t134, t91, -t6, -t7, 0, 0;];
tauJ_reg = t5;
