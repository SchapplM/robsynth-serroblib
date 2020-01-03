% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S4RPRR6
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
% Datum: 2019-12-31 16:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S4RPRR6_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR6_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR6_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR6_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR6_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR6_invdynJ_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:52:41
% EndTime: 2019-12-31 16:52:46
% DurationCPUTime: 1.67s
% Computational Cost: add. (4333->239), mult. (10738->345), div. (0->0), fcn. (7816->8), ass. (0->152)
t124 = sin(qJ(3));
t121 = sin(pkin(7));
t122 = cos(pkin(7));
t127 = cos(qJ(3));
t105 = (-t121 * t124 + t122 * t127) * qJD(1);
t136 = t121 * t127 + t122 * t124;
t107 = t136 * qJD(1);
t159 = t107 * t105;
t180 = qJDD(3) + t159;
t185 = t124 * t180;
t184 = t127 * t180;
t129 = qJD(1) ^ 2;
t125 = sin(qJ(1));
t175 = cos(qJ(1));
t135 = t175 * g(1) + t125 * g(2);
t181 = -t129 * pkin(1) + qJDD(1) * qJ(2) + 0.2e1 * qJD(1) * qJD(2) - t135;
t118 = t121 ^ 2;
t119 = t122 ^ 2;
t155 = t118 + t119;
t123 = sin(qJ(4));
t117 = qJDD(3) + qJDD(4);
t126 = cos(qJ(4));
t78 = -t126 * t105 + t123 * t107;
t80 = t123 * t105 + t126 * t107;
t61 = t80 * t78;
t179 = -t61 + t117;
t183 = t123 * t179;
t182 = t126 * t179;
t173 = t122 * g(3);
t174 = t122 * pkin(2);
t178 = -t173 + (-pkin(5) * qJDD(1) + t129 * t174 - t181) * t121;
t145 = t125 * g(1) - t175 * g(2);
t139 = -qJDD(2) + t145;
t156 = t129 * qJ(2);
t160 = qJDD(1) * pkin(1);
t101 = t139 + t156 + t160;
t177 = t155 * t156 - t101 - t160;
t154 = t105 * qJD(3);
t140 = -t121 * g(3) + t181 * t122;
t149 = t122 * qJDD(1);
t75 = -t119 * t129 * pkin(2) + pkin(5) * t149 + t140;
t52 = t124 * t75 - t127 * t178;
t104 = t136 * qJDD(1);
t91 = t104 + t154;
t176 = -t52 + (-t91 + t154) * pkin(6);
t76 = t78 ^ 2;
t77 = t80 ^ 2;
t102 = t105 ^ 2;
t103 = t107 ^ 2;
t120 = qJD(3) + qJD(4);
t116 = t120 ^ 2;
t132 = pkin(3) * t180 + t176;
t138 = qJD(3) * pkin(3) - t107 * pkin(6);
t53 = t178 * t124 + t127 * t75;
t153 = t107 * qJD(3);
t150 = t121 * qJDD(1);
t68 = -t124 * t150 + t127 * t149;
t89 = t68 - t153;
t30 = -t102 * pkin(3) + t89 * pkin(6) - qJD(3) * t138 + t53;
t13 = t123 * t30 - t126 * t132;
t165 = t126 * t30;
t14 = t123 * t132 + t165;
t6 = t123 * t14 - t126 * t13;
t172 = t124 * t6;
t171 = t127 * t6;
t25 = t124 * t53 - t127 * t52;
t170 = t121 * t25;
t84 = (pkin(1) + t174) * qJDD(1) + (t155 * pkin(5) + qJ(2)) * t129 + t139;
t47 = t89 * pkin(3) + t102 * pkin(6) - t107 * t138 + t84;
t169 = t123 * t47;
t58 = t61 + t117;
t168 = t123 * t58;
t167 = t124 * t84;
t86 = qJDD(3) - t159;
t166 = t124 * t86;
t164 = t126 * t47;
t163 = t126 * t58;
t162 = t127 * t84;
t161 = t127 * t86;
t158 = t120 * t123;
t157 = t120 * t126;
t152 = qJD(4) + t120;
t7 = t123 * t13 + t126 * t14;
t143 = t121 * (t181 * t121 + t173) + t122 * t140;
t142 = t123 * t91 - t126 * t89;
t26 = t124 * t52 + t127 * t53;
t137 = t123 * t89 + t126 * t91;
t134 = (-qJD(4) + t120) * t80 - t142;
t50 = -t78 * qJD(4) + t137;
t128 = qJD(3) ^ 2;
t114 = t119 * qJDD(1);
t113 = t118 * qJDD(1);
t109 = t155 * t129;
t97 = -t103 - t128;
t96 = -t103 + t128;
t95 = t102 - t128;
t90 = t104 + 0.2e1 * t154;
t88 = -t68 + 0.2e1 * t153;
t83 = -t128 - t102;
t74 = t120 * t78;
t73 = -t77 + t116;
t72 = t76 - t116;
t71 = -t77 - t116;
t69 = -t102 - t103;
t67 = -t124 * t97 - t161;
t66 = t127 * t97 - t166;
t65 = t124 * t104 + t127 * t68;
t64 = -t127 * t104 + t124 * t68;
t63 = t127 * t83 - t185;
t62 = t124 * t83 + t184;
t60 = t77 - t76;
t56 = -t116 - t76;
t55 = (t123 * t80 - t126 * t78) * t120;
t54 = (-t123 * t78 - t126 * t80) * t120;
t49 = -t80 * qJD(4) - t142;
t48 = -t76 - t77;
t46 = t126 * t72 - t168;
t45 = -t123 * t73 + t182;
t44 = t123 * t72 + t163;
t43 = t126 * t73 + t183;
t42 = -t123 * t71 - t163;
t41 = t126 * t71 - t168;
t40 = t50 + t74;
t39 = t50 - t74;
t38 = -t152 * t78 + t137;
t35 = t152 * t80 + t142;
t34 = t126 * t50 - t80 * t158;
t33 = t123 * t50 + t80 * t157;
t32 = -t123 * t49 + t78 * t157;
t31 = t126 * t49 + t78 * t158;
t29 = t126 * t56 - t183;
t28 = t123 * t56 + t182;
t24 = -pkin(6) * t41 - t164;
t23 = -t124 * t41 + t127 * t42;
t22 = t124 * t42 + t127 * t41;
t21 = -pkin(6) * t28 - t169;
t20 = t123 * t40 + t126 * t134;
t19 = -t123 * t39 - t126 * t35;
t18 = t123 * t134 - t126 * t40;
t17 = -t123 * t35 + t126 * t39;
t16 = -t124 * t28 + t127 * t29;
t15 = t124 * t29 + t127 * t28;
t11 = -pkin(3) * t38 + pkin(6) * t42 - t169;
t10 = -pkin(3) * t35 + pkin(6) * t29 + t164;
t9 = -t124 * t18 + t127 * t20;
t8 = t124 * t20 + t127 * t18;
t5 = pkin(3) * t47 + pkin(6) * t7;
t4 = -pkin(6) * t18 - t6;
t3 = -pkin(3) * t48 + pkin(6) * t20 + t7;
t2 = t127 * t7 - t172;
t1 = t124 * t7 + t171;
t12 = [0, 0, 0, 0, 0, qJDD(1), t145, t135, 0, 0, t113, 0.2e1 * t121 * t149, 0, t114, 0, 0, -t177 * t122, t177 * t121, pkin(1) * t109 + qJ(2) * (t114 + t113) + t143, pkin(1) * t101 + qJ(2) * t143, t121 * (-t124 * t153 + t127 * t91) + t122 * (t124 * t91 + t127 * t153), t121 * (-t124 * t90 - t127 * t88) + t122 * (-t124 * t88 + t127 * t90), t121 * (-t124 * t96 + t184) + t122 * (t127 * t96 + t185), t121 * (-t124 * t89 - t127 * t154) + t122 * (-t124 * t154 + t127 * t89), t121 * (t127 * t95 - t166) + t122 * (t124 * t95 + t161), (t121 * (t105 * t127 + t107 * t124) + t122 * (t105 * t124 - t107 * t127)) * qJD(3), t121 * (-pkin(5) * t62 - t167) + t122 * (-pkin(2) * t88 + pkin(5) * t63 + t162) - pkin(1) * t88 + qJ(2) * (-t121 * t62 + t122 * t63), t121 * (-pkin(5) * t66 - t162) + t122 * (-pkin(2) * t90 + pkin(5) * t67 - t167) - pkin(1) * t90 + qJ(2) * (-t121 * t66 + t122 * t67), t121 * (-pkin(5) * t64 - t25) + t122 * (-pkin(2) * t69 + pkin(5) * t65 + t26) - pkin(1) * t69 + qJ(2) * (-t121 * t64 + t122 * t65), -pkin(5) * t170 + t122 * (pkin(2) * t84 + pkin(5) * t26) + pkin(1) * t84 + qJ(2) * (t122 * t26 - t170), t121 * (-t124 * t33 + t127 * t34) + t122 * (t124 * t34 + t127 * t33), t121 * (-t124 * t17 + t127 * t19) + t122 * (t124 * t19 + t127 * t17), t121 * (-t124 * t43 + t127 * t45) + t122 * (t124 * t45 + t127 * t43), t121 * (-t124 * t31 + t127 * t32) + t122 * (t124 * t32 + t127 * t31), t121 * (-t124 * t44 + t127 * t46) + t122 * (t124 * t46 + t127 * t44), t121 * (-t124 * t54 + t127 * t55) + t122 * (t124 * t55 + t127 * t54), t121 * (-pkin(5) * t15 - t124 * t10 + t127 * t21) + t122 * (-pkin(2) * t35 + pkin(5) * t16 + t127 * t10 + t124 * t21) - pkin(1) * t35 + qJ(2) * (-t121 * t15 + t122 * t16), t121 * (-pkin(5) * t22 - t124 * t11 + t127 * t24) + t122 * (-pkin(2) * t38 + pkin(5) * t23 + t127 * t11 + t124 * t24) - pkin(1) * t38 + qJ(2) * (-t121 * t22 + t122 * t23), t121 * (-pkin(5) * t8 - t124 * t3 + t127 * t4) + t122 * (-pkin(2) * t48 + pkin(5) * t9 + t124 * t4 + t127 * t3) - pkin(1) * t48 + qJ(2) * (-t121 * t8 + t122 * t9), t121 * (-pkin(5) * t1 - pkin(6) * t171 - t124 * t5) + t122 * (pkin(2) * t47 + pkin(5) * t2 - pkin(6) * t172 + t127 * t5) + pkin(1) * t47 + qJ(2) * (-t121 * t1 + t122 * t2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t149, t150, -t109, -t101, 0, 0, 0, 0, 0, 0, t88, t90, t69, -t84, 0, 0, 0, 0, 0, 0, t35, t38, t48, -t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t159, t103 - t102, t104, t159, t68, qJDD(3), -t52, -t53, 0, 0, t61, t60, t40, -t61, t134, t117, pkin(3) * t28 - t13, -t165 - t123 * t176 + (-t123 * t180 + t41) * pkin(3), pkin(3) * t18, pkin(3) * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, t60, t40, -t61, t134, t117, -t13, -t14, 0, 0;];
tauJ_reg = t12;
