% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPPRPR6
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 14:27
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPPRPR6_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR6_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR6_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR6_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR6_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRPR6_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:27:23
% EndTime: 2019-05-05 14:27:27
% DurationCPUTime: 1.33s
% Computational Cost: add. (2835->223), mult. (5613->236), div. (0->0), fcn. (2728->6), ass. (0->148)
t174 = qJ(2) - pkin(7);
t112 = sin(qJ(4));
t115 = cos(qJ(4));
t118 = qJD(1) ^ 2;
t157 = t115 * t118;
t92 = t112 * t157;
t86 = qJDD(4) - t92;
t162 = t115 * t86;
t117 = qJD(4) ^ 2;
t108 = t112 ^ 2;
t159 = t108 * t118;
t88 = -t117 - t159;
t49 = t112 * t88 + t162;
t190 = t174 * t49;
t175 = pkin(1) + qJ(3);
t85 = qJDD(4) + t92;
t169 = t112 * t85;
t109 = t115 ^ 2;
t100 = t109 * t118;
t90 = -t100 - t117;
t51 = -t115 * t90 + t169;
t189 = t174 * t51;
t111 = sin(qJ(6));
t114 = cos(qJ(6));
t155 = qJD(1) * t112;
t71 = qJD(4) * t111 - t114 * t155;
t73 = qJD(4) * t114 + t111 * t155;
t47 = t73 * t71;
t149 = t115 * qJDD(1);
t152 = qJD(1) * qJD(4);
t98 = t112 * t152;
t78 = -t98 + t149;
t67 = qJDD(6) + t78;
t186 = -t47 + t67;
t188 = t111 * t186;
t187 = t114 * t186;
t144 = t115 * t152;
t150 = t112 * qJDD(1);
t77 = t144 + t150;
t37 = -t71 * qJD(6) + t114 * qJDD(4) + t111 * t77;
t154 = qJD(1) * t115;
t94 = qJD(6) + t154;
t58 = t94 * t71;
t185 = -t58 + t37;
t158 = t115 * qJ(5);
t136 = pkin(4) * t112 - t158;
t75 = t136 * qJD(1);
t184 = qJDD(4) * qJ(5) - t75 * t155;
t113 = sin(qJ(1));
t116 = cos(qJ(1));
t143 = t113 * g(1) - t116 * g(2);
t139 = -qJDD(2) + t143;
t142 = t175 * qJDD(1);
t110 = t118 * pkin(7);
t181 = t77 * pkin(4) - t78 * qJ(5) - t110;
t121 = t142 + t139 + t181;
t161 = qJ(5) * t112;
t137 = pkin(4) * t115 + t161;
t153 = qJD(5) * t115;
t156 = t118 * qJ(2);
t179 = 2 * qJD(3);
t22 = (t137 * qJD(4) - 0.2e1 * t153 + t179) * qJD(1) + t121 + t156;
t128 = qJD(1) * t179 + t139;
t54 = t142 + t128 + t156;
t183 = t162 - t112 * (-t100 + t117);
t182 = t75 * t154 + qJDD(5);
t65 = t71 ^ 2;
t66 = t73 ^ 2;
t91 = t94 ^ 2;
t180 = 0.2e1 * qJD(1);
t178 = pkin(4) + pkin(8);
t177 = t112 * g(3);
t176 = t117 * pkin(4);
t141 = t111 * qJDD(4) - t114 * t77;
t127 = (-qJD(6) + t94) * t73 - t141;
t31 = t58 + t37;
t11 = t111 * t127 - t114 * t31;
t173 = t11 * t115;
t107 = qJDD(1) * qJ(2);
t140 = t116 * g(1) + t113 * g(2);
t134 = qJD(2) * t180 - t140;
t131 = qJDD(3) + t134;
t55 = -t175 * t118 + t107 + t131;
t46 = -qJDD(1) * pkin(7) + t55;
t43 = t112 * t46;
t39 = t115 * g(3) - t43;
t126 = -t39 + t184;
t122 = -t126 + t176;
t146 = pkin(5) * t154;
t151 = qJD(5) * qJD(4);
t15 = pkin(8) * t159 + t77 * pkin(5) - 0.2e1 * t151 - qJD(4) * (-qJD(4) * pkin(8) + t146) + t122;
t172 = t111 * t15;
t35 = t47 + t67;
t171 = t111 * t35;
t170 = t111 * t94;
t166 = t114 * t15;
t165 = t114 * t35;
t164 = t114 * t94;
t163 = t115 * t46;
t160 = qJDD(1) * pkin(1);
t148 = -t66 - t91;
t147 = t115 * t47;
t13 = t77 * pkin(8) + (-pkin(5) * t108 + qJ(2)) * t118 + (qJD(4) * t161 + t179 + (t178 * qJD(4) - 0.2e1 * qJD(5) - t146) * t115) * qJD(1) + t121;
t130 = -t117 * qJ(5) - t163 + t182;
t16 = t78 * pkin(5) - t178 * qJDD(4) + (pkin(5) * t152 + pkin(8) * t157 - g(3)) * t112 + t130;
t5 = t111 * t13 - t114 * t16;
t6 = t111 * t16 + t114 * t13;
t2 = t111 * t6 - t114 * t5;
t3 = t111 * t5 + t114 * t6;
t38 = t163 + t177;
t138 = t112 * t15 + t115 * t2;
t19 = -t112 * t39 + t115 * t38;
t135 = -t115 * (-t117 + t159) + t169;
t82 = (-t108 - t109) * qJDD(1);
t83 = t100 + t159;
t133 = t174 * t82 - t175 * t83;
t132 = t136 + t175;
t125 = t112 * t178 - t158 + t175;
t124 = -qJDD(4) * pkin(4) + t130;
t102 = 0.2e1 * t151;
t101 = 0.2e1 * t107;
t84 = t100 - t159;
t79 = -0.2e1 * t98 + t149;
t76 = 0.2e1 * t144 + t150;
t59 = t139 + t156 + t160;
t57 = -t66 + t91;
t56 = t65 - t91;
t53 = (t78 - t98) * t115;
t48 = (t77 + t144) * t112;
t45 = -t110 + t54;
t44 = t66 - t65;
t42 = -t112 * t79 - t115 * t76;
t40 = -t91 - t65;
t36 = -qJD(6) * t73 - t141;
t33 = -t65 - t66;
t26 = (qJD(6) + t94) * t73 + t141;
t24 = -t124 + t177;
t23 = t102 - t122;
t21 = -t111 * t148 - t165;
t20 = t114 * t148 - t171;
t18 = t114 * t40 - t188;
t17 = t111 * t40 + t187;
t12 = t111 * t31 + t114 * t127;
t10 = t112 * t23 + t115 * t24;
t9 = t112 * t185 - t115 * t20;
t8 = t112 * t26 - t115 * t17;
t7 = t112 * t33 - t173;
t1 = [0, 0, 0, 0, 0, qJDD(1), t143, t140, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, -t139 - 0.2e1 * t160, t101 + t134, qJ(2) * (-pkin(1) * t118 + t107 + t134) + pkin(1) * t59, qJDD(1), 0, 0, 0, 0, 0, 0, t101 + t131, t128 + 0.2e1 * t142, qJ(2) * t55 + t175 * t54, t53, t42, t183, t48, -t135, 0, t112 * t45 + t175 * t76 + t190, t115 * t45 + t175 * t79 - t189, t133 - t19, t174 * t19 + t175 * t45, 0, -t183, t135, t53, t42, t48, t115 * (qJ(5) * t83 + t124) - t112 * (pkin(4) * t83 + t102 - t176 + t184 + t43) + t133, -t112 * t22 - t132 * t76 - t190, t115 * (-t137 * t152 + t153 * t180 - t181 - t54) + t189 - t132 * t79, t174 * t10 + t132 * t22, t147 - t112 * (-t111 * t37 - t73 * t164), t115 * t44 - t112 * (t111 * t26 - t114 * t185), t115 * t31 - t112 * (-t114 * t57 - t188), -t147 - t112 * (-t114 * t36 - t71 * t170), t115 * t127 - t112 * (-t111 * t56 - t165), t115 * t67 - t112 * (t111 * t71 + t114 * t73) * t94, t115 * (pkin(5) * t17 - t5) - t112 * (pkin(5) * t26 - t166) + t174 * t8 + t125 * t18, t115 * (pkin(5) * t20 - t6) - t112 * (pkin(5) * t185 + t172) + t174 * t9 + t125 * t21, pkin(5) * t173 - t112 * (pkin(5) * t33 - t3) + t174 * t7 + t125 * t12, t125 * t3 + (pkin(5) - t174) * t138; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t118, -t59, 0, 0, 0, 0, 0, 0, 0, -t118, -qJDD(1), -t54, 0, 0, 0, 0, 0, 0, -t76, -t79, t83, -t45, 0, 0, 0, 0, 0, 0, t83, t76, t79, -t22, 0, 0, 0, 0, 0, 0, -t18, -t21, -t12, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t118, t55, 0, 0, 0, 0, 0, 0, t49, -t51, t82, t19, 0, 0, 0, 0, 0, 0, t82, -t49, t51, t10, 0, 0, 0, 0, 0, 0, t8, t9, t7, -t138; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t92, t84, t149, -t92, -t150, qJDD(4), t38, t39, 0, 0, qJDD(4), -t149, t150, t92, t84, -t92, -t137 * qJDD(1), (-t117 - t88) * qJ(5) + (-qJDD(4) - t86) * pkin(4) - t38 + t182, qJ(5) * t85 + t102 + (-t117 - t90) * pkin(4) + t126, pkin(4) * t24 + qJ(5) * t23, t114 * t37 - t73 * t170, -t111 * t185 - t114 * t26, -t111 * t57 + t187, -t111 * t36 + t71 * t164, t114 * t56 - t171, (t111 * t73 - t114 * t71) * t94, qJ(5) * t26 - t178 * t17 - t172, qJ(5) * t185 - t178 * t20 - t166, qJ(5) * t33 - t178 * t11 - t2, -qJ(5) * t15 - t178 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t149, t86, t90, -t24, 0, 0, 0, 0, 0, 0, t17, t20, t11, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, t44, t31, -t47, t127, t67, -t5, -t6, 0, 0;];
tauJ_reg  = t1;
