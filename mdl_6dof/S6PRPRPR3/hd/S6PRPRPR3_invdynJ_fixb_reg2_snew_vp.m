% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6PRPRPR3
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 22:37
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6PRPRPR3_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR3_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR3_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR3_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_invdynJ_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 22:36:21
% EndTime: 2019-05-04 22:36:28
% DurationCPUTime: 2.49s
% Computational Cost: add. (5251->261), mult. (10130->354), div. (0->0), fcn. (6819->12), ass. (0->165)
t156 = cos(qJ(4));
t187 = qJD(2) * qJD(4);
t178 = t156 * t187;
t153 = sin(qJ(4));
t185 = t153 * qJDD(2);
t106 = 0.2e1 * t178 + t185;
t144 = sin(pkin(11));
t146 = sin(pkin(6));
t147 = cos(pkin(11));
t149 = cos(pkin(6));
t154 = sin(qJ(2));
t157 = cos(qJ(2));
t139 = t153 ^ 2;
t159 = qJD(2) ^ 2;
t135 = t139 * t159;
t158 = qJD(4) ^ 2;
t122 = -t135 - t158;
t190 = t156 * t159;
t180 = t153 * t190;
t119 = qJDD(4) - t180;
t191 = t156 * t119;
t85 = t153 * t122 + t191;
t58 = t147 * t106 + t144 * t85;
t81 = t153 * t119 - t156 * t122;
t231 = t149 * t81 + (t154 * (-t144 * t106 + t147 * t85) + t157 * t58) * t146;
t230 = pkin(2) * t58 + pkin(8) * t85;
t179 = t153 * t187;
t184 = t156 * qJDD(2);
t109 = -0.2e1 * t179 + t184;
t140 = t156 ^ 2;
t136 = t140 * t159;
t124 = -t136 - t158;
t118 = qJDD(4) + t180;
t194 = t153 * t118;
t84 = -t156 * t124 + t194;
t57 = -t147 * t109 + t144 * t84;
t79 = t156 * t118 + t153 * t124;
t228 = (t154 * (t144 * t109 + t147 * t84) + t157 * t57) * t146 - t149 * t79;
t226 = -pkin(2) * t57 - pkin(8) * t84;
t152 = sin(qJ(6));
t155 = cos(qJ(6));
t189 = qJD(2) * t156;
t101 = t152 * qJD(4) + t155 * t189;
t103 = t155 * qJD(4) - t152 * t189;
t75 = t103 * t101;
t107 = t178 + t185;
t96 = qJDD(6) + t107;
t214 = -t75 + t96;
t218 = t152 * t214;
t217 = t155 * t214;
t216 = t107 + t178;
t195 = t153 * qJ(5);
t172 = -t156 * pkin(4) - t195;
t145 = sin(pkin(10));
t148 = cos(pkin(10));
t117 = -t148 * g(1) - t145 * g(2);
t141 = -g(3) + qJDD(1);
t116 = t145 * g(1) - t148 * g(2);
t196 = t149 * t116;
t170 = t141 * t146 + t196;
t65 = -t154 * t117 + t170 * t157;
t161 = qJDD(2) * pkin(2) + t65;
t66 = t157 * t117 + t170 * t154;
t61 = -t159 * pkin(2) + t66;
t40 = t144 * t161 + t147 * t61;
t38 = -t159 * pkin(3) + qJDD(2) * pkin(8) + t40;
t175 = t159 * t172 + t38;
t213 = -t158 * pkin(4) + t175 * t156;
t212 = t191 + t153 * (t136 - t158);
t94 = t101 ^ 2;
t95 = t103 ^ 2;
t188 = t153 * qJD(2);
t128 = qJD(6) + t188;
t125 = t128 ^ 2;
t211 = 2 * qJD(5);
t210 = -pkin(4) - pkin(9);
t108 = -t179 + t184;
t120 = pkin(5) * t188 - qJD(4) * pkin(9);
t186 = qJDD(4) * qJ(5);
t132 = t149 * t141;
t87 = -t146 * t116 + qJDD(3) + t132;
t206 = t153 * t87;
t20 = t186 + t206 - pkin(9) * t136 + t108 * pkin(5) + (t211 + t120) * qJD(4) + t213;
t209 = t152 * t20;
t64 = t75 + t96;
t208 = t152 * t64;
t174 = t152 * qJDD(4) + t155 * t108;
t165 = (-qJD(6) + t128) * t103 - t174;
t203 = t128 * t101;
t68 = -t101 * qJD(6) + t155 * qJDD(4) - t152 * t108;
t52 = t68 + t203;
t35 = t152 * t165 - t155 * t52;
t207 = t153 * t35;
t205 = t155 * t20;
t204 = t155 * t64;
t76 = t156 * t87;
t202 = t128 * t152;
t201 = t128 * t155;
t183 = -t95 - t125;
t113 = (t139 + t140) * qJDD(2);
t114 = t135 + t136;
t73 = t144 * t113 + t147 * t114;
t182 = pkin(2) * t73 + pkin(3) * t114 + pkin(8) * t113;
t181 = t153 * t75;
t177 = t144 * t61 - t147 * t161;
t173 = -qJDD(4) * pkin(4) - t158 * qJ(5) + qJDD(5) - t76;
t21 = -qJDD(4) * pkin(9) + (t107 - t178) * pkin(5) + (-pkin(9) * t190 + t175) * t153 + t173;
t37 = -qJDD(2) * pkin(3) - t159 * pkin(8) + t177;
t164 = -t108 * pkin(4) - t216 * qJ(5) + t37;
t176 = pkin(4) * qJD(4) - (2 * qJD(5));
t22 = -pkin(5) * t136 - t108 * pkin(9) + (-t120 + t176) * t188 + t164;
t9 = t152 * t22 - t155 * t21;
t33 = t153 * t38 - t76;
t34 = t156 * t38 + t206;
t15 = t153 * t33 + t156 * t34;
t10 = t152 * t21 + t155 * t22;
t5 = t155 * t10 + t152 * t9;
t4 = t152 * t10 - t155 * t9;
t3 = t153 * t4 + t156 * t20;
t169 = t156 * (-t135 + t158) + t194;
t168 = t68 - t203;
t167 = pkin(3) - t172;
t166 = t156 * t210 - pkin(3) - t195;
t163 = qJD(4) * t211 + t213;
t26 = t175 * t153 + t173;
t160 = t163 + t186;
t27 = t176 * t188 + t164;
t115 = t135 - t136;
t112 = -t144 * qJDD(2) - t147 * t159;
t111 = t147 * qJDD(2) - t144 * t159;
t89 = -t95 + t125;
t88 = t94 - t125;
t78 = t216 * t153;
t77 = (t108 - t179) * t156;
t74 = t95 - t94;
t71 = t156 * t106 + t153 * t109;
t69 = -t125 - t94;
t67 = -t103 * qJD(6) - t174;
t62 = -t94 - t95;
t47 = (qJD(6) + t128) * t103 + t174;
t45 = (t154 * (t147 * t113 - t144 * t114) + t157 * t73) * t146;
t44 = -t152 * t183 - t204;
t43 = t155 * t183 - t208;
t42 = t155 * t69 - t218;
t41 = t152 * t69 + t217;
t36 = t152 * t52 + t155 * t165;
t31 = t153 * t43 + t156 * t168;
t30 = t153 * t168 - t156 * t43;
t29 = t153 * t41 + t156 * t47;
t28 = t153 * t47 - t156 * t41;
t25 = t160 + t206;
t24 = t156 * t62 + t207;
t23 = t153 * t62 - t156 * t35;
t18 = t144 * t40 - t147 * t177;
t17 = t144 * t31 - t147 * t44;
t16 = t144 * t29 - t147 * t42;
t14 = t153 * t34 - t156 * t33;
t13 = t144 * t24 - t147 * t36;
t12 = t153 * t26 + t156 * t25;
t11 = t153 * t25 - t156 * t26;
t8 = t144 * t15 - t147 * t37;
t6 = t144 * t12 - t147 * t27;
t2 = t153 * t20 - t156 * t4;
t1 = t144 * t3 - t147 * t5;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t141, 0, 0, 0, 0, 0, 0, (qJDD(2) * t157 - t154 * t159) * t146, (-qJDD(2) * t154 - t157 * t159) * t146, 0, t149 * t132 + (t154 * t66 + t157 * t65 - t196) * t146, 0, 0, 0, 0, 0, 0, (t111 * t157 + t112 * t154) * t146, (-t111 * t154 + t112 * t157) * t146, 0, t149 * t87 + (t154 * (t144 * t177 + t147 * t40) + t157 * t18) * t146, 0, 0, 0, 0, 0, 0, -t228, -t231, t45, t149 * t14 + (t154 * (t144 * t37 + t147 * t15) + t157 * t8) * t146, 0, 0, 0, 0, 0, 0, t45, t228, t231, t149 * t11 + (t154 * (t147 * t12 + t144 * t27) + t157 * t6) * t146, 0, 0, 0, 0, 0, 0, t149 * t28 + (t154 * (t144 * t42 + t147 * t29) + t157 * t16) * t146, t149 * t30 + (t154 * (t144 * t44 + t147 * t31) + t157 * t17) * t146, t149 * t23 + (t154 * (t144 * t36 + t147 * t24) + t157 * t13) * t146, t149 * t2 + (t154 * (t144 * t5 + t147 * t3) + t157 * t1) * t146; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t65, -t66, 0, 0, 0, 0, 0, 0, 0, qJDD(2), pkin(2) * t111 - t177, pkin(2) * t112 - t40, 0, pkin(2) * t18, t78, t71, t169, t77, t212, 0, pkin(3) * t109 - t156 * t37 + t226, -pkin(3) * t106 + t153 * t37 - t230, t15 + t182, pkin(2) * t8 - pkin(3) * t37 + pkin(8) * t15, 0, -t169, -t212, t78, t71, t77, (pkin(4) * t114 + t160) * t156 + (qJ(5) * t114 + t26 + t76) * t153 + t182, -t167 * t109 + t156 * t27 - t226, t153 * (-pkin(4) * t179 + t188 * t211 - t164) + t167 * t106 + t230, pkin(2) * t6 + pkin(8) * t12 - t167 * t27, t181 + t156 * (-t103 * t201 - t152 * t68), t153 * t74 + t156 * (t152 * t47 - t155 * t168), t153 * t52 + t156 * (-t155 * t89 - t218), -t181 + t156 * (-t101 * t202 - t155 * t67), t153 * t165 + t156 * (-t152 * t88 - t204), t153 * t96 + t156 * (t101 * t152 + t103 * t155) * t128, t153 * (pkin(5) * t41 - t9) + t156 * (pkin(5) * t47 + t205) + pkin(8) * t29 + pkin(2) * t16 + t166 * t42, t153 * (pkin(5) * t43 - t10) + t156 * (pkin(5) * t168 - t209) + pkin(8) * t31 + pkin(2) * t17 + t166 * t44, pkin(5) * t207 + t156 * (pkin(5) * t62 - t5) + pkin(8) * t24 + pkin(2) * t13 + t166 * t36, pkin(2) * t1 + t166 * t5 + (pkin(5) + pkin(8)) * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, 0, 0, 0, 0, 0, 0, t79, -t81, 0, t14, 0, 0, 0, 0, 0, 0, 0, -t79, t81, t11, 0, 0, 0, 0, 0, 0, t28, t30, t23, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t180, t115, t185, t180, t184, qJDD(4), -t33, -t34, 0, 0, qJDD(4), -t185, -t184, -t180, t115, t180, (-pkin(4) * t153 + qJ(5) * t156) * qJDD(2), -pkin(4) * t118 - qJ(5) * t124 + t26, -pkin(4) * t122 + t206 + (qJDD(4) + t119) * qJ(5) + t163, -pkin(4) * t26 + qJ(5) * t25, -t103 * t202 + t155 * t68, -t152 * t168 - t155 * t47, -t152 * t89 + t217, t101 * t201 - t152 * t67, t155 * t88 - t208, (-t101 * t155 + t103 * t152) * t128, qJ(5) * t47 + t210 * t41 + t209, qJ(5) * t168 + t210 * t43 + t205, qJ(5) * t62 + t210 * t35 - t4, qJ(5) * t20 + t210 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t185, t118, t122, t26, 0, 0, 0, 0, 0, 0, t41, t43, t35, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, t74, t52, -t75, t165, t96, -t9, -t10, 0, 0;];
tauJ_reg  = t7;
