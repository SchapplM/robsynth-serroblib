% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6PRPRRP6
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 00:03
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6PRPRRP6_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP6_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP6_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP6_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP6_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP6_invdynJ_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 00:02:09
% EndTime: 2019-05-05 00:02:15
% DurationCPUTime: 2.50s
% Computational Cost: add. (4989->264), mult. (9404->336), div. (0->0), fcn. (6331->10), ass. (0->183)
t148 = cos(qJ(5));
t145 = sin(qJ(5));
t149 = cos(qJ(4));
t188 = qJD(2) * qJD(4);
t182 = t149 * t188;
t146 = sin(qJ(4));
t187 = t146 * qJDD(2);
t117 = -t182 - t187;
t110 = qJDD(5) - t117;
t190 = qJD(2) * t149;
t112 = -qJD(4) * t148 + t145 * t190;
t114 = qJD(4) * t145 + t148 * t190;
t199 = t114 * t112;
t226 = t110 + t199;
t212 = t145 * t226;
t109 = t114 ^ 2;
t131 = qJD(2) * t146 + qJD(5);
t220 = t131 ^ 2;
t230 = -t109 - t220;
t34 = -t148 * t230 + t212;
t271 = pkin(4) * t34;
t270 = pkin(9) * t34;
t206 = t148 * t226;
t36 = t145 * t230 + t206;
t269 = pkin(9) * t36;
t268 = qJ(3) * t34;
t267 = t146 * t36;
t147 = sin(qJ(2));
t266 = t147 * t34;
t265 = t149 * t36;
t133 = t149 * qJDD(2);
t183 = t146 * t188;
t118 = t133 - t183;
t178 = qJDD(4) * t148 - t145 * t118;
t164 = qJD(5) * t114 - t178;
t96 = t114 * t131;
t48 = t164 - t96;
t221 = t112 ^ 2;
t90 = t221 - t220;
t264 = t146 * t48 + t149 * (-t148 * t90 + t212);
t263 = t145 * t90 + t206;
t167 = -t145 * qJDD(4) - t148 * t118;
t160 = -qJD(5) * t112 - t167;
t200 = t112 * t131;
t225 = t160 - t200;
t214 = t145 * t225;
t229 = t109 - t221;
t231 = t164 + t96;
t261 = -t146 * t229 + t149 * (t148 * t231 + t214);
t218 = -pkin(8) - pkin(2);
t223 = -t220 - t221;
t227 = t110 - t199;
t61 = t148 * t227;
t234 = t145 * t223 + t61;
t211 = t145 * t227;
t233 = t148 * t223 - t211;
t250 = t146 * t233 - t149 * t231;
t260 = qJ(3) * t234 + t218 * t250;
t142 = sin(pkin(6));
t143 = cos(pkin(6));
t150 = cos(qJ(2));
t259 = t143 * (t146 * t231 + t149 * t233) + (t147 * t234 - t150 * t250) * t142;
t258 = pkin(4) * t234;
t257 = pkin(9) * t233;
t256 = pkin(9) * t234;
t224 = t160 + t200;
t91 = -t109 + t220;
t249 = t146 * t224 + t149 * (-t145 * t91 + t61);
t228 = t109 + t221;
t247 = pkin(4) * t228;
t245 = t148 * t91 + t211;
t244 = qJ(6) * t225;
t241 = t146 * t228;
t237 = t149 * t228;
t202 = sin(pkin(10));
t203 = cos(pkin(10));
t161 = g(1) * t202 - g(2) * t203;
t192 = -g(3) + qJDD(1);
t235 = t142 * t192 + t143 * t161;
t232 = -t145 * t231 + t148 * t225;
t173 = pkin(4) * t146 - pkin(9) * t149;
t115 = t173 * qJD(2);
t141 = qJDD(2) * pkin(2);
t151 = qJD(2) ^ 2;
t162 = -g(1) * t203 - g(2) * t202;
t69 = -t147 * t162 + t150 * t235;
t154 = qJDD(3) - t69;
t59 = -qJ(3) * t151 - t141 + t154;
t152 = -qJDD(2) * pkin(8) + t59;
t94 = -t142 * t161 + t143 * t192;
t205 = t149 * t94;
t219 = qJD(4) ^ 2;
t30 = -t219 * pkin(4) + qJDD(4) * pkin(9) + t205 + (-qJD(2) * t115 + t152) * t146;
t169 = -t118 + t183;
t170 = -t117 + t182;
t136 = qJDD(2) * qJ(3);
t70 = t147 * t235 + t150 * t162;
t177 = 0.2e1 * qJD(3) * qJD(2) + t70;
t168 = t136 + t177;
t57 = t151 * t218 + t168;
t33 = pkin(4) * t170 + pkin(9) * t169 + t57;
t16 = t145 * t33 + t148 * t30;
t77 = pkin(5) * t112 - qJ(6) * t114;
t179 = -qJ(6) * t110 + t112 * t77 - t16;
t222 = -pkin(5) * (t220 + t230) + qJ(6) * t226 - t179;
t217 = pkin(5) * t148;
t38 = t146 * t94 - t149 * t152;
t29 = -qJDD(4) * pkin(4) - pkin(9) * t219 + t115 * t190 + t38;
t216 = t145 * t29;
t213 = t145 * t224;
t208 = t148 * t29;
t201 = qJ(6) * t148;
t198 = t131 * t145;
t197 = t131 * t148;
t138 = t146 ^ 2;
t196 = t138 * t151;
t139 = t149 ^ 2;
t195 = t139 * t151;
t184 = t149 * t151 * t146;
t122 = qJDD(4) + t184;
t194 = t146 * t122;
t123 = qJDD(4) - t184;
t193 = t149 * t123;
t191 = t138 + t139;
t189 = qJD(6) * t131;
t186 = t146 * t199;
t185 = t112 * t197;
t181 = -qJ(6) * t145 - pkin(4);
t15 = t145 * t30 - t148 * t33;
t6 = t145 * t15 + t148 * t16;
t89 = t114 * t198;
t176 = t149 * (t148 * t160 - t89) + t186;
t175 = t112 * t198 - t148 * t164;
t124 = 0.2e1 * t189;
t174 = t124 - t179;
t11 = -pkin(5) * t220 + t174;
t12 = -pkin(5) * t110 - qJ(6) * t220 + t114 * t77 + qJDD(6) + t15;
t172 = -pkin(5) * t12 + qJ(6) * t11;
t171 = -pkin(5) * t224 - qJ(6) * t48;
t5 = t145 * t16 - t148 * t15;
t39 = t146 * t152 + t205;
t19 = t146 * t39 - t149 * t38;
t166 = qJ(3) + t173;
t163 = (-t112 * t145 - t114 * t148) * t131;
t158 = t149 * (t89 - t185) + t146 * t110;
t157 = pkin(5) * t164 - t244 + t29;
t156 = t149 * (t145 * t164 + t185) - t186;
t155 = 0.2e1 * qJD(6) * t114 - t157;
t153 = pkin(5) * t227 + qJ(6) * t223 - t12;
t129 = -t195 - t219;
t128 = -t196 - t219;
t121 = t191 * t151;
t120 = t191 * qJDD(2);
t119 = t133 - 0.2e1 * t183;
t116 = 0.2e1 * t182 + t187;
t100 = (-qJDD(2) * t150 + t147 * t151) * t142;
t99 = (qJDD(2) * t147 + t150 * t151) * t142;
t88 = t129 * t149 - t194;
t87 = t128 * t146 + t193;
t86 = t143 * t94;
t58 = -t151 * pkin(2) + t168;
t55 = (qJD(5) + t131) * t112 + t167;
t50 = (-qJD(5) + t131) * t114 + t178;
t46 = t148 * t224;
t45 = t114 * t197 + t145 * t160;
t27 = t148 * t50 + t213;
t26 = -t148 * t48 + t213;
t25 = t145 * t50 - t46;
t24 = -t145 * t48 - t46;
t22 = t149 * t55 - t267;
t20 = t149 * t225 + t267;
t18 = t146 * t27 + t237;
t17 = t146 * t26 + t237;
t13 = (pkin(5) * t131 - 0.2e1 * qJD(6)) * t114 + t157;
t10 = (-t231 - t96) * pkin(5) + t155;
t9 = -pkin(5) * t96 + t155 + t244;
t8 = qJ(6) * t228 + t12;
t7 = (-t220 + t228) * pkin(5) + t174;
t4 = t146 * t6 - t149 * t29;
t3 = t11 * t148 + t12 * t145;
t2 = t11 * t145 - t12 * t148;
t1 = -t13 * t149 + t146 * t3;
t14 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t192, 0, 0, 0, 0, 0, 0, -t100, -t99, 0, t86 + (t147 * t70 + t150 * t69) * t142, 0, 0, 0, 0, 0, 0, 0, t100, t99, t86 + (t147 * t58 - t150 * t59) * t142, 0, 0, 0, 0, 0, 0, t143 * (-t123 * t146 + t128 * t149) + (t116 * t147 - t150 * t87) * t142, t143 * (-t122 * t149 - t129 * t146) + (t119 * t147 - t150 * t88) * t142, (t120 * t150 - t121 * t147) * t142, t143 * (t146 * t38 + t149 * t39) + (t147 * t57 - t150 * t19) * t142, 0, 0, 0, 0, 0, 0, t259, t143 * (-t146 * t55 - t265) + (-t150 * t22 - t266) * t142, t143 * (t149 * t27 - t241) + (t147 * t25 - t150 * t18) * t142, t143 * (t146 * t29 + t149 * t6) + (t147 * t5 - t150 * t4) * t142, 0, 0, 0, 0, 0, 0, t259, t143 * (t149 * t26 - t241) + (t147 * t24 - t150 * t17) * t142, t143 * (-t146 * t225 + t265) + (-t150 * t20 + t266) * t142, t143 * (t13 * t146 + t149 * t3) + (-t1 * t150 + t147 * t2) * t142; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t69, -t70, 0, 0, qJDD(2), 0, 0, 0, 0, 0, 0, -0.2e1 * t141 + t154, 0.2e1 * t136 + t177, -pkin(2) * t59 + qJ(3) * t58, -t169 * t149, -t116 * t149 - t119 * t146, t193 - t146 * (-t195 + t219), t170 * t146, t149 * (t196 - t219) - t194, 0, qJ(3) * t116 + t146 * t57 + t218 * t87, qJ(3) * t119 + t149 * t57 + t218 * t88, -qJ(3) * t121 - t120 * t218 - t19, qJ(3) * t57 + t19 * t218, t176, -t261, t249, t156, -t264, t158, t149 * (t216 - t256) - t146 * (t15 - t258) + t260, t149 * (t208 + t270) - t146 * (t16 + t271) - t268 + t218 * t22, -t149 * t5 + t166 * t25 + t18 * t218, t166 * t5 + t218 * t4, t176, t249, t261, t158, t264, t156, t149 * (-t10 * t145 - t201 * t231 - t256) - t146 * (-t153 - t258) + t260, t149 * (-pkin(9) * t24 - t145 * t7 + t148 * t8) - t146 * (-pkin(4) * t24 - t171) + qJ(3) * t24 + t218 * t17, t149 * (-pkin(5) * t214 + t148 * t9 - t270) - t146 * (-0.2e1 * t189 - t222 - t271) + t268 + t218 * t20, t149 * (-pkin(9) * t2 + (pkin(5) * t145 - t201) * t13) - t146 * (-pkin(4) * t2 - t172) + qJ(3) * t2 + t218 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t151, t59, 0, 0, 0, 0, 0, 0, t87, t88, -t120, t19, 0, 0, 0, 0, 0, 0, t250, t22, t18, t4, 0, 0, 0, 0, 0, 0, t250, t17, t20, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t184, (-t138 + t139) * t151, t133, -t184, -t187, qJDD(4), -t38, -t39, 0, 0, t45, t232, t245, t175, t263, t163, -pkin(4) * t231 - t208 + t257, pkin(4) * t55 + t216 - t269, pkin(9) * t27 + t247 + t6, -pkin(4) * t29 + pkin(9) * t6, t45, t245, -t232, t163, -t263, t175, t148 * t10 + t181 * t231 + t257, pkin(9) * t26 + t145 * t8 + t148 * t7 + t247, t269 + t145 * t9 + (pkin(4) + t217) * t225, pkin(9) * t3 + (t181 - t217) * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t199, t229, t224, -t199, -t48, t110, -t15, -t16, 0, 0, t199, t224, -t229, t110, t48, -t199, t153, t171, t124 + t222, t172; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t227, t224, t230, t12;];
tauJ_reg  = t14;
