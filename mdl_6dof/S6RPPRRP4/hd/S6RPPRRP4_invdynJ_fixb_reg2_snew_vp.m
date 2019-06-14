% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPPRRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 14:56
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPPRRP4_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP4_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP4_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP4_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP4_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP4_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:55:03
% EndTime: 2019-05-05 14:55:11
% DurationCPUTime: 2.90s
% Computational Cost: add. (6750->290), mult. (12277->352), div. (0->0), fcn. (6600->8), ass. (0->192)
t162 = cos(qJ(5));
t159 = sin(qJ(5));
t160 = sin(qJ(4));
t201 = qJD(1) * qJD(4);
t142 = t160 * t201;
t163 = cos(qJ(4));
t199 = t163 * qJDD(1);
t127 = t142 - t199;
t119 = qJDD(5) - t127;
t203 = qJD(1) * t160;
t120 = qJD(4) * t162 + t159 * t203;
t122 = -qJD(4) * t159 + t162 * t203;
t212 = t120 * t122;
t237 = t119 + t212;
t223 = t159 * t237;
t118 = t122 ^ 2;
t140 = qJD(1) * t163 + qJD(5);
t231 = t140 ^ 2;
t241 = -t118 - t231;
t44 = -t162 * t241 + t223;
t283 = pkin(3) * t44;
t282 = pkin(4) * t44;
t281 = pkin(8) * t44;
t218 = t162 * t237;
t46 = t159 * t241 + t218;
t280 = pkin(8) * t46;
t155 = sin(pkin(9));
t279 = t155 * t44;
t156 = cos(pkin(9));
t278 = t156 * t44;
t277 = t160 * t46;
t276 = t163 * t46;
t232 = t120 ^ 2;
t102 = t232 - t231;
t109 = t140 * t122;
t196 = t163 * t201;
t200 = t160 * qJDD(1);
t125 = -t196 - t200;
t190 = qJDD(4) * t162 - t159 * t125;
t82 = qJD(5) * t122 + t190;
t60 = t109 - t82;
t275 = t160 * (-t102 * t162 + t223) - t163 * t60;
t229 = pkin(1) + pkin(2);
t238 = t119 - t212;
t217 = t162 * t238;
t234 = -t231 - t232;
t245 = t159 * t234 + t217;
t242 = -t109 - t82;
t222 = t159 * t238;
t244 = t162 * t234 - t222;
t258 = t160 * t242 + t163 * t244;
t271 = t155 * t258 - t156 * t245;
t274 = -t229 * t271 + pkin(3) * t245 + qJ(2) * (t155 * t245 + t156 * t258) - pkin(7) * t258;
t273 = t102 * t159 + t218;
t181 = -t159 * qJDD(4) - t162 * t125;
t172 = qJD(5) * t120 - t181;
t211 = t140 * t120;
t235 = t211 + t172;
t225 = t159 * t235;
t240 = t118 - t232;
t270 = t160 * (t162 * t242 + t225) + t163 * t240;
t267 = pkin(4) * t245;
t266 = pkin(8) * t244;
t265 = pkin(8) * t245;
t259 = t160 * t244 - t163 * t242;
t103 = -t118 + t231;
t236 = -t211 + t172;
t257 = t163 * t236 - t160 * (-t103 * t159 + t217);
t239 = t118 + t232;
t256 = pkin(4) * t239;
t255 = t103 * t162 + t222;
t254 = qJ(6) * t235;
t252 = t160 * t239;
t248 = t163 * t239;
t243 = -t159 * t242 + t162 * t235;
t182 = -t127 - t142;
t183 = -t125 + t196;
t165 = qJD(1) ^ 2;
t151 = qJDD(1) * qJ(2);
t161 = sin(qJ(1));
t164 = cos(qJ(1));
t184 = t164 * g(1) + t161 * g(2);
t180 = 0.2e1 * qJD(2) * qJD(1) - t184;
t175 = t151 + t180;
t101 = -t165 * t229 + t175;
t204 = g(1) * t161 - g(2) * t164;
t195 = -qJDD(2) + t204;
t179 = -t165 * qJ(2) - t195;
t169 = -qJDD(1) * t229 + t179;
t192 = t155 * t101 - t156 * t169;
t56 = qJDD(1) * pkin(3) - pkin(7) * t165 + t192;
t40 = pkin(4) * t182 + pkin(8) * t183 + t56;
t187 = pkin(4) * t163 + pkin(8) * t160;
t176 = t165 * t187;
t230 = qJD(4) ^ 2;
t205 = g(3) + qJDD(3);
t68 = t101 * t156 + t155 * t169;
t57 = -pkin(3) * t165 - qJDD(1) * pkin(7) + t68;
t53 = t160 * t205 + t163 * t57;
t43 = -pkin(4) * t230 + qJDD(4) * pkin(8) - t163 * t176 + t53;
t25 = t159 * t40 + t162 * t43;
t85 = -pkin(5) * t120 + qJ(6) * t122;
t191 = qJ(6) * t119 + t120 * t85 + t25;
t233 = -(t231 + t241) * pkin(5) + qJ(6) * t237 + t191;
t228 = pkin(5) * t162;
t52 = t160 * t57 - t163 * t205;
t42 = -qJDD(4) * pkin(4) - pkin(8) * t230 - t160 * t176 + t52;
t227 = t159 * t42;
t224 = t159 * t236;
t221 = t162 * t42;
t219 = t162 * t236;
t214 = qJ(6) * t162;
t213 = qJDD(1) * pkin(1);
t210 = t140 * t159;
t209 = t140 * t162;
t139 = t163 * t165 * t160;
t132 = qJDD(4) + t139;
t207 = t160 * t132;
t133 = qJDD(4) - t139;
t206 = t163 * t133;
t202 = qJD(6) * t140;
t198 = t120 * t209;
t197 = t163 * t212;
t194 = -qJ(6) * t159 - pkin(4);
t193 = qJ(2) * t156 - pkin(7);
t24 = t159 * t43 - t162 * t40;
t8 = t159 * t24 + t162 * t25;
t189 = -t120 * t210 + t162 * t82;
t134 = 0.2e1 * t202;
t188 = t134 + t191;
t15 = -pkin(5) * t231 + t188;
t16 = -pkin(5) * t119 - qJ(6) * t231 - t122 * t85 + qJDD(6) + t24;
t186 = -pkin(5) * t16 + qJ(6) * t15;
t185 = -pkin(5) * t236 - qJ(6) * t60;
t7 = t159 * t25 - t162 * t24;
t30 = t160 * t52 + t163 * t53;
t100 = t122 * t210;
t178 = -t160 * (t162 * t172 + t100) + t197;
t174 = (t120 * t159 + t122 * t162) * t140;
t173 = qJ(2) * t155 + pkin(3) + t187;
t171 = -t82 * pkin(5) - t254 + t42;
t170 = t163 * t119 - t160 * (-t100 + t198);
t168 = -0.2e1 * qJD(6) * t122 - t171;
t167 = -t197 - t160 * (-t159 * t82 - t198);
t166 = pkin(5) * t238 + qJ(6) * t234 - t16;
t153 = t163 ^ 2;
t152 = t160 ^ 2;
t147 = t153 * t165;
t146 = t152 * t165;
t137 = -t147 - t230;
t136 = -t146 - t230;
t131 = t146 + t147;
t130 = (-t152 - t153) * qJDD(1);
t129 = qJDD(1) * t156 + t155 * t165;
t128 = -qJDD(1) * t155 + t156 * t165;
t126 = 0.2e1 * t142 - t199;
t124 = 0.2e1 * t196 + t200;
t110 = -t179 + t213;
t97 = -t136 * t160 - t206;
t96 = t137 * t163 - t207;
t86 = t130 * t155 + t131 * t156;
t70 = t124 * t156 + t155 * t97;
t69 = t126 * t156 + t155 * t96;
t66 = (-qJD(5) - t140) * t120 + t181;
t61 = (qJD(5) - t140) * t122 + t190;
t54 = -t122 * t209 + t159 * t172;
t37 = t155 * t68 - t156 * t192;
t36 = -t162 * t60 + t224;
t35 = t162 * t61 + t224;
t34 = -t159 * t60 - t219;
t33 = t159 * t61 - t219;
t31 = -t160 * t66 - t276;
t28 = -t160 * t235 + t276;
t27 = t163 * t36 - t252;
t26 = t163 * t35 - t252;
t22 = t155 * t30 - t156 * t56;
t21 = (-pkin(5) * t140 + 0.2e1 * qJD(6)) * t122 + t171;
t19 = t155 * t31 + t278;
t17 = t155 * t28 - t278;
t14 = (-t242 + t109) * pkin(5) + t168;
t13 = pkin(5) * t109 + t168 + t254;
t12 = t155 * t27 - t156 * t34;
t11 = t155 * t26 - t156 * t33;
t10 = qJ(6) * t239 + t16;
t9 = (-t231 + t239) * pkin(5) + t188;
t6 = t160 * t42 + t163 * t8;
t5 = t15 * t162 + t159 * t16;
t4 = t15 * t159 - t16 * t162;
t3 = t160 * t21 + t163 * t5;
t2 = t155 * t6 - t156 * t7;
t1 = t155 * t3 - t156 * t4;
t18 = [0, 0, 0, 0, 0, qJDD(1), t204, t184, 0, 0, 0, 0, 0, qJDD(1), 0, 0, t195 + 0.2e1 * t213, 0, 0.2e1 * t151 + t180, qJ(2) * (-pkin(1) * t165 + t175) + pkin(1) * t110, 0, 0, 0, 0, 0, qJDD(1), -qJ(2) * t128 + t129 * t229 + t192, qJ(2) * t129 + t128 * t229 + t68, 0, qJ(2) * (t155 * t192 + t156 * t68) - t229 * t37, t183 * t160, t124 * t163 - t126 * t160, -t163 * (-t146 + t230) - t207, t182 * t163, -t206 - t160 * (t147 - t230), 0, qJ(2) * (-t126 * t155 + t156 * t96) + t163 * t56 - pkin(3) * t126 - pkin(7) * t96 - t229 * t69, qJ(2) * (-t124 * t155 + t156 * t97) - pkin(3) * t124 - pkin(7) * t97 - t160 * t56 - t229 * t70, qJ(2) * (t130 * t156 - t131 * t155) - pkin(3) * t131 - pkin(7) * t130 - t229 * t86 - t30, qJ(2) * (t155 * t56 + t156 * t30) + pkin(3) * t56 - pkin(7) * t30 - t229 * t22, t178, t270, t257, t167, t275, t170, -t163 * (t24 - t267) - t160 * (t227 - t265) + t274, qJ(2) * (t156 * t31 - t279) - t163 * (t25 + t282) - t283 - pkin(7) * t31 - t160 * (t221 + t281) - t229 * t19, -t11 * t229 + t160 * t7 + t173 * t33 + t193 * t26, t173 * t7 + t193 * t6 - t2 * t229, t178, t257, -t270, t170, -t275, t167, -t163 * (-t166 - t267) - t160 * (-t14 * t159 - t214 * t242 - t265) + t274, qJ(2) * (t155 * t34 + t156 * t27) - t163 * (-pkin(4) * t34 - t185) + pkin(3) * t34 - pkin(7) * t27 - t160 * (-pkin(8) * t34 + t10 * t162 - t159 * t9) - t229 * t12, qJ(2) * (t156 * t28 + t279) - t163 * (-0.2e1 * t202 - t233 - t282) + t283 - pkin(7) * t28 - t160 * (-pkin(5) * t225 + t13 * t162 - t281) - t229 * t17, qJ(2) * (t155 * t4 + t156 * t3) - t163 * (-pkin(4) * t4 - t186) + pkin(3) * t4 - pkin(7) * t3 - t160 * (-pkin(8) * t4 + (pkin(5) * t159 - t214) * t21) - t229 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), 0, -t165, -t110, 0, 0, 0, 0, 0, 0, -t129, -t128, 0, t37, 0, 0, 0, 0, 0, 0, t69, t70, t86, t22, 0, 0, 0, 0, 0, 0, t271, t19, t11, t2, 0, 0, 0, 0, 0, 0, t271, t12, t17, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t205, 0, 0, 0, 0, 0, 0, t132 * t163 + t137 * t160, -t133 * t160 + t136 * t163, 0, t160 * t53 - t163 * t52, 0, 0, 0, 0, 0, 0, t259, t163 * t66 - t277, t160 * t35 + t248, t160 * t8 - t163 * t42, 0, 0, 0, 0, 0, 0, t259, t160 * t36 + t248, t163 * t235 + t277, t160 * t5 - t163 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t139, t146 - t147, -t200, t139, -t199, qJDD(4), -t52, -t53, 0, 0, t54, t243, t255, t189, t273, t174, -pkin(4) * t242 - t221 + t266, pkin(4) * t66 + t227 - t280, pkin(8) * t35 + t256 + t8, -pkin(4) * t42 + pkin(8) * t8, t54, t255, -t243, t174, -t273, t189, t162 * t14 + t194 * t242 + t266, pkin(8) * t36 + t10 * t159 + t162 * t9 + t256, t280 + t159 * t13 + (pkin(4) + t228) * t235, pkin(8) * t5 + (t194 - t228) * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t212, t240, t236, -t212, -t60, t119, -t24, -t25, 0, 0, t212, t236, -t240, t119, t60, -t212, t166, t185, t134 + t233, t186; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t238, t236, t241, t16;];
tauJ_reg  = t18;
