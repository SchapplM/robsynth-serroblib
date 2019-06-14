% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPRPRP10
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 18:12
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPRPRP10_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP10_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP10_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP10_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP10_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPRP10_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 18:11:21
% EndTime: 2019-05-05 18:11:29
% DurationCPUTime: 2.81s
% Computational Cost: add. (5218->282), mult. (10496->295), div. (0->0), fcn. (5688->6), ass. (0->201)
t158 = cos(qJ(5));
t155 = sin(qJ(5));
t156 = sin(qJ(3));
t213 = qJD(1) * qJD(3);
t145 = t156 * t213;
t159 = cos(qJ(3));
t146 = t159 * qJDD(1);
t123 = t146 - t145;
t114 = qJDD(5) + t123;
t216 = qJD(1) * t156;
t115 = t155 * qJD(3) - t158 * t216;
t117 = t158 * qJD(3) + t155 * t216;
t223 = t117 * t115;
t68 = -t223 - t114;
t239 = t155 * t68;
t113 = t117 ^ 2;
t215 = t159 * qJD(1);
t141 = qJD(5) + t215;
t252 = t141 ^ 2;
t258 = -t113 - t252;
t38 = t158 * t258 + t239;
t293 = pkin(4) * t38;
t292 = t159 * t38;
t249 = pkin(3) + pkin(8);
t291 = t249 * t38;
t219 = t159 * qJ(4);
t181 = t156 * t249 + qJ(2) - t219;
t231 = t158 * t68;
t290 = t181 * (-t155 * t258 + t231);
t102 = t141 * t117;
t203 = t159 * t213;
t210 = t156 * qJDD(1);
t122 = t203 + t210;
t200 = t155 * qJDD(3) - t158 * t122;
t182 = t117 * qJD(5) + t200;
t50 = -t102 + t182;
t253 = t115 ^ 2;
t96 = t253 - t252;
t289 = t156 * (-t155 * t96 + t231) + t159 * t50;
t286 = -t158 * t96 - t239;
t224 = t115 * t141;
t189 = t158 * qJDD(3) + t155 * t122;
t71 = -t115 * qJD(5) + t189;
t259 = t71 - t224;
t233 = t158 * t259;
t49 = t102 + t182;
t81 = t113 - t253;
t285 = t156 * (-t155 * t49 + t233) + t159 * t81;
t62 = -t253 - t113;
t284 = pkin(4) * t62;
t255 = -t223 + t114;
t230 = t158 * t255;
t254 = -t252 - t253;
t264 = t155 * t254 + t230;
t283 = pkin(4) * t264;
t282 = qJ(4) * t62;
t280 = t156 * t62;
t277 = t159 * t264;
t276 = t249 * t264;
t97 = -t113 + t252;
t275 = -t155 * t97 + t230;
t238 = t155 * t255;
t260 = t224 + t71;
t273 = t159 * t260 - t156 * (-t158 * t97 - t238);
t250 = pkin(1) + pkin(7);
t153 = t156 ^ 2;
t162 = qJD(1) ^ 2;
t147 = t153 * t162;
t161 = qJD(3) ^ 2;
t134 = -t147 - t161;
t217 = t159 * t162;
t207 = t156 * t217;
t129 = -qJDD(3) + t207;
t218 = t159 * t129;
t88 = -t156 * t134 + t218;
t272 = t250 * t88;
t154 = t159 ^ 2;
t148 = t154 * t162;
t136 = -t148 - t161;
t130 = qJDD(3) + t207;
t220 = t156 * t130;
t89 = -t159 * t136 + t220;
t271 = t250 * t89;
t267 = t259 * qJ(6);
t266 = t123 - t145;
t152 = qJDD(1) * qJ(2);
t157 = sin(qJ(1));
t160 = cos(qJ(1));
t196 = t160 * g(1) + t157 * g(2);
t188 = -t152 + t196;
t175 = -pkin(3) * t203 + t266 * qJ(4) + t188;
t209 = t250 * t162;
t265 = t122 * pkin(3) - t175 - t209;
t263 = t155 * t259 + t158 * t49;
t262 = t181 * (t158 * t254 - t238);
t261 = qJ(6) * t155 + pkin(4);
t212 = qJD(2) * qJD(1);
t150 = 0.2e1 * t212;
t257 = -0.2e1 * qJD(4) * t215 + t150;
t193 = t156 * pkin(3) - t219;
t118 = t193 * qJD(1);
t256 = qJDD(3) * qJ(4) - t118 * t216;
t251 = 2 * qJD(6);
t248 = t182 * pkin(5);
t247 = pkin(5) * t155;
t246 = t156 * g(3);
t245 = t161 * pkin(3);
t195 = t157 * g(1) - t160 * g(2);
t185 = qJDD(2) - t195;
t178 = -t162 * qJ(2) + t185;
t172 = -t250 * qJDD(1) + t178;
t170 = t159 * t172;
t169 = t118 * t215 + qJDD(4) - t170;
t167 = t161 * qJ(4) - t169;
t163 = t123 * pkin(4) - t249 * qJDD(3) + (pkin(4) * t213 + pkin(8) * t217 - g(3)) * t156 - t167;
t187 = pkin(4) * t215 - qJD(3) * pkin(8);
t30 = -t187 * t215 + t249 * t122 + (-t153 * pkin(4) - t250) * t162 - t175 + t257;
t17 = t155 * t163 + t158 * t30;
t76 = t115 * pkin(5) - t117 * qJ(6);
t197 = t114 * qJ(6) - t115 * t76 + t141 * t251 + t17;
t12 = -pkin(5) * t252 + t197;
t16 = t155 * t30 - t158 * t163;
t14 = -t114 * pkin(5) - qJ(6) * t252 + t117 * t76 + qJDD(6) + t16;
t244 = -pkin(5) * t14 + qJ(6) * t12;
t243 = -pkin(5) * t260 - qJ(6) * t50;
t92 = t156 * t172;
t78 = t159 * g(3) - t92;
t179 = -t78 + t256;
t176 = -t179 + t245;
t211 = qJD(4) * qJD(3);
t33 = t122 * pkin(4) + pkin(8) * t147 - qJD(3) * t187 + t176 - 0.2e1 * t211;
t242 = t155 * t33;
t240 = t155 * t260;
t165 = -pkin(5) * t102 + t117 * t251 + t33;
t164 = t165 + t267;
t15 = t164 - t248;
t236 = t156 * t15;
t235 = t158 * t33;
t232 = t158 * t260;
t51 = (-qJD(5) + t141) * t117 - t200;
t25 = t155 * t51 - t232;
t228 = t159 * t25;
t225 = qJDD(1) * pkin(1);
t222 = t141 * t155;
t221 = t141 * t158;
t214 = qJD(5) + t141;
t208 = t115 * t222;
t206 = t115 * t221;
t205 = t159 * t223;
t201 = -qJ(6) * t158 + qJ(4);
t126 = (t153 + t154) * qJDD(1);
t128 = -t148 - t147;
t199 = qJ(2) * t128 + t250 * t126;
t93 = t117 * t222;
t198 = t93 - t206;
t5 = t155 * t17 - t158 * t16;
t194 = t156 * t33 + t159 * t5;
t192 = t155 * t16 + t158 * t17;
t77 = t170 + t246;
t43 = -t156 * t78 + t159 * t77;
t191 = t156 * (-t148 + t161) + t218;
t190 = -t159 * (t147 - t161) + t220;
t186 = qJ(2) + t193;
t94 = t117 * t221;
t184 = -t156 * (-t155 * t71 - t94) + t205;
t183 = t155 * t182 + t206;
t180 = t159 * t114 - t156 * (t94 + t208);
t177 = -pkin(5) * t258 - qJ(6) * t68 + t12;
t174 = pkin(5) * t255 + qJ(6) * t254 - t14;
t173 = -t205 - t156 * (t158 * t182 - t208);
t168 = t257 + t265;
t166 = -qJDD(3) * pkin(3) - t167;
t149 = 0.2e1 * t211;
t127 = t148 - t147;
t124 = t146 - 0.2e1 * t145;
t121 = 0.2e1 * t203 + t210;
t103 = -t178 + t225;
t95 = t188 + t209 - 0.2e1 * t212;
t91 = t266 * t159;
t86 = (t122 + t203) * t156;
t74 = -t159 * t121 - t156 * t124;
t58 = -t166 + t246;
t57 = t149 - t176;
t54 = -t115 * t214 + t189;
t48 = t117 * t214 + t200;
t45 = t158 * t71 - t93;
t28 = t156 * t57 + t159 * t58;
t24 = -t155 * t50 - t232;
t23 = t156 * t54 - t292;
t22 = t156 * t48 - t277;
t21 = -t156 * t259 + t292;
t20 = t156 * t49 - t277;
t19 = -t159 * t24 + t280;
t18 = -t228 + t280;
t10 = t164 + (-t48 - t182) * pkin(5);
t9 = t165 - t248 + 0.2e1 * t267;
t8 = -qJ(6) * t62 + t14;
t7 = (-t252 - t62) * pkin(5) + t197;
t2 = t155 * t12 - t158 * t14;
t1 = -t159 * t2 - t236;
t3 = [0, 0, 0, 0, 0, qJDD(1), t195, t196, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t185 - 0.2e1 * t225, t150 + 0.2e1 * t152 - t196, pkin(1) * t103 + qJ(2) * (-t162 * pkin(1) + t150 - t188), t91, t74, -t191, t86, -t190, 0, qJ(2) * t121 - t156 * t95 + t272, qJ(2) * t124 - t159 * t95 + t271, t199 - t43, -qJ(2) * t95 - t250 * t43, 0, t191, t190, t91, t74, t86, t159 * (-qJ(4) * t128 + t166) - t156 * (-pkin(3) * t128 + t149 - t245 + t256 + t92) + t199, -t186 * t121 - t156 * t168 - t272, t159 * (0.2e1 * (qJD(4) * t159 - qJD(2)) * qJD(1) - t265) - t271 - t186 * t124, t186 * t168 - t250 * t28, t184, t285, t273, t173, -t289, t180, t159 * (-t16 + t283) - t156 * (pkin(4) * t49 - t235) + t262 - t250 * t20, t159 * (-t17 + t293) - t156 * (pkin(4) * t54 + t242) + t290 - t250 * t23, pkin(4) * t228 - t156 * (-t192 + t284) + t181 * (t158 * t51 + t240) - t250 * t18, t181 * t192 + (pkin(4) + t250) * t194, t184, t273, -t285, t180, t289, t173, t159 * (t174 + t283) - t156 * (-t158 * t10 + t261 * t48) + t262 - t250 * t22, t159 * (pkin(4) * t24 + t243) - t156 * (-t155 * t8 - t158 * t7 + t284) + t181 * (-t158 * t50 + t240) - t250 * t19, t159 * (t177 - t293) - t156 * (-pkin(4) * t259 - pkin(5) * t233 - t155 * t9) - t290 - t250 * t21, t159 * (pkin(4) * t2 + t244) + t181 * (t158 * t12 + t155 * t14) - (-pkin(5) * t158 - t261) * t236 - t250 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t162, -t103, 0, 0, 0, 0, 0, 0, -t88, -t89, -t126, t43, 0, 0, 0, 0, 0, 0, -t126, t88, t89, t28, 0, 0, 0, 0, 0, 0, t20, t23, t18, -t194, 0, 0, 0, 0, 0, 0, t22, t19, t21, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t207, t127, t146, -t207, -t210, qJDD(3), t77, t78, 0, 0, qJDD(3), -t146, t210, t207, t127, -t207, (-pkin(3) * t159 - qJ(4) * t156) * qJDD(1), -t246 + (-t161 - t134) * qJ(4) + (-qJDD(3) + t129) * pkin(3) + t169, qJ(4) * t130 + t149 + (-t136 - t161) * pkin(3) + t179, pkin(3) * t58 + qJ(4) * t57, t45, -t263, t275, t183, -t286, t198, qJ(4) * t49 - t242 - t276, qJ(4) * t54 - t235 - t291, -t249 * t25 + t282 - t5, -qJ(4) * t33 - t249 * t5, t45, t275, t263, t198, t286, t183, -t155 * t10 + t201 * t48 - t276, -t155 * t7 + t158 * t8 - t249 * t24 + t282, t158 * t9 + (-qJ(4) - t247) * t259 + t291, -t249 * t2 + (-t201 - t247) * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t146, -t129, t136, -t58, 0, 0, 0, 0, 0, 0, t264, t38, t25, t5, 0, 0, 0, 0, 0, 0, t264, t24, -t38, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t223, t81, t260, -t223, -t50, t114, -t16, -t17, 0, 0, t223, t260, -t81, t114, t50, -t223, t174, t243, t177, t244; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t255, t260, t258, t14;];
tauJ_reg  = t3;
