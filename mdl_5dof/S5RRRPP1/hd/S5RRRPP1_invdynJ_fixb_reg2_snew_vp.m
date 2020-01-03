% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRRPP1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRRPP1_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP1_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP1_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP1_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP1_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP1_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:49:58
% EndTime: 2019-12-31 20:50:05
% DurationCPUTime: 2.50s
% Computational Cost: add. (7818->233), mult. (10601->289), div. (0->0), fcn. (6639->8), ass. (0->158)
t186 = sin(qJ(2));
t188 = cos(qJ(2));
t187 = cos(qJ(3));
t177 = qJD(1) + qJD(2);
t221 = qJD(3) * t177;
t212 = t187 * t221;
t176 = qJDD(1) + qJDD(2);
t185 = sin(qJ(3));
t227 = t185 * t176;
t156 = t212 + t227;
t183 = sin(pkin(8));
t184 = cos(pkin(8));
t213 = t185 * t221;
t224 = t187 * t176;
t199 = -t213 + t224;
t119 = t184 * t156 + t183 * t199;
t234 = t177 * t185;
t147 = -t184 * t187 * t177 + t183 * t234;
t223 = qJD(3) * t147;
t263 = t223 - t119;
t226 = t185 * t184;
t149 = (t187 * t183 + t226) * t177;
t236 = t149 * t147;
t252 = qJDD(3) + t236;
t233 = t183 * t252;
t146 = t149 ^ 2;
t189 = qJD(3) ^ 2;
t257 = -t146 - t189;
t65 = -t184 * t257 + t233;
t230 = t184 * t252;
t67 = t183 * t257 + t230;
t52 = t185 * t65 - t187 * t67;
t296 = pkin(1) * (t186 * t52 + t188 * t263);
t295 = pkin(2) * t263 + pkin(7) * t52;
t293 = pkin(3) * t65;
t292 = qJ(4) * t65;
t291 = qJ(4) * t67;
t251 = t147 ^ 2;
t133 = t251 - t189;
t48 = t185 * (-t184 * t133 + t233) - t187 * (t183 * t133 + t230);
t290 = t263 * qJ(5);
t253 = qJDD(3) - t236;
t232 = t183 * t253;
t254 = -t251 - t189;
t259 = t184 * t254 - t232;
t98 = t184 * t253;
t261 = t183 * t254 + t98;
t271 = -t185 * t261 + t187 * t259;
t203 = t183 * t156 - t184 * t199;
t222 = qJD(3) * t149;
t88 = t203 + t222;
t287 = -pkin(2) * t88 + pkin(7) * t271;
t90 = -t203 + t222;
t93 = t223 + t119;
t260 = t183 * t93 + t184 * t90;
t262 = t183 * t90 - t184 * t93;
t270 = -t185 * t262 + t187 * t260;
t86 = -t251 - t146;
t286 = -pkin(2) * t86 + pkin(7) * t270;
t243 = t183 * t263;
t28 = t185 * (t184 * t88 - t243) - t187 * (-t183 * t88 - t184 * t263);
t285 = pkin(1) * (t186 * t271 - t188 * t88);
t284 = pkin(1) * (t186 * t270 - t188 * t86);
t282 = pkin(3) * t262;
t280 = qJ(4) * t259;
t279 = qJ(4) * t261;
t278 = qJ(4) * t262;
t217 = qJD(4) * t149;
t273 = pkin(3) * t261 - 0.2e1 * t217;
t272 = -pkin(3) * t86 + qJ(4) * t260;
t134 = -t146 + t189;
t269 = t185 * (-t183 * t134 + t98) + t187 * (t184 * t134 + t232);
t248 = sin(qJ(1));
t249 = cos(qJ(1));
t197 = t249 * g(1) + t248 * g(2);
t162 = -qJD(1) ^ 2 * pkin(1) - t197;
t196 = t248 * g(1) - t249 * g(2);
t195 = qJDD(1) * pkin(1) + t196;
t122 = t188 * t162 + t186 * t195;
t175 = t177 ^ 2;
t112 = -t175 * pkin(2) + t176 * pkin(7) + t122;
t229 = t185 * t112;
t96 = t187 * g(3) + t229;
t97 = -t185 * g(3) + t187 * t112;
t63 = t185 * t96 + t187 * t97;
t256 = t146 - t251;
t235 = t175 * t185;
t190 = qJDD(3) * pkin(3) - t156 * qJ(4) - t229 + (pkin(3) * t235 + qJ(4) * t221 - g(3)) * t187;
t163 = qJD(3) * pkin(3) - qJ(4) * t234;
t181 = t187 ^ 2;
t171 = t181 * t175;
t71 = -pkin(3) * t171 + t199 * qJ(4) - qJD(3) * t163 + t97;
t38 = -0.2e1 * qJD(4) * t147 + t183 * t190 + t184 * t71;
t250 = 2 * qJD(5);
t247 = pkin(4) * t184;
t121 = -t186 * t162 + t188 * t195;
t111 = -t176 * pkin(2) - t175 * pkin(7) - t121;
t73 = -t199 * pkin(3) - qJ(4) * t171 + t163 * t234 + qJDD(4) + t111;
t245 = t183 * t73;
t241 = t184 * t73;
t208 = t183 * t71 - t184 * t190;
t37 = t208 + 0.2e1 * t217;
t18 = t183 * t38 - t184 * t37;
t239 = t185 * t18;
t238 = -pkin(2) * t111 + pkin(7) * t63;
t237 = qJ(5) * t184;
t167 = t187 * t235;
t228 = t185 * (qJDD(3) + t167);
t225 = t187 * (qJDD(3) - t167);
t220 = qJD(3) * t183;
t219 = qJD(3) * t184;
t180 = t185 ^ 2;
t170 = t180 * t175;
t130 = -t225 - t185 * (-t170 - t189);
t155 = 0.2e1 * t212 + t227;
t216 = -pkin(2) * t155 + pkin(7) * t130 + t185 * t111;
t129 = t187 * (-t171 - t189) - t228;
t157 = -0.2e1 * t213 + t224;
t215 = pkin(2) * t157 + pkin(7) * t129 - t187 * t111;
t101 = t147 * pkin(4) - t149 * qJ(5);
t201 = qJDD(3) * qJ(5) + qJD(3) * t250 - t147 * t101 + t38;
t24 = -t189 * pkin(4) + t201;
t198 = -qJDD(3) * pkin(4) - t189 * qJ(5) + qJDD(5) + t208;
t25 = (0.2e1 * qJD(4) + t101) * t149 + t198;
t10 = t183 * t24 - t184 * t25;
t11 = t183 * t25 + t184 * t24;
t209 = -qJ(5) * t183 - pkin(3);
t3 = -t185 * t10 + t187 * t11;
t193 = t203 * pkin(4) + t290 + t73;
t34 = (pkin(4) * qJD(3) - (2 * qJD(5))) * t149 + t193;
t214 = t185 * (-qJ(4) * t10 + (pkin(4) * t183 - t237) * t34) + t187 * (qJ(4) * t11 + (t209 - t247) * t34) - pkin(2) * t34 + pkin(7) * t3;
t20 = (-t189 - t86) * pkin(4) + t201;
t21 = -qJ(5) * t86 + t25;
t211 = t185 * (-t183 * t20 + t184 * t21 - t278) + t187 * (t183 * t21 + t184 * t20 + t272) + t286;
t19 = t183 * t37 + t184 * t38;
t210 = t185 * (-t18 - t278) + t187 * (t19 + t272) + t286;
t191 = t149 * t250 - t193;
t22 = -pkin(4) * t222 + t191 - t290;
t207 = t185 * (pkin(4) * t243 + t184 * t22 - t292) + t187 * (t291 + t183 * t22 - (pkin(3) + t247) * t263) - t295;
t23 = t191 + (-t88 - t222) * pkin(4);
t206 = t185 * (-t183 * t23 - t88 * t237 - t279) + t187 * (t184 * t23 + t209 * t88 + t280) + t287;
t205 = t185 * (t245 - t279) + t187 * (-pkin(3) * t88 - t241 + t280) + t287;
t204 = t185 * (t241 + t292) + t187 * (pkin(3) * t263 + t245 - t291) + t295;
t160 = (t180 + t181) * t176;
t161 = t170 + t171;
t202 = pkin(2) * t161 + pkin(7) * t160 + t63;
t7 = t187 * t19 - t239;
t200 = pkin(7) * t7 - qJ(4) * t239 + t187 * (-pkin(3) * t73 + qJ(4) * t19) - pkin(2) * t73;
t194 = t185 * (t147 * t219 + t183 * t203) + t187 * (t147 * t220 - t184 * t203);
t132 = t149 * t220;
t192 = t185 * t132 + (-t147 * t226 + t187 * (-t147 * t183 - t149 * t184)) * qJD(3);
t128 = t228 + t187 * (-t170 + t189);
t127 = t185 * (t171 - t189) + t225;
t124 = (t156 + t212) * t185;
t123 = t157 * t187;
t120 = t187 * t155 + t185 * t157;
t55 = t185 * (t184 * t119 - t132) + t187 * (t183 * t119 + t149 * t219);
t1 = [0, 0, 0, 0, 0, qJDD(1), t196, t197, 0, 0, 0, 0, 0, 0, 0, t176, pkin(1) * (-t186 * t175 + t188 * t176) + t121, pkin(1) * (-t188 * t175 - t186 * t176) - t122, 0, pkin(1) * (t188 * t121 + t186 * t122), t124, t120, t128, t123, t127, 0, pkin(1) * (t186 * t129 + t188 * t157) + t215, pkin(1) * (t186 * t130 - t188 * t155) + t216, pkin(1) * (t186 * t160 + t188 * t161) + t202, pkin(1) * (-t188 * t111 + t186 * t63) + t238, t55, -t28, t269, t194, -t48, t192, t285 + t205, t204 + t296, t284 + t210, pkin(1) * (t186 * t7 - t188 * t73) + t200, t55, t269, t28, t192, t48, t194, t285 + t206, t284 + t211, t207 - t296, pkin(1) * (t186 * t3 - t188 * t34) + t214; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t176, t121, -t122, 0, 0, t124, t120, t128, t123, t127, 0, t215, t216, t202, t238, t55, -t28, t269, t194, -t48, t192, t205, t204, t210, t200, t55, t269, t28, t192, t48, t194, t206, t211, t207, t214; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t167, t170 - t171, t227, t167, t224, qJDD(3), -t96, -t97, 0, 0, t236, t256, t93, -t236, t90, qJDD(3), -t208 + t273, -t38 - t293, t282, pkin(3) * t18, t236, t93, -t256, qJDD(3), -t90, -t236, pkin(4) * t253 + qJ(5) * t254 - t149 * t101 - t198 + t273, -pkin(4) * t93 + qJ(5) * t90 + t282, t293 + qJ(5) * t252 + (-t257 - t189) * pkin(4) + t201, pkin(3) * t10 - pkin(4) * t25 + qJ(5) * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88, -t263, t86, t73, 0, 0, 0, 0, 0, 0, t88, t86, t263, t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t253, t93, t257, t25;];
tauJ_reg = t1;
