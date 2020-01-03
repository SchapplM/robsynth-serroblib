% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRPR7_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR7_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR7_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR7_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR7_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR7_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR7_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:16:07
% EndTime: 2019-12-31 21:16:20
% DurationCPUTime: 6.03s
% Computational Cost: add. (7011->491), mult. (17728->703), div. (0->0), fcn. (12540->8), ass. (0->235)
t243 = sin(qJ(3));
t246 = cos(qJ(2));
t316 = cos(qJ(3));
t266 = t316 * t246;
t244 = sin(qJ(2));
t278 = qJD(1) * t244;
t203 = -qJD(1) * t266 + t243 * t278;
t240 = sin(pkin(9));
t241 = cos(pkin(9));
t242 = sin(qJ(5));
t245 = cos(qJ(5));
t250 = t240 * t242 - t241 * t245;
t140 = t250 * t203;
t199 = t250 * qJD(5);
t354 = t140 + t199;
t215 = t243 * t246 + t244 * t316;
t204 = t215 * qJD(1);
t239 = qJD(2) + qJD(3);
t185 = t204 * t241 + t239 * t240;
t236 = -pkin(2) * t246 - pkin(1);
t226 = qJD(1) * t236;
t141 = t203 * pkin(3) - t204 * qJ(4) + t226;
t329 = -pkin(7) - pkin(6);
t228 = t329 * t246;
t218 = qJD(1) * t228;
t206 = t316 * t218;
t227 = t329 * t244;
t217 = qJD(1) * t227;
t210 = qJD(2) * pkin(2) + t217;
t175 = t243 * t210 - t206;
t152 = t239 * qJ(4) + t175;
t76 = t241 * t141 - t152 * t240;
t49 = pkin(4) * t203 - pkin(8) * t185 + t76;
t261 = -t204 * t240 + t241 * t239;
t77 = t240 * t141 + t241 * t152;
t54 = pkin(8) * t261 + t77;
t13 = -t242 * t54 + t245 * t49;
t353 = t13 * mrSges(6,1);
t14 = t242 * t49 + t245 * t54;
t352 = t14 * mrSges(6,2);
t213 = t240 * t245 + t241 * t242;
t139 = t213 * t203;
t200 = t213 * qJD(5);
t351 = t139 + t200;
t350 = -t185 * t242 + t245 * t261;
t115 = t185 * t245 + t242 * t261;
t248 = -t243 * t244 + t266;
t178 = t239 * t248;
t167 = t178 * qJD(1);
t43 = qJD(5) * t350 - t167 * t250;
t333 = t43 / 0.2e1;
t44 = -qJD(5) * t115 - t167 * t213;
t332 = t44 / 0.2e1;
t179 = t239 * t215;
t168 = t179 * qJD(1);
t324 = t168 / 0.2e1;
t315 = pkin(2) * t243;
t232 = qJ(4) + t315;
t208 = (-pkin(8) - t232) * t240;
t238 = t241 * pkin(8);
t281 = t232 * t241;
t209 = t238 + t281;
t162 = t208 * t245 - t209 * t242;
t265 = qJD(3) * t316;
t260 = pkin(2) * t265;
t229 = t260 + qJD(4);
t283 = t203 * t241;
t257 = t204 * pkin(4) + pkin(8) * t283;
t169 = pkin(3) * t204 + qJ(4) * t203;
t271 = pkin(2) * t278;
t146 = t169 + t271;
t205 = t243 * t218;
t177 = t217 * t316 + t205;
t89 = t241 * t146 - t177 * t240;
t57 = t257 + t89;
t284 = t203 * t240;
t273 = pkin(8) * t284;
t90 = t240 * t146 + t241 * t177;
t66 = t273 + t90;
t349 = qJD(5) * t162 - t229 * t250 - t242 * t57 - t245 * t66;
t163 = t208 * t242 + t209 * t245;
t348 = -qJD(5) * t163 - t213 * t229 + t242 * t66 - t245 * t57;
t222 = (-pkin(8) - qJ(4)) * t240;
t289 = qJ(4) * t241;
t223 = t238 + t289;
t182 = t222 * t242 + t223 * t245;
t174 = t210 * t316 + t205;
t94 = t241 * t169 - t174 * t240;
t60 = t257 + t94;
t95 = t240 * t169 + t241 * t174;
t70 = t273 + t95;
t347 = -qJD(4) * t213 - qJD(5) * t182 + t242 * t70 - t245 * t60;
t181 = t222 * t245 - t223 * t242;
t346 = -qJD(4) * t250 + qJD(5) * t181 - t242 * t60 - t245 * t70;
t145 = -t239 * pkin(3) + qJD(4) - t174;
t256 = mrSges(5,1) * t240 + mrSges(5,2) * t241;
t345 = t145 * t256;
t298 = t204 * mrSges(4,3);
t342 = mrSges(4,1) * t239 + mrSges(5,1) * t261 - mrSges(5,2) * t185 - t298;
t341 = t316 * t227 + t243 * t228;
t277 = qJD(1) * t246;
t290 = Ifges(3,6) * qJD(2);
t309 = Ifges(3,4) * t244;
t340 = pkin(6) * (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t277) + t290 / 0.2e1 + (t246 * Ifges(3,2) + t309) * qJD(1) / 0.2e1;
t274 = qJD(1) * qJD(2);
t264 = t244 * t274;
t65 = pkin(2) * t264 + pkin(3) * t168 - qJ(4) * t167 - qJD(4) * t204;
t267 = qJD(2) * t329;
t259 = qJD(1) * t267;
t211 = t244 * t259;
t249 = t246 * t259;
t275 = qJD(3) * t243;
t101 = t210 * t265 + t316 * t211 + t218 * t275 + t243 * t249;
t91 = qJD(4) * t239 + t101;
t28 = -t240 * t91 + t241 * t65;
t286 = t167 * t241;
t12 = pkin(4) * t168 - pkin(8) * t286 + t28;
t287 = t167 * t240;
t29 = t240 * t65 + t241 * t91;
t15 = -pkin(8) * t287 + t29;
t2 = qJD(5) * t13 + t12 * t242 + t15 * t245;
t3 = -qJD(5) * t14 + t12 * t245 - t15 * t242;
t339 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t43 + Ifges(6,6) * t44;
t102 = t175 * qJD(3) + t243 * t211 - t316 * t249;
t103 = -pkin(4) * t261 + t145;
t294 = t239 * Ifges(4,6);
t296 = t204 * Ifges(4,4);
t150 = -t203 * Ifges(4,2) + t294 + t296;
t196 = Ifges(4,4) * t203;
t295 = t239 * Ifges(4,5);
t297 = t204 * Ifges(4,1);
t151 = -t196 + t295 + t297;
t253 = Ifges(5,5) * t241 - Ifges(5,6) * t240;
t305 = Ifges(5,2) * t240;
t307 = Ifges(5,4) * t241;
t254 = -t305 + t307;
t308 = Ifges(5,4) * t240;
t255 = Ifges(5,1) * t241 - t308;
t317 = t241 / 0.2e1;
t268 = (t185 * Ifges(5,1) + Ifges(5,4) * t261 + t203 * Ifges(5,5)) * t317;
t318 = -t240 / 0.2e1;
t269 = (t185 * Ifges(5,4) + Ifges(5,2) * t261 + t203 * Ifges(5,6)) * t318;
t292 = t241 * t29;
t310 = mrSges(4,3) * t203;
t319 = -t204 / 0.2e1;
t320 = t203 / 0.2e1;
t321 = -t203 / 0.2e1;
t198 = qJD(5) + t203;
t322 = t198 / 0.2e1;
t323 = -t198 / 0.2e1;
t325 = t115 / 0.2e1;
t326 = -t115 / 0.2e1;
t327 = t350 / 0.2e1;
t328 = -t350 / 0.2e1;
t111 = Ifges(6,4) * t350;
t47 = Ifges(6,1) * t115 + Ifges(6,5) * t198 + t111;
t330 = t47 / 0.2e1;
t306 = Ifges(6,4) * t115;
t46 = Ifges(6,2) * t350 + Ifges(6,6) * t198 + t306;
t331 = t46 / 0.2e1;
t334 = Ifges(6,1) * t333 + Ifges(6,4) * t332 + Ifges(6,5) * t324;
t335 = Ifges(6,4) * t333 + Ifges(6,2) * t332 + Ifges(6,6) * t324;
t299 = t198 * Ifges(6,3);
t303 = t115 * Ifges(6,5);
t304 = t350 * Ifges(6,6);
t45 = t299 + t303 + t304;
t55 = t168 * Ifges(5,6) + t167 * t254;
t56 = t168 * Ifges(5,5) + t167 * t255;
t61 = pkin(4) * t287 + t102;
t300 = t185 * Ifges(5,5);
t301 = t261 * Ifges(5,6);
t97 = t203 * Ifges(5,3) + t300 + t301;
t338 = (Ifges(5,5) * t240 + Ifges(6,5) * t213 + Ifges(5,6) * t241 - Ifges(6,6) * t250) * t324 + t61 * (mrSges(6,1) * t250 + mrSges(6,2) * t213) + (Ifges(6,4) * t213 - Ifges(6,2) * t250) * t332 + (Ifges(6,1) * t213 - Ifges(6,4) * t250) * t333 - t250 * t335 + (-Ifges(4,2) * t204 + t151 - t196) * t320 - t101 * mrSges(4,2) + (-Ifges(6,5) * t199 - Ifges(6,6) * t200) * t322 + (-Ifges(6,1) * t199 - Ifges(6,4) * t200) * t325 + (-Ifges(6,4) * t199 - Ifges(6,2) * t200) * t327 - t77 * (-mrSges(5,2) * t204 + mrSges(5,3) * t284) - t76 * (mrSges(5,1) * t204 + mrSges(5,3) * t283) - t226 * (mrSges(4,1) * t204 - mrSges(4,2) * t203) - t261 * (Ifges(5,6) * t204 - t203 * t254) / 0.2e1 + t204 * t352 - t174 * t310 + mrSges(5,3) * t292 + (-mrSges(5,1) * t241 + mrSges(5,2) * t240 - mrSges(4,1)) * t102 - t204 * t353 + t203 * t268 + t203 * t269 - t185 * (Ifges(5,5) * t204 - t203 * t255) / 0.2e1 + (-Ifges(4,1) * t203 - t296 + t45 + t97) * t319 + (t13 * t354 - t14 * t351 - t2 * t250 - t3 * t213) * mrSges(6,3) + (t351 * mrSges(6,1) - mrSges(6,2) * t354) * t103 + t55 * t317 + (Ifges(5,3) * t204 - t203 * t253) * t321 + (Ifges(6,5) * t140 + Ifges(6,6) * t139 + Ifges(6,3) * t204) * t323 - t239 * (-Ifges(4,5) * t203 - Ifges(4,6) * t204) / 0.2e1 + t240 * t56 / 0.2e1 + t204 * t150 / 0.2e1 + Ifges(4,5) * t167 - Ifges(4,6) * t168 - t139 * t46 / 0.2e1 - t140 * t47 / 0.2e1 + t203 * t345 + (Ifges(6,1) * t140 + Ifges(6,4) * t139 + Ifges(6,5) * t204) * t326 + (Ifges(6,4) * t140 + Ifges(6,2) * t139 + Ifges(6,6) * t204) * t328 - t199 * t330 - t200 * t331 + t213 * t334 + (Ifges(5,1) * t240 + t307) * t286 / 0.2e1 - (Ifges(5,2) * t241 + t308) * t287 / 0.2e1;
t337 = -0.2e1 * pkin(1);
t314 = pkin(6) * (qJD(2) * mrSges(3,1) - mrSges(3,3) * t278);
t312 = t241 * pkin(4);
t293 = t240 * t28;
t219 = t244 * t267;
t220 = t246 * t267;
t117 = qJD(3) * t341 + t316 * t219 + t243 * t220;
t276 = qJD(2) * t244;
t83 = pkin(2) * t276 + pkin(3) * t179 - qJ(4) * t178 - qJD(4) * t215;
t39 = t241 * t117 + t240 * t83;
t291 = Ifges(3,5) * qJD(2);
t288 = t102 * t341;
t285 = t178 * t240;
t282 = t215 * t240;
t88 = mrSges(5,1) * t287 + mrSges(5,2) * t286;
t173 = -pkin(3) * t248 - t215 * qJ(4) + t236;
t188 = t243 * t227 - t228 * t316;
t108 = t240 * t173 + t241 * t188;
t272 = t316 * pkin(2);
t10 = -t44 * mrSges(6,1) + t43 * mrSges(6,2);
t263 = t246 * t274;
t38 = -t117 * t240 + t241 * t83;
t107 = t241 * t173 - t188 * t240;
t176 = t217 * t243 - t206;
t235 = -t272 - pkin(3);
t252 = t292 - t293;
t251 = -t240 * t76 + t241 * t77;
t69 = -pkin(4) * t248 - t215 * t238 + t107;
t82 = -pkin(8) * t282 + t108;
t21 = -t242 * t82 + t245 * t69;
t22 = t242 * t69 + t245 * t82;
t118 = qJD(3) * t188 + t243 * t219 - t316 * t220;
t237 = Ifges(3,4) * t277;
t233 = -pkin(3) - t312;
t221 = t235 - t312;
t202 = Ifges(3,1) * t278 + t237 + t291;
t191 = pkin(4) * t284;
t189 = -mrSges(4,2) * t239 - t310;
t171 = mrSges(4,1) * t203 + mrSges(4,2) * t204;
t164 = Ifges(6,3) * t168;
t156 = t250 * t215;
t155 = t213 * t215;
t143 = pkin(4) * t282 - t341;
t134 = mrSges(5,1) * t203 - mrSges(5,3) * t185;
t133 = -mrSges(5,2) * t203 + mrSges(5,3) * t261;
t127 = t176 - t191;
t123 = -t191 + t175;
t93 = mrSges(5,1) * t168 - mrSges(5,3) * t286;
t92 = -mrSges(5,2) * t168 - mrSges(5,3) * t287;
t86 = mrSges(6,1) * t198 - mrSges(6,3) * t115;
t85 = -mrSges(6,2) * t198 + mrSges(6,3) * t350;
t73 = pkin(4) * t285 + t118;
t59 = -t178 * t213 + t199 * t215;
t58 = -t178 * t250 - t200 * t215;
t51 = -mrSges(6,1) * t350 + mrSges(6,2) * t115;
t30 = -pkin(8) * t285 + t39;
t27 = -mrSges(6,2) * t168 + mrSges(6,3) * t44;
t26 = mrSges(6,1) * t168 - mrSges(6,3) * t43;
t20 = pkin(4) * t179 - t178 * t238 + t38;
t5 = -qJD(5) * t22 + t20 * t245 - t242 * t30;
t4 = qJD(5) * t21 + t20 * t242 + t245 * t30;
t1 = [-(t28 * mrSges(5,1) - t29 * mrSges(5,2) + t164 / 0.2e1 + t253 * t167 + (Ifges(4,2) + Ifges(5,3) + Ifges(6,3) / 0.2e1) * t168 + t339) * t248 + (-t290 / 0.2e1 + (mrSges(3,1) * t337 - 0.3e1 / 0.2e1 * t309 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t246) * qJD(1) + (qJD(1) * (-mrSges(4,1) * t248 + mrSges(4,2) * t215) + 0.2e1 * m(4) * t226 + t171) * pkin(2) - t340) * t276 + (t167 * t248 - t168 * t215 + t178 * t321 + t179 * t319) * Ifges(4,4) + m(5) * (t107 * t28 + t108 * t29 + t118 * t145 + t38 * t76 + t39 * t77 - t288) + m(4) * (t101 * t188 + t117 * t175 - t118 * t174 - t288) + (t101 * t248 + t102 * t215 - t167 * t341 - t168 * t188 - t174 * t178 - t175 * t179) * mrSges(4,3) - t341 * t88 - t342 * t118 + t107 * t93 + t108 * t92 + t103 * (-mrSges(6,1) * t59 + mrSges(6,2) * t58) + t4 * t85 + t5 * t86 + t73 * t51 + m(6) * (t103 * t73 + t13 * t5 + t14 * t4 + t143 * t61 + t2 * t22 + t21 * t3) + (-t13 * t58 + t14 * t59 - t155 * t2 + t156 * t3) * mrSges(6,3) + t61 * (mrSges(6,1) * t155 - mrSges(6,2) * t156) + (-Ifges(6,5) * t156 - Ifges(6,6) * t155) * t324 + (-Ifges(6,4) * t156 - Ifges(6,2) * t155) * t332 + (-Ifges(6,1) * t156 - Ifges(6,4) * t155) * t333 + (t297 / 0.2e1 + t253 * t320 + t261 * t254 / 0.2e1 + t185 * t255 / 0.2e1 + t345 + t295 / 0.2e1 + t226 * mrSges(4,2) + t151 / 0.2e1 + t269 + t268 + (-t240 * t77 - t241 * t76) * mrSges(5,3)) * t178 + (t301 / 0.2e1 + t300 / 0.2e1 + t76 * mrSges(5,1) - t294 / 0.2e1 - t77 * mrSges(5,2) + t226 * mrSges(4,1) + t299 / 0.2e1 + t97 / 0.2e1 + t353 - t352 + t45 / 0.2e1 + t304 / 0.2e1 + t303 / 0.2e1 - t150 / 0.2e1 + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t203) * t179 + (Ifges(6,5) * t58 + Ifges(6,6) * t59) * t322 + t236 * (mrSges(4,1) * t168 + mrSges(4,2) * t167) + t117 * t189 + t143 * t10 + t39 * t133 + t38 * t134 + t21 * t26 + t22 * t27 + (Ifges(6,1) * t58 + Ifges(6,4) * t59) * t325 + (Ifges(6,4) * t58 + Ifges(6,2) * t59) * t327 + t58 * t330 + t59 * t331 - t156 * t334 - t155 * t335 + (t202 / 0.2e1 - t314 + t291 / 0.2e1 + (mrSges(3,2) * t337 + 0.3e1 / 0.2e1 * Ifges(3,4) * t246) * qJD(1)) * t246 * qJD(2) + (t102 * t256 + t253 * t324 + t55 * t318 + t56 * t317 + (-t240 * t29 - t241 * t28) * mrSges(5,3) + (Ifges(4,1) + Ifges(5,1) * t241 ^ 2 / 0.2e1 + (-t307 + t305 / 0.2e1) * t240) * t167) * t215; t338 + (m(5) * t145 + m(6) * t103 - t342 + t51) * pkin(2) * t275 + t340 * t278 + (t102 * t235 - t145 * t176 + t229 * t251 + t232 * t252 - t76 * t89 - t77 * t90) * m(5) + (t260 - t177) * t189 + (-mrSges(3,1) * t263 + mrSges(3,2) * t264) * pkin(6) + (-t167 * t272 - t168 * t315) * mrSges(4,3) + ((-t316 * t102 + t101 * t243 + (-t174 * t243 + t175 * t316) * qJD(3)) * pkin(2) + t174 * t176 - t175 * t177 - t226 * t271) * m(4) + t342 * t176 - (-Ifges(3,2) * t278 + t202 + t237) * t277 / 0.2e1 + (-t244 * (Ifges(3,1) * t246 - t309) / 0.2e1 + pkin(1) * (mrSges(3,1) * t244 + mrSges(3,2) * t246)) * qJD(1) ^ 2 + t92 * t281 + t175 * t298 + t348 * t86 + (-t103 * t127 + t13 * t348 + t14 * t349 + t162 * t3 + t163 * t2 + t221 * t61) * m(6) + t349 * t85 + (-t240 * t229 - t89) * t134 + (t229 * t241 - t90) * t133 + t277 * t314 + t235 * t88 + t221 * t10 + t163 * t27 + t162 * t26 - t127 * t51 - Ifges(3,6) * t264 - t240 * t232 * t93 - t171 * t271 - (Ifges(3,5) * t246 - Ifges(3,6) * t244) * t274 / 0.2e1 + Ifges(3,5) * t263 - mrSges(5,3) * t293; -pkin(3) * t88 + t233 * t10 + (-t28 * mrSges(5,3) - qJ(4) * t93 - qJD(4) * t134) * t240 + (t342 + t298) * t175 + t181 * t26 + t182 * t27 - t174 * t189 - t123 * t51 - t94 * t134 + t347 * t86 + t346 * t85 + t92 * t289 + (qJD(4) * t241 - t95) * t133 + (-t103 * t123 + t13 * t347 + t14 * t346 + t181 * t3 + t182 * t2 + t233 * t61) * m(6) + (-pkin(3) * t102 + qJ(4) * t252 + qJD(4) * t251 - t145 * t175 - t76 * t94 - t77 * t95) * m(5) + t338; t115 * t86 - t350 * t85 - t261 * t133 + t185 * t134 + t10 + t88 + (t115 * t13 - t14 * t350 + t61) * m(6) + (t185 * t76 - t261 * t77 + t102) * m(5); t164 - t103 * (mrSges(6,1) * t115 + mrSges(6,2) * t350) + (Ifges(6,1) * t350 - t306) * t326 + t46 * t325 + (Ifges(6,5) * t350 - Ifges(6,6) * t115) * t323 - t13 * t85 + t14 * t86 + (t115 * t14 + t13 * t350) * mrSges(6,3) + (-Ifges(6,2) * t115 + t111 + t47) * t328 + t339;];
tauc = t1(:);
