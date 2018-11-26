% Calculate vector of centrifugal and coriolis load on the joints for
% S6RPRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 16:19
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RPRRPR7_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR7_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR7_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR7_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR7_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR7_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR7_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:19:43
% EndTime: 2018-11-23 16:19:48
% DurationCPUTime: 4.97s
% Computational Cost: add. (10016->487), mult. (22459->675), div. (0->0), fcn. (15498->8), ass. (0->241)
t220 = qJD(3) + qJD(4);
t225 = sin(qJ(4));
t226 = sin(qJ(3));
t279 = qJD(4) * t225;
t280 = qJD(3) * t226;
t228 = cos(qJ(4));
t229 = cos(qJ(3));
t283 = t228 * t229;
t147 = t220 * t283 - t225 * t280 - t226 * t279;
t191 = -t225 * t229 - t228 * t226;
t148 = t220 * t191;
t223 = sin(pkin(10));
t295 = cos(pkin(10));
t242 = t147 * t295 + t223 * t148;
t230 = -pkin(1) - pkin(7);
t202 = qJD(1) * t230 + qJD(2);
t282 = qJD(1) * t226;
t170 = -pkin(8) * t282 + t202 * t226;
t160 = t228 * t170;
t281 = qJD(1) * t229;
t171 = -pkin(8) * t281 + t229 * t202;
t315 = qJD(3) * pkin(3);
t162 = t171 + t315;
t121 = t162 * t225 + t160;
t183 = t191 * qJD(1);
t292 = qJ(5) * t183;
t105 = t121 + t292;
t286 = t223 * t105;
t159 = t225 * t170;
t120 = t228 * t162 - t159;
t184 = -t225 * t282 + t228 * t281;
t174 = t184 * qJ(5);
t104 = t120 - t174;
t96 = pkin(4) * t220 + t104;
t48 = t295 * t96 - t286;
t97 = t295 * t105;
t49 = t223 * t96 + t97;
t92 = t147 * t223 - t295 * t148;
t360 = t242 * t49 - t48 * t92;
t224 = sin(qJ(6));
t227 = cos(qJ(6));
t241 = t223 * t183 + t184 * t295;
t114 = t220 * t227 - t224 * t241;
t268 = t295 * t183 - t184 * t223;
t126 = qJD(6) - t268;
t263 = mrSges(7,1) * t224 + mrSges(7,2) * t227;
t46 = -t220 * pkin(5) - t48;
t244 = t46 * t263;
t258 = Ifges(7,5) * t227 - Ifges(7,6) * t224;
t317 = Ifges(7,4) * t227;
t260 = -Ifges(7,2) * t224 + t317;
t318 = Ifges(7,4) * t224;
t262 = Ifges(7,1) * t227 - t318;
t333 = t227 / 0.2e1;
t335 = -t224 / 0.2e1;
t115 = t220 * t224 + t227 * t241;
t342 = t115 / 0.2e1;
t319 = Ifges(7,4) * t115;
t44 = Ifges(7,2) * t114 + Ifges(7,6) * t126 + t319;
t113 = Ifges(7,4) * t114;
t45 = t115 * Ifges(7,1) + t126 * Ifges(7,5) + t113;
t359 = t126 * t258 / 0.2e1 + t262 * t342 + t114 * t260 / 0.2e1 + t45 * t333 + t44 * t335 + t244;
t356 = mrSges(6,3) * t268;
t307 = t241 * Ifges(6,4);
t310 = t268 * Ifges(6,4);
t124 = t228 * t171 - t159;
t107 = -t174 + t124;
t123 = -t171 * t225 - t160;
t246 = t123 - t292;
t270 = t295 * t225;
t316 = pkin(3) * qJD(4);
t354 = -t107 * t223 + t246 * t295 + (t223 * t228 + t270) * t316;
t285 = t223 * t225;
t169 = (t228 * t295 - t285) * t316;
t57 = t107 * t295 + t223 * t246;
t353 = t169 - t57;
t296 = -mrSges(6,1) * t220 - mrSges(7,1) * t114 + mrSges(7,2) * t115 + mrSges(6,3) * t241;
t352 = (qJ(2) * (m(3) + m(4)));
t323 = pkin(8) - t230;
t196 = t323 * t226;
t197 = t323 * t229;
t150 = -t228 * t196 - t225 * t197;
t137 = t148 * qJD(1);
t249 = t225 * t226 - t283;
t138 = t220 * t249 * qJD(1);
t90 = t137 * t295 + t223 * t138;
t52 = qJD(6) * t114 + t227 * t90;
t89 = t137 * t223 - t138 * t295;
t23 = mrSges(7,1) * t89 - mrSges(7,3) * t52;
t53 = -qJD(6) * t115 - t224 * t90;
t24 = -mrSges(7,2) * t89 + mrSges(7,3) * t53;
t351 = -t224 * t23 + t227 * t24;
t47 = pkin(9) * t220 + t49;
t200 = pkin(3) * t282 + qJD(1) * qJ(2);
t151 = -pkin(4) * t183 + qJD(5) + t200;
t64 = -pkin(5) * t268 - pkin(9) * t241 + t151;
t15 = -t224 * t47 + t227 * t64;
t16 = t224 * t64 + t227 * t47;
t350 = -t15 * t224 + t16 * t227;
t269 = pkin(8) * qJD(1) - t202;
t251 = t269 * t226;
t252 = t269 * t229;
t234 = (t225 * t252 + t228 * t251) * qJD(3);
t248 = -t137 * qJ(5) - t184 * qJD(5);
t278 = qJD(4) * t228;
t69 = t162 * t278 - t170 * t279 + (t225 * t251 - t228 * t252) * qJD(3);
t40 = qJ(5) * t138 + qJD(5) * t183 + t69;
t13 = t295 * t40 + (-t162 * t279 - t170 * t278 + t234 + t248) * t223;
t217 = pkin(3) * t281;
t195 = qJD(1) * qJD(2) + qJD(3) * t217;
t116 = -pkin(4) * t138 + t195;
t25 = pkin(5) * t89 - pkin(9) * t90 + t116;
t2 = qJD(6) * t15 + t13 * t227 + t224 * t25;
t291 = qJD(6) * t16;
t3 = -t13 * t224 + t227 * t25 - t291;
t349 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,5) * t52 + Ifges(7,6) * t53;
t348 = t52 / 0.2e1;
t347 = t53 / 0.2e1;
t346 = t89 / 0.2e1;
t344 = -t114 / 0.2e1;
t343 = -t115 / 0.2e1;
t341 = -t126 / 0.2e1;
t340 = t148 / 0.2e1;
t339 = -t183 / 0.2e1;
t338 = -t184 / 0.2e1;
t337 = t184 / 0.2e1;
t334 = t226 / 0.2e1;
t332 = -t229 / 0.2e1;
t331 = pkin(4) * t184;
t330 = pkin(4) * t223;
t70 = -qJD(4) * t121 + t234;
t12 = t223 * t40 - t295 * (t248 + t70);
t122 = qJ(5) * t191 + t150;
t149 = t196 * t225 - t228 * t197;
t245 = qJ(5) * t249 + t149;
t71 = t122 * t223 - t245 * t295;
t329 = t12 * t71;
t328 = t15 * mrSges(7,3);
t327 = t2 * t227;
t326 = t224 * t3;
t325 = t89 * mrSges(6,3);
t324 = t90 * mrSges(6,3);
t322 = mrSges(5,3) * t183;
t321 = Ifges(4,4) * t226;
t320 = Ifges(4,4) * t229;
t314 = t114 * Ifges(7,6);
t313 = t115 * Ifges(7,5);
t146 = t191 * t223 - t249 * t295;
t312 = t12 * t146;
t311 = t126 * Ifges(7,3);
t309 = t268 * Ifges(6,2);
t308 = t241 * Ifges(6,1);
t304 = t184 * mrSges(5,3);
t303 = t184 * Ifges(5,4);
t302 = t220 * Ifges(6,5);
t301 = t220 * Ifges(6,6);
t298 = t49 * t241;
t294 = Ifges(4,5) * qJD(3);
t293 = Ifges(4,6) * qJD(3);
t290 = t268 * t224;
t289 = t268 * t227;
t288 = t137 * t249;
t287 = t138 * t191;
t211 = t226 * pkin(3) + qJ(2);
t214 = pkin(3) * t228 + pkin(4);
t173 = pkin(3) * t270 + t223 * t214;
t277 = qJD(6) * t224;
t276 = qJD(6) * t227;
t203 = t229 * t315 + qJD(2);
t273 = t295 * pkin(4);
t272 = t89 * mrSges(6,1) + t90 * mrSges(6,2);
t140 = -mrSges(5,1) * t183 + mrSges(5,2) * t184;
t271 = -m(5) * t200 - t140;
t163 = -pkin(4) * t191 + t211;
t129 = pkin(4) * t147 + t203;
t266 = -t2 * t224 - t227 * t3;
t265 = mrSges(4,1) * t226 + mrSges(4,2) * t229;
t264 = mrSges(7,1) * t227 - mrSges(7,2) * t224;
t261 = Ifges(7,1) * t224 + t317;
t259 = Ifges(7,2) * t227 + t318;
t257 = Ifges(7,5) * t224 + Ifges(7,6) * t227;
t256 = t15 * t227 + t16 * t224;
t72 = t122 * t295 + t223 * t245;
t240 = -t191 * t295 - t223 * t249;
t76 = pkin(5) * t240 - pkin(9) * t146 + t163;
t27 = t224 * t76 + t227 * t72;
t26 = -t224 * t72 + t227 * t76;
t73 = -mrSges(7,2) * t126 + mrSges(7,3) * t114;
t74 = mrSges(7,1) * t126 - mrSges(7,3) * t115;
t254 = -t224 * t74 + t227 * t73;
t253 = -t224 * t73 - t227 * t74;
t198 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t282;
t199 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t281;
t250 = t229 * t198 - t226 * t199;
t118 = -mrSges(6,2) * t220 + t356;
t247 = -t118 - t254;
t243 = qJ(2) * (mrSges(4,1) * t229 - mrSges(4,2) * t226);
t189 = t323 * t280;
t190 = qJD(3) * t197;
t99 = t225 * t189 - t228 * t190 + t196 * t279 - t197 * t278;
t67 = pkin(5) * t241 - pkin(9) * t268 + t331;
t172 = -pkin(3) * t285 + t214 * t295;
t236 = t120 * t148 + t121 * t147 - t191 * t69 - t249 * t70;
t100 = -qJD(4) * t150 + t228 * t189 + t190 * t225;
t235 = -qJD(6) * t256 - t326 + t327;
t233 = -qJ(5) * t148 + qJD(5) * t249 + t100;
t127 = t183 * Ifges(5,2) + t220 * Ifges(5,6) + t303;
t175 = Ifges(5,4) * t183;
t128 = t184 * Ifges(5,1) + t220 * Ifges(5,5) + t175;
t43 = t311 + t313 + t314;
t78 = t301 + t307 + t309;
t79 = t302 + t308 + t310;
t8 = t52 * Ifges(7,4) + t53 * Ifges(7,2) + t89 * Ifges(7,6);
t9 = t52 * Ifges(7,1) + t53 * Ifges(7,4) + t89 * Ifges(7,5);
t232 = t359 * qJD(6) + (-Ifges(5,2) * t184 + t128 + t175) * t339 + t48 * t356 - (Ifges(5,5) * t183 + Ifges(6,5) * t268 - Ifges(5,6) * t184 - Ifges(6,6) * t241) * t220 / 0.2e1 - (-Ifges(6,2) * t241 + t310 + t79) * t268 / 0.2e1 - (Ifges(6,1) * t268 - t307 + t43) * t241 / 0.2e1 + t241 * t78 / 0.2e1 + (Ifges(7,3) * t241 + t258 * t268) * t341 + (Ifges(7,5) * t241 + t262 * t268) * t343 + (Ifges(7,6) * t241 + t260 * t268) * t344 - t151 * (mrSges(6,1) * t241 + mrSges(6,2) * t268) - t200 * (mrSges(5,1) * t184 + mrSges(5,2) * t183) - t16 * (-mrSges(7,2) * t241 - mrSges(7,3) * t290) - t15 * (mrSges(7,1) * t241 - mrSges(7,3) * t289) + t224 * t9 / 0.2e1 + t257 * t346 + t259 * t347 + t261 * t348 + t8 * t333 + t127 * t337 + (Ifges(5,1) * t183 - t303) * t338 + mrSges(7,3) * t327 + t120 * t322 + (-t264 - mrSges(6,1)) * t12 + t121 * t304 - t45 * t289 / 0.2e1 + t44 * t290 / 0.2e1 - t268 * t244 - t13 * mrSges(6,2) - t69 * mrSges(5,2) + t70 * mrSges(5,1) - Ifges(6,6) * t89 + Ifges(6,5) * t90 + Ifges(5,5) * t137 + Ifges(5,6) * t138;
t210 = -t273 - pkin(5);
t193 = t265 * qJD(1);
t182 = t294 + (Ifges(4,1) * t229 - t321) * qJD(1);
t181 = t293 + (-Ifges(4,2) * t226 + t320) * qJD(1);
t167 = pkin(9) + t173;
t166 = -pkin(5) - t172;
t158 = mrSges(5,1) * t220 - t304;
t157 = -mrSges(5,2) * t220 + t322;
t156 = t217 + t331;
t86 = Ifges(7,3) * t89;
t84 = -mrSges(6,1) * t268 + mrSges(6,2) * t241;
t65 = t217 + t67;
t63 = -qJ(5) * t147 + qJD(5) * t191 + t99;
t55 = t104 * t295 - t286;
t54 = t104 * t223 + t97;
t34 = pkin(5) * t242 + pkin(9) * t92 + t129;
t22 = t223 * t233 + t295 * t63;
t21 = t223 * t63 - t233 * t295;
t20 = t224 * t65 + t227 * t57;
t19 = -t224 * t57 + t227 * t65;
t18 = t224 * t67 + t227 * t55;
t17 = -t224 * t55 + t227 * t67;
t14 = -mrSges(7,1) * t53 + mrSges(7,2) * t52;
t5 = -qJD(6) * t27 - t22 * t224 + t227 * t34;
t4 = qJD(6) * t26 + t22 * t227 + t224 * t34;
t1 = [(t71 * t90 - t72 * t89 - t360) * mrSges(6,3) - (t302 / 0.2e1 + t79 / 0.2e1 + t310 / 0.2e1 + t308 / 0.2e1 + t151 * mrSges(6,2) - t256 * mrSges(7,3) + t359) * t92 + (t262 * t348 + t260 * t347 + t258 * t346 + t8 * t335 + t9 * t333 + Ifges(6,1) * t90 - Ifges(6,4) * t89 + t116 * mrSges(6,2) + (mrSges(6,3) + t263) * t12 + t266 * mrSges(7,3) + (t46 * t264 + t259 * t344 + t261 * t343 + t257 * t341 + t45 * t335 - t227 * t44 / 0.2e1 - t350 * mrSges(7,3)) * qJD(6)) * t146 + (t148 * t337 - t288) * Ifges(5,1) + m(5) * (t100 * t120 + t121 * t99 + t149 * t70 + t150 * t69 + t195 * t211 + t200 * t203) + (-t137 * t149 + t138 * t150 - t236) * mrSges(5,3) + t195 * (-mrSges(5,1) * t191 - mrSges(5,2) * t249) + (t137 * t191 - t138 * t249 + t147 * t338 + t183 * t340) * Ifges(5,4) + (t15 * mrSges(7,1) - t16 * mrSges(7,2) + t314 / 0.2e1 + t313 / 0.2e1 + t311 / 0.2e1 - t301 / 0.2e1 + t43 / 0.2e1 - t78 / 0.2e1 - t309 / 0.2e1 - t307 / 0.2e1 + t151 * mrSges(6,1)) * t242 + t220 * (Ifges(5,5) * t148 - Ifges(5,6) * t147) / 0.2e1 + t200 * (mrSges(5,1) * t147 + mrSges(5,2) * t148) + t203 * t140 + t211 * (-mrSges(5,1) * t138 + mrSges(5,2) * t137) + qJD(2) * t193 + t128 * t340 + t163 * t272 + ((-t181 / 0.2e1 + t230 * t198 - t293 / 0.2e1) * t229 + (-t182 / 0.2e1 - t230 * t199 - t294 / 0.2e1) * t226) * qJD(3) + (-t13 * mrSges(6,3) + t86 / 0.2e1 - Ifges(6,4) * t90 + t116 * mrSges(6,1) + (Ifges(7,3) / 0.2e1 + Ifges(6,2)) * t89 + t349) * t240 + t296 * t21 + m(7) * (t15 * t5 + t16 * t4 + t2 * t27 + t21 * t46 + t26 * t3 + t329) + m(6) * (t116 * t163 + t129 * t151 + t13 * t72 - t21 * t48 + t22 * t49 + t329) + (t147 * t339 + t287) * Ifges(5,2) + t26 * t23 + t27 * t24 + (((2 * mrSges(3,3)) + t265 + (2 * t352)) * qJD(2) + (0.2e1 * t243 + (0.3e1 / 0.2e1 * Ifges(4,2) - 0.3e1 / 0.2e1 * Ifges(4,1)) * t226 * t229 + (0.3e1 / 0.2e1 * t226 ^ 2 - 0.3e1 / 0.2e1 * t229 ^ 2) * Ifges(4,4)) * qJD(3)) * qJD(1) + t71 * t14 + t4 * t73 + t5 * t74 + t22 * t118 + t129 * t84 - t147 * t127 / 0.2e1 + t99 * t157 + t100 * t158; t147 * t157 + t148 * t158 + t296 * t92 - (t14 + t324) * t146 + t250 * qJD(3) + (-t287 + t288) * mrSges(5,3) - t247 * t242 + (qJD(6) * t253 - t325 + t351) * t240 + m(6) * (t13 * t240 - t312 + t360) + m(5) * t236 + m(7) * (t235 * t240 + t242 * t350 + t46 * t92 - t312) + (-m(6) * t151 - m(7) * t256 - t193 + t253 + t271 - t84 + (-mrSges(3,3) - t352) * qJD(1)) * qJD(1); t232 + ((t157 * t228 - t158 * t225) * qJD(4) + (-t137 * t228 + t138 * t225) * mrSges(5,3)) * pkin(3) + (-qJD(3) * t265 - t250) * t202 + (t167 * t24 + t169 * t73 + (-t167 * t74 - t328) * qJD(6)) * t227 + (t182 * t334 + t229 * t181 / 0.2e1 + (-t243 + (-Ifges(4,1) * t226 - t320) * t332 + (-Ifges(4,2) * t229 - t321) * t334) * qJD(1) + t271 * t229 * pkin(3) + (-Ifges(4,5) * t226 / 0.2e1 + Ifges(4,6) * t332) * qJD(3)) * qJD(1) + (-t172 * t90 - t173 * t89 + t298) * mrSges(6,3) + t353 * t118 + (-t169 * t74 + (-qJD(6) * t73 - t23) * t167 + (-t3 - t291) * mrSges(7,3)) * t224 + t166 * t14 - t20 * t73 - t19 * t74 - t156 * t84 - t124 * t157 - t123 * t158 + t296 * t354 + (-t12 * t172 + t13 * t173 - t151 * t156 + t353 * t49 - t354 * t48) * m(6) + ((t225 * t69 + t228 * t70 + (-t120 * t225 + t121 * t228) * qJD(4)) * pkin(3) - t120 * t123 - t121 * t124) * m(5) + (t12 * t166 - t15 * t19 - t16 * t20 + t235 * t167 + t169 * t350 + t354 * t46) * m(7); mrSges(6,3) * t298 - t55 * t118 - t120 * t157 + t121 * t158 + t210 * t14 - t17 * t74 - t18 * t73 - t273 * t324 - t276 * t328 - t325 * t330 - t84 * t331 + t232 - t296 * t54 + (-t16 * t277 - t326) * mrSges(7,3) + (t12 * t210 - t15 * t17 - t16 * t18 - t46 * t54) * m(7) + ((-t12 * t295 + t13 * t223) * pkin(4) - t151 * t331 + t48 * t54 - t49 * t55) * m(6) + (m(7) * t235 - t276 * t74 - t277 * t73 + t351) * (pkin(9) + t330); t254 * qJD(6) - t296 * t241 + t247 * t268 + t224 * t24 + t227 * t23 + t272 + (t126 * t350 - t241 * t46 - t266) * m(7) + (t241 * t48 - t268 * t49 + t116) * m(6); t86 - t46 * (mrSges(7,1) * t115 + mrSges(7,2) * t114) + (Ifges(7,1) * t114 - t319) * t343 + t44 * t342 + (Ifges(7,5) * t114 - Ifges(7,6) * t115) * t341 - t15 * t73 + t16 * t74 + (t114 * t15 + t115 * t16) * mrSges(7,3) + (-Ifges(7,2) * t115 + t113 + t45) * t344 + t349;];
tauc  = t1(:);
