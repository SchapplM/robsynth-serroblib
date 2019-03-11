% Calculate time derivative of joint inertia matrix for
% S6RRRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRP4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP4_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP4_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP4_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP4_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP4_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP4_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:12:18
% EndTime: 2019-03-10 01:12:32
% DurationCPUTime: 6.28s
% Computational Cost: add. (10204->512), mult. (22697->703), div. (0->0), fcn. (21105->8), ass. (0->221)
t235 = sin(qJ(3));
t236 = sin(qJ(2));
t239 = cos(qJ(3));
t240 = cos(qJ(2));
t201 = t235 * t240 + t239 * t236;
t348 = qJD(2) + qJD(3);
t165 = t348 * t201;
t360 = Ifges(7,4) + Ifges(6,5);
t368 = t360 * t165;
t367 = Ifges(6,1) + Ifges(7,1);
t366 = -Ifges(6,4) + Ifges(7,5);
t358 = -Ifges(6,6) + Ifges(7,6);
t365 = -mrSges(6,1) - mrSges(7,1);
t364 = -mrSges(6,2) + mrSges(7,3);
t335 = -pkin(8) - pkin(7);
t214 = t335 * t236;
t216 = t335 * t240;
t182 = t214 * t235 - t216 * t239;
t279 = qJD(2) * t335;
t209 = t236 * t279;
t264 = t240 * t279;
t121 = qJD(3) * t182 + t209 * t235 - t239 * t264;
t234 = sin(qJ(4));
t289 = qJD(4) * t234;
t277 = t201 * t289;
t199 = t235 * t236 - t239 * t240;
t164 = t348 * t199;
t238 = cos(qJ(4));
t303 = t164 * t238;
t253 = t277 + t303;
t288 = qJD(4) * t238;
t304 = t164 * t234;
t254 = t201 * t288 - t304;
t64 = mrSges(5,1) * t254 - mrSges(5,2) * t253;
t363 = -m(5) * t121 - t64;
t233 = sin(qJ(5));
t237 = cos(qJ(5));
t198 = t233 * t234 - t237 * t238;
t347 = qJD(4) + qJD(5);
t162 = t347 * t198;
t295 = t233 * t238;
t200 = t234 * t237 + t295;
t163 = t347 * t200;
t265 = -t162 * t360 + t163 * t358;
t362 = Ifges(5,5) * t288 + t265;
t349 = (t234 ^ 2 + t238 ^ 2) * t239;
t352 = (mrSges(6,3) + mrSges(7,2));
t361 = 2 * t352;
t359 = Ifges(7,2) + Ifges(6,3);
t137 = t200 * t201;
t48 = -t137 * t347 + t198 * t164;
t287 = qJD(5) * t233;
t300 = t201 * t238;
t301 = t201 * t234;
t49 = -t164 * t295 - t233 * t277 - t287 * t301 + (t300 * t347 - t304) * t237;
t357 = t366 * t49 + t367 * t48 + t368;
t138 = t198 * t201;
t356 = t366 * t137 - t138 * t367 + t360 * t199;
t355 = -t162 * t367 + t366 * t163;
t354 = t366 * t198 + t200 * t367;
t350 = t239 * t214 + t216 * t235;
t120 = qJD(3) * t350 + t239 * t209 + t235 * t264;
t225 = -pkin(2) * t240 - pkin(1);
t149 = t199 * pkin(3) - t201 * pkin(9) + t225;
t98 = pkin(2) * qJD(2) * t236 + pkin(3) * t165 + pkin(9) * t164;
t28 = t238 * t120 + t149 * t288 - t182 * t289 + t234 * t98;
t174 = t238 * t182;
t101 = t234 * t149 + t174;
t267 = -t120 * t234 + t238 * t98;
t29 = -qJD(4) * t101 + t267;
t353 = -t29 * t234 + t238 * t28;
t351 = -Ifges(5,5) * t303 + Ifges(5,3) * t165;
t221 = pkin(2) * t235 + pkin(9);
t229 = t238 * pkin(10);
t197 = t221 * t238 + t229;
t321 = -pkin(10) - t221;
t273 = t321 * t234;
t142 = t233 * t197 - t237 * t273;
t269 = qJD(4) * t321;
t312 = pkin(2) * qJD(3);
t284 = t239 * t312;
t183 = t234 * t269 + t238 * t284;
t251 = -t234 * t284 + t238 * t269;
t74 = -qJD(5) * t142 + t237 * t183 + t233 * t251;
t143 = t237 * t197 + t233 * t273;
t75 = qJD(5) * t143 + t233 * t183 - t237 * t251;
t346 = t364 * t74 + t365 * t75;
t22 = pkin(10) * t303 + pkin(4) * t165 + (-t174 + (pkin(10) * t201 - t149) * t234) * qJD(4) + t267;
t24 = -pkin(10) * t254 + t28;
t100 = t238 * t149 - t182 * t234;
t71 = pkin(4) * t199 - pkin(10) * t300 + t100;
t87 = -pkin(10) * t301 + t101;
t318 = t233 * t71 + t237 * t87;
t8 = -qJD(5) * t318 + t22 * t237 - t233 * t24;
t215 = pkin(9) * t238 + t229;
t334 = -pkin(10) - pkin(9);
t280 = t334 * t234;
t179 = t233 * t215 - t237 * t280;
t278 = qJD(4) * t334;
t207 = t234 * t278;
t263 = t238 * t278;
t118 = -qJD(5) * t179 + t237 * t207 + t233 * t263;
t181 = t237 * t215 + t233 * t280;
t119 = qJD(5) * t181 + t233 * t207 - t237 * t263;
t345 = t364 * t118 + t119 * t365;
t344 = 0.2e1 * m(5);
t343 = 2 * m(6);
t342 = 2 * m(7);
t102 = mrSges(7,1) * t163 + mrSges(7,3) * t162;
t341 = 0.2e1 * t102;
t103 = mrSges(6,1) * t163 - mrSges(6,2) * t162;
t340 = 0.2e1 * t103;
t339 = 0.2e1 * t121;
t168 = mrSges(7,1) * t198 - mrSges(7,3) * t200;
t338 = 0.2e1 * t168;
t337 = 0.2e1 * t225;
t326 = pkin(2) * t239;
t169 = mrSges(6,1) * t198 + mrSges(6,2) * t200;
t325 = pkin(4) * t169;
t324 = pkin(4) * t233;
t34 = mrSges(6,1) * t165 - mrSges(6,3) * t48;
t35 = -t165 * mrSges(7,1) + t48 * mrSges(7,2);
t320 = t35 - t34;
t36 = -mrSges(6,2) * t165 - mrSges(6,3) * t49;
t37 = -mrSges(7,2) * t49 + mrSges(7,3) * t165;
t319 = t36 + t37;
t317 = mrSges(6,3) * t163;
t316 = mrSges(6,3) * t198;
t315 = Ifges(5,4) * t234;
t314 = Ifges(5,4) * t238;
t313 = Ifges(5,6) * t234;
t152 = t162 * mrSges(7,2);
t309 = t235 * mrSges(4,1);
t307 = t239 * mrSges(4,2);
t305 = t121 * t350;
t227 = pkin(4) * t289;
t228 = t235 * t312;
t208 = t228 + t227;
t299 = t208 * t169;
t211 = -mrSges(5,1) * t238 + mrSges(5,2) * t234;
t294 = t235 * t211;
t125 = mrSges(6,1) * t199 + mrSges(6,3) * t138;
t126 = -mrSges(7,1) * t199 - mrSges(7,2) * t138;
t293 = -t125 + t126;
t124 = -mrSges(6,2) * t199 - mrSges(6,3) * t137;
t127 = -mrSges(7,2) * t137 + mrSges(7,3) * t199;
t292 = t127 + t124;
t286 = qJD(5) * t237;
t285 = 0.2e1 * t240;
t283 = pkin(4) * t287;
t282 = pkin(4) * t286;
t281 = t237 * t162 * mrSges(6,3);
t224 = -pkin(4) * t238 - pkin(3);
t275 = t142 * t287;
t274 = t179 * t287;
t272 = -t289 / 0.2e1;
t271 = -(2 * Ifges(4,4)) - t313;
t270 = t142 * t75 + t143 * t74;
t268 = t181 * t118 + t119 * t179;
t266 = t200 * t283;
t262 = mrSges(5,3) * t349;
t132 = pkin(4) * t301 - t350;
t261 = mrSges(5,1) * t234 + mrSges(5,2) * t238;
t260 = Ifges(5,1) * t238 - t315;
t259 = -Ifges(5,2) * t234 + t314;
t258 = Ifges(5,5) * t234 + Ifges(5,6) * t238;
t32 = -t233 * t87 + t237 * t71;
t255 = t165 * t359 + t358 * t49 + t360 * t48;
t7 = t233 * t22 + t237 * t24 + t71 * t286 - t287 * t87;
t252 = t118 * t143 + t119 * t142 + t179 * t75 + t181 * t74;
t146 = pkin(5) * t198 - qJ(6) * t200 + t224;
t86 = pkin(5) * t163 + qJ(6) * t162 - qJD(6) * t200 + t227;
t104 = -Ifges(7,5) * t162 + Ifges(7,3) * t163;
t105 = -Ifges(6,4) * t162 - Ifges(6,2) * t163;
t170 = Ifges(7,5) * t200 + Ifges(7,3) * t198;
t171 = Ifges(6,4) * t200 - Ifges(6,2) * t198;
t205 = t259 * qJD(4);
t206 = t260 * qJD(4);
t213 = Ifges(5,1) * t234 + t314;
t250 = t238 * t205 + t234 * t206 + t213 * t288 + t355 * t200 + (t104 - t105) * t198 + (t170 - t171) * t163 - t354 * t162;
t2 = qJ(6) * t165 + qJD(6) * t199 + t7;
t4 = -pkin(5) * t165 - t8;
t249 = t8 * mrSges(6,1) - t4 * mrSges(7,1) - t7 * mrSges(6,2) + t2 * mrSges(7,3) + t255;
t62 = pkin(4) * t254 + t121;
t246 = (pkin(5) * t162 - qJ(6) * t163 - qJD(6) * t198) * mrSges(7,2) + t265;
t217 = qJD(6) + t282;
t245 = -mrSges(6,2) * t282 + t217 * mrSges(7,3) + t283 * t365;
t218 = qJ(6) + t324;
t222 = -pkin(4) * t237 - pkin(5);
t243 = mrSges(6,3) * t266 - t222 * t152 - t282 * t316 - t317 * t324 + (-t163 * t218 - t198 * t217 + t266) * mrSges(7,2) + t362;
t150 = -mrSges(5,2) * t199 - mrSges(5,3) * t301;
t151 = mrSges(5,1) * t199 - mrSges(5,3) * t300;
t83 = mrSges(5,1) * t165 + mrSges(5,3) * t253;
t84 = -mrSges(5,2) * t165 - mrSges(5,3) * t254;
t242 = m(5) * (-t100 * t288 - t101 * t289 + t353) + t238 * t84 - t234 * t83 - t151 * t288 - t150 * t289;
t10 = pkin(5) * t49 - qJ(6) * t48 + qJD(6) * t138 + t62;
t128 = Ifges(5,6) * t199 + t201 * t259;
t129 = Ifges(5,5) * t199 + t201 * t260;
t15 = Ifges(7,5) * t48 + Ifges(7,6) * t165 + Ifges(7,3) * t49;
t16 = Ifges(6,4) * t48 - Ifges(6,2) * t49 + Ifges(6,6) * t165;
t204 = t261 * qJD(4);
t212 = Ifges(5,2) * t238 + t315;
t30 = qJ(6) * t199 + t318;
t31 = -pkin(5) * t199 - t32;
t53 = -Ifges(5,4) * t253 - Ifges(5,2) * t254 + Ifges(5,6) * t165;
t54 = -Ifges(5,1) * t253 - Ifges(5,4) * t254 + Ifges(5,5) * t165;
t60 = pkin(5) * t137 + qJ(6) * t138 + t132;
t77 = -Ifges(7,5) * t138 + Ifges(7,6) * t199 + Ifges(7,3) * t137;
t78 = -Ifges(6,4) * t138 - Ifges(6,2) * t137 + Ifges(6,6) * t199;
t241 = t128 * t272 - t318 * t317 + (t162 * t32 - t200 * t8) * mrSges(6,3) + (t170 / 0.2e1 - t171 / 0.2e1) * t49 - t350 * t204 + (-Ifges(5,6) * t289 + t362) * t199 / 0.2e1 + (t77 / 0.2e1 - t78 / 0.2e1) * t163 + (t104 / 0.2e1 - t105 / 0.2e1) * t137 + (t211 - mrSges(4,1)) * t121 + (-t163 * t30 - t198 * t2 + t200 * t4) * mrSges(7,2) + (t15 / 0.2e1 - t16 / 0.2e1) * t198 + t238 * t53 / 0.2e1 + t234 * t54 / 0.2e1 - Ifges(4,6) * t165 + t10 * t168 + t62 * t169 - Ifges(4,5) * t164 + t132 * t103 - t120 * mrSges(4,2) + t60 * t102 + ((-t100 * t238 - t101 * t234) * qJD(4) + t353) * mrSges(5,3) - t254 * t212 / 0.2e1 + t354 * t48 / 0.2e1 - t355 * t138 / 0.2e1 + (t201 * t272 - t303 / 0.2e1) * t213 - t356 * t162 / 0.2e1 + t357 * t200 / 0.2e1 + (t358 * t198 + t360 * t200 + t258) * t165 / 0.2e1 - t31 * t152 - t7 * t316 + t129 * t288 / 0.2e1 + t206 * t300 / 0.2e1 - t205 * t301 / 0.2e1;
t232 = qJD(6) * mrSges(7,3);
t223 = -pkin(3) - t326;
t210 = t224 - t326;
t145 = t261 * t201;
t135 = t146 - t326;
t89 = mrSges(6,1) * t137 - mrSges(6,2) * t138;
t88 = mrSges(7,1) * t137 + mrSges(7,3) * t138;
t81 = t228 + t86;
t20 = mrSges(6,1) * t49 + mrSges(6,2) * t48;
t19 = mrSges(7,1) * t49 - mrSges(7,3) * t48;
t1 = [(t77 - t78) * t49 + (t132 * t62 + t318 * t7 + t32 * t8) * t343 + 0.2e1 * t318 * t36 - (t357 + t368) * t138 - 0.2e1 * t350 * t64 + 0.2e1 * (t164 * t350 - t165 * t182) * mrSges(4,3) + t145 * t339 + (t10 * t60 + t2 * t30 + t31 * t4) * t342 + 0.2e1 * t28 * t150 + 0.2e1 * t29 * t151 + 0.2e1 * t8 * t125 + 0.2e1 * t4 * t126 + 0.2e1 * t2 * t127 + 0.2e1 * t132 * t20 + 0.2e1 * t7 * t124 + 0.2e1 * t100 * t83 + 0.2e1 * t101 * t84 + 0.2e1 * t10 * t88 + 0.2e1 * t62 * t89 + 0.2e1 * t60 * t19 + 0.2e1 * t30 * t37 + 0.2e1 * t32 * t34 + 0.2e1 * t31 * t35 + (mrSges(4,1) * t165 - mrSges(4,2) * t164) * t337 + (mrSges(4,3) * t339 - 0.2e1 * Ifges(4,1) * t164 - t234 * t53 + t238 * t54 + (Ifges(5,5) * t238 + t271) * t165 + (-t238 * t128 - t234 * t129 - t199 * t258) * qJD(4)) * t201 + (t100 * t29 + t101 * t28 - t305) * t344 + 0.2e1 * m(4) * (t120 * t182 - t305) + t356 * t48 + (t358 * t165 + t15 - t16) * t137 + (-0.2e1 * t120 * mrSges(4,3) - t271 * t164 + ((2 * Ifges(4,2)) + Ifges(5,3) + t359) * t165 + t255 + t351) * t199 - t129 * t303 + t128 * t304 + ((-mrSges(3,2) * pkin(1) + Ifges(3,4) * t240) * t285 + (m(4) * pkin(2) * t337 + 0.2e1 * pkin(2) * (mrSges(4,1) * t199 + mrSges(4,2) * t201) - 0.2e1 * pkin(1) * mrSges(3,1) - 0.2e1 * Ifges(3,4) * t236 + (Ifges(3,1) - Ifges(3,2)) * t285) * t236) * qJD(2); (m(4) * (t120 * t235 - t121 * t239) + (t164 * t239 - t165 * t235) * mrSges(4,3) + ((-t199 * mrSges(4,3) + t238 * t150 - t234 * t151 + m(5) * (-t100 * t234 + t101 * t238) + m(4) * t182) * t239 + (t201 * mrSges(4,3) + t145 - (m(5) + m(4)) * t350) * t235) * qJD(3)) * pkin(2) + t242 * t221 + t241 + (Ifges(3,5) * t240 - Ifges(3,6) * t236 + (-mrSges(3,1) * t240 + mrSges(3,2) * t236) * pkin(7)) * qJD(2) + t293 * t75 + t292 * t74 + t319 * t143 + t320 * t142 + t210 * t20 + t208 * t89 + t135 * t19 + t81 * t88 + m(6) * (t132 * t208 - t142 * t8 + t143 * t7 + t210 * t62 + t318 * t74 - t32 * t75) + m(7) * (t10 * t135 + t142 * t4 + t143 * t2 + t30 * t74 + t31 * t75 + t60 * t81) - t363 * t223; t250 + ((t221 * t349 + t223 * t235) * t344 - 0.2e1 * t307 + 0.2e1 * t294 - 0.2e1 * t309 + 0.2e1 * t262) * t312 - t212 * t289 + 0.2e1 * t223 * t204 + t210 * t340 + 0.2e1 * t299 + t81 * t338 + t135 * t341 + (t208 * t210 + t270) * t343 + (t135 * t81 + t270) * t342 + (-t142 * t162 - t143 * t163 - t74 * t198 + t75 * t200) * t361; t242 * pkin(9) + t319 * t181 + t320 * t179 + t293 * t119 + t292 * t118 + m(6) * (t118 * t318 - t119 * t32 + t132 * t227 - t179 * t8 + t181 * t7 + t224 * t62) + m(7) * (t10 * t146 + t118 * t30 + t119 * t31 + t179 * t4 + t181 * t2 + t60 * t86) + t241 + t89 * t227 + t224 * t20 + t146 * t19 + t86 * t88 + t363 * pkin(3); t250 + (m(5) * (-pkin(3) * t235 + pkin(9) * t349) - t307 + t294 - t309 + t262) * t312 + t299 + m(6) * (t208 * t224 + t210 * t227 + t252) + m(7) * (t135 * t86 + t146 * t81 + t252) + (-t212 + t325) * t289 + (t223 - pkin(3)) * t204 + (t146 + t135) * t102 + (t224 + t210) * t103 + (t86 + t81) * t168 + t352 * ((t119 + t75) * t200 + (-t118 - t74) * t198 + (-t143 - t181) * t163 - (t142 + t179) * t162); t250 + t224 * t340 - 0.2e1 * pkin(3) * t204 + t86 * t338 + t146 * t341 + (-t212 + 0.2e1 * t325) * t289 + (t224 * t227 + t268) * t343 + (t146 * t86 + t268) * t342 + (-t118 * t198 + t119 * t200 - t162 * t179 - t163 * t181) * t361; m(7) * (t2 * t218 + t217 * t30 + t222 * t4) + ((m(6) * t8 + t34 + (m(6) * t318 + t124) * qJD(5)) * t237 + (m(6) * t7 + t36 + (-m(6) * t32 + m(7) * t31 + t293) * qJD(5)) * t233) * pkin(4) - t254 * Ifges(5,6) - Ifges(5,5) * t277 + t249 + t217 * t127 + t218 * t37 + t222 * t35 - t28 * mrSges(5,2) + t29 * mrSges(5,1) + t351; m(7) * (t143 * t217 + t218 * t74 + t222 * t75) + (t281 + m(7) * t275 + m(6) * (t143 * t286 + t233 * t74 - t237 * t75 + t275)) * pkin(4) + t243 + (t211 * t221 - t313) * qJD(4) - t261 * t284 + t346; m(7) * (t118 * t218 + t119 * t222 + t181 * t217) + (t281 + m(7) * t274 + m(6) * (t118 * t233 - t119 * t237 + t181 * t286 + t274)) * pkin(4) + t243 + (pkin(9) * t211 - t313) * qJD(4) + t345; 0.2e1 * m(7) * (t217 * t218 + t222 * t283) + 0.2e1 * t245; m(7) * (-pkin(5) * t4 + qJ(6) * t2 + qJD(6) * t30) + t249 + qJD(6) * t127 + qJ(6) * t37 - pkin(5) * t35; m(7) * (-pkin(5) * t75 + qJ(6) * t74 + qJD(6) * t143) + t246 + t346; m(7) * (-pkin(5) * t119 + qJ(6) * t118 + qJD(6) * t181) + t246 + t345; m(7) * (-pkin(5) * t283 + qJ(6) * t217 + qJD(6) * t218) + t232 + t245; 0.2e1 * m(7) * qJ(6) * qJD(6) + 0.2e1 * t232; m(7) * t4 + t35; m(7) * t75 - t152; m(7) * t119 - t152; m(7) * t283; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
