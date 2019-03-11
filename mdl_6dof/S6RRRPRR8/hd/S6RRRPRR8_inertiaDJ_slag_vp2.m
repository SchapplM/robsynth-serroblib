% Calculate time derivative of joint inertia matrix for
% S6RRRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 18:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR8_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR8_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR8_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR8_inertiaDJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR8_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR8_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR8_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:44:12
% EndTime: 2019-03-09 18:44:32
% DurationCPUTime: 8.54s
% Computational Cost: add. (17328->738), mult. (44262->1096), div. (0->0), fcn. (44715->12), ass. (0->307)
t384 = Ifges(4,3) + Ifges(5,3);
t266 = sin(pkin(12));
t275 = cos(qJ(3));
t271 = sin(qJ(3));
t328 = cos(pkin(12));
t291 = t328 * t271;
t234 = t266 * t275 + t291;
t270 = sin(qJ(5));
t313 = qJD(5) * t270;
t322 = t266 * t271;
t278 = t275 * t328 - t322;
t225 = t278 * qJD(3);
t274 = cos(qJ(5));
t318 = t274 * t225;
t279 = t234 * t313 - t318;
t269 = sin(qJ(6));
t273 = cos(qJ(6));
t286 = t269 * t270 - t273 * t274;
t357 = -t286 / 0.2e1;
t237 = t269 * t274 + t270 * t273;
t356 = t237 / 0.2e1;
t349 = t270 / 0.2e1;
t348 = t274 / 0.2e1;
t295 = -t313 / 0.2e1;
t165 = t286 * t234;
t263 = -pkin(3) * t275 - pkin(2);
t175 = -pkin(4) * t278 - pkin(10) * t234 + t263;
t344 = -qJ(4) - pkin(9);
t247 = t344 * t275;
t197 = -t247 * t328 + t322 * t344;
t187 = t274 * t197;
t120 = t270 * t175 + t187;
t268 = cos(pkin(6));
t272 = sin(qJ(2));
t267 = sin(pkin(6));
t276 = cos(qJ(2));
t320 = t267 * t276;
t229 = t268 * t272 * pkin(1) + pkin(8) * t320;
t210 = pkin(9) * t268 + t229;
t211 = (-pkin(2) * t276 - pkin(9) * t272 - pkin(1)) * t267;
t149 = t275 * t210 + t271 * t211;
t321 = t267 * t272;
t256 = pkin(8) * t321;
t347 = pkin(1) * t276;
t228 = t268 * t347 - t256;
t312 = qJD(5) * t274;
t280 = t270 * t225 + t234 * t312;
t293 = qJD(3) * t344;
t223 = qJD(4) * t275 + t271 * t293;
t277 = -qJD(4) * t271 + t275 * t293;
t157 = t223 * t328 + t266 * t277;
t224 = t234 * qJD(3);
t315 = qJD(3) * t271;
t307 = pkin(3) * t315;
t158 = pkin(4) * t224 - pkin(10) * t225 + t307;
t55 = t274 * t157 + t270 * t158 + t175 * t312 - t197 * t313;
t290 = -t157 * t270 + t274 * t158;
t56 = -qJD(5) * t120 + t290;
t383 = -t270 * t56 + t274 * t55;
t316 = qJD(2) * t272;
t301 = t267 * t316;
t317 = qJD(2) * t267;
t213 = (pkin(2) * t272 - pkin(9) * t276) * t317;
t214 = t228 * qJD(2);
t102 = -qJD(3) * t149 + t275 * t213 - t214 * t271;
t226 = t268 * t275 - t271 * t321;
t300 = t276 * t317;
t195 = qJD(3) * t226 + t275 * t300;
t227 = t268 * t271 + t275 * t321;
t66 = pkin(3) * t301 - qJ(4) * t195 - qJD(4) * t227 + t102;
t314 = qJD(3) * t275;
t101 = -t210 * t315 + t211 * t314 + t271 * t213 + t275 * t214;
t194 = -qJD(3) * t227 - t271 * t300;
t75 = qJ(4) * t194 + qJD(4) * t226 + t101;
t33 = t266 * t66 + t328 * t75;
t31 = pkin(10) * t301 + t33;
t128 = -t194 * t328 + t195 * t266;
t129 = t266 * t194 + t195 * t328;
t215 = t229 * qJD(2);
t147 = -t194 * pkin(3) + t215;
t54 = t128 * pkin(4) - t129 * pkin(10) + t147;
t148 = -t271 * t210 + t275 * t211;
t116 = -pkin(3) * t320 - t227 * qJ(4) + t148;
t130 = qJ(4) * t226 + t149;
t74 = t266 * t116 + t328 * t130;
t65 = -pkin(10) * t320 + t74;
t162 = -t226 * t328 + t227 * t266;
t163 = t266 * t226 + t227 * t328;
t209 = t256 + (-pkin(2) - t347) * t268;
t168 = -t226 * pkin(3) + t209;
t88 = t162 * pkin(4) - t163 * pkin(10) + t168;
t10 = t270 * t54 + t274 * t31 + t88 * t312 - t313 * t65;
t35 = t270 * t88 + t274 * t65;
t11 = -qJD(5) * t35 - t270 * t31 + t274 * t54;
t382 = t10 * t274 - t11 * t270;
t381 = Ifges(6,5) * t349 + Ifges(7,5) * t356 + Ifges(6,6) * t348 + Ifges(7,6) * t357;
t379 = qJD(5) + qJD(6);
t188 = t379 * t286;
t189 = t379 * t237;
t132 = -Ifges(7,5) * t188 - Ifges(7,6) * t189;
t264 = Ifges(6,5) * t312;
t380 = Ifges(6,6) * t295 + t264 / 0.2e1 + t132 / 0.2e1;
t289 = mrSges(6,1) * t270 + mrSges(6,2) * t274;
t238 = t289 * qJD(5);
t378 = 2 * m(5);
t377 = 2 * m(6);
t376 = 2 * m(7);
t375 = -2 * mrSges(3,3);
t374 = -2 * mrSges(5,3);
t156 = t223 * t266 - t328 * t277;
t373 = 0.2e1 * t156;
t196 = -t247 * t266 - t344 * t291;
t372 = 0.2e1 * t196;
t371 = 0.2e1 * t215;
t370 = m(5) * pkin(3);
t281 = -t274 * t163 + t270 * t320;
t80 = qJD(5) * t281 - t270 * t129 + t274 * t301;
t369 = t80 / 0.2e1;
t140 = -t270 * t163 - t274 * t320;
t89 = t140 * t273 + t269 * t281;
t368 = t89 / 0.2e1;
t90 = t140 * t269 - t273 * t281;
t367 = t90 / 0.2e1;
t365 = t140 / 0.2e1;
t164 = t237 * t234;
t364 = -t164 / 0.2e1;
t363 = -t165 / 0.2e1;
t362 = -t188 / 0.2e1;
t361 = -t189 / 0.2e1;
t192 = Ifges(7,4) * t237 - Ifges(7,2) * t286;
t359 = t192 / 0.2e1;
t193 = Ifges(7,1) * t237 - Ifges(7,4) * t286;
t358 = t193 / 0.2e1;
t340 = Ifges(6,4) * t270;
t288 = Ifges(6,1) * t274 - t340;
t243 = t288 * qJD(5);
t354 = t243 / 0.2e1;
t339 = Ifges(6,4) * t274;
t251 = Ifges(6,1) * t270 + t339;
t351 = t251 / 0.2e1;
t350 = -t270 / 0.2e1;
t346 = pkin(3) * t266;
t261 = pkin(10) + t346;
t345 = pkin(11) + t261;
t343 = mrSges(7,3) * t237;
t342 = Ifges(4,4) * t271;
t341 = Ifges(4,4) * t275;
t338 = Ifges(6,6) * t270;
t335 = t188 * mrSges(7,3);
t334 = t189 * mrSges(7,3);
t333 = t214 * mrSges(3,2);
t332 = t215 * mrSges(3,1);
t331 = t286 * mrSges(7,3);
t327 = t156 * t196;
t326 = t234 * t270;
t325 = t234 * t274;
t324 = t261 * t270;
t323 = t261 * t274;
t311 = qJD(6) * t269;
t310 = qJD(6) * t273;
t309 = 0.2e1 * t267;
t79 = qJD(5) * t140 + t274 * t129 + t270 * t301;
t22 = qJD(6) * t89 + t269 * t80 + t273 * t79;
t23 = -qJD(6) * t90 - t269 * t79 + t273 * t80;
t6 = Ifges(7,5) * t22 + Ifges(7,6) * t23 + Ifges(7,3) * t128;
t25 = Ifges(6,5) * t79 + Ifges(6,6) * t80 + Ifges(6,3) * t128;
t83 = -t189 * t234 - t225 * t286;
t84 = t165 * t379 - t237 * t225;
t36 = Ifges(7,5) * t83 + Ifges(7,6) * t84 + Ifges(7,3) * t224;
t306 = pkin(5) * t313;
t305 = Ifges(4,6) * t320;
t304 = mrSges(6,3) * t313;
t303 = mrSges(6,3) * t312;
t302 = t328 * pkin(3);
t299 = t261 * t313;
t298 = t261 * t312;
t294 = t312 / 0.2e1;
t34 = -t270 * t65 + t274 * t88;
t292 = qJD(5) * t345;
t78 = t128 * mrSges(5,1) + t129 * mrSges(5,2);
t169 = t224 * mrSges(5,1) + t225 * mrSges(5,2);
t131 = t189 * mrSges(7,1) - t188 * mrSges(7,2);
t119 = t274 * t175 - t197 * t270;
t262 = -t302 - pkin(4);
t32 = -t266 * t75 + t328 * t66;
t246 = -mrSges(6,1) * t274 + mrSges(6,2) * t270;
t287 = -Ifges(6,2) * t270 + t339;
t24 = pkin(5) * t162 + pkin(11) * t281 + t34;
t28 = pkin(11) * t140 + t35;
t12 = t24 * t273 - t269 * t28;
t13 = t24 * t269 + t273 * t28;
t73 = t116 * t328 - t266 * t130;
t105 = -pkin(11) * t326 + t120;
t92 = -pkin(5) * t278 - pkin(11) * t325 + t119;
t51 = t105 * t273 + t269 * t92;
t50 = -t105 * t269 + t273 * t92;
t230 = t345 * t270;
t231 = t345 * t274;
t173 = -t230 * t273 - t231 * t269;
t174 = -t230 * t269 + t231 * t273;
t221 = t270 * t292;
t222 = t274 * t292;
t111 = qJD(6) * t173 - t221 * t273 - t222 * t269;
t112 = -qJD(6) * t174 + t221 * t269 - t222 * t273;
t285 = t112 * mrSges(7,1) - t111 * mrSges(7,2) + t132;
t4 = pkin(5) * t128 - pkin(11) * t79 + t11;
t5 = pkin(11) * t80 + t10;
t2 = qJD(6) * t12 + t269 * t4 + t273 * t5;
t3 = -qJD(6) * t13 - t269 * t5 + t273 * t4;
t284 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + t6;
t45 = -pkin(11) * t318 + pkin(5) * t224 + (-t187 + (pkin(11) * t234 - t175) * t270) * qJD(5) + t290;
t49 = -pkin(11) * t280 + t55;
t15 = qJD(6) * t50 + t269 * t45 + t273 * t49;
t16 = -qJD(6) * t51 - t269 * t49 + t273 * t45;
t283 = t16 * mrSges(7,1) - t15 * mrSges(7,2) + t36;
t64 = pkin(4) * t320 - t73;
t282 = Ifges(4,5) * t195 + Ifges(5,5) * t129 + Ifges(4,6) * t194 - Ifges(5,6) * t128 + t301 * t384;
t30 = -pkin(4) * t301 - t32;
t93 = -Ifges(6,5) * t279 - Ifges(6,6) * t280 + Ifges(6,3) * t224;
t265 = Ifges(4,5) * t314;
t255 = Ifges(3,5) * t300;
t252 = Ifges(4,1) * t271 + t341;
t250 = Ifges(4,2) * t275 + t342;
t249 = Ifges(6,2) * t274 + t340;
t245 = -t274 * pkin(5) + t262;
t244 = (Ifges(4,1) * t275 - t342) * qJD(3);
t242 = (-Ifges(4,2) * t271 + t341) * qJD(3);
t241 = t287 * qJD(5);
t239 = (mrSges(4,1) * t271 + mrSges(4,2) * t275) * qJD(3);
t232 = (-mrSges(7,1) * t269 - mrSges(7,2) * t273) * qJD(6) * pkin(5);
t220 = Ifges(5,5) * t225;
t219 = Ifges(5,6) * t224;
t199 = -mrSges(4,1) * t320 - t227 * mrSges(4,3);
t198 = mrSges(4,2) * t320 + t226 * mrSges(4,3);
t190 = mrSges(7,1) * t286 + mrSges(7,2) * t237;
t184 = Ifges(5,1) * t234 + Ifges(5,4) * t278;
t183 = Ifges(5,4) * t234 + Ifges(5,2) * t278;
t182 = -mrSges(5,1) * t278 + mrSges(5,2) * t234;
t177 = -mrSges(6,1) * t278 - mrSges(6,3) * t325;
t176 = mrSges(6,2) * t278 - mrSges(6,3) * t326;
t172 = t289 * t234;
t171 = Ifges(5,1) * t225 - Ifges(5,4) * t224;
t170 = Ifges(5,4) * t225 - Ifges(5,2) * t224;
t160 = mrSges(4,1) * t301 - mrSges(4,3) * t195;
t159 = -mrSges(4,2) * t301 + mrSges(4,3) * t194;
t155 = Ifges(4,1) * t227 + Ifges(4,4) * t226 - Ifges(4,5) * t320;
t154 = Ifges(4,4) * t227 + Ifges(4,2) * t226 - t305;
t153 = pkin(5) * t326 + t196;
t146 = -mrSges(5,1) * t320 - t163 * mrSges(5,3);
t145 = mrSges(5,2) * t320 - t162 * mrSges(5,3);
t144 = -Ifges(6,5) * t278 + t234 * t288;
t143 = -Ifges(6,6) * t278 + t234 * t287;
t142 = -Ifges(6,3) * t278 + (Ifges(6,5) * t274 - t338) * t234;
t139 = -mrSges(7,1) * t278 + mrSges(7,3) * t165;
t138 = mrSges(7,2) * t278 - mrSges(7,3) * t164;
t137 = -mrSges(6,2) * t224 - mrSges(6,3) * t280;
t136 = mrSges(6,1) * t224 + mrSges(6,3) * t279;
t135 = -mrSges(4,1) * t194 + mrSges(4,2) * t195;
t134 = -Ifges(7,1) * t188 - Ifges(7,4) * t189;
t133 = -Ifges(7,4) * t188 - Ifges(7,2) * t189;
t118 = Ifges(4,1) * t195 + Ifges(4,4) * t194 + Ifges(4,5) * t301;
t117 = Ifges(4,4) * t195 + Ifges(4,2) * t194 + Ifges(4,6) * t301;
t114 = mrSges(5,1) * t301 - mrSges(5,3) * t129;
t113 = -mrSges(5,2) * t301 - mrSges(5,3) * t128;
t110 = mrSges(6,1) * t280 - mrSges(6,2) * t279;
t109 = pkin(5) * t280 + t156;
t107 = mrSges(7,1) * t164 - mrSges(7,2) * t165;
t106 = mrSges(5,1) * t162 + mrSges(5,2) * t163;
t104 = Ifges(5,1) * t163 - Ifges(5,4) * t162 - Ifges(5,5) * t320;
t103 = Ifges(5,4) * t163 - Ifges(5,2) * t162 - Ifges(5,6) * t320;
t100 = mrSges(6,1) * t162 + mrSges(6,3) * t281;
t99 = -mrSges(6,2) * t162 + mrSges(6,3) * t140;
t98 = -Ifges(7,1) * t165 - Ifges(7,4) * t164 - Ifges(7,5) * t278;
t97 = -Ifges(7,4) * t165 - Ifges(7,2) * t164 - Ifges(7,6) * t278;
t96 = -Ifges(7,5) * t165 - Ifges(7,6) * t164 - Ifges(7,3) * t278;
t95 = -Ifges(6,1) * t279 - Ifges(6,4) * t280 + Ifges(6,5) * t224;
t94 = -Ifges(6,4) * t279 - Ifges(6,2) * t280 + Ifges(6,6) * t224;
t91 = -mrSges(6,1) * t140 - mrSges(6,2) * t281;
t70 = Ifges(5,1) * t129 - Ifges(5,4) * t128 + Ifges(5,5) * t301;
t69 = Ifges(5,4) * t129 - Ifges(5,2) * t128 + Ifges(5,6) * t301;
t68 = -mrSges(7,2) * t224 + mrSges(7,3) * t84;
t67 = mrSges(7,1) * t224 - mrSges(7,3) * t83;
t63 = -Ifges(6,1) * t281 + Ifges(6,4) * t140 + Ifges(6,5) * t162;
t62 = -Ifges(6,4) * t281 + Ifges(6,2) * t140 + Ifges(6,6) * t162;
t61 = -Ifges(6,5) * t281 + Ifges(6,6) * t140 + Ifges(6,3) * t162;
t58 = mrSges(7,1) * t162 - mrSges(7,3) * t90;
t57 = -mrSges(7,2) * t162 + mrSges(7,3) * t89;
t48 = -t140 * pkin(5) + t64;
t47 = -mrSges(6,2) * t128 + mrSges(6,3) * t80;
t46 = mrSges(6,1) * t128 - mrSges(6,3) * t79;
t44 = -mrSges(7,1) * t89 + mrSges(7,2) * t90;
t43 = -mrSges(7,1) * t84 + mrSges(7,2) * t83;
t42 = Ifges(7,1) * t90 + Ifges(7,4) * t89 + Ifges(7,5) * t162;
t41 = Ifges(7,4) * t90 + Ifges(7,2) * t89 + Ifges(7,6) * t162;
t40 = Ifges(7,5) * t90 + Ifges(7,6) * t89 + Ifges(7,3) * t162;
t39 = -mrSges(6,1) * t80 + mrSges(6,2) * t79;
t38 = Ifges(7,1) * t83 + Ifges(7,4) * t84 + Ifges(7,5) * t224;
t37 = Ifges(7,4) * t83 + Ifges(7,2) * t84 + Ifges(7,6) * t224;
t27 = Ifges(6,1) * t79 + Ifges(6,4) * t80 + Ifges(6,5) * t128;
t26 = Ifges(6,4) * t79 + Ifges(6,2) * t80 + Ifges(6,6) * t128;
t19 = -t80 * pkin(5) + t30;
t18 = -mrSges(7,2) * t128 + mrSges(7,3) * t23;
t17 = mrSges(7,1) * t128 - mrSges(7,3) * t22;
t9 = -mrSges(7,1) * t23 + mrSges(7,2) * t22;
t8 = Ifges(7,1) * t22 + Ifges(7,4) * t23 + Ifges(7,5) * t128;
t7 = Ifges(7,4) * t22 + Ifges(7,2) * t23 + Ifges(7,6) * t128;
t1 = [(mrSges(3,3) * t272 * t371 + (0.2e1 * mrSges(3,3) * t214 - t282) * t276 + ((t228 * t375 + Ifges(3,5) * t268 + (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t276) * t309) * t276 + (t229 * t375 + Ifges(4,5) * t227 + Ifges(5,5) * t163 - 0.2e1 * Ifges(3,6) * t268 + Ifges(4,6) * t226 - Ifges(5,6) * t162 + (-mrSges(3,1) * pkin(1) - Ifges(3,4) * t272) * t309 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) - t384) * t320) * t272) * qJD(2)) * t267 + (t6 + t25 - t69) * t162 + (t255 - 0.2e1 * t332 - 0.2e1 * t333) * t268 + 0.2e1 * m(4) * (t101 * t149 + t102 * t148 + t209 * t215) + 0.2e1 * m(3) * (t214 * t229 - t215 * t228) - t281 * t27 + (t40 + t61 - t103) * t128 + t227 * t118 + t226 * t117 + 0.2e1 * t209 * t135 + t194 * t154 + t195 * t155 + 0.2e1 * t101 * t198 + 0.2e1 * t102 * t199 + 0.2e1 * t168 * t78 + 0.2e1 * t149 * t159 + 0.2e1 * t148 * t160 + t163 * t70 + t140 * t26 + 0.2e1 * t33 * t145 + 0.2e1 * t32 * t146 + 0.2e1 * t147 * t106 + t129 * t104 + (-mrSges(4,1) * t226 + mrSges(4,2) * t227) * t371 + (t12 * t3 + t13 * t2 + t19 * t48) * t376 + (t10 * t35 + t11 * t34 + t30 * t64) * t377 + (t147 * t168 + t32 * t73 + t33 * t74) * t378 + 0.2e1 * t12 * t17 + 0.2e1 * t13 * t18 + t23 * t41 + t22 * t42 + 0.2e1 * t19 * t44 + 0.2e1 * t34 * t46 + 0.2e1 * t35 * t47 + 0.2e1 * t48 * t9 + 0.2e1 * t2 * t57 + 0.2e1 * t3 * t58 + 0.2e1 * t64 * t39 + t79 * t63 + t80 * t62 + t89 * t7 + t90 * t8 + 0.2e1 * t30 * t91 + 0.2e1 * t10 * t99 + 0.2e1 * t11 * t100 + 0.2e1 * t74 * t113 + 0.2e1 * t73 * t114; -t332 - t333 + (-t170 / 0.2e1 + t36 / 0.2e1 + t93 / 0.2e1) * t162 + (-t215 * mrSges(4,1) + t117 / 0.2e1 + t101 * mrSges(4,3) + pkin(9) * t159) * t275 + (t215 * mrSges(4,2) + t118 / 0.2e1 - t102 * mrSges(4,3) - pkin(9) * t160) * t271 + (t40 / 0.2e1 + t61 / 0.2e1 - t103 / 0.2e1 - t74 * mrSges(5,3)) * t224 + (t104 / 0.2e1 - t73 * mrSges(5,3) + t62 * t350 + t63 * t348) * t225 + (t70 / 0.2e1 - t32 * mrSges(5,3) + t26 * t350 + t27 * t348 + (-t274 * t62 / 0.2e1 + t63 * t350) * qJD(5)) * t234 - (-t69 / 0.2e1 + t6 / 0.2e1 + t25 / 0.2e1 - t33 * mrSges(5,3)) * t278 + ((-t265 / 0.2e1 - t220 / 0.2e1 + t219 / 0.2e1) * t276 + (Ifges(4,5) * t271 / 0.2e1 + Ifges(4,6) * t275 / 0.2e1 - Ifges(3,6) + Ifges(5,5) * t234 / 0.2e1 + Ifges(5,6) * t278 / 0.2e1) * t316) * t267 - t281 * t95 / 0.2e1 + m(5) * (t147 * t263 - t156 * t73 + t157 * t74 + t168 * t307 - t196 * t32 + t197 * t33) + ((-pkin(9) * t199 - t148 * mrSges(4,3) + t155 / 0.2e1) * t275 + (-pkin(9) * t198 - t149 * mrSges(4,3) + pkin(3) * t106 - t154 / 0.2e1 + t305 / 0.2e1) * t271) * qJD(3) + (t39 - t114) * t196 + (t91 - t146) * t156 + ((t101 * t275 - t102 * t271 + (-t148 * t275 - t149 * t271) * qJD(3)) * pkin(9) - pkin(2) * t215) * m(4) + m(7) * (t109 * t48 + t12 * t16 + t13 * t15 + t153 * t19 + t2 * t51 + t3 * t50) + m(6) * (t10 * t120 + t11 * t119 + t156 * t64 + t196 * t30 + t34 * t56 + t35 * t55) + t263 * t78 + t226 * t242 / 0.2e1 + t227 * t244 / 0.2e1 + t194 * t250 / 0.2e1 + t195 * t252 / 0.2e1 + t209 * t239 + t197 * t113 + t11 * t177 + t147 * t182 + t129 * t184 / 0.2e1 + t168 * t169 + t163 * t171 / 0.2e1 + t30 * t172 + t10 * t176 + t153 * t9 + t157 * t145 + t79 * t144 / 0.2e1 - pkin(2) * t135 + t34 * t136 + t35 * t137 + t2 * t138 + t3 * t139 + t8 * t363 + t7 * t364 + t94 * t365 + t38 * t367 + t37 * t368 + t143 * t369 + (-t183 / 0.2e1 + t142 / 0.2e1 + t96 / 0.2e1) * t128 + t255 + t48 * t43 + t50 * t17 + t51 * t18 + t15 * t57 + t16 * t58 + t12 * t67 + t13 * t68 + t83 * t42 / 0.2e1 + t84 * t41 / 0.2e1 + t23 * t97 / 0.2e1 + t22 * t98 / 0.2e1 + t55 * t99 + t56 * t100 + t19 * t107 + t109 * t44 + t64 * t110 + t119 * t46 + t120 * t47; (mrSges(5,3) * t372 - t143 * t270 + t144 * t274 + t184) * t225 + (mrSges(5,3) * t373 - t270 * t94 + t274 * t95 + t171 + (-t143 * t274 - t144 * t270) * qJD(5)) * t234 + (t197 * t374 + t142 - t183 + t96) * t224 - (t157 * t374 - t170 + t36 + t93) * t278 + t275 * t242 + t271 * t244 + 0.2e1 * t263 * t169 - 0.2e1 * pkin(2) * t239 + 0.2e1 * t56 * t177 - t164 * t37 - t165 * t38 + 0.2e1 * t55 * t176 + 0.2e1 * t153 * t43 + 0.2e1 * t119 * t136 + 0.2e1 * t120 * t137 + 0.2e1 * t15 * t138 + 0.2e1 * t16 * t139 + t110 * t372 + t172 * t373 + (t109 * t153 + t15 * t51 + t16 * t50) * t376 + (t119 * t56 + t120 * t55 + t327) * t377 + (t157 * t197 + t263 * t307 + t327) * t378 + (t275 * t252 + (0.2e1 * pkin(3) * t182 - t250) * t271) * qJD(3) + 0.2e1 * t50 * t67 + 0.2e1 * t51 * t68 + t84 * t97 + t83 * t98 + 0.2e1 * t109 * t107; t282 + t380 * t162 + t381 * t128 + m(6) * (t262 * t30 + ((-t270 * t35 - t274 * t34) * qJD(5) + t382) * t261) + t382 * mrSges(6,3) + t44 * t306 - t3 * t343 + t12 * t335 - t46 * t324 - t281 * t354 + m(7) * (t111 * t13 + t112 * t12 + t173 * t3 + t174 * t2 + t19 * t245 + t306 * t48) - t34 * t303 - t35 * t304 - t100 * t298 - t99 * t299 + t262 * t39 + t245 * t9 + t30 * t246 + t64 * t238 + t19 * t190 + t173 * t17 + t174 * t18 + t48 * t131 + t114 * t302 + t63 * t294 + t62 * t295 + t47 * t323 - t2 * t331 - t13 * t334 + t113 * t346 + t26 * t348 + t27 * t349 + t79 * t351 + t8 * t356 + t7 * t357 + t22 * t358 + t23 * t359 + t41 * t361 + t42 * t362 + t241 * t365 + t134 * t367 + t133 * t368 + t249 * t369 + (t266 * t33 + t32 * t328) * t370 + t32 * mrSges(5,1) - t33 * mrSges(5,2) - t101 * mrSges(4,2) + t102 * mrSges(4,1) + t111 * t57 + t112 * t58; (t234 * t251 + t143) * t295 + (t266 * t370 - mrSges(5,2)) * t157 + (m(6) * t262 - t328 * t370 - mrSges(5,1) + t246) * t156 + (-mrSges(4,1) * t314 + mrSges(4,2) * t315) * pkin(9) + (-mrSges(5,3) * t346 + t381) * t224 + t107 * t306 + t383 * mrSges(6,3) - t280 * t249 / 0.2e1 - t16 * t343 + t50 * t335 - t136 * t324 - t241 * t326 / 0.2e1 - Ifges(4,6) * t315 - t380 * t278 + m(7) * (t109 * t245 + t111 * t51 + t112 * t50 + t15 * t174 + t153 * t306 + t16 * t173) - t119 * t303 - t120 * t304 - t177 * t298 - t176 * t299 + t262 * t110 + t245 * t43 + t196 * t238 + t109 * t190 + t173 * t67 + t174 * t68 + t153 * t131 + t111 * t138 + t112 * t139 + m(6) * ((-t119 * t274 - t120 * t270) * qJD(5) + t383) * t261 - t219 + t220 + t144 * t294 + t137 * t323 - t15 * t331 - t51 * t334 + t94 * t348 + t95 * t349 + t318 * t351 + t325 * t354 + t38 * t356 + t37 * t357 + t83 * t358 + t84 * t359 + t97 * t361 + t98 * t362 + t134 * t363 + t133 * t364 - t225 * mrSges(5,3) * t302 + t265; -t188 * t193 + t237 * t134 - t189 * t192 - t286 * t133 + (t111 * t174 + t112 * t173 + t245 * t306) * t376 + 0.2e1 * t190 * t306 + 0.2e1 * t245 * t131 + 0.2e1 * t262 * t238 - t249 * t313 + t270 * t243 + (qJD(5) * t251 + t241) * t274 + 0.2e1 * (-t111 * t286 - t112 * t237 + t173 * t188 - t174 * t189) * mrSges(7,3); -t286 * t17 + t237 * t18 - t188 * t57 - t189 * t58 + t270 * t47 + t274 * t46 + (-t100 * t270 + t274 * t99) * qJD(5) + m(7) * (-t12 * t189 - t13 * t188 + t2 * t237 - t286 * t3) + m(6) * (t10 * t270 + t11 * t274 + (-t270 * t34 + t274 * t35) * qJD(5)) + m(5) * t147 + t78; m(5) * t307 + t274 * t136 + t270 * t137 - t188 * t138 - t189 * t139 - t286 * t67 + t237 * t68 + (t176 * t274 - t177 * t270) * qJD(5) + m(7) * (t15 * t237 - t16 * t286 - t188 * t51 - t189 * t50) + m(6) * (t270 * t55 + t274 * t56 + (-t119 * t270 + t120 * t274) * qJD(5)) + t169; m(7) * (t111 * t237 - t112 * t286 - t173 * t189 - t174 * t188); (-t188 * t237 + t189 * t286) * t376; t11 * mrSges(6,1) - t10 * mrSges(6,2) + (m(7) * (-t12 * t311 + t13 * t310 + t2 * t269 + t273 * t3) + t57 * t310 + t269 * t18 - t58 * t311 + t273 * t17) * pkin(5) + t284 + t25; t56 * mrSges(6,1) - t55 * mrSges(6,2) + (m(7) * (t15 * t269 + t16 * t273 + t310 * t51 - t311 * t50) + t138 * t310 + t269 * t68 - t139 * t311 + t273 * t67) * pkin(5) + t93 + t283; t264 + (t246 * t261 - t338) * qJD(5) + (m(7) * (t111 * t269 + t112 * t273 + (-t173 * t269 + t174 * t273) * qJD(6)) + (t273 * t188 - t269 * t189 + (t237 * t269 - t273 * t286) * qJD(6)) * mrSges(7,3)) * pkin(5) + t285; -t238 + m(7) * (-t188 * t269 - t189 * t273 + (t237 * t273 + t269 * t286) * qJD(6)) * pkin(5) - t131; 0.2e1 * t232; t284; t283; t285; -t131; t232; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
