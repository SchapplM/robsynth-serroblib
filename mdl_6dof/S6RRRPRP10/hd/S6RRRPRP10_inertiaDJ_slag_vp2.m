% Calculate time derivative of joint inertia matrix for
% S6RRRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-03-09 17:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRP10_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP10_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP10_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP10_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP10_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP10_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP10_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:29:33
% EndTime: 2019-03-09 17:29:48
% DurationCPUTime: 6.06s
% Computational Cost: add. (9499->720), mult. (25626->1020), div. (0->0), fcn. (24441->10), ass. (0->277)
t278 = cos(pkin(6));
t280 = sin(qJ(3));
t283 = cos(qJ(3));
t276 = sin(pkin(6));
t281 = sin(qJ(2));
t331 = t276 * t281;
t226 = t278 * t280 + t283 * t331;
t284 = cos(qJ(2));
t326 = qJD(2) * t276;
t306 = t284 * t326;
t177 = qJD(3) * t226 + t280 * t306;
t225 = -t278 * t283 + t280 * t331;
t178 = -qJD(3) * t225 + t283 * t306;
t275 = sin(pkin(11));
t277 = cos(pkin(11));
t325 = qJD(2) * t281;
t307 = t276 * t325;
t128 = -t178 * t275 + t277 * t307;
t129 = t178 * t277 + t275 * t307;
t279 = sin(qJ(5));
t282 = cos(qJ(5));
t330 = t276 * t284;
t175 = -t226 * t275 - t277 * t330;
t176 = t226 * t277 - t275 * t330;
t286 = t282 * t175 - t176 * t279;
t41 = qJD(5) * t286 + t128 * t279 + t129 * t282;
t97 = t175 * t279 + t176 * t282;
t42 = qJD(5) * t97 - t282 * t128 + t129 * t279;
t7 = Ifges(6,5) * t41 - Ifges(6,6) * t42 + Ifges(6,3) * t177;
t8 = Ifges(7,4) * t41 + Ifges(7,2) * t177 + Ifges(7,6) * t42;
t363 = t7 + t8;
t348 = t277 / 0.2e1;
t334 = t275 * t279;
t236 = -t282 * t277 + t334;
t237 = t275 * t282 + t277 * t279;
t321 = qJD(5) * t280;
t323 = qJD(3) * t283;
t139 = -t236 * t323 - t237 * t321;
t320 = qJD(5) * t282;
t329 = t277 * t280;
t140 = t237 * t323 + t320 * t329 - t321 * t334;
t324 = qJD(3) * t280;
t69 = Ifges(6,5) * t139 - Ifges(6,6) * t140 + Ifges(6,3) * t324;
t70 = Ifges(7,4) * t139 + Ifges(7,2) * t324 + Ifges(7,6) * t140;
t361 = -t69 - t70;
t244 = -pkin(3) * t283 - qJ(4) * t280 - pkin(2);
t234 = t277 * t244;
t158 = -pkin(10) * t329 + t234 + (-pkin(9) * t275 - pkin(4)) * t283;
t328 = t277 * t283;
t198 = pkin(9) * t328 + t275 * t244;
t333 = t275 * t280;
t179 = -pkin(10) * t333 + t198;
t360 = t279 * t158 + t282 * t179;
t260 = pkin(8) * t331;
t346 = pkin(1) * t284;
t230 = t278 * t346 - t260;
t359 = m(7) * qJ(6) + mrSges(7,3);
t222 = -qJD(4) * t280 + (pkin(3) * t280 - qJ(4) * t283) * qJD(3);
t315 = pkin(9) * t324;
t173 = t277 * t222 + t275 * t315;
t127 = (pkin(4) * t280 - pkin(10) * t328) * qJD(3) + t173;
t202 = t275 * t222;
t332 = t275 * t283;
t146 = t202 + (-pkin(9) * t329 - pkin(10) * t332) * qJD(3);
t44 = -qJD(5) * t360 + t127 * t282 - t146 * t279;
t231 = t278 * t281 * pkin(1) + pkin(8) * t330;
t208 = pkin(9) * t278 + t231;
t209 = (-pkin(2) * t284 - pkin(9) * t281 - pkin(1)) * t276;
t213 = (pkin(2) * t281 - pkin(9) * t284) * t326;
t214 = t230 * qJD(2);
t75 = -t208 * t324 + t209 * t323 + t280 * t213 + t283 * t214;
t63 = (qJ(4) * t325 - qJD(4) * t284) * t276 + t75;
t215 = t231 * qJD(2);
t67 = t177 * pkin(3) - t178 * qJ(4) - t226 * qJD(4) + t215;
t23 = -t275 * t63 + t277 * t67;
t19 = pkin(4) * t177 - pkin(10) * t129 + t23;
t24 = t275 * t67 + t277 * t63;
t21 = pkin(10) * t128 + t24;
t207 = t260 + (-pkin(2) - t346) * t278;
t103 = t225 * pkin(3) - t226 * qJ(4) + t207;
t117 = t283 * t208 + t280 * t209;
t104 = -qJ(4) * t330 + t117;
t60 = t277 * t103 - t104 * t275;
t40 = pkin(4) * t225 - pkin(10) * t176 + t60;
t61 = t275 * t103 + t277 * t104;
t52 = pkin(10) * t175 + t61;
t344 = t279 * t40 + t282 * t52;
t4 = -qJD(5) * t344 + t19 * t282 - t21 * t279;
t358 = 2 * m(5);
t357 = 2 * m(6);
t356 = 2 * m(7);
t355 = 0.2e1 * pkin(9);
t354 = -2 * mrSges(3,3);
t353 = t128 / 0.2e1;
t352 = t129 / 0.2e1;
t339 = Ifges(5,4) * t277;
t292 = -Ifges(5,2) * t275 + t339;
t192 = (Ifges(5,6) * t280 + t283 * t292) * qJD(3);
t351 = t192 / 0.2e1;
t340 = Ifges(5,4) * t275;
t293 = Ifges(5,1) * t277 - t340;
t193 = (Ifges(5,5) * t280 + t283 * t293) * qJD(3);
t350 = t193 / 0.2e1;
t349 = -t275 / 0.2e1;
t347 = m(5) * t283;
t272 = t280 * pkin(9);
t345 = pkin(10) + qJ(4);
t343 = mrSges(5,2) * t277;
t342 = Ifges(4,4) * t280;
t341 = Ifges(4,4) * t283;
t338 = t214 * mrSges(3,2);
t337 = t215 * mrSges(3,1);
t336 = t215 * mrSges(4,1);
t335 = t215 * mrSges(4,2);
t223 = t236 * qJD(5);
t224 = t237 * qJD(5);
t152 = -Ifges(6,5) * t223 - Ifges(6,6) * t224;
t153 = -Ifges(7,4) * t223 + Ifges(7,6) * t224;
t305 = t275 * t323;
t212 = mrSges(5,1) * t305 + t323 * t343;
t271 = pkin(9) * t323;
t232 = pkin(4) * t305 + t271;
t243 = pkin(4) * t333 + t272;
t322 = qJD(5) * t279;
t319 = t279 * qJD(4);
t318 = t282 * qJD(4);
t6 = Ifges(7,5) * t41 + Ifges(7,6) * t177 + Ifges(7,3) * t42;
t9 = Ifges(6,4) * t41 - Ifges(6,2) * t42 + Ifges(6,6) * t177;
t317 = t6 / 0.2e1 - t9 / 0.2e1;
t314 = Ifges(4,6) * t330;
t10 = Ifges(7,1) * t41 + Ifges(7,4) * t177 + Ifges(7,5) * t42;
t11 = Ifges(6,1) * t41 - Ifges(6,4) * t42 + Ifges(6,5) * t177;
t313 = t10 / 0.2e1 + t11 / 0.2e1;
t45 = Ifges(7,5) * t97 + Ifges(7,6) * t225 - Ifges(7,3) * t286;
t48 = Ifges(6,4) * t97 + Ifges(6,2) * t286 + Ifges(6,6) * t225;
t312 = -t48 / 0.2e1 + t45 / 0.2e1;
t49 = Ifges(7,1) * t97 + Ifges(7,4) * t225 - Ifges(7,5) * t286;
t50 = Ifges(6,1) * t97 + Ifges(6,4) * t286 + Ifges(6,5) * t225;
t311 = t49 / 0.2e1 + t50 / 0.2e1;
t68 = Ifges(7,5) * t139 + Ifges(7,6) * t324 + Ifges(7,3) * t140;
t71 = Ifges(6,4) * t139 - Ifges(6,2) * t140 + Ifges(6,6) * t324;
t310 = t68 / 0.2e1 - t71 / 0.2e1;
t72 = Ifges(7,1) * t139 + Ifges(7,4) * t324 + Ifges(7,5) * t140;
t73 = Ifges(6,1) * t139 - Ifges(6,4) * t140 + Ifges(6,5) * t324;
t309 = t73 / 0.2e1 + t72 / 0.2e1;
t308 = Ifges(4,5) * t178 - Ifges(4,6) * t177 + Ifges(4,3) * t307;
t267 = -pkin(4) * t277 - pkin(3);
t210 = t237 * t280;
t211 = t236 * t280;
t121 = -Ifges(7,5) * t211 - Ifges(7,6) * t283 + Ifges(7,3) * t210;
t124 = -Ifges(6,4) * t211 - Ifges(6,2) * t210 - Ifges(6,6) * t283;
t304 = -t124 / 0.2e1 + t121 / 0.2e1;
t125 = -Ifges(7,1) * t211 - Ifges(7,4) * t283 + Ifges(7,5) * t210;
t126 = -Ifges(6,1) * t211 - Ifges(6,4) * t210 - Ifges(6,5) * t283;
t303 = t125 / 0.2e1 + t126 / 0.2e1;
t151 = -Ifges(7,5) * t223 + Ifges(7,3) * t224;
t154 = -Ifges(6,4) * t223 - Ifges(6,2) * t224;
t302 = t151 / 0.2e1 - t154 / 0.2e1;
t301 = t152 / 0.2e1 + t153 / 0.2e1;
t155 = -Ifges(7,1) * t223 + Ifges(7,5) * t224;
t156 = -Ifges(6,1) * t223 - Ifges(6,4) * t224;
t300 = t155 / 0.2e1 + t156 / 0.2e1;
t163 = Ifges(7,5) * t237 + Ifges(7,3) * t236;
t166 = Ifges(6,4) * t237 - Ifges(6,2) * t236;
t299 = t163 / 0.2e1 - t166 / 0.2e1;
t167 = Ifges(7,1) * t237 + Ifges(7,5) * t236;
t168 = Ifges(6,1) * t237 - Ifges(6,4) * t236;
t298 = t167 / 0.2e1 + t168 / 0.2e1;
t15 = t42 * mrSges(6,1) + t41 * mrSges(6,2);
t14 = t42 * mrSges(7,1) - t41 * mrSges(7,3);
t297 = t345 * t275;
t26 = -t177 * mrSges(7,1) + t41 * mrSges(7,2);
t296 = qJD(5) * t345;
t74 = -t128 * mrSges(5,1) + t129 * mrSges(5,2);
t78 = t140 * mrSges(6,1) + t139 * mrSges(6,2);
t150 = t224 * mrSges(6,1) - t223 * mrSges(6,2);
t77 = t140 * mrSges(7,1) - t139 * mrSges(7,3);
t149 = t224 * mrSges(7,1) + t223 * mrSges(7,3);
t247 = t345 * t277;
t118 = t277 * t318 - t247 * t322 + (-t282 * t296 - t319) * t275;
t119 = t277 * t319 + t247 * t320 + (-t279 * t296 + t318) * t275;
t180 = t247 * t279 + t282 * t297;
t181 = t282 * t247 - t279 * t297;
t295 = t181 * t118 + t119 * t180;
t116 = -t280 * t208 + t209 * t283;
t110 = -mrSges(7,1) * t324 + t139 * mrSges(7,2);
t105 = pkin(3) * t330 - t116;
t294 = Ifges(5,5) * t275 / 0.2e1 + Ifges(5,6) * t348 + (Ifges(6,5) + Ifges(7,4)) * t237 / 0.2e1 + (-Ifges(6,6) / 0.2e1 + Ifges(7,6) / 0.2e1) * t236;
t291 = Ifges(5,5) * t277 - Ifges(5,6) * t275;
t16 = -t279 * t52 + t282 * t40;
t94 = t158 * t282 - t179 * t279;
t76 = -t208 * t323 - t209 * t324 + t213 * t283 - t280 * t214;
t3 = t279 * t19 + t282 * t21 + t40 * t320 - t322 * t52;
t43 = t279 * t127 + t282 * t146 + t158 * t320 - t179 * t322;
t83 = -pkin(4) * t175 + t105;
t66 = -pkin(3) * t307 - t76;
t53 = -pkin(4) * t128 + t66;
t270 = Ifges(4,5) * t323;
t256 = Ifges(3,5) * t306;
t252 = Ifges(4,1) * t280 + t341;
t251 = Ifges(4,2) * t283 + t342;
t250 = Ifges(5,1) * t275 + t339;
t249 = Ifges(5,2) * t277 + t340;
t246 = -mrSges(5,1) * t277 + mrSges(5,2) * t275;
t242 = (Ifges(4,1) * t283 - t342) * qJD(3);
t241 = (-Ifges(4,2) * t280 + t341) * qJD(3);
t240 = (mrSges(4,1) * t280 + mrSges(4,2) * t283) * qJD(3);
t239 = -mrSges(5,1) * t283 - mrSges(5,3) * t329;
t238 = mrSges(5,2) * t283 - mrSges(5,3) * t333;
t229 = (mrSges(5,1) * t280 - mrSges(5,3) * t328) * qJD(3);
t228 = (-mrSges(5,2) * t280 - mrSges(5,3) * t332) * qJD(3);
t227 = (mrSges(5,1) * t275 + t343) * t280;
t206 = -Ifges(5,5) * t283 + t280 * t293;
t205 = -Ifges(5,6) * t283 + t280 * t292;
t204 = -Ifges(5,3) * t283 + t280 * t291;
t197 = -pkin(9) * t332 + t234;
t191 = (Ifges(5,3) * t280 + t283 * t291) * qJD(3);
t187 = mrSges(7,1) * t283 - mrSges(7,2) * t211;
t186 = -mrSges(6,1) * t283 + mrSges(6,3) * t211;
t185 = mrSges(6,2) * t283 - mrSges(6,3) * t210;
t184 = -mrSges(7,2) * t210 - mrSges(7,3) * t283;
t183 = -mrSges(4,1) * t330 - t226 * mrSges(4,3);
t182 = mrSges(4,2) * t330 - t225 * mrSges(4,3);
t174 = -t277 * t315 + t202;
t162 = mrSges(6,1) * t236 + mrSges(6,2) * t237;
t161 = mrSges(7,1) * t236 - mrSges(7,3) * t237;
t157 = pkin(5) * t236 - qJ(6) * t237 + t267;
t145 = mrSges(4,1) * t307 - mrSges(4,3) * t178;
t144 = -mrSges(4,2) * t307 - mrSges(4,3) * t177;
t143 = mrSges(6,1) * t210 - mrSges(6,2) * t211;
t142 = mrSges(7,1) * t210 + mrSges(7,3) * t211;
t138 = Ifges(4,1) * t226 - Ifges(4,4) * t225 - Ifges(4,5) * t330;
t137 = Ifges(4,4) * t226 - Ifges(4,2) * t225 - t314;
t123 = -Ifges(7,4) * t211 - Ifges(7,2) * t283 + Ifges(7,6) * t210;
t122 = -Ifges(6,5) * t211 - Ifges(6,6) * t210 - Ifges(6,3) * t283;
t113 = mrSges(5,1) * t225 - mrSges(5,3) * t176;
t112 = -mrSges(5,2) * t225 + mrSges(5,3) * t175;
t111 = -mrSges(6,2) * t324 - mrSges(6,3) * t140;
t109 = mrSges(6,1) * t324 - mrSges(6,3) * t139;
t108 = -mrSges(7,2) * t140 + mrSges(7,3) * t324;
t107 = pkin(5) * t224 + qJ(6) * t223 - qJD(6) * t237;
t106 = pkin(5) * t210 + qJ(6) * t211 + t243;
t99 = mrSges(4,1) * t177 + mrSges(4,2) * t178;
t98 = -mrSges(5,1) * t175 + mrSges(5,2) * t176;
t93 = Ifges(4,1) * t178 - Ifges(4,4) * t177 + Ifges(4,5) * t307;
t92 = Ifges(4,4) * t178 - Ifges(4,2) * t177 + Ifges(4,6) * t307;
t91 = pkin(5) * t283 - t94;
t89 = -qJ(6) * t283 + t360;
t88 = mrSges(5,1) * t177 - mrSges(5,3) * t129;
t87 = -mrSges(5,2) * t177 + mrSges(5,3) * t128;
t86 = Ifges(5,1) * t176 + Ifges(5,4) * t175 + Ifges(5,5) * t225;
t85 = Ifges(5,4) * t176 + Ifges(5,2) * t175 + Ifges(5,6) * t225;
t84 = Ifges(5,5) * t176 + Ifges(5,6) * t175 + Ifges(5,3) * t225;
t82 = -mrSges(7,1) * t225 + mrSges(7,2) * t97;
t81 = mrSges(6,1) * t225 - mrSges(6,3) * t97;
t80 = -mrSges(6,2) * t225 + mrSges(6,3) * t286;
t79 = mrSges(7,2) * t286 + mrSges(7,3) * t225;
t59 = pkin(5) * t140 - qJ(6) * t139 + qJD(6) * t211 + t232;
t58 = Ifges(5,1) * t129 + Ifges(5,4) * t128 + Ifges(5,5) * t177;
t57 = Ifges(5,4) * t129 + Ifges(5,2) * t128 + Ifges(5,6) * t177;
t56 = Ifges(5,5) * t129 + Ifges(5,6) * t128 + Ifges(5,3) * t177;
t55 = -mrSges(6,1) * t286 + mrSges(6,2) * t97;
t54 = -mrSges(7,1) * t286 - mrSges(7,3) * t97;
t47 = Ifges(7,4) * t97 + Ifges(7,2) * t225 - Ifges(7,6) * t286;
t46 = Ifges(6,5) * t97 + Ifges(6,6) * t286 + Ifges(6,3) * t225;
t32 = -pkin(5) * t324 - t44;
t29 = qJ(6) * t324 - qJD(6) * t283 + t43;
t28 = -mrSges(7,2) * t42 + mrSges(7,3) * t177;
t27 = -mrSges(6,2) * t177 - mrSges(6,3) * t42;
t25 = mrSges(6,1) * t177 - mrSges(6,3) * t41;
t22 = -pkin(5) * t286 - qJ(6) * t97 + t83;
t13 = -pkin(5) * t225 - t16;
t12 = qJ(6) * t225 + t344;
t5 = pkin(5) * t42 - qJ(6) * t41 - qJD(6) * t97 + t53;
t2 = -pkin(5) * t177 - t4;
t1 = qJ(6) * t177 + qJD(6) * t225 + t3;
t17 = [(-t284 * t308 + 0.2e1 * (t214 * t284 + t215 * t281) * mrSges(3,3) + ((t230 * t354 + Ifges(3,5) * t278 + 0.2e1 * (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t284) * t276) * t284 + (t231 * t354 + Ifges(4,5) * t226 - 0.2e1 * Ifges(3,6) * t278 - Ifges(4,6) * t225 + (-0.2e1 * pkin(1) * mrSges(3,1) - 0.2e1 * Ifges(3,4) * t281 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) - Ifges(4,3)) * t284) * t276) * t281) * qJD(2)) * t276 + 0.2e1 * m(3) * (t214 * t231 - t215 * t230) + 0.2e1 * m(4) * (t116 * t76 + t117 * t75 + t207 * t215) + (t10 + t11) * t97 + (-t48 + t45) * t42 + (t49 + t50) * t41 + (t105 * t66 + t23 * t60 + t24 * t61) * t358 + (t1 * t12 + t13 * t2 + t22 * t5) * t356 - (t6 - t9) * t286 + (t16 * t4 + t3 * t344 + t53 * t83) * t357 + 0.2e1 * t344 * t27 + (-t137 + t84 + t46 + t47) * t177 + 0.2e1 * t207 * t99 + (t93 + 0.2e1 * t335) * t226 + (t256 - 0.2e1 * t337 - 0.2e1 * t338) * t278 + t178 * t138 + 0.2e1 * t75 * t182 + 0.2e1 * t76 * t183 + t175 * t57 + t176 * t58 + 0.2e1 * t117 * t144 + 0.2e1 * t116 * t145 + t128 * t85 + t129 * t86 + 0.2e1 * t24 * t112 + 0.2e1 * t23 * t113 + 0.2e1 * t105 * t74 + 0.2e1 * t66 * t98 + 0.2e1 * t1 * t79 + 0.2e1 * t3 * t80 + 0.2e1 * t4 * t81 + 0.2e1 * t2 * t82 + 0.2e1 * t83 * t15 + 0.2e1 * t61 * t87 + 0.2e1 * t60 * t88 + 0.2e1 * t53 * t55 + 0.2e1 * t5 * t54 + 0.2e1 * t16 * t25 + 0.2e1 * t13 * t26 + 0.2e1 * t12 * t28 + 0.2e1 * t22 * t14 + (t56 - t92 + 0.2e1 * t336 + t363) * t225; t256 - t337 - t338 + (t283 * t144 + (-t145 + t74) * t280) * pkin(9) + t176 * t350 + t175 * t351 + t206 * t352 + t205 * t353 + m(6) * (t16 * t44 + t232 * t83 + t243 * t53 + t3 * t360 + t344 * t43 + t4 * t94) + t360 * t27 - t310 * t286 + t344 * t111 + m(7) * (t1 * t89 + t106 * t5 + t12 * t29 + t13 * t32 + t2 * t91 + t22 * t59) + (t58 * t348 + t57 * t349 - t76 * mrSges(4,3) + t335 + t93 / 0.2e1) * t280 + t178 * t252 / 0.2e1 + t24 * t238 + t23 * t239 + t207 * t240 + t226 * t242 / 0.2e1 + t243 * t15 + ((t138 / 0.2e1 + t86 * t348 + t85 * t349 - t116 * mrSges(4,3)) * t283 + (t314 / 0.2e1 - t137 / 0.2e1 + t84 / 0.2e1 + t46 / 0.2e1 + t47 / 0.2e1 - t117 * mrSges(4,3)) * t280 + (-t280 * t182 + (-t183 + t98) * t283 + m(4) * (-t116 * t283 - t117 * t280) + t105 * t347) * pkin(9)) * qJD(3) + t66 * t227 + t61 * t228 + t60 * t229 + t232 * t55 + m(4) * (pkin(9) * t283 * t75 - pkin(2) * t215 - t272 * t76) + m(5) * (t173 * t60 + t174 * t61 + t197 * t23 + t198 * t24 + t272 * t66) + t105 * t212 + t4 * t186 + t2 * t187 + t197 * t88 + t198 * t87 + (t75 * mrSges(4,3) - t336 + t92 / 0.2e1 - t56 / 0.2e1 - t7 / 0.2e1 - t8 / 0.2e1) * t283 + t1 * t184 + t3 * t185 + t173 * t113 + t174 * t112 + t5 * t142 + t53 * t143 + (-t284 * t270 / 0.2e1 + (Ifges(4,5) * t280 / 0.2e1 + Ifges(4,6) * t283 / 0.2e1 - Ifges(3,6)) * t325) * t276 + t106 * t14 + t12 * t108 + t16 * t109 + t13 * t110 - pkin(2) * t99 + t89 * t28 + t91 * t26 + t94 * t25 + t29 * t79 + t43 * t80 + t44 * t81 + t32 * t82 + t83 * t78 + t22 * t77 + t317 * t210 + t311 * t139 + t312 * t140 - t313 * t211 + t309 * t97 + t59 * t54 + t303 * t41 + t304 * t42 + (-t251 / 0.2e1 + t204 / 0.2e1 + t122 / 0.2e1 + t123 / 0.2e1) * t177 + (-t241 / 0.2e1 + t191 / 0.2e1 + t69 / 0.2e1 + t70 / 0.2e1) * t225; -(t73 + t72) * t211 + (-t71 + t68) * t210 + (t121 - t124) * t140 + (t125 + t126) * t139 + (t173 * t197 + t174 * t198) * t358 + (t106 * t59 + t29 * t89 + t32 * t91) * t356 + (t232 * t243 + t360 * t43 + t44 * t94) * t357 + 0.2e1 * t360 * t111 + ((-t205 * t275 + t206 * t277 + t227 * t355 + t252) * t283 + (0.2e1 * pkin(9) ^ 2 * t347 + t122 + t123 + t204 - t251) * t280) * qJD(3) + (-t192 * t275 + t193 * t277 + t212 * t355 + t242) * t280 + 0.2e1 * t174 * t238 + 0.2e1 * t173 * t239 - 0.2e1 * pkin(2) * t240 + 0.2e1 * t243 * t78 + 0.2e1 * t198 * t228 + 0.2e1 * t197 * t229 + 0.2e1 * t232 * t143 + 0.2e1 * t44 * t186 + 0.2e1 * t32 * t187 + 0.2e1 * t29 * t184 + 0.2e1 * t43 * t185 + 0.2e1 * t59 * t142 + 0.2e1 * t106 * t77 + 0.2e1 * t89 * t108 + 0.2e1 * t94 * t109 + 0.2e1 * t91 * t110 + (t241 - t191 + t361) * t283; m(5) * (-pkin(3) * t66 + (-t275 * t60 + t277 * t61) * qJD(4) + (-t23 * t275 + t24 * t277) * qJ(4)) + (qJD(4) * t112 + qJ(4) * t87 + t24 * mrSges(5,3) + t57 / 0.2e1) * t277 + t308 + (-qJD(4) * t113 - qJ(4) * t88 - t23 * mrSges(5,3) + t58 / 0.2e1) * t275 + (t27 + t28) * t181 + (t26 - t25) * t180 + (t82 - t81) * t119 + (t79 + t80) * t118 + t250 * t352 + t249 * t353 - t302 * t286 + m(6) * (t118 * t344 - t119 * t16 - t180 * t4 + t181 * t3 + t267 * t53) + t267 * t15 + t66 * t246 + t5 * t161 + t53 * t162 + t157 * t14 + t22 * t149 + t83 * t150 + t107 * t54 - pkin(3) * t74 - t75 * mrSges(4,2) + t76 * mrSges(4,1) + t317 * t236 - t311 * t223 + t312 * t224 + t313 * t237 + t298 * t41 + t299 * t42 + t300 * t97 + t301 * t225 + t294 * t177 + m(7) * (t1 * t181 + t107 * t22 + t118 * t12 + t119 * t13 + t157 * t5 + t180 * t2) + (-t1 * t236 - t12 * t224 - t13 * t223 + t2 * t237) * mrSges(7,2) + (t16 * t223 - t224 * t344 - t236 * t3 - t237 * t4) * mrSges(6,3); t270 + (m(5) * (qJ(4) * t174 + qJD(4) * t198) + t174 * mrSges(5,3) + qJ(4) * t228 + qJD(4) * t238 + t351) * t277 + (m(5) * (-qJ(4) * t173 - qJD(4) * t197) - t173 * mrSges(5,3) - qJ(4) * t229 - qJD(4) * t239 + t350) * t275 + (t108 + t111) * t181 + m(7) * (t106 * t107 + t118 * t89 + t119 * t91 + t157 * t59 + t180 * t32 + t181 * t29) + m(6) * (t118 * t360 - t119 * t94 - t180 * t44 + t181 * t43 + t232 * t267) + (t223 * t94 - t224 * t360 - t236 * t43 - t237 * t44) * mrSges(6,3) + (t110 - t109) * t180 + (-t186 + t187) * t119 + (t184 + t185) * t118 + t267 * t78 + t243 * t150 + ((t250 * t348 + t249 * t349 + (-m(5) * pkin(3) - mrSges(4,1) + t246) * pkin(9)) * t283 + (pkin(9) * mrSges(4,2) - Ifges(4,6) + t294) * t280) * qJD(3) + t232 * t162 - pkin(3) * t212 + t59 * t161 + t157 * t77 + t107 * t142 + t106 * t149 + t310 * t236 + t309 * t237 + t298 * t139 + t299 * t140 - t300 * t211 - t301 * t283 + t302 * t210 - t303 * t223 + t304 * t224 + (-t223 * t91 - t224 * t89 - t236 * t29 + t237 * t32) * mrSges(7,2); 0.2e1 * t107 * t161 + 0.2e1 * t157 * t149 + 0.2e1 * t267 * t150 + (t156 + t155) * t237 + (t151 - t154) * t236 + (t163 - t166) * t224 - (t168 + t167) * t223 + t295 * t357 + (t107 * t157 + t295) * t356 + (qJ(4) * t358 + 0.2e1 * mrSges(5,3)) * qJD(4) * (t275 ^ 2 + t277 ^ 2) + 0.2e1 * (mrSges(6,3) + mrSges(7,2)) * (-t118 * t236 + t119 * t237 - t180 * t223 - t181 * t224); m(5) * t66 + m(6) * t53 + m(7) * t5 + t14 + t15 + t74; m(5) * t271 + m(6) * t232 + m(7) * t59 + t212 + t77 + t78; m(7) * t107 + t149 + t150; 0; m(7) * (-pkin(5) * t2 + qJ(6) * t1 + qJD(6) * t12) - pkin(5) * t26 + t1 * mrSges(7,3) + qJD(6) * t79 + qJ(6) * t28 - t3 * mrSges(6,2) + t4 * mrSges(6,1) - t2 * mrSges(7,1) + t363; -t43 * mrSges(6,2) + t44 * mrSges(6,1) - t32 * mrSges(7,1) - pkin(5) * t110 + t29 * mrSges(7,3) + qJD(6) * t184 + qJ(6) * t108 + m(7) * (-pkin(5) * t32 + qJ(6) * t29 + qJD(6) * t89) - t361; m(7) * qJD(6) * t181 + (pkin(5) * t223 - qJ(6) * t224 - qJD(6) * t236) * mrSges(7,2) + t153 + t152 + (-m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1)) * t119 + (-mrSges(6,2) + t359) * t118; 0; 0.2e1 * t359 * qJD(6); m(7) * t2 + t26; m(7) * t32 + t110; m(7) * t119 - t223 * mrSges(7,2); 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t17(1) t17(2) t17(4) t17(7) t17(11) t17(16); t17(2) t17(3) t17(5) t17(8) t17(12) t17(17); t17(4) t17(5) t17(6) t17(9) t17(13) t17(18); t17(7) t17(8) t17(9) t17(10) t17(14) t17(19); t17(11) t17(12) t17(13) t17(14) t17(15) t17(20); t17(16) t17(17) t17(18) t17(19) t17(20) t17(21);];
Mq  = res;
