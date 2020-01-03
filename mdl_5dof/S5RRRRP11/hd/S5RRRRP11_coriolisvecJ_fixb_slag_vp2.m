% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 22:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRRP11_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP11_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP11_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP11_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP11_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP11_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP11_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:14:26
% EndTime: 2019-12-31 22:15:00
% DurationCPUTime: 14.50s
% Computational Cost: add. (7398->629), mult. (19979->851), div. (0->0), fcn. (14680->8), ass. (0->267)
t361 = Ifges(5,1) + Ifges(6,1);
t360 = Ifges(6,4) + Ifges(5,5);
t212 = sin(qJ(3));
t216 = cos(qJ(2));
t209 = sin(pkin(5));
t282 = qJD(1) * t209;
t267 = t216 * t282;
t279 = qJD(3) * t212;
t373 = t212 * t267 - t279;
t359 = Ifges(6,5) - Ifges(5,4);
t370 = Ifges(5,6) - Ifges(6,6);
t215 = cos(qJ(3));
t259 = pkin(3) * t212 - pkin(9) * t215;
t191 = t259 * qJD(3);
t199 = -pkin(3) * t215 - pkin(9) * t212 - pkin(2);
t211 = sin(qJ(4));
t214 = cos(qJ(4));
t275 = qJD(4) * t215;
t276 = qJD(4) * t214;
t108 = t211 * t191 + t199 * t276 + (-t211 * t275 - t214 * t279) * pkin(8);
t213 = sin(qJ(2));
t268 = t213 * t282;
t210 = cos(pkin(5));
t320 = pkin(1) * t216;
t274 = t210 * t320;
t176 = -pkin(7) * t268 + qJD(1) * t274;
t231 = t209 * (pkin(2) * t213 - pkin(8) * t216);
t177 = qJD(1) * t231;
t115 = t215 * t176 + t212 * t177;
t102 = pkin(9) * t268 + t115;
t206 = t210 * t213 * pkin(1);
t287 = t209 * t216;
t118 = (t206 + (pkin(7) + t259) * t287) * qJD(1);
t56 = t214 * t102 + t211 * t118;
t372 = t108 - t56;
t277 = qJD(4) * t211;
t109 = pkin(8) * (t211 * t279 - t214 * t275) + t191 * t214 - t199 * t277;
t55 = -t211 * t102 + t214 * t118;
t371 = t109 - t55;
t202 = qJD(1) * t210 + qJD(2);
t158 = t202 * t212 + t215 * t268;
t239 = -qJD(3) + t267;
t119 = t211 * t158 + t214 * t239;
t117 = Ifges(5,4) * t119;
t120 = t214 * t158 - t211 * t239;
t157 = t202 * t215 - t212 * t268;
t153 = qJD(4) - t157;
t300 = Ifges(6,5) * t119;
t354 = t120 * t361 + t360 * t153 - t117 + t300;
t369 = -qJ(5) * t373 - qJD(5) * t215 + t372;
t368 = pkin(4) * t373 - t371;
t114 = -t212 * t176 + t177 * t215;
t101 = -pkin(3) * t268 - t114;
t286 = t211 * t215;
t146 = -t214 * t268 + t267 * t286;
t285 = t214 * t216;
t147 = (t211 * t213 + t215 * t285) * t282;
t240 = pkin(4) * t211 - qJ(5) * t214;
t232 = pkin(8) + t240;
t241 = pkin(4) * t214 + qJ(5) * t211;
t278 = qJD(3) * t215;
t367 = -pkin(4) * t146 + qJ(5) * t147 - t101 + (qJD(4) * t241 - qJD(5) * t214) * t212 + t232 * t278;
t139 = -t202 * pkin(2) - t176;
t79 = -t157 * pkin(3) - t158 * pkin(9) + t139;
t283 = pkin(7) * t287 + t206;
t179 = t283 * qJD(1);
t140 = t202 * pkin(8) + t179;
t172 = (-pkin(2) * t216 - pkin(8) * t213 - pkin(1)) * t209;
t151 = qJD(1) * t172;
t94 = t215 * t140 + t212 * t151;
t81 = -pkin(9) * t239 + t94;
t23 = t211 * t79 + t214 * t81;
t20 = qJ(5) * t153 + t23;
t93 = -t212 * t140 + t215 * t151;
t80 = pkin(3) * t239 - t93;
t21 = t119 * pkin(4) - t120 * qJ(5) + t80;
t116 = Ifges(6,5) * t120;
t45 = Ifges(6,6) * t153 + Ifges(6,3) * t119 + t116;
t304 = Ifges(5,4) * t120;
t48 = -Ifges(5,2) * t119 + Ifges(5,6) * t153 + t304;
t366 = t20 * mrSges(6,2) + t23 * mrSges(5,3) - t45 / 0.2e1 + t48 / 0.2e1 - t21 * mrSges(6,1) - t80 * mrSges(5,1);
t280 = qJD(2) * t216;
t266 = t212 * t280;
t126 = t202 * t279 + (t213 * t278 + t266) * t282;
t265 = t215 * t280;
t125 = t202 * t278 + (-t213 * t279 + t265) * t282;
t281 = qJD(2) * t209;
t262 = qJD(1) * t281;
t260 = t213 * t262;
t64 = -qJD(4) * t119 + t214 * t125 + t211 * t260;
t65 = qJD(4) * t120 + t211 * t125 - t214 * t260;
t355 = t360 * t126 + t359 * t65 + t361 * t64;
t365 = -t211 * t370 + t214 * t360;
t299 = Ifges(6,5) * t211;
t303 = Ifges(5,4) * t211;
t364 = t214 * t361 + t299 - t303;
t22 = -t211 * t81 + t214 * t79;
t352 = qJD(5) - t22;
t19 = -pkin(4) * t153 + t352;
t271 = -Ifges(5,3) / 0.2e1 - Ifges(6,2) / 0.2e1;
t272 = Ifges(5,6) / 0.2e1 - Ifges(6,6) / 0.2e1;
t273 = Ifges(5,5) / 0.2e1 + Ifges(6,4) / 0.2e1;
t305 = Ifges(4,4) * t158;
t46 = Ifges(5,5) * t120 - Ifges(5,6) * t119 + Ifges(5,3) * t153;
t47 = Ifges(6,4) * t120 + Ifges(6,2) * t153 + Ifges(6,6) * t119;
t228 = Ifges(4,6) * t239;
t297 = Ifges(4,2) * t157;
t89 = -t228 + t297 + t305;
t219 = t272 * t119 - t273 * t120 + t271 * t153 + t19 * mrSges(6,1) + t23 * mrSges(5,2) + t94 * mrSges(4,3) - t46 / 0.2e1 - t47 / 0.2e1 + t89 / 0.2e1 + t305 / 0.2e1 - t139 * mrSges(4,1) - t20 * mrSges(6,3) - t22 * mrSges(5,1);
t363 = t297 / 0.2e1 + t219;
t362 = -Ifges(3,6) * t202 / 0.2e1;
t358 = Ifges(6,2) + Ifges(5,3);
t12 = Ifges(5,5) * t64 - Ifges(5,6) * t65 + Ifges(5,3) * t126;
t13 = Ifges(6,4) * t64 + Ifges(6,2) * t126 + Ifges(6,6) * t65;
t356 = t13 + t12;
t353 = -qJD(5) * t211 + t153 * t240 - t94;
t288 = t209 * t213;
t203 = pkin(7) * t288;
t186 = -t203 + t274;
t152 = Ifges(4,4) * t157;
t229 = Ifges(4,5) * t239;
t307 = Ifges(4,1) * t158;
t90 = t152 - t229 + t307;
t351 = t93 * mrSges(4,3) - t90 / 0.2e1 - t139 * mrSges(4,2) - t152 / 0.2e1;
t350 = t211 * t360 + t214 * t370;
t298 = Ifges(6,5) * t214;
t302 = Ifges(5,4) * t214;
t349 = t211 * t361 - t298 + t302;
t178 = qJD(2) * t231;
t167 = qJD(1) * t178;
t180 = t186 * qJD(2);
t168 = qJD(1) * t180;
t43 = -t140 * t279 + t151 * t278 + t212 * t167 + t215 * t168;
t33 = pkin(9) * t260 + t43;
t181 = t283 * qJD(2);
t169 = qJD(1) * t181;
t58 = t126 * pkin(3) - t125 * pkin(9) + t169;
t3 = t211 * t58 + t214 * t33 + t79 * t276 - t277 * t81;
t4 = -qJD(4) * t23 - t211 * t33 + t214 * t58;
t348 = -t211 * t4 + t214 * t3;
t1 = qJ(5) * t126 + qJD(5) * t153 + t3;
t2 = -pkin(4) * t126 - t4;
t347 = t1 * t214 + t2 * t211;
t170 = t203 + (-pkin(2) - t320) * t210;
t183 = -t210 * t215 + t212 * t288;
t184 = t210 * t212 + t215 * t288;
t97 = t183 * pkin(3) - t184 * pkin(9) + t170;
t171 = pkin(8) * t210 + t283;
t112 = t215 * t171 + t212 * t172;
t99 = -pkin(9) * t287 + t112;
t310 = t211 * t97 + t214 * t99;
t264 = t213 * t281;
t66 = -t171 * t279 + t172 * t278 + t212 * t178 + t215 * t180;
t53 = pkin(9) * t264 + t66;
t131 = qJD(3) * t184 + t209 * t266;
t132 = -qJD(3) * t183 + t209 * t265;
t74 = t131 * pkin(3) - t132 * pkin(9) + t181;
t9 = -qJD(4) * t310 - t211 * t53 + t214 * t74;
t237 = t211 * t23 + t214 * t22;
t238 = t19 * t214 - t20 * t211;
t243 = Ifges(6,3) * t211 + t298;
t250 = -Ifges(5,2) * t211 + t302;
t255 = mrSges(6,1) * t211 - mrSges(6,3) * t214;
t257 = mrSges(5,1) * t211 + mrSges(5,2) * t214;
t321 = t214 / 0.2e1;
t324 = t211 / 0.2e1;
t325 = -t211 / 0.2e1;
t330 = t153 / 0.2e1;
t337 = t120 / 0.2e1;
t339 = t119 / 0.2e1;
t340 = -t119 / 0.2e1;
t217 = t238 * mrSges(6,2) - t237 * mrSges(5,3) + t21 * t255 + t243 * t339 + t250 * t340 + t257 * t80 + t354 * t321 + t324 * t45 + t325 * t48 + t330 * t365 + t337 * t364;
t346 = t64 / 0.2e1;
t345 = -t65 / 0.2e1;
t344 = t65 / 0.2e1;
t342 = pkin(1) * mrSges(3,1);
t341 = pkin(1) * mrSges(3,2);
t338 = -t120 / 0.2e1;
t336 = t126 / 0.2e1;
t331 = -t153 / 0.2e1;
t329 = -t183 / 0.2e1;
t327 = t184 / 0.2e1;
t326 = t210 / 0.2e1;
t322 = -t214 / 0.2e1;
t315 = t43 * mrSges(4,2);
t44 = -t140 * t278 - t151 * t279 + t167 * t215 - t212 * t168;
t314 = t44 * mrSges(4,1);
t313 = qJD(3) / 0.2e1;
t82 = -mrSges(6,2) * t119 + mrSges(6,3) * t153;
t309 = mrSges(5,3) * t119;
t83 = -mrSges(5,2) * t153 - t309;
t312 = -t82 - t83;
t308 = mrSges(5,3) * t120;
t84 = mrSges(5,1) * t153 - t308;
t85 = -mrSges(6,1) * t153 + mrSges(6,2) * t120;
t311 = -t84 + t85;
t306 = Ifges(3,4) * t213;
t294 = t125 * Ifges(4,1);
t293 = t125 * Ifges(4,4);
t292 = t126 * Ifges(4,4);
t128 = -mrSges(4,1) * t239 - t158 * mrSges(4,3);
t73 = mrSges(5,1) * t119 + mrSges(5,2) * t120;
t290 = t128 - t73;
t110 = pkin(3) * t158 - pkin(9) * t157;
t39 = t211 * t110 + t214 * t93;
t289 = t199 * t214;
t284 = -mrSges(3,1) * t202 - mrSges(4,1) * t157 + mrSges(4,2) * t158 + mrSges(3,3) * t268;
t166 = t214 * t215 * pkin(8) + t211 * t199;
t270 = t211 * t287;
t269 = Ifges(4,5) * t125 - Ifges(4,6) * t126 + Ifges(4,3) * t260;
t263 = -t287 / 0.2e1;
t29 = -t126 * mrSges(6,1) + t64 * mrSges(6,2);
t111 = -t212 * t171 + t172 * t215;
t98 = pkin(3) * t287 - t111;
t258 = mrSges(5,1) * t214 - mrSges(5,2) * t211;
t256 = mrSges(6,1) * t214 + mrSges(6,3) * t211;
t249 = Ifges(5,2) * t214 + t303;
t246 = Ifges(4,5) * t132 - Ifges(4,6) * t131;
t242 = -Ifges(6,3) * t214 + t299;
t35 = -t211 * t99 + t214 * t97;
t38 = t110 * t214 - t211 * t93;
t233 = qJD(1) * t263 + qJD(3);
t67 = -t171 * t278 - t172 * t279 + t178 * t215 - t212 * t180;
t133 = t184 * t211 + t209 * t285;
t8 = t211 * t74 + t214 * t53 + t97 * t276 - t277 * t99;
t200 = Ifges(3,4) * t267;
t225 = Ifges(3,1) * t268 / 0.2e1 + t200 / 0.2e1 + t202 * Ifges(3,5) - t176 * mrSges(3,3);
t223 = -t4 * mrSges(5,1) + t2 * mrSges(6,1) + t3 * mrSges(5,2) - t1 * mrSges(6,3);
t54 = -pkin(3) * t264 - t67;
t34 = -pkin(3) * t260 - t44;
t222 = -t307 / 0.2e1 + t351;
t221 = t93 * mrSges(4,1) + t362 - (Ifges(3,2) * t216 + t306) * t282 / 0.2e1 + Ifges(4,6) * t157 + Ifges(4,5) * t158 - t179 * mrSges(3,3) - t94 * mrSges(4,2) + (-t239 / 0.2e1 + t313) * Ifges(4,3);
t198 = Ifges(3,5) * t216 * t262;
t193 = -pkin(3) - t241;
t175 = -t202 * mrSges(3,2) + mrSges(3,3) * t267;
t173 = t232 * t212;
t165 = -pkin(8) * t286 + t289;
t142 = -t289 + (pkin(8) * t211 + pkin(4)) * t215;
t141 = -qJ(5) * t215 + t166;
t134 = t184 * t214 - t270;
t127 = mrSges(4,2) * t239 + t157 * mrSges(4,3);
t105 = -mrSges(4,2) * t260 - mrSges(4,3) * t126;
t104 = mrSges(4,1) * t260 - mrSges(4,3) * t125;
t78 = -qJD(4) * t133 + t132 * t214 + t211 * t264;
t77 = -qJD(4) * t270 + t132 * t211 + t184 * t276 - t214 * t264;
t76 = mrSges(4,1) * t126 + mrSges(4,2) * t125;
t72 = mrSges(6,1) * t119 - mrSges(6,3) * t120;
t71 = pkin(4) * t120 + qJ(5) * t119;
t69 = Ifges(4,5) * t260 - t292 + t294;
t68 = -t126 * Ifges(4,2) + Ifges(4,6) * t260 + t293;
t37 = pkin(4) * t133 - qJ(5) * t134 + t98;
t31 = -mrSges(6,2) * t65 + mrSges(6,3) * t126;
t30 = -mrSges(5,2) * t126 - mrSges(5,3) * t65;
t28 = mrSges(5,1) * t126 - mrSges(5,3) * t64;
t27 = -pkin(4) * t183 - t35;
t26 = qJ(5) * t183 + t310;
t25 = -pkin(4) * t158 - t38;
t24 = qJ(5) * t158 + t39;
t18 = mrSges(5,1) * t65 + mrSges(5,2) * t64;
t17 = mrSges(6,1) * t65 - mrSges(6,3) * t64;
t14 = Ifges(5,4) * t64 - Ifges(5,2) * t65 + Ifges(5,6) * t126;
t11 = Ifges(6,5) * t64 + Ifges(6,6) * t126 + Ifges(6,3) * t65;
t10 = pkin(4) * t77 - qJ(5) * t78 - qJD(5) * t134 + t54;
t7 = -pkin(4) * t131 - t9;
t6 = pkin(4) * t65 - qJ(5) * t64 - qJD(5) * t120 + t34;
t5 = qJ(5) * t131 + qJD(5) * t183 + t8;
t15 = [t180 * t175 + t2 * (-mrSges(6,1) * t183 + mrSges(6,2) * t134) + t4 * (mrSges(5,1) * t183 - mrSges(5,3) * t134) + t170 * t76 + t139 * (mrSges(4,1) * t131 + mrSges(4,2) * t132) + t22 * (mrSges(5,1) * t131 - mrSges(5,3) * t78) + t19 * (-mrSges(6,1) * t131 + mrSges(6,2) * t78) - t131 * t89 / 0.2e1 + t132 * t90 / 0.2e1 + t66 * t127 + t67 * t128 + t158 * (Ifges(4,1) * t132 - Ifges(4,4) * t131) / 0.2e1 + (-t131 * t94 - t132 * t93 - t183 * t43 - t184 * t44) * mrSges(4,3) + (Ifges(5,4) * t78 + Ifges(5,6) * t131) * t340 + (Ifges(5,4) * t134 + Ifges(5,6) * t183) * t345 + t198 * t326 + t69 * t327 + t68 * t329 + t246 * t313 + t287 * t315 + t157 * (Ifges(4,4) * t132 - Ifges(4,2) * t131) / 0.2e1 - t287 * t314 - t126 * (Ifges(4,4) * t184 - Ifges(4,2) * t183 - Ifges(4,6) * t287) / 0.2e1 + t125 * (Ifges(4,1) * t184 - Ifges(4,4) * t183 - Ifges(4,5) * t287) / 0.2e1 + t168 * (-t210 * mrSges(3,2) + mrSges(3,3) * t287) + (-mrSges(3,1) * t210 + mrSges(4,1) * t183 + mrSges(4,2) * t184 + mrSges(3,3) * t288) * t169 + (t47 + t46) * t131 / 0.2e1 + m(5) * (t22 * t9 + t23 * t8 + t3 * t310 + t34 * t98 + t35 * t4 + t54 * t80) + t310 * t30 + (-t216 * t246 / 0.2e1 + ((-t186 * mrSges(3,3) + Ifges(3,5) * t326 + (-0.2e1 * t341 + 0.3e1 / 0.2e1 * Ifges(3,4) * t216) * t209) * t216 + (-t283 * mrSges(3,3) - Ifges(3,6) * t210 + Ifges(4,5) * t327 + Ifges(4,6) * t329 + (-0.2e1 * t342 - 0.3e1 / 0.2e1 * t306 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2) - Ifges(4,3)) * t216) * t209) * t213) * qJD(2)) * t282 + m(3) * (t168 * t283 - t169 * t186 - t176 * t181 + t179 * t180) + t284 * t181 + (Ifges(6,5) * t78 + Ifges(6,6) * t131) * t339 + (Ifges(6,5) * t134 + Ifges(6,6) * t183) * t344 + m(6) * (t1 * t26 + t10 * t21 + t19 * t7 + t2 * t27 + t20 * t5 + t37 * t6) + m(4) * (t111 * t44 + t112 * t43 + t139 * t181 + t169 * t170 + t66 * t94 + t67 * t93) + t111 * t104 + t112 * t105 + t98 * t18 + t5 * t82 + t8 * t83 + t9 * t84 + t7 * t85 + t10 * t72 + t54 * t73 + (t1 * t183 + t131 * t20 - t134 * t6 - t21 * t78) * mrSges(6,3) + t35 * t28 + t37 * t17 + t27 * t29 + t26 * t31 + (-Ifges(5,2) * t340 + Ifges(6,3) * t339 - t330 * t370 + t359 * t337 - t366) * t77 + (-t3 * mrSges(5,3) - t1 * mrSges(6,2) + t6 * mrSges(6,1) + t34 * mrSges(5,1) + t11 / 0.2e1 - t14 / 0.2e1 + Ifges(6,3) * t344 - Ifges(5,2) * t345 + t359 * t346 - t370 * t336) * t133 + t269 * t263 + (-t131 * t23 + t134 * t34 - t183 * t3 + t78 * t80) * mrSges(5,2) + t354 * t78 / 0.2e1 + t355 * t134 / 0.2e1 + t356 * t183 / 0.2e1 + (t131 * t358 + t360 * t78) * t330 + (t134 * t360 + t183 * t358) * t336 + (t131 * t360 + t361 * t78) * t337 + (t134 * t361 + t183 * t360) * t346 + (t225 * t216 + (t362 + t221) * t213) * t281; t165 * t28 + t166 * t30 - t168 * mrSges(3,2) - t169 * mrSges(3,1) + t173 * t17 - t176 * t175 + t142 * t29 + t141 * t31 - t115 * t127 - t114 * t128 + (-m(4) * t169 - t76) * pkin(2) + t367 * t72 + t368 * t85 + m(5) * (t108 * t23 + t109 * t22 + t165 * t4 + t166 * t3) + (t68 / 0.2e1 - t13 / 0.2e1 - t12 / 0.2e1 - t169 * mrSges(4,1) + t293 / 0.2e1 + t43 * mrSges(4,3) + t272 * t65 - t273 * t64 + (m(4) * t43 + t105) * pkin(8) + (-Ifges(4,2) / 0.2e1 + t271) * t126 + t223) * t215 + (-t354 / 0.2e1 - mrSges(5,2) * t80 - t19 * mrSges(6,2) + t22 * mrSges(5,3) + mrSges(6,3) * t21 + Ifges(5,4) * t339 + Ifges(6,5) * t340 + t360 * t331 + t361 * t338) * t147 - t284 * t179 - m(5) * (t101 * t80 + t22 * t55 + t23 * t56) - m(4) * (t114 * t93 + t115 * t94 + t139 * t179) + t371 * t84 + t372 * t83 + (((-m(4) * t94 - t127) * pkin(8) - qJD(3) * Ifges(4,6) / 0.2e1 - t363) * t212 + ((-m(4) * t93 + m(5) * t80 - t290) * pkin(8) + Ifges(4,5) * t313 + t217 - t222) * t215) * qJD(3) + ((qJD(2) * (Ifges(4,5) * t212 + Ifges(4,6) * t215) / 0.2e1 + (t342 + t306 / 0.2e1) * t282 + (t202 / 0.2e1 - qJD(2)) * Ifges(3,6) - t221) * t213 + ((t341 - Ifges(3,1) * t213 / 0.2e1 + (Ifges(3,2) + Ifges(4,3)) * t213 / 0.2e1) * t282 - t200 / 0.2e1 + (-Ifges(4,5) * t233 + t222) * t215 + (t233 * Ifges(4,6) + t363) * t212 - t225) * t216) * t282 + (t1 * t141 + t142 * t2 + t173 * t6 + t19 * t368 + t20 * t369 + t21 * t367) * m(6) + t369 * t82 + (-Ifges(5,2) * t339 + Ifges(6,3) * t340 - t331 * t370 + t338 * t359 + t366) * t146 - t101 * t73 + (t169 * mrSges(4,2) - t292 / 0.2e1 + t294 / 0.2e1 - t44 * mrSges(4,3) + t243 * t344 + t250 * t345 + t11 * t324 + t14 * t325 + t34 * t257 + t6 * t255 + t69 / 0.2e1 + (-t211 * t3 - t214 * t4) * mrSges(5,3) + (-t1 * t211 + t2 * t214) * mrSges(6,2) + (-m(4) * t44 + m(5) * t34 - t104 + t18) * pkin(8) + (t21 * t256 + t80 * t258 + t242 * t340 + t249 * t339 + t48 * t322 + (t211 * t22 - t214 * t23) * mrSges(5,3) + (-t19 * t211 - t20 * t214) * mrSges(6,2) + t349 * t338 + t350 * t331 + t354 * t325) * qJD(4) + t364 * t346 + t365 * t336 + (qJD(4) * t45 + t355) * t321) * t212 + t198; -t6 * t256 - t34 * t258 + t193 * t17 - t93 * t127 + t314 - t315 + t14 * t321 + t11 * t322 + t290 * t94 + (t219 - t228 / 0.2e1) * t158 + t269 - t24 * t82 - t39 * t83 - t38 * t84 - t25 * t85 - pkin(3) * t18 + t217 * qJD(4) + t242 * t344 + t249 * t345 + (-pkin(3) * t34 - t22 * t38 - t23 * t39 - t80 * t94) * m(5) + t347 * mrSges(6,2) + t348 * mrSges(5,3) + ((-m(5) * t237 + m(6) * t238 + t211 * t312 + t214 * t311) * qJD(4) + m(6) * t347 + m(5) * t348 + (t30 + t31) * t214 + (-t28 + t29) * t211) * pkin(9) + t349 * t346 + t350 * t336 + ((-Ifges(4,1) / 0.2e1 + Ifges(4,2) / 0.2e1) * t158 + t229 / 0.2e1 - t217 + t351) * t157 + (-t19 * t25 + t193 * t6 - t20 * t24 + t21 * t353) * m(6) + t353 * t72 + t355 * t324; -t80 * (mrSges(5,1) * t120 - mrSges(5,2) * t119) + t48 * t337 + (Ifges(6,3) * t120 - t300) * t340 - t21 * (mrSges(6,1) * t120 + mrSges(6,3) * t119) - t223 + (t308 - t311) * t23 + (-t309 + t312) * t22 + qJD(5) * t82 - t71 * t72 - pkin(4) * t29 + qJ(5) * t31 + (t119 * t19 + t120 * t20) * mrSges(6,2) + (-t119 * t360 - t120 * t370) * t331 + (-pkin(4) * t2 + qJ(5) * t1 - t19 * t23 + t20 * t352 - t21 * t71) * m(6) + (-Ifges(5,2) * t120 - t117 + t354) * t339 + (-t119 * t361 + t116 - t304 + t45) * t338 + t356; t120 * t72 - t153 * t82 + 0.2e1 * (t2 / 0.2e1 + t21 * t337 + t20 * t331) * m(6) + t29;];
tauc = t15(:);
