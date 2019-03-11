% Calculate time derivative of joint inertia matrix for
% S6RRRPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha4,d1,d2,d3,theta4]';
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
% Datum: 2019-03-09 15:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPPP1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPP1_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPP1_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPP1_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPP1_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPP1_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPP1_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:14:30
% EndTime: 2019-03-09 15:14:43
% DurationCPUTime: 5.70s
% Computational Cost: add. (5530->679), mult. (16262->950), div. (0->0), fcn. (14080->8), ass. (0->285)
t346 = Ifges(5,5) + Ifges(7,5);
t345 = Ifges(6,5) + Ifges(7,4);
t236 = sin(pkin(6));
t344 = -t236 / 0.2e1;
t343 = t236 / 0.2e1;
t238 = cos(pkin(6));
t342 = -t238 / 0.2e1;
t341 = t238 / 0.2e1;
t239 = sin(qJ(3));
t241 = cos(qJ(3));
t255 = qJ(4) * t236 * t239 + pkin(3) * t241;
t195 = -pkin(2) - t255;
t322 = qJ(4) * t238;
t277 = pkin(9) + t322;
t206 = t277 * t239;
t237 = cos(pkin(10));
t340 = (t195 * t236 - t206 * t238) * t237;
t240 = sin(qJ(2));
t242 = cos(qJ(2));
t219 = -pkin(2) * t242 - t240 * pkin(9) - pkin(1);
t306 = t241 * t242;
t232 = pkin(8) * t306;
t165 = t239 * t219 + t232;
t310 = t239 * t240;
t253 = t236 * t242 + t238 * t310;
t298 = qJD(3) * t239;
t317 = t236 * t241;
t256 = pkin(3) * t239 - qJ(4) * t317;
t294 = qJD(4) * t239;
t150 = t256 * qJD(3) - t236 * t294;
t259 = qJD(3) * t277;
t293 = qJD(4) * t241;
t280 = t238 * t293;
t151 = -t239 * t259 + t280;
t281 = t238 * t294;
t152 = -t241 * t259 - t281;
t235 = sin(pkin(10));
t319 = t235 * t238;
t320 = t235 * t236;
t49 = t150 * t320 + t237 * t151 + t152 * t319;
t34 = -(qJ(5) * t298 - qJD(5) * t241) * t236 - t49;
t295 = qJD(4) * t236;
t296 = qJD(3) * t241;
t300 = qJD(2) * t241;
t213 = (pkin(2) * t240 - pkin(9) * t242) * qJD(2);
t304 = t239 * t213 + t219 * t296;
t37 = (-t295 + (-pkin(8) * qJD(3) - qJD(2) * t322) * t239) * t242 + (-pkin(8) * t300 - t281 + (qJD(2) * t236 - t238 * t296) * qJ(4)) * t240 + t304;
t312 = t238 * t241;
t290 = qJ(4) * t312;
t301 = qJD(2) * t240;
t330 = pkin(8) * t239;
t303 = t241 * t213 + t301 * t330;
t313 = t238 * t240;
t50 = -t240 * t280 + (pkin(3) * t240 - t242 * t290) * qJD(2) + (-t232 + (qJ(4) * t313 - t219) * t239) * qJD(3) + t303;
t249 = pkin(8) + t256;
t299 = qJD(2) * t242;
t76 = (t255 * qJD(3) - t236 * t293) * t240 + t249 * t299;
t7 = -t235 * t37 + (t236 * t76 + t238 * t50) * t237;
t114 = -t253 * qJ(4) + t165;
t199 = t241 * t219;
t118 = -t240 * t290 + t199 + (-pkin(3) - t330) * t242;
t153 = t249 * t240;
t32 = -t235 * t114 + (t118 * t238 + t153 * t236) * t237;
t339 = 2 * m(4);
t338 = 2 * m(5);
t337 = 2 * m(6);
t336 = 2 * m(7);
t335 = -0.2e1 * pkin(1);
t334 = 0.2e1 * pkin(8);
t22 = -t236 * t50 + t238 * t76;
t333 = m(5) * t22;
t79 = t238 * t150 - t152 * t236;
t332 = m(5) * t79;
t331 = -t240 / 0.2e1;
t329 = -mrSges(6,2) + mrSges(5,1);
t328 = pkin(4) + qJ(6);
t327 = m(5) * qJD(4);
t326 = Ifges(4,4) * t239;
t325 = Ifges(4,4) * t241;
t324 = Ifges(4,6) * t239;
t323 = t242 * Ifges(4,6);
t321 = t150 * t237;
t318 = t236 * t237;
t315 = t237 * t238;
t314 = t237 * t239;
t220 = Ifges(4,2) * t241 + t326;
t311 = t239 * t220;
t309 = t239 * t242;
t308 = t240 * t241;
t221 = Ifges(4,1) * t239 + t325;
t307 = t241 * t221;
t207 = t277 * t241;
t187 = t235 * t207;
t305 = pkin(4) * t317 + t187;
t286 = t241 * t299;
t302 = -Ifges(4,5) * t286 - Ifges(4,3) * t301;
t186 = pkin(3) * t319 + qJ(4) * t318;
t297 = qJD(3) * t240;
t292 = qJD(5) * t235;
t291 = 0.2e1 * t235;
t8 = t237 * t37 + t50 * t319 + t76 * t320;
t33 = t237 * t114 + t118 * t319 + t153 * t320;
t78 = t195 * t320 - t206 * t319 + t237 * t207;
t288 = -pkin(3) * t237 - pkin(4);
t287 = t239 * t299;
t285 = t236 * t298;
t284 = t239 * t297;
t283 = t240 * t296;
t282 = t235 * t295;
t254 = t235 * t241 + t238 * t314;
t276 = t238 * t283;
t74 = -t235 * t284 + t237 * t276 + t254 * t299 - t301 * t318;
t75 = -t235 * t276 - t287 * t319 - t237 * t284 + (t237 * t306 + t240 * t320) * qJD(2);
t28 = -t75 * mrSges(7,2) + t74 * mrSges(7,3);
t134 = t236 * t283 + (t236 * t309 + t313) * qJD(2);
t45 = t75 * mrSges(6,1) + t134 * mrSges(6,2);
t44 = -t74 * mrSges(7,1) + t134 * mrSges(7,2);
t42 = t75 * mrSges(7,1) - t134 * mrSges(7,3);
t279 = (2 * Ifges(3,4)) + t324;
t278 = -qJ(5) * t235 - pkin(3);
t176 = t237 * t296 - t298 * t319;
t129 = t176 * mrSges(6,1) + mrSges(6,2) * t285;
t175 = t254 * qJD(3);
t128 = -t175 * mrSges(7,1) + mrSges(7,2) * t285;
t106 = -t176 * mrSges(7,2) + t175 * mrSges(7,3);
t65 = -t118 * t236 + t238 * t153;
t120 = t238 * t195 + t206 * t236;
t13 = Ifges(7,5) * t134 + Ifges(7,6) * t74 + Ifges(7,3) * t75;
t17 = Ifges(6,4) * t134 - Ifges(6,2) * t75 + Ifges(6,6) * t74;
t21 = Ifges(5,1) * t75 - Ifges(5,4) * t74 + Ifges(5,5) * t134;
t275 = t13 / 0.2e1 + t21 / 0.2e1 - t17 / 0.2e1;
t14 = Ifges(6,5) * t134 - Ifges(6,6) * t75 + Ifges(6,3) * t74;
t16 = Ifges(7,4) * t134 + Ifges(7,2) * t74 + Ifges(7,6) * t75;
t18 = Ifges(5,4) * t75 - Ifges(5,2) * t74 + Ifges(5,6) * t134;
t274 = t14 / 0.2e1 + t16 / 0.2e1 - t18 / 0.2e1;
t15 = Ifges(5,5) * t75 - Ifges(5,6) * t74 + Ifges(5,3) * t134;
t19 = Ifges(7,1) * t134 + Ifges(7,4) * t74 + Ifges(7,5) * t75;
t20 = Ifges(6,1) * t134 - Ifges(6,4) * t75 + Ifges(6,5) * t74;
t273 = -t19 / 0.2e1 - t20 / 0.2e1 - t15 / 0.2e1;
t80 = Ifges(7,5) * t285 + Ifges(7,6) * t175 + Ifges(7,3) * t176;
t83 = Ifges(6,4) * t285 - Ifges(6,2) * t176 + Ifges(6,6) * t175;
t88 = Ifges(5,1) * t176 - Ifges(5,4) * t175 + Ifges(5,5) * t285;
t272 = t80 / 0.2e1 + t88 / 0.2e1 - t83 / 0.2e1;
t81 = Ifges(6,5) * t285 - Ifges(6,6) * t176 + Ifges(6,3) * t175;
t82 = Ifges(7,4) * t285 + Ifges(7,2) * t175 + Ifges(7,6) * t176;
t87 = Ifges(5,4) * t176 - Ifges(5,2) * t175 + Ifges(5,6) * t285;
t271 = t81 / 0.2e1 + t82 / 0.2e1 - t87 / 0.2e1;
t84 = Ifges(7,1) * t285 + Ifges(7,4) * t175 + Ifges(7,5) * t176;
t85 = Ifges(6,1) * t285 - Ifges(6,4) * t176 + Ifges(6,5) * t175;
t86 = Ifges(5,5) * t176 - Ifges(5,6) * t175 + Ifges(5,3) * t285;
t270 = t84 / 0.2e1 + t85 / 0.2e1 + t86 / 0.2e1;
t138 = t235 * t151;
t269 = -t152 * t315 + t138;
t154 = -qJ(5) * t238 - t186;
t268 = ((Ifges(5,6) - t345) * t237 + (-Ifges(6,4) + t346) * t235) * t343 + (Ifges(5,3) + Ifges(7,1) + Ifges(6,1)) * t341;
t267 = Ifges(6,4) * t342 + (-Ifges(6,2) * t235 - Ifges(6,6) * t237) * t344 + ((Ifges(5,4) - Ifges(7,6)) * t237 + (Ifges(5,1) + Ifges(7,3)) * t235) * t343 + t346 * t341;
t266 = Ifges(5,6) * t342 + (Ifges(5,4) * t235 + Ifges(5,2) * t237) * t344 + ((-Ifges(7,2) - Ifges(6,3)) * t237 + (-Ifges(6,6) + Ifges(7,6)) * t235) * t343 + t345 * t341;
t265 = -mrSges(4,1) * t241 + mrSges(4,2) * t239;
t264 = mrSges(4,1) * t239 + mrSges(4,2) * t241;
t263 = Ifges(4,1) * t241 - t326;
t262 = -Ifges(4,2) * t239 + t325;
t261 = Ifges(4,5) * t239 + Ifges(4,6) * t241;
t126 = t176 * mrSges(7,1) - mrSges(7,3) * t285;
t189 = t236 * t310 - t238 * t242;
t24 = -qJ(5) * t189 - t33;
t123 = -t253 * t235 + t237 * t308;
t252 = -qJ(5) * t123 + t65;
t181 = t235 * t312 + t314;
t251 = -qJ(5) * t181 + t120;
t66 = qJ(5) * t317 - t78;
t248 = -t284 + t286;
t247 = t283 + t287;
t4 = -qJ(5) * t134 - qJD(5) * t189 - t8;
t244 = -qJ(5) * t75 - qJD(5) * t123 + t22;
t243 = -qJ(5) * t176 - qJD(5) * t181 + t79;
t234 = Ifges(4,5) * t296;
t228 = qJ(4) * t320;
t212 = -mrSges(4,1) * t242 - mrSges(4,3) * t308;
t211 = mrSges(4,2) * t242 - mrSges(4,3) * t310;
t210 = t263 * qJD(3);
t209 = t262 * qJD(3);
t208 = t264 * qJD(3);
t205 = mrSges(6,1) * t320 + mrSges(6,2) * t238;
t204 = mrSges(7,1) * t318 + mrSges(7,2) * t238;
t203 = -mrSges(6,1) * t318 - mrSges(6,3) * t238;
t202 = mrSges(7,1) * t320 - mrSges(7,3) * t238;
t201 = -mrSges(5,2) * t238 + mrSges(5,3) * t318;
t200 = mrSges(5,1) * t238 - mrSges(5,3) * t320;
t194 = qJD(5) * t238 + t237 * t295;
t193 = -qJD(6) * t238 + t282;
t185 = (-mrSges(7,2) * t235 - mrSges(7,3) * t237) * t236;
t184 = (-mrSges(5,1) * t237 + mrSges(5,2) * t235) * t236;
t183 = (mrSges(6,2) * t237 - mrSges(6,3) * t235) * t236;
t182 = pkin(3) * t315 - t228;
t180 = t235 * t239 - t237 * t312;
t179 = (-qJD(6) * t237 - t292) * t236;
t174 = -Ifges(4,5) * t242 + t263 * t240;
t173 = t262 * t240 - t323;
t171 = (-pkin(4) * t237 + t278) * t236;
t170 = t288 * t238 + t228;
t169 = t176 * mrSges(6,3);
t167 = t176 * mrSges(5,2);
t164 = -pkin(8) * t309 + t199;
t149 = -mrSges(4,2) * t301 - t247 * mrSges(4,3);
t148 = mrSges(4,1) * t301 - t248 * mrSges(4,3);
t146 = -mrSges(5,1) * t317 - mrSges(5,3) * t181;
t145 = mrSges(5,2) * t317 - mrSges(5,3) * t180;
t144 = mrSges(6,1) * t181 - mrSges(6,2) * t317;
t143 = -mrSges(7,1) * t180 - mrSges(7,2) * t317;
t142 = mrSges(6,1) * t180 + mrSges(6,3) * t317;
t141 = mrSges(7,1) * t181 + mrSges(7,3) * t317;
t136 = (-t328 * t237 + t278) * t236;
t135 = pkin(5) * t318 - t154;
t131 = mrSges(5,1) * t285 - mrSges(5,3) * t176;
t130 = -mrSges(5,2) * t285 - mrSges(5,3) * t175;
t127 = mrSges(6,1) * t175 - mrSges(6,3) * t285;
t122 = t235 * t308 + t253 * t237;
t121 = pkin(5) * t320 + t228 + (-qJ(6) + t288) * t238;
t119 = t247 * mrSges(4,1) + t248 * mrSges(4,2);
t117 = -mrSges(6,2) * t180 - mrSges(6,3) * t181;
t116 = mrSges(5,1) * t180 + mrSges(5,2) * t181;
t115 = -mrSges(7,2) * t181 + mrSges(7,3) * t180;
t113 = -t221 * t297 + (Ifges(4,5) * t240 + t263 * t242) * qJD(2);
t112 = -t220 * t297 + (Ifges(4,6) * t240 + t262 * t242) * qJD(2);
t108 = -t175 * mrSges(6,2) - t169;
t107 = t175 * mrSges(5,1) + t167;
t105 = Ifges(5,1) * t181 - Ifges(5,4) * t180 - Ifges(5,5) * t317;
t104 = Ifges(5,4) * t181 - Ifges(5,2) * t180 - Ifges(5,6) * t317;
t103 = Ifges(5,5) * t181 - Ifges(5,6) * t180 - Ifges(5,3) * t317;
t102 = -Ifges(6,1) * t317 - Ifges(6,4) * t181 + Ifges(6,5) * t180;
t101 = -Ifges(7,1) * t317 + Ifges(7,4) * t180 + Ifges(7,5) * t181;
t100 = -Ifges(6,4) * t317 - Ifges(6,2) * t181 + Ifges(6,6) * t180;
t99 = -Ifges(7,4) * t317 + Ifges(7,2) * t180 + Ifges(7,6) * t181;
t98 = -Ifges(6,5) * t317 - Ifges(6,6) * t181 + Ifges(6,3) * t180;
t97 = -Ifges(7,5) * t317 + Ifges(7,6) * t180 + Ifges(7,3) * t181;
t96 = -t165 * qJD(3) + t303;
t95 = (-t240 * t300 - t242 * t298) * pkin(8) + t304;
t94 = mrSges(6,1) * t123 + mrSges(6,2) * t189;
t93 = -mrSges(7,1) * t122 + mrSges(7,2) * t189;
t92 = mrSges(6,1) * t122 - mrSges(6,3) * t189;
t91 = mrSges(7,1) * t123 - mrSges(7,3) * t189;
t90 = mrSges(5,1) * t189 - mrSges(5,3) * t123;
t89 = -mrSges(5,2) * t189 - mrSges(5,3) * t122;
t77 = -t187 + t340;
t72 = t75 * mrSges(6,3);
t70 = t75 * mrSges(5,2);
t67 = t305 - t340;
t64 = -mrSges(6,2) * t122 - mrSges(6,3) * t123;
t63 = mrSges(5,1) * t122 + mrSges(5,2) * t123;
t62 = -mrSges(7,2) * t123 + mrSges(7,3) * t122;
t61 = pkin(4) * t180 + t251;
t60 = -pkin(5) * t180 - t66;
t59 = Ifges(5,1) * t123 - Ifges(5,4) * t122 + Ifges(5,5) * t189;
t58 = Ifges(6,1) * t189 - Ifges(6,4) * t123 + Ifges(6,5) * t122;
t57 = Ifges(7,1) * t189 + Ifges(7,4) * t122 + Ifges(7,5) * t123;
t56 = Ifges(5,4) * t123 - Ifges(5,2) * t122 + Ifges(5,6) * t189;
t55 = Ifges(6,4) * t189 - Ifges(6,2) * t123 + Ifges(6,6) * t122;
t54 = Ifges(7,4) * t189 + Ifges(7,2) * t122 + Ifges(7,6) * t123;
t53 = Ifges(5,5) * t123 - Ifges(5,6) * t122 + Ifges(5,3) * t189;
t52 = Ifges(6,5) * t189 - Ifges(6,6) * t123 + Ifges(6,3) * t122;
t51 = Ifges(7,5) * t189 + Ifges(7,6) * t122 + Ifges(7,3) * t123;
t48 = -t138 + (t150 * t236 + t152 * t238) * t237;
t47 = t328 * t180 + t251;
t46 = t206 * t315 + pkin(5) * t181 + (qJ(6) * t241 - t195 * t237) * t236 + t305;
t43 = mrSges(6,1) * t74 - mrSges(6,3) * t134;
t41 = mrSges(5,1) * t134 - mrSges(5,3) * t75;
t40 = -mrSges(5,2) * t134 - mrSges(5,3) * t74;
t38 = (-pkin(4) * t298 - t321) * t236 + t269;
t31 = pkin(4) * t175 + t243;
t30 = -t74 * mrSges(6,2) - t72;
t29 = t74 * mrSges(5,1) + t70;
t27 = -pkin(5) * t175 - t34;
t26 = pkin(4) * t122 + t252;
t25 = -pkin(4) * t189 - t32;
t23 = pkin(5) * t176 + (qJD(6) * t241 - t328 * t298 - t321) * t236 + t269;
t12 = qJD(6) * t180 + t328 * t175 + t243;
t11 = t328 * t122 + t252;
t10 = -pkin(5) * t122 - t24;
t9 = pkin(5) * t123 - t328 * t189 - t32;
t6 = -pkin(4) * t134 - t7;
t5 = pkin(4) * t74 + t244;
t3 = -pkin(5) * t74 - t4;
t2 = pkin(5) * t75 - qJD(6) * t189 - t328 * t134 - t7;
t1 = qJD(6) * t122 + t328 * t74 + t244;
t35 = [(t15 + t19 + t20) * t189 + (t53 + t57 + t58) * t134 + (t13 + t21 - t17) * t123 + (t51 + t59 - t55) * t75 + (t52 + t54 - t56) * t74 + (t14 + t16 - t18) * t122 + (t1 * t11 + t10 * t3 + t2 * t9) * t336 + (t24 * t4 + t25 * t6 + t26 * t5) * t337 + (t22 * t65 + t32 * t7 + t33 * t8) * t338 + (t164 * t96 + t165 * t95) * t339 + (t119 * t334 - t239 * t112 + t241 * t113 + (-t241 * t173 - t239 * t174 + t242 * t261) * qJD(3) + (mrSges(3,1) * t335 + (Ifges(4,5) * t241 - t279) * t240 + (pkin(8) ^ 2 * t339 + t264 * t334 + (2 * Ifges(3,1)) - (2 * Ifges(3,2)) - Ifges(4,3)) * t242) * qJD(2)) * t240 + ((mrSges(3,2) * t335 - t239 * t173 + t241 * t174 + t279 * t242) * qJD(2) + t302) * t242 + 0.2e1 * t11 * t28 + 0.2e1 * t26 * t30 + 0.2e1 * t33 * t40 + 0.2e1 * t32 * t41 + 0.2e1 * t9 * t42 + 0.2e1 * t24 * t43 + 0.2e1 * t10 * t44 + 0.2e1 * t25 * t45 + 0.2e1 * t1 * t62 + 0.2e1 * t22 * t63 + 0.2e1 * t5 * t64 + 0.2e1 * t65 * t29 + 0.2e1 * t8 * t89 + 0.2e1 * t7 * t90 + 0.2e1 * t2 * t91 + 0.2e1 * t4 * t92 + 0.2e1 * t3 * t93 + 0.2e1 * t6 * t94 + 0.2e1 * t164 * t148 + 0.2e1 * t165 * t149 + 0.2e1 * t95 * t211 + 0.2e1 * t96 * t212; m(4) * (-pkin(2) * pkin(8) * t299 + (-t239 * t96 + t241 * t95) * pkin(9)) + (t97 / 0.2e1 - t100 / 0.2e1 + t105 / 0.2e1) * t75 + (-pkin(9) * t148 - t96 * mrSges(4,3) + t113 / 0.2e1) * t239 + (t98 / 0.2e1 + t99 / 0.2e1 - t104 / 0.2e1) * t74 + (t241 * t210 / 0.2e1 - t239 * t209 / 0.2e1 + pkin(8) * t208 + (t261 / 0.2e1 - Ifges(3,6) + pkin(8) * mrSges(3,2)) * qJD(2)) * t240 + (t51 / 0.2e1 + t59 / 0.2e1 - t55 / 0.2e1) * t176 + (t52 / 0.2e1 + t54 / 0.2e1 - t56 / 0.2e1) * t175 + ((t220 * t331 - t164 * mrSges(4,3) + t174 / 0.2e1 + (-m(4) * t164 - t212) * pkin(9)) * t241 + (t221 * t331 - t165 * mrSges(4,3) - t173 / 0.2e1 + t323 / 0.2e1 + (-m(4) * t165 - t211) * pkin(9) + (t53 / 0.2e1 + t57 / 0.2e1 + t58 / 0.2e1) * t236) * t239) * qJD(3) + (-t234 / 0.2e1 + (t307 / 0.2e1 - t311 / 0.2e1 + Ifges(3,5) + (-mrSges(3,1) + t265) * pkin(8)) * qJD(2)) * t242 + (t101 / 0.2e1 + t102 / 0.2e1 + t103 / 0.2e1) * t134 + m(5) * (t120 * t22 + t32 * t48 + t33 * t49 + t65 * t79 + t7 * t77 + t78 * t8) + m(7) * (t1 * t47 + t10 * t27 + t11 * t12 + t2 * t46 + t23 * t9 + t3 * t60) + m(6) * (t24 * t34 + t25 * t38 + t26 * t31 + t4 * t66 + t5 * t61 + t6 * t67) + t274 * t180 + t275 * t181 + t270 * t189 + t271 * t122 + t272 * t123 + (pkin(9) * t149 + t95 * mrSges(4,3) + t112 / 0.2e1 + t273 * t236) * t241 + t46 * t42 + t47 * t28 + t60 * t44 + t61 * t30 + t12 * t62 + t31 * t64 + t66 * t43 + t67 * t45 + t77 * t41 + t78 * t40 + t79 * t63 + t49 * t89 + t48 * t90 + t23 * t91 + t34 * t92 + t27 * t93 + t38 * t94 + t11 * t106 + t65 * t107 + t26 * t108 + t1 * t115 + t22 * t116 + t5 * t117 - pkin(2) * t119 + t120 * t29 + t9 * t126 + t24 * t127 + t10 * t128 + t25 * t129 + t33 * t130 + t32 * t131 + t2 * t141 + t4 * t142 + t3 * t143 + t6 * t144 + t8 * t145 + t7 * t146; (t12 * t47 + t23 * t46 + t27 * t60) * t336 + (t31 * t61 + t34 * t66 + t38 * t67) * t337 + (t120 * t79 + t48 * t77 + t49 * t78) * t338 + (t307 - t311) * qJD(3) + ((-t84 - t85 - t86) * t241 + (t101 + t102 + t103) * t298) * t236 + (t80 + t88 - t83) * t181 + (t81 + t82 - t87) * t180 + (t97 + t105 - t100) * t176 + (t98 + t99 - t104) * t175 + 0.2e1 * t47 * t106 + 0.2e1 * t61 * t108 + 0.2e1 * t12 * t115 + 0.2e1 * t79 * t116 + 0.2e1 * t31 * t117 + 0.2e1 * t120 * t107 + 0.2e1 * t46 * t126 + 0.2e1 * t66 * t127 + 0.2e1 * t60 * t128 + 0.2e1 * t67 * t129 + 0.2e1 * t78 * t130 + 0.2e1 * t77 * t131 + 0.2e1 * t23 * t141 + 0.2e1 * t34 * t142 + 0.2e1 * t27 * t143 + 0.2e1 * t38 * t144 + 0.2e1 * t49 * t145 + 0.2e1 * t48 * t146 - 0.2e1 * pkin(2) * t208 + t239 * t210 + t241 * t209; (t93 - t92) * t194 + m(5) * (t182 * t7 + t186 * t8) + m(6) * (t154 * t4 + t170 * t6 + t171 * t5 - t194 * t24) + m(7) * (t1 * t136 + t10 * t194 + t11 * t179 + t121 * t2 + t135 * t3 + t193 * t9) + ((-t29 - t333) * pkin(3) + ((m(5) * t33 + t89) * qJD(4) - t274) * t237 + (-qJD(5) * t64 + (-t90 + t94) * qJD(4) - t32 * t327 + m(6) * (qJD(4) * t25 - qJD(5) * t26) + t275) * t235) * t236 - Ifges(4,5) * t284 - t273 * t238 + t266 * t74 + t267 * t75 + t268 * t134 - t302 - t95 * mrSges(4,2) + t96 * mrSges(4,1) + t121 * t42 + t135 * t44 + t136 * t28 + t154 * t43 + t170 * t45 + t171 * t30 + t179 * t62 + t182 * t41 + t5 * t183 + t22 * t184 + t1 * t185 + t186 * t40 + t193 * t91 + t7 * t200 + t8 * t201 + t2 * t202 + t4 * t203 + t3 * t204 + t6 * t205 - t247 * Ifges(4,6); (t143 - t142) * t194 + m(6) * (t154 * t34 + t170 * t38 + t171 * t31 - t194 * t66) + m(5) * (t182 * t48 + t186 * t49) + m(7) * (t12 * t136 + t121 * t23 + t135 * t27 + t179 * t47 + t193 * t46 + t194 * t60) + t234 + ((-t107 - t332) * pkin(3) + ((m(5) * t78 + t145) * qJD(4) - t271) * t237 + t268 * t298 + (-qJD(5) * t117 + (t144 - t146) * qJD(4) - t77 * t327 + m(6) * (qJD(4) * t67 - qJD(5) * t61) + t272) * t235) * t236 + (t265 * pkin(9) - t324) * qJD(3) + t270 * t238 + t266 * t175 + t267 * t176 + t121 * t126 + t135 * t128 + t136 * t106 + t154 * t127 + t170 * t129 + t171 * t108 + t179 * t115 + t182 * t131 + t31 * t183 + t79 * t184 + t12 * t185 + t186 * t130 + t193 * t141 + t48 * t200 + t49 * t201 + t23 * t202 + t34 * t203 + t27 * t204 + t38 * t205; 0.2e1 * t179 * t185 + 0.2e1 * t193 * t202 + (t121 * t193 + t136 * t179) * t336 + ((-m(6) * t171 - t183) * qJD(5) * t291 + (0.2e1 * t237 * t201 + (-t200 + t205) * t291 + t170 * t235 * t337 + (-t182 * t235 + t186 * t237) * t338) * qJD(4)) * t236 + (-0.2e1 * m(6) * t154 + t135 * t336 - 0.2e1 * t203 + 0.2e1 * t204) * t194; m(6) * t5 + m(7) * t1 + t329 * t74 + t28 + t333 + t70 - t72; m(6) * t31 + m(7) * t12 + t329 * t175 + t106 + t167 - t169 + t332; -m(6) * t236 * t292 + m(7) * t179; 0; m(6) * t6 + m(7) * t2 + t42 + t45; m(6) * t38 + m(7) * t23 + t126 + t129; m(6) * t282 + m(7) * t193; 0; 0; m(7) * t3 + t44; m(7) * t27 + t128; m(7) * t194; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t35(1) t35(2) t35(4) t35(7) t35(11) t35(16); t35(2) t35(3) t35(5) t35(8) t35(12) t35(17); t35(4) t35(5) t35(6) t35(9) t35(13) t35(18); t35(7) t35(8) t35(9) t35(10) t35(14) t35(19); t35(11) t35(12) t35(13) t35(14) t35(15) t35(20); t35(16) t35(17) t35(18) t35(19) t35(20) t35(21);];
Mq  = res;
