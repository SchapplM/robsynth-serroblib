% Calculate time derivative of joint inertia matrix for
% S6PRRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-03-09 01:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRRR5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR5_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR5_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR5_inertiaDJ_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR5_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR5_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR5_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:04:02
% EndTime: 2019-03-09 01:04:18
% DurationCPUTime: 6.71s
% Computational Cost: add. (10462->715), mult. (29295->1057), div. (0->0), fcn. (29131->14), ass. (0->282)
t240 = sin(qJ(5));
t245 = cos(qJ(5));
t210 = -mrSges(6,1) * t245 + mrSges(6,2) * t240;
t347 = -m(6) * pkin(4) + t210;
t241 = sin(qJ(4));
t283 = qJD(5) * t241;
t246 = cos(qJ(4));
t285 = qJD(4) * t246;
t251 = -t240 * t283 + t245 * t285;
t209 = -pkin(4) * t246 - pkin(11) * t241 - pkin(3);
t293 = t245 * t246;
t228 = pkin(10) * t293;
t164 = t240 * t209 + t228;
t346 = qJD(5) * t164;
t239 = sin(qJ(6));
t244 = cos(qJ(6));
t258 = t239 * t240 - t244 * t245;
t321 = -t258 / 0.2e1;
t194 = t239 * t245 + t240 * t244;
t320 = t194 / 0.2e1;
t345 = t240 / 0.2e1;
t316 = t245 / 0.2e1;
t177 = t258 * t241;
t235 = sin(pkin(7));
t242 = sin(qJ(3));
t300 = t235 * t242;
t224 = pkin(9) * t300;
t237 = cos(pkin(7));
t247 = cos(qJ(3));
t315 = pkin(2) * t247;
t183 = t237 * t315 - t224;
t248 = cos(qJ(2));
t292 = t247 * t248;
t243 = sin(qJ(2));
t296 = t242 * t243;
t344 = t237 * t292 - t296;
t343 = -mrSges(5,1) + t347;
t342 = -m(5) * pkin(3) - mrSges(5,1) * t246 + mrSges(5,2) * t241 - mrSges(4,1);
t341 = qJD(5) + qJD(6);
t340 = 0.2e1 * m(6);
t339 = 2 * m(7);
t338 = 0.2e1 * pkin(10);
t337 = -2 * mrSges(4,3);
t179 = -t246 * t237 + t241 * t300;
t288 = qJD(3) * t235;
t268 = t247 * t288;
t140 = -qJD(4) * t179 + t246 * t268;
t180 = t237 * t241 + t246 * t300;
t299 = t235 * t247;
t253 = -t180 * t245 + t240 * t299;
t287 = qJD(3) * t242;
t269 = t235 * t287;
t71 = qJD(5) * t253 - t140 * t240 + t245 * t269;
t334 = t71 / 0.2e1;
t142 = -t180 * t240 - t245 * t299;
t79 = t142 * t244 + t239 * t253;
t333 = t79 / 0.2e1;
t80 = t142 * t239 - t244 * t253;
t332 = t80 / 0.2e1;
t331 = -pkin(12) - pkin(11);
t134 = t341 * t258;
t330 = -t134 / 0.2e1;
t135 = t341 * t194;
t329 = -t135 / 0.2e1;
t138 = Ifges(7,4) * t194 - Ifges(7,2) * t258;
t328 = t138 / 0.2e1;
t139 = Ifges(7,1) * t194 - Ifges(7,4) * t258;
t327 = t139 / 0.2e1;
t326 = t142 / 0.2e1;
t325 = -t253 / 0.2e1;
t311 = Ifges(6,4) * t240;
t261 = Ifges(6,1) * t245 - t311;
t171 = -Ifges(6,5) * t246 + t241 * t261;
t324 = t171 / 0.2e1;
t176 = t194 * t241;
t323 = -t176 / 0.2e1;
t322 = -t177 / 0.2e1;
t310 = Ifges(6,4) * t245;
t215 = Ifges(6,1) * t240 + t310;
t319 = t215 / 0.2e1;
t318 = -t240 / 0.2e1;
t317 = -t245 / 0.2e1;
t314 = pkin(10) * t240;
t313 = Ifges(5,4) * t241;
t312 = Ifges(5,4) * t246;
t309 = Ifges(6,6) * t240;
t236 = sin(pkin(6));
t238 = cos(pkin(6));
t294 = t243 * t247;
t295 = t242 * t248;
t252 = t237 * t295 + t294;
t129 = t236 * t252 + t238 * t300;
t178 = -t235 * t236 * t248 + t238 * t237;
t101 = t129 * t241 - t178 * t246;
t102 = t129 * t246 + t178 * t241;
t289 = qJD(2) * t236;
t270 = t243 * t289;
t263 = t235 * t270;
t97 = t238 * t268 + (t344 * qJD(3) + (-t237 * t296 + t292) * qJD(2)) * t236;
t42 = qJD(4) * t102 + t241 * t97 - t246 * t263;
t30 = t101 * t42;
t128 = -t344 * t236 - t238 * t299;
t96 = t238 * t269 + (t252 * qJD(3) + (t237 * t294 + t295) * qJD(2)) * t236;
t308 = t128 * t96;
t307 = t42 * t241;
t43 = -qJD(4) * t101 + t241 * t263 + t246 * t97;
t306 = t43 * t246;
t166 = t224 + (-pkin(3) - t315) * t237;
t103 = pkin(4) * t179 - pkin(11) * t180 + t166;
t184 = t237 * t242 * pkin(2) + pkin(9) * t299;
t167 = pkin(10) * t237 + t184;
t168 = (-pkin(3) * t247 - pkin(10) * t242 - pkin(2)) * t235;
t111 = t246 * t167 + t241 * t168;
t105 = -pkin(11) * t299 + t111;
t49 = t240 * t103 + t245 * t105;
t117 = mrSges(5,1) * t269 - mrSges(5,3) * t140;
t70 = qJD(5) * t142 + t140 * t245 + t240 * t269;
t36 = -mrSges(6,1) * t71 + mrSges(6,2) * t70;
t304 = -t117 + t36;
t146 = -mrSges(5,1) * t299 - mrSges(5,3) * t180;
t89 = -mrSges(6,1) * t142 - mrSges(6,2) * t253;
t303 = -t146 + t89;
t175 = t184 * qJD(3);
t301 = t128 * t175;
t298 = t240 * t241;
t297 = t241 * t245;
t82 = -Ifges(7,5) * t134 - Ifges(7,6) * t135;
t206 = (pkin(4) * t241 - pkin(11) * t246) * qJD(4);
t286 = qJD(4) * t241;
t291 = t245 * t206 + t286 * t314;
t290 = -mrSges(4,1) * t237 + mrSges(5,1) * t179 + mrSges(5,2) * t180 + mrSges(4,3) * t300;
t284 = qJD(5) * t240;
t282 = qJD(5) * t245;
t281 = qJD(6) * t239;
t280 = qJD(6) * t244;
t141 = qJD(4) * t180 + t241 * t268;
t23 = qJD(6) * t79 + t239 * t71 + t244 * t70;
t24 = -qJD(6) * t80 - t239 * t70 + t244 * t71;
t7 = Ifges(7,5) * t23 + Ifges(7,6) * t24 + Ifges(7,3) * t141;
t31 = Ifges(6,5) * t70 + Ifges(6,6) * t71 + Ifges(6,3) * t141;
t90 = -t135 * t241 - t258 * t285;
t91 = t341 * t177 - t194 * t285;
t44 = Ifges(7,5) * t90 + Ifges(7,6) * t91 + Ifges(7,3) * t286;
t278 = pkin(5) * t284;
t277 = Ifges(5,6) * t299;
t56 = t102 * t245 + t128 * t240;
t17 = -qJD(5) * t56 - t240 * t43 + t245 * t96;
t55 = -t102 * t240 + t128 * t245;
t18 = qJD(5) * t55 + t240 * t96 + t245 * t43;
t28 = -t239 * t56 + t244 * t55;
t5 = qJD(6) * t28 + t17 * t239 + t18 * t244;
t29 = t239 * t55 + t244 * t56;
t6 = -qJD(6) * t29 + t17 * t244 - t18 * t239;
t276 = t6 * mrSges(7,1) - t5 * mrSges(7,2);
t66 = -Ifges(6,1) * t253 + Ifges(6,4) * t142 + Ifges(6,5) * t179;
t274 = t66 * t316;
t233 = Ifges(6,5) * t282;
t273 = -Ifges(6,6) * t284 / 0.2e1 + t233 / 0.2e1 + t82 / 0.2e1;
t272 = Ifges(5,5) * t140 - Ifges(5,6) * t141 + Ifges(5,3) * t269;
t271 = qJD(5) * t331;
t265 = Ifges(6,5) * t345 + Ifges(7,5) * t320 + Ifges(6,6) * t316 + Ifges(7,6) * t321;
t48 = t245 * t103 - t105 * t240;
t110 = -t241 * t167 + t168 * t246;
t108 = t240 * t206 + t209 * t282 + (-t245 * t286 - t246 * t284) * pkin(10);
t190 = t245 * t209;
t163 = -t246 * t314 + t190;
t264 = -qJD(5) * t163 + t108;
t104 = pkin(4) * t299 - t110;
t262 = mrSges(6,1) * t240 + mrSges(6,2) * t245;
t260 = -Ifges(6,2) * t240 + t310;
t213 = Ifges(6,2) * t245 + t311;
t173 = (pkin(3) * t242 - pkin(10) * t247) * t288;
t174 = t183 * qJD(3);
t57 = -t167 * t286 + t168 * t285 + t241 * t173 + t246 * t174;
t53 = pkin(11) * t269 + t57;
t67 = pkin(4) * t141 - pkin(11) * t140 + t175;
t15 = t103 * t282 - t105 * t284 + t240 * t67 + t245 * t53;
t16 = -qJD(5) * t49 - t240 * t53 + t245 * t67;
t259 = t15 * t245 - t16 * t240;
t35 = pkin(5) * t179 + pkin(12) * t253 + t48;
t40 = pkin(12) * t142 + t49;
t13 = -t239 * t40 + t244 * t35;
t14 = t239 * t35 + t244 * t40;
t125 = -pkin(12) * t297 + t190 + (-pkin(5) - t314) * t246;
t144 = -pkin(12) * t298 + t164;
t77 = t125 * t244 - t144 * t239;
t78 = t125 * t239 + t144 * t244;
t218 = t331 * t240;
t219 = t331 * t245;
t147 = t218 * t244 + t219 * t239;
t148 = t218 * t239 - t219 * t244;
t204 = t240 * t271;
t205 = t245 * t271;
t94 = qJD(6) * t147 + t204 * t244 + t205 * t239;
t95 = -qJD(6) * t148 - t204 * t239 + t205 * t244;
t257 = t95 * mrSges(7,1) - t94 * mrSges(7,2) + t82;
t58 = -t167 * t285 - t168 * t286 + t173 * t246 - t241 * t174;
t11 = pkin(5) * t141 - pkin(12) * t70 + t16;
t12 = pkin(12) * t71 + t15;
t2 = qJD(6) * t13 + t11 * t239 + t12 * t244;
t3 = -qJD(6) * t14 + t11 * t244 - t12 * t239;
t256 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + t7;
t74 = (pkin(5) * t241 - pkin(12) * t293) * qJD(4) + (-t228 + (pkin(12) * t241 - t209) * t240) * qJD(5) + t291;
t250 = t240 * t285 + t241 * t282;
t85 = -pkin(12) * t250 + t108;
t26 = qJD(6) * t77 + t239 * t74 + t244 * t85;
t27 = -qJD(6) * t78 - t239 * t85 + t244 * t74;
t255 = t27 * mrSges(7,1) - t26 * mrSges(7,2) + t44;
t254 = t101 * t285 + t307;
t54 = -pkin(4) * t269 - t58;
t119 = t251 * Ifges(6,5) - Ifges(6,6) * t250 + Ifges(6,3) * t286;
t234 = Ifges(5,5) * t285;
t230 = -pkin(5) * t245 - pkin(4);
t221 = Ifges(4,5) * t268;
t216 = Ifges(5,1) * t241 + t312;
t214 = Ifges(5,2) * t246 + t313;
t208 = (pkin(5) * t240 + pkin(10)) * t241;
t203 = -mrSges(6,1) * t246 - mrSges(6,3) * t297;
t202 = mrSges(6,2) * t246 - mrSges(6,3) * t298;
t201 = (Ifges(5,1) * t246 - t313) * qJD(4);
t200 = t261 * qJD(5);
t199 = (-Ifges(5,2) * t241 + t312) * qJD(4);
t198 = t260 * qJD(5);
t196 = (mrSges(5,1) * t241 + mrSges(5,2) * t246) * qJD(4);
t195 = t262 * qJD(5);
t192 = -mrSges(4,2) * t237 + mrSges(4,3) * t299;
t187 = (-mrSges(7,1) * t239 - mrSges(7,2) * t244) * qJD(6) * pkin(5);
t185 = t262 * t241;
t172 = (mrSges(4,1) * t242 + mrSges(4,2) * t247) * t288;
t170 = -Ifges(6,6) * t246 + t241 * t260;
t169 = -Ifges(6,3) * t246 + (Ifges(6,5) * t245 - t309) * t241;
t159 = pkin(5) * t250 + pkin(10) * t285;
t155 = -mrSges(6,2) * t286 - mrSges(6,3) * t250;
t154 = mrSges(6,1) * t286 - mrSges(6,3) * t251;
t150 = -mrSges(7,1) * t246 + mrSges(7,3) * t177;
t149 = mrSges(7,2) * t246 - mrSges(7,3) * t176;
t145 = mrSges(5,2) * t299 - mrSges(5,3) * t179;
t136 = mrSges(7,1) * t258 + mrSges(7,2) * t194;
t124 = mrSges(6,1) * t250 + mrSges(6,2) * t251;
t122 = mrSges(7,1) * t176 - mrSges(7,2) * t177;
t121 = -t215 * t283 + (Ifges(6,5) * t241 + t246 * t261) * qJD(4);
t120 = -t213 * t283 + (Ifges(6,6) * t241 + t246 * t260) * qJD(4);
t118 = -mrSges(5,2) * t269 - mrSges(5,3) * t141;
t116 = Ifges(5,1) * t180 - Ifges(5,4) * t179 - Ifges(5,5) * t299;
t115 = Ifges(5,4) * t180 - Ifges(5,2) * t179 - t277;
t114 = -Ifges(7,1) * t177 - Ifges(7,4) * t176 - Ifges(7,5) * t246;
t113 = -Ifges(7,4) * t177 - Ifges(7,2) * t176 - Ifges(7,6) * t246;
t112 = -Ifges(7,5) * t177 - Ifges(7,6) * t176 - Ifges(7,3) * t246;
t109 = t291 - t346;
t107 = mrSges(6,1) * t179 + mrSges(6,3) * t253;
t106 = -mrSges(6,2) * t179 + mrSges(6,3) * t142;
t86 = mrSges(5,1) * t141 + mrSges(5,2) * t140;
t84 = -Ifges(7,1) * t134 - Ifges(7,4) * t135;
t83 = -Ifges(7,4) * t134 - Ifges(7,2) * t135;
t81 = mrSges(7,1) * t135 - mrSges(7,2) * t134;
t76 = -mrSges(7,2) * t286 + mrSges(7,3) * t91;
t75 = mrSges(7,1) * t286 - mrSges(7,3) * t90;
t73 = Ifges(5,1) * t140 - Ifges(5,4) * t141 + Ifges(5,5) * t269;
t72 = Ifges(5,4) * t140 - Ifges(5,2) * t141 + Ifges(5,6) * t269;
t65 = -Ifges(6,4) * t253 + Ifges(6,2) * t142 + Ifges(6,6) * t179;
t64 = -Ifges(6,5) * t253 + Ifges(6,6) * t142 + Ifges(6,3) * t179;
t61 = -pkin(5) * t142 + t104;
t60 = mrSges(7,1) * t179 - mrSges(7,3) * t80;
t59 = -mrSges(7,2) * t179 + mrSges(7,3) * t79;
t51 = -mrSges(6,2) * t141 + mrSges(6,3) * t71;
t50 = mrSges(6,1) * t141 - mrSges(6,3) * t70;
t47 = -mrSges(7,1) * t91 + mrSges(7,2) * t90;
t46 = Ifges(7,1) * t90 + Ifges(7,4) * t91 + Ifges(7,5) * t286;
t45 = Ifges(7,4) * t90 + Ifges(7,2) * t91 + Ifges(7,6) * t286;
t41 = -mrSges(7,1) * t79 + mrSges(7,2) * t80;
t39 = Ifges(7,1) * t80 + Ifges(7,4) * t79 + Ifges(7,5) * t179;
t38 = Ifges(7,4) * t80 + Ifges(7,2) * t79 + Ifges(7,6) * t179;
t37 = Ifges(7,5) * t80 + Ifges(7,6) * t79 + Ifges(7,3) * t179;
t34 = -pkin(5) * t71 + t54;
t33 = Ifges(6,1) * t70 + Ifges(6,4) * t71 + Ifges(6,5) * t141;
t32 = Ifges(6,4) * t70 + Ifges(6,2) * t71 + Ifges(6,6) * t141;
t20 = -mrSges(7,2) * t141 + mrSges(7,3) * t24;
t19 = mrSges(7,1) * t141 - mrSges(7,3) * t23;
t10 = -mrSges(7,1) * t24 + mrSges(7,2) * t23;
t9 = Ifges(7,1) * t23 + Ifges(7,4) * t24 + Ifges(7,5) * t141;
t8 = Ifges(7,4) * t23 + Ifges(7,2) * t24 + Ifges(7,6) * t141;
t1 = [0.2e1 * m(7) * (t28 * t6 + t29 * t5 + t30) + 0.2e1 * m(6) * (t17 * t55 + t18 * t56 + t30) + 0.2e1 * m(5) * (t102 * t43 + t30 + t308) + 0.2e1 * m(4) * (t129 * t97 + t178 * t263 + t308); t102 * t118 + t18 * t106 + t17 * t107 + t128 * t86 + t43 * t145 + t178 * t172 + t28 * t19 + t97 * t192 + t29 * t20 + t5 * t59 + t55 * t50 + t56 * t51 + t6 * t60 + t290 * t96 + (-mrSges(3,1) * t243 - mrSges(3,2) * t248) * t289 + (t41 + t303) * t42 + (t10 + t304) * t101 + ((-mrSges(4,1) * t247 + mrSges(4,2) * t242) * t263 + (t128 * t247 - t129 * t242) * qJD(3) * mrSges(4,3)) * t235 + m(4) * (-pkin(2) * t235 ^ 2 * t270 + t129 * t174 - t183 * t96 + t184 * t97 + t301) + m(5) * (-t101 * t58 + t102 * t57 - t110 * t42 + t111 * t43 + t166 * t96 + t301) + m(7) * (t101 * t34 + t13 * t6 + t14 * t5 + t2 * t29 + t28 * t3 + t42 * t61) + m(6) * (t101 * t54 + t104 * t42 + t15 * t56 + t16 * t55 + t17 * t48 + t18 * t49); 0.2e1 * t290 * t175 - t253 * t33 + t180 * t73 + 0.2e1 * t174 * t192 + 0.2e1 * t166 * t86 + t142 * t32 + 0.2e1 * t57 * t145 + 0.2e1 * t58 * t146 + t140 * t116 + 0.2e1 * t110 * t117 + 0.2e1 * t111 * t118 + 0.2e1 * t104 * t36 + 0.2e1 * t15 * t106 + 0.2e1 * t16 * t107 + 0.2e1 * t54 * t89 + t79 * t8 + t80 * t9 + t70 * t66 + t71 * t65 + 0.2e1 * t2 * t59 + 0.2e1 * t3 * t60 + 0.2e1 * t61 * t10 + 0.2e1 * t48 * t50 + 0.2e1 * t49 * t51 + 0.2e1 * t34 * t41 + t24 * t38 + t23 * t39 + 0.2e1 * t14 * t20 + 0.2e1 * t13 * t19 + (t13 * t3 + t14 * t2 + t34 * t61) * t339 + (t104 * t54 + t15 * t49 + t16 * t48) * t340 + (t31 + t7 - t72) * t179 + 0.2e1 * m(4) * (t174 * t184 - t175 * t183) + 0.2e1 * m(5) * (t110 * t58 + t111 * t57 + t166 * t175) + (-t115 + t64 + t37) * t141 + t237 * t221 + (-t247 * t272 - 0.2e1 * pkin(2) * t172 + ((0.2e1 * Ifges(4,4) * t299 + Ifges(4,5) * t237 + t183 * t337) * t247 + (-0.2e1 * Ifges(4,4) * t300 + t184 * t337 + Ifges(5,5) * t180 - 0.2e1 * Ifges(4,6) * t237 - Ifges(5,6) * t179 + ((2 * Ifges(4,1)) - (2 * Ifges(4,2)) - Ifges(5,3)) * t299) * t242) * qJD(3)) * t235; -t97 * mrSges(4,2) + t128 * t196 + t5 * t149 + t6 * t150 + t55 * t154 + t56 * t155 + t17 * t203 + t18 * t202 + t28 * t75 + t29 * t76 + (t122 + t185) * t42 + (t124 + t47) * t101 + m(7) * (t101 * t159 + t208 * t42 + t26 * t29 + t27 * t28 + t5 * t78 + t6 * t77) + m(6) * (t108 * t56 + t109 * t55 + t163 * t17 + t164 * t18) + (m(6) * t254 / 0.2e1 + m(5) * (-t102 * t286 + t254 + t306) / 0.2e1) * t338 + (t307 + t306 + (t101 * t246 - t102 * t241) * qJD(4)) * mrSges(5,3) + t342 * t96; t221 + t140 * t216 / 0.2e1 + t166 * t196 + t180 * t201 / 0.2e1 + t15 * t202 + t16 * t203 + t208 * t10 + t54 * t185 - t174 * mrSges(4,2) + t2 * t149 + t3 * t150 + t48 * t154 + t49 * t155 + t159 * t41 + t163 * t50 + t164 * t51 + t24 * t113 / 0.2e1 + t23 * t114 / 0.2e1 + t34 * t122 + t104 * t124 + t108 * t106 + t109 * t107 + t90 * t39 / 0.2e1 + t91 * t38 / 0.2e1 - pkin(3) * t86 + t13 * t75 + t14 * t76 + t77 * t19 + t78 * t20 + t26 * t59 + t27 * t60 + t61 * t47 + m(6) * (t108 * t49 + t109 * t48 + t15 * t164 + t16 * t163) + t9 * t322 + t8 * t323 + t70 * t324 + t121 * t325 + t120 * t326 + t46 * t332 + t45 * t333 + t170 * t334 + m(7) * (t13 * t27 + t14 * t26 + t159 * t61 + t2 * t78 + t208 * t34 + t3 * t77) + (-t199 / 0.2e1 + t119 / 0.2e1 + t44 / 0.2e1) * t179 + (-t214 / 0.2e1 + t169 / 0.2e1 + t112 / 0.2e1) * t141 + t342 * t175 + (t57 * mrSges(5,3) + t72 / 0.2e1 - t31 / 0.2e1 - t7 / 0.2e1) * t246 + (-t247 * t234 / 0.2e1 + (Ifges(5,5) * t241 / 0.2e1 + Ifges(5,6) * t246 / 0.2e1 - Ifges(4,6)) * t287) * t235 + (t246 * t118 + t304 * t241 + (-t241 * t145 + t246 * t303) * qJD(4) + m(6) * (t104 * t285 + t241 * t54) + m(5) * (-t110 * t285 - t111 * t286 - t58 * t241 + t57 * t246)) * pkin(10) + ((t116 / 0.2e1 + t274 + t65 * t318 - t110 * mrSges(5,3)) * t246 + (t277 / 0.2e1 - t115 / 0.2e1 + t64 / 0.2e1 + t37 / 0.2e1 - t111 * mrSges(5,3)) * t241) * qJD(4) + (t33 * t316 + t32 * t318 - t58 * mrSges(5,3) + t73 / 0.2e1 + (t317 * t65 + t318 * t66) * qJD(5)) * t241; -0.2e1 * pkin(3) * t196 + 0.2e1 * t108 * t202 + 0.2e1 * t109 * t203 + t91 * t113 + t90 * t114 + 0.2e1 * t159 * t122 + 0.2e1 * t26 * t149 + 0.2e1 * t27 * t150 + 0.2e1 * t163 * t154 + 0.2e1 * t164 * t155 - t176 * t45 - t177 * t46 + 0.2e1 * t208 * t47 + 0.2e1 * t77 * t75 + 0.2e1 * t78 * t76 + (t108 * t164 + t109 * t163) * t340 + (t159 * t208 + t26 * t78 + t27 * t77) * t339 + (-t119 + t199 - t44 + (-t240 * t170 + t171 * t245 + t185 * t338 + t216) * qJD(4)) * t246 + (t124 * t338 - t240 * t120 + t245 * t121 + t201 + (-t170 * t245 - t171 * t240) * qJD(5) + (pkin(10) ^ 2 * t246 * t340 + t112 + t169 - t214) * qJD(4)) * t241; -t43 * mrSges(5,2) + (t195 + t81) * t101 + m(7) * (t101 * t278 + t147 * t6 + t148 * t5 + t28 * t95 + t29 * t94) + (t134 * t28 - t135 * t29 - t194 * t6 - t258 * t5) * mrSges(7,3) + (m(7) * t230 + t136 + t343) * t42 + (m(6) * pkin(11) + mrSges(6,3)) * (-t17 * t240 + t18 * t245 + (-t240 * t56 - t245 * t55) * qJD(5)); t347 * t54 + (t13 * t134 - t135 * t14 - t194 * t3 - t2 * t258) * mrSges(7,3) + t230 * t10 + t104 * t195 + t147 * t19 + t148 * t20 + t34 * t136 + t95 * t60 + t94 * t59 + t61 * t81 + t58 * mrSges(5,1) - t57 * mrSges(5,2) - pkin(4) * t36 + t33 * t345 + t70 * t319 + t9 * t320 + t8 * t321 + t200 * t325 + t198 * t326 + t23 * t327 + t24 * t328 + t38 * t329 + t39 * t330 + t84 * t332 + t83 * t333 + t213 * t334 + t32 * t316 + t272 + ((-t240 * t49 - t245 * t48) * qJD(5) + t259) * mrSges(6,3) + t265 * t141 + t273 * t179 + (t274 + (-t65 / 0.2e1 + pkin(5) * t41) * t240) * qJD(5) + m(7) * (t13 * t95 + t14 * t94 + t147 * t3 + t148 * t2 + t230 * t34 + t278 * t61) + (m(6) * (-t282 * t48 - t284 * t49 + t259) + t245 * t51 - t240 * t50 - t107 * t282 - t106 * t284) * pkin(11); t234 + m(7) * (t147 * t27 + t148 * t26 + t159 * t230 + t77 * t95 + t78 * t94) + t230 * t47 + t208 * t81 + t45 * t321 + t46 * t320 + t84 * t322 + t83 * t323 + t94 * t149 + t95 * t150 + t159 * t136 + t147 * t75 + t148 * t76 + t114 * t330 + t113 * t329 + t91 * t328 + t90 * t327 - pkin(4) * t124 + (t343 * qJD(4) * pkin(10) - t273) * t246 + (-t109 * mrSges(6,3) + t121 / 0.2e1 - t213 * t285 / 0.2e1 + (-t170 / 0.2e1 - t164 * mrSges(6,3) + (m(7) * t208 + t122) * pkin(5)) * qJD(5) + (m(6) * (-t109 - t346) - t154 - qJD(5) * t202) * pkin(11)) * t240 + (t134 * t77 - t135 * t78 - t194 * t27 - t258 * t26) * mrSges(7,3) + (qJD(5) * t324 + t120 / 0.2e1 + t285 * t319 + t264 * mrSges(6,3) + (m(6) * t264 - qJD(5) * t203 + t155) * pkin(11)) * t245 + (t200 * t316 + t198 * t318 + pkin(10) * t195 + (t213 * t317 + t215 * t318) * qJD(5) + (pkin(10) * mrSges(5,2) - Ifges(5,6) + t265) * qJD(4)) * t241; -t134 * t139 + t194 * t84 - t135 * t138 - t258 * t83 + 0.2e1 * t136 * t278 + 0.2e1 * t230 * t81 + (t147 * t95 + t148 * t94 + t230 * t278) * t339 + t240 * t200 - t213 * t284 - 0.2e1 * pkin(4) * t195 + (qJD(5) * t215 + t198) * t245 + 0.2e1 * (t134 * t147 - t135 * t148 - t194 * t95 - t258 * t94) * mrSges(7,3); t17 * mrSges(6,1) - t18 * mrSges(6,2) + m(7) * (t239 * t5 + t244 * t6 + (-t239 * t28 + t244 * t29) * qJD(6)) * pkin(5) + t276; t16 * mrSges(6,1) - t15 * mrSges(6,2) + (m(7) * (-t13 * t281 + t14 * t280 + t2 * t239 + t244 * t3) + t59 * t280 + t239 * t20 - t60 * t281 + t244 * t19) * pkin(5) + t256 + t31; t109 * mrSges(6,1) - t108 * mrSges(6,2) + (m(7) * (t239 * t26 + t244 * t27 + t280 * t78 - t281 * t77) + t149 * t280 + t239 * t76 - t150 * t281 + t244 * t75) * pkin(5) + t119 + t255; t233 + (pkin(11) * t210 - t309) * qJD(5) + (m(7) * (t239 * t94 + t244 * t95 + (-t147 * t239 + t148 * t244) * qJD(6)) + (t244 * t134 - t239 * t135 + (t194 * t239 - t244 * t258) * qJD(6)) * mrSges(7,3)) * pkin(5) + t257; 0.2e1 * t187; t276; t256; t255; t257; t187; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
