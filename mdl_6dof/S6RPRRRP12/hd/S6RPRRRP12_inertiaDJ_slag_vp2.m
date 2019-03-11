% Calculate time derivative of joint inertia matrix for
% S6RPRRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 06:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRP12_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP12_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP12_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP12_inertiaDJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP12_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP12_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP12_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:44:13
% EndTime: 2019-03-09 06:44:30
% DurationCPUTime: 7.37s
% Computational Cost: add. (13450->662), mult. (39317->929), div. (0->0), fcn. (41206->12), ass. (0->268)
t217 = sin(qJ(5));
t343 = (Ifges(6,5) + Ifges(7,4)) * t217;
t212 = sin(pkin(7));
t215 = cos(pkin(7));
t216 = cos(pkin(6));
t211 = sin(pkin(12));
t213 = sin(pkin(6));
t214 = cos(pkin(12));
t295 = t213 * t214;
t327 = pkin(1) * t216;
t282 = qJ(2) * t295 + t211 * t327;
t126 = (t212 * t216 + t215 * t295) * pkin(9) + t282;
t219 = sin(qJ(3));
t222 = cos(qJ(3));
t202 = t214 * t327;
t299 = t211 * t213;
t129 = pkin(2) * t216 + t202 + (-pkin(9) * t215 - qJ(2)) * t299;
t142 = (-pkin(9) * t211 * t212 - pkin(2) * t214 - pkin(1)) * t213;
t232 = t129 * t215 + t142 * t212;
t78 = -t126 * t219 + t222 * t232;
t218 = sin(qJ(4));
t221 = cos(qJ(4));
t277 = qJD(4) * t221;
t257 = t217 * t277;
t220 = cos(qJ(5));
t274 = qJD(5) * t220;
t226 = t218 * t274 + t257;
t293 = t215 * t222;
t296 = t212 * t222;
t341 = t213 * (-t211 * t219 + t214 * t293) + t216 * t296;
t340 = -2 * Ifges(4,4);
t294 = t215 * t219;
t297 = t212 * t219;
t128 = t216 * t297 + (t211 * t222 + t214 * t294) * t213;
t122 = t128 * qJD(3);
t121 = t341 * qJD(3);
t154 = -t212 * t295 + t215 * t216;
t99 = t128 * t218 - t154 * t221;
t86 = -qJD(4) * t99 + t121 * t221;
t100 = t128 * t221 + t154 * t218;
t88 = t100 * t220 - t217 * t341;
t49 = qJD(5) * t88 - t122 * t220 + t217 * t86;
t87 = t100 * t217 + t220 * t341;
t50 = -qJD(5) * t87 + t122 * t217 + t220 * t86;
t85 = qJD(4) * t100 + t121 * t218;
t10 = Ifges(7,4) * t50 + Ifges(7,2) * t85 + Ifges(7,6) * t49;
t9 = Ifges(6,5) * t50 - Ifges(6,6) * t49 + Ifges(6,3) * t85;
t338 = t10 + t9;
t336 = m(7) * qJ(6) + mrSges(7,3);
t178 = (pkin(4) * t218 - pkin(11) * t221) * qJD(4);
t182 = -pkin(4) * t221 - pkin(11) * t218 - pkin(3);
t273 = qJD(5) * t221;
t276 = qJD(5) * t217;
t278 = qJD(4) * t218;
t106 = pkin(10) * (t217 * t278 - t220 * t273) + t178 * t220 - t182 * t276;
t281 = qJD(2) * t213;
t66 = (-t211 * t294 + t214 * t222) * t281 + t78 * qJD(3);
t97 = -t129 * t212 + t142 * t215;
t71 = -pkin(3) * t341 - pkin(10) * t128 + t97;
t117 = t222 * t126;
t79 = t129 * t294 + t142 * t297 + t117;
t76 = pkin(10) * t154 + t79;
t262 = t211 * t281;
t244 = t212 * t262;
t92 = pkin(3) * t122 - pkin(10) * t121 + t244;
t19 = t218 * t92 + t221 * t66 + t277 * t71 - t278 * t76;
t17 = pkin(11) * t122 + t19;
t67 = (t211 * t293 + t214 * t219) * t281 + (t219 * t232 + t117) * qJD(3);
t29 = t85 * pkin(4) - t86 * pkin(11) + t67;
t34 = t218 * t71 + t221 * t76;
t32 = -pkin(11) * t341 + t34;
t75 = -pkin(3) * t154 - t78;
t48 = pkin(4) * t99 - pkin(11) * t100 + t75;
t322 = t217 * t48 + t220 * t32;
t4 = -qJD(5) * t322 - t17 * t217 + t220 * t29;
t236 = pkin(5) * t220 + qJ(6) * t217;
t272 = qJD(6) * t220;
t335 = qJD(5) * t236 - t272;
t334 = 2 * m(6);
t333 = 2 * m(7);
t332 = 0.2e1 * pkin(10);
t331 = -2 * mrSges(4,3);
t330 = m(6) / 0.2e1;
t329 = m(5) * pkin(3);
t20 = -t218 * t66 + t221 * t92 - t277 * t76 - t278 * t71;
t18 = -pkin(4) * t122 - t20;
t328 = m(6) * t18;
t326 = pkin(10) * t221;
t23 = mrSges(6,1) * t49 + mrSges(6,2) * t50;
t63 = mrSges(5,1) * t122 - mrSges(5,3) * t86;
t325 = t23 - t63;
t24 = -mrSges(6,2) * t85 - mrSges(6,3) * t49;
t27 = -mrSges(7,2) * t49 + mrSges(7,3) * t85;
t324 = t24 + t27;
t25 = mrSges(6,1) * t85 - mrSges(6,3) * t50;
t26 = -mrSges(7,1) * t85 + mrSges(7,2) * t50;
t323 = -t25 + t26;
t55 = mrSges(6,1) * t87 + mrSges(6,2) * t88;
t90 = -mrSges(5,1) * t341 - mrSges(5,3) * t100;
t321 = t55 - t90;
t56 = -mrSges(7,2) * t87 + mrSges(7,3) * t99;
t57 = -mrSges(6,2) * t99 - mrSges(6,3) * t87;
t320 = t56 + t57;
t58 = mrSges(6,1) * t99 - mrSges(6,3) * t88;
t59 = -mrSges(7,1) * t99 + mrSges(7,2) * t88;
t319 = -t58 + t59;
t318 = Ifges(5,4) * t218;
t317 = Ifges(5,4) * t221;
t316 = Ifges(6,4) * t217;
t315 = Ifges(6,4) * t220;
t314 = Ifges(7,5) * t217;
t313 = Ifges(7,5) * t220;
t312 = Ifges(6,6) * t220;
t311 = t122 * Ifges(5,5);
t310 = t122 * Ifges(5,6);
t309 = t341 * Ifges(5,6);
t308 = t218 * mrSges(5,2);
t307 = t222 * t67;
t155 = -t215 * t221 + t218 * t297;
t279 = qJD(3) * t222;
t260 = t212 * t279;
t131 = -qJD(4) * t155 + t221 * t260;
t156 = t215 * t218 + t221 * t297;
t133 = t156 * t217 + t220 * t296;
t280 = qJD(3) * t219;
t261 = t212 * t280;
t95 = -qJD(5) * t133 + t131 * t220 + t217 * t261;
t306 = t95 * t220;
t265 = t217 * t296;
t96 = -qJD(5) * t265 + t131 * t217 + t156 * t274 - t220 * t261;
t305 = t96 * t217;
t304 = -mrSges(5,1) * t221 - mrSges(4,1) + t308;
t184 = -t220 * mrSges(6,1) + t217 * mrSges(6,2);
t303 = t184 - mrSges(5,1);
t302 = t131 * t221;
t132 = qJD(4) * t156 + t218 * t260;
t301 = t132 * t218;
t107 = t155 * t132;
t300 = t182 * t220;
t292 = t217 * t218;
t291 = t218 * t220;
t290 = Ifges(4,5) * t121 - Ifges(4,6) * t122;
t256 = t220 * t277;
t275 = qJD(5) * t218;
t259 = t217 * t275;
t227 = t256 - t259;
t138 = mrSges(6,1) * t278 - mrSges(6,3) * t227;
t139 = mrSges(7,2) * t256 + (-mrSges(7,1) * qJD(4) - mrSges(7,2) * t276) * t218;
t289 = -t138 + t139;
t140 = -mrSges(6,2) * t278 - mrSges(6,3) * t226;
t141 = -mrSges(7,2) * t226 + mrSges(7,3) * t278;
t288 = t140 + t141;
t237 = Ifges(7,3) * t217 + t313;
t146 = -Ifges(7,6) * t221 + t218 * t237;
t238 = -Ifges(6,2) * t217 + t315;
t149 = -Ifges(6,6) * t221 + t218 * t238;
t287 = t146 - t149;
t239 = Ifges(7,1) * t220 + t314;
t150 = -Ifges(7,4) * t221 + t218 * t239;
t240 = Ifges(6,1) * t220 - t316;
t151 = -Ifges(6,5) * t221 + t218 * t240;
t286 = t150 + t151;
t174 = mrSges(6,2) * t221 - mrSges(6,3) * t292;
t177 = -mrSges(7,2) * t292 - mrSges(7,3) * t221;
t285 = t174 + t177;
t175 = -mrSges(6,1) * t221 - mrSges(6,3) * t291;
t176 = mrSges(7,1) * t221 + mrSges(7,2) * t291;
t284 = -t175 + t176;
t283 = Ifges(6,5) * t256 + Ifges(6,3) * t278;
t144 = t182 * t217 + t220 * t326;
t167 = Ifges(7,4) * t274 + Ifges(7,6) * t276;
t11 = Ifges(6,4) * t50 - Ifges(6,2) * t49 + Ifges(6,6) * t85;
t8 = Ifges(7,5) * t50 + Ifges(7,6) * t85 + Ifges(7,3) * t49;
t270 = t8 / 0.2e1 - t11 / 0.2e1;
t269 = Ifges(5,5) * t86 - Ifges(5,6) * t85 + Ifges(5,3) * t122;
t12 = Ifges(7,1) * t50 + Ifges(7,4) * t85 + Ifges(7,5) * t49;
t13 = Ifges(6,1) * t50 - Ifges(6,4) * t49 + Ifges(6,5) * t85;
t268 = t12 / 0.2e1 + t13 / 0.2e1;
t35 = Ifges(7,5) * t88 + Ifges(7,6) * t99 + Ifges(7,3) * t87;
t38 = Ifges(6,4) * t88 - Ifges(6,2) * t87 + Ifges(6,6) * t99;
t267 = t35 / 0.2e1 - t38 / 0.2e1;
t39 = Ifges(7,1) * t88 + Ifges(7,4) * t99 + Ifges(7,5) * t87;
t40 = Ifges(6,1) * t88 - Ifges(6,4) * t87 + Ifges(6,5) * t99;
t266 = t39 / 0.2e1 + t40 / 0.2e1;
t186 = -Ifges(7,3) * t220 + t314;
t108 = -t186 * t275 + (Ifges(7,6) * t218 + t221 * t237) * qJD(4);
t189 = Ifges(6,2) * t220 + t316;
t111 = -t189 * t275 + (Ifges(6,6) * t218 + t221 * t238) * qJD(4);
t255 = t108 / 0.2e1 - t111 / 0.2e1;
t191 = Ifges(7,1) * t217 - t313;
t112 = -t191 * t275 + (Ifges(7,4) * t218 + t221 * t239) * qJD(4);
t192 = Ifges(6,1) * t217 + t315;
t113 = -t192 * t275 + (Ifges(6,5) * t218 + t221 * t240) * qJD(4);
t254 = t112 / 0.2e1 + t113 / 0.2e1;
t253 = t146 / 0.2e1 - t149 / 0.2e1;
t252 = t150 / 0.2e1 + t151 / 0.2e1;
t165 = t237 * qJD(5);
t168 = t238 * qJD(5);
t251 = t165 / 0.2e1 - t168 / 0.2e1;
t166 = Ifges(6,5) * t274 - Ifges(6,6) * t276;
t250 = t166 / 0.2e1 + t167 / 0.2e1;
t170 = t239 * qJD(5);
t171 = t240 * qJD(5);
t249 = t170 / 0.2e1 + t171 / 0.2e1;
t248 = t186 / 0.2e1 - t189 / 0.2e1;
t247 = t312 / 0.2e1 - Ifges(7,6) * t220 / 0.2e1 + t343 / 0.2e1;
t246 = t191 / 0.2e1 + t192 / 0.2e1;
t33 = -t218 * t76 + t221 * t71;
t245 = Ifges(7,4) * t256 + Ifges(7,2) * t278 + Ifges(7,6) * t226;
t242 = mrSges(6,1) * t217 + mrSges(6,2) * t220;
t183 = -t220 * mrSges(7,1) - t217 * mrSges(7,3);
t241 = mrSges(7,1) * t217 - mrSges(7,3) * t220;
t235 = pkin(5) * t217 - qJ(6) * t220;
t14 = -t217 * t32 + t220 * t48;
t31 = pkin(4) * t341 - t33;
t231 = pkin(10) + t235;
t3 = t17 * t220 + t217 * t29 + t274 * t48 - t276 * t32;
t230 = t155 * t277 + t301;
t105 = t217 * t178 + t182 * t274 + (-t217 * t273 - t220 * t278) * pkin(10);
t209 = Ifges(5,5) * t277;
t193 = Ifges(5,1) * t218 + t317;
t190 = Ifges(5,2) * t221 + t318;
t179 = -pkin(4) - t236;
t172 = (Ifges(5,1) * t221 - t318) * qJD(4);
t169 = (-Ifges(5,2) * t218 + t317) * qJD(4);
t164 = (mrSges(5,1) * t218 + mrSges(5,2) * t221) * qJD(4);
t163 = t242 * qJD(5);
t162 = t241 * qJD(5);
t159 = t242 * t218;
t158 = t241 * t218;
t153 = qJD(5) * t235 - qJD(6) * t217;
t152 = t231 * t218;
t148 = -Ifges(7,2) * t221 + (Ifges(7,4) * t220 + Ifges(7,6) * t217) * t218;
t147 = -Ifges(6,3) * t221 + (Ifges(6,5) * t220 - Ifges(6,6) * t217) * t218;
t143 = -t217 * t326 + t300;
t137 = -t300 + (pkin(10) * t217 + pkin(5)) * t221;
t136 = -qJ(6) * t221 + t144;
t134 = t156 * t220 - t265;
t125 = mrSges(6,1) * t226 + mrSges(6,2) * t227;
t124 = mrSges(7,1) * t226 - mrSges(7,3) * t227;
t110 = -Ifges(7,4) * t259 + t245;
t109 = -Ifges(6,5) * t259 - Ifges(6,6) * t226 + t283;
t104 = t218 * t335 + t231 * t277;
t103 = mrSges(4,1) * t154 - mrSges(4,3) * t128;
t102 = -mrSges(4,2) * t154 + mrSges(4,3) * t341;
t101 = -pkin(5) * t278 - t106;
t98 = qJ(6) * t278 - qJD(6) * t221 + t105;
t94 = mrSges(4,1) * t122 + mrSges(4,2) * t121;
t93 = pkin(11) * t306;
t89 = mrSges(5,2) * t341 - mrSges(5,3) * t99;
t77 = mrSges(5,1) * t99 + mrSges(5,2) * t100;
t62 = -mrSges(5,2) * t122 - mrSges(5,3) * t85;
t61 = Ifges(5,1) * t100 - Ifges(5,4) * t99 - Ifges(5,5) * t341;
t60 = Ifges(5,4) * t100 - Ifges(5,2) * t99 - t309;
t54 = mrSges(7,1) * t87 - mrSges(7,3) * t88;
t53 = mrSges(5,1) * t85 + mrSges(5,2) * t86;
t52 = Ifges(5,1) * t86 - Ifges(5,4) * t85 + t311;
t51 = Ifges(5,4) * t86 - Ifges(5,2) * t85 + t310;
t37 = Ifges(7,4) * t88 + Ifges(7,2) * t99 + Ifges(7,6) * t87;
t36 = Ifges(6,5) * t88 - Ifges(6,6) * t87 + Ifges(6,3) * t99;
t22 = mrSges(7,1) * t49 - mrSges(7,3) * t50;
t21 = pkin(5) * t87 - qJ(6) * t88 + t31;
t7 = -pkin(5) * t99 - t14;
t6 = qJ(6) * t99 + t322;
t5 = pkin(5) * t49 - qJ(6) * t50 - qJD(6) * t88 + t18;
t2 = -pkin(5) * t85 - t4;
t1 = qJ(6) * t85 + qJD(6) * t99 + t3;
t15 = [(t12 + t13) * t88 + 0.2e1 * m(5) * (t19 * t34 + t20 * t33 + t67 * t75) + 0.2e1 * (t214 * (-mrSges(3,2) * t216 + mrSges(3,3) * t295) + m(3) * (t282 * t214 + (qJ(2) * t299 - t202) * t211)) * t281 + (t14 * t4 + t18 * t31 + t3 * t322) * t334 + 0.2e1 * t322 * t24 + (t39 + t40) * t50 + (t35 - t38) * t49 + (-t60 + t37 + t36) * t85 + (t1 * t6 + t2 * t7 + t21 * t5) * t333 + (-t11 + t8) * t87 + 0.2e1 * (-mrSges(4,1) * t341 + mrSges(4,2) * t128) * t244 + (t128 * t340 + Ifges(5,5) * t100 - Ifges(4,6) * t154 - Ifges(5,6) * t99 + t331 * t79 - ((2 * Ifges(4,2)) + Ifges(5,3)) * t341) * t122 + (0.2e1 * Ifges(4,1) * t128 + Ifges(4,5) * t154 + t331 * t78 - t340 * t341) * t121 - t341 * t269 - 0.2e1 * t67 * t103 + t100 * t52 + 0.2e1 * t66 * t102 + 0.2e1 * t97 * t94 + 0.2e1 * t19 * t89 + 0.2e1 * t20 * t90 + t86 * t61 + 0.2e1 * t75 * t53 + 0.2e1 * t67 * t77 + 0.2e1 * t1 * t56 + 0.2e1 * t3 * t57 + 0.2e1 * t4 * t58 + 0.2e1 * t2 * t59 + 0.2e1 * t34 * t62 + 0.2e1 * t33 * t63 + 0.2e1 * t5 * t54 + 0.2e1 * t18 * t55 + 0.2e1 * t31 * t23 + 0.2e1 * t14 * t25 + 0.2e1 * t7 * t26 + 0.2e1 * t6 * t27 + 0.2e1 * t21 * t22 + 0.2e1 * m(4) * (t244 * t97 + t66 * t79 - t67 * t78) + (-t51 + t338) * t99 + t154 * t290 - 0.2e1 * (mrSges(3,1) * t216 - mrSges(3,3) * t299) * t262; t131 * t89 + t156 * t62 + t215 * t94 + t319 * t96 + t320 * t95 + t324 * t134 + t323 * t133 + (t22 + t325) * t155 + (t54 + t321) * t132 + m(5) * (t131 * t34 - t132 * t33 - t155 * t20 + t156 * t19) + m(7) * (t1 * t134 + t132 * t21 + t133 * t2 + t155 * t5 + t6 * t95 + t7 * t96) + m(6) * (t132 * t31 - t133 * t4 + t134 * t3 - t14 * t96 + t155 * t18 + t322 * t95) + (-t222 * t53 + (-t121 * t222 - t122 * t219) * mrSges(4,3) + (t222 * t102 + (-t103 + t77) * t219) * qJD(3) + m(4) * (t215 * t262 + t219 * t66 + t279 * t79 - t280 * t78 - t307) + m(5) * (t280 * t75 - t307)) * t212; 0.2e1 * m(5) * (-t212 ^ 2 * t219 * t279 + t131 * t156 + t107) + 0.2e1 * (m(6) + m(7)) * (t133 * t96 + t134 * t95 + t107); m(6) * (t105 * t322 + t106 * t14 + t143 * t4 + t144 * t3) + t322 * t140 + m(7) * (t1 * t136 + t101 * t7 + t104 * t21 + t137 * t2 + t152 * t5 + t6 * t98) + (-t169 / 0.2e1 + t109 / 0.2e1 + t110 / 0.2e1) * t99 + (-t190 / 0.2e1 + t147 / 0.2e1 + t148 / 0.2e1) * t85 - t341 * t209 / 0.2e1 + (t311 / 0.2e1 + t52 / 0.2e1 - t20 * mrSges(5,3) + t268 * t220 + t270 * t217 + (-t217 * t266 + t220 * t267) * qJD(5) + (t309 / 0.2e1 + t37 / 0.2e1 - t60 / 0.2e1 + t36 / 0.2e1 - t34 * mrSges(5,3)) * qJD(4) + (-qJD(4) * t89 + t328 + m(5) * (-qJD(4) * t34 - t20) + t325) * pkin(10)) * t218 + t86 * t193 / 0.2e1 + t100 * t172 / 0.2e1 + t3 * t174 + t4 * t175 + t2 * t176 + t1 * t177 + t5 * t158 + t18 * t159 + t75 * t164 + t152 * t22 + t136 * t27 + t137 * t26 + t14 * t138 + t7 * t139 + t6 * t141 + t143 * t25 + t144 * t24 + t21 * t124 + t31 * t125 + t104 * t54 + t105 * t57 + t106 * t58 + t101 * t59 + t98 * t56 - t66 * mrSges(4,2) - pkin(3) * t53 + t290 + (t304 - t329) * t67 + t252 * t50 + t253 * t49 + t254 * t88 + t255 * t87 + (t310 / 0.2e1 + t51 / 0.2e1 - t9 / 0.2e1 - t10 / 0.2e1 + t19 * mrSges(5,3) + (m(5) * t19 + t62) * pkin(10) + (t61 / 0.2e1 - t33 * mrSges(5,3) + t266 * t220 + t267 * t217 + (-m(5) * t33 + m(6) * t31 + t321) * pkin(10)) * qJD(4)) * t221; t284 * t96 + t285 * t95 + (t124 + t125) * t155 + t288 * t134 + t289 * t133 + (t158 + t159) * t132 + (-t222 * t164 + (-mrSges(4,2) * t222 + t219 * t304) * qJD(3)) * t212 + m(6) * (t105 * t134 - t106 * t133 - t143 * t96 + t144 * t95) + m(7) * (t101 * t133 + t104 * t155 + t132 * t152 + t134 * t98 + t136 * t95 + t137 * t96) - t261 * t329 + (t230 * t330 + m(5) * (-t156 * t278 + t230 + t302) / 0.2e1) * t332 + (t302 + t301 + (t155 * t221 - t156 * t218) * qJD(4)) * mrSges(5,3); -0.2e1 * pkin(3) * t164 + 0.2e1 * t101 * t176 + 0.2e1 * t104 * t158 + 0.2e1 * t105 * t174 + 0.2e1 * t106 * t175 + 0.2e1 * t152 * t124 + 0.2e1 * t136 * t141 + 0.2e1 * t137 * t139 + 0.2e1 * t143 * t138 + 0.2e1 * t144 * t140 + 0.2e1 * t98 * t177 + (t105 * t144 + t106 * t143) * t334 + (t101 * t137 + t104 * t152 + t136 * t98) * t333 + (-t109 - t110 + t169 + (t159 * t332 + t217 * t287 + t220 * t286 + t193) * qJD(4)) * t221 + (t125 * t332 + t172 + (t112 + t113) * t220 + (t108 - t111) * t217 + (-t217 * t286 + t220 * t287) * qJD(5) + (pkin(10) ^ 2 * t221 * t334 + t147 + t148 - t190) * qJD(4)) * t218; m(7) * (t153 * t21 + t179 * t5) + (t1 * mrSges(7,2) + t3 * mrSges(6,3) - t270) * t220 + (t324 * t220 + t323 * t217 + (-t217 * t320 + t220 * t319) * qJD(5) + m(7) * (t1 * t220 + t2 * t217 + t274 * t7 - t276 * t6) + m(6) * (-t14 * t274 - t217 * t4 + t220 * t3 - t276 * t322)) * pkin(11) + t18 * t184 + t179 * t22 + t5 * t183 + t21 * t162 + t31 * t163 + t153 * t54 + t20 * mrSges(5,1) - t19 * mrSges(5,2) + t269 + t246 * t50 + t247 * t85 + t248 * t49 + t249 * t88 + t250 * t99 + t251 * t87 + ((mrSges(7,2) * t7 - mrSges(6,3) * t14 + t266) * t220 + (-mrSges(7,2) * t6 - mrSges(6,3) * t322 + t267) * t217) * qJD(5) + (t2 * mrSges(7,2) - t4 * mrSges(6,3) + t268) * t217 + (-t23 - t328) * pkin(4); -t131 * mrSges(5,2) + (t162 + t163) * t155 + (t183 + t303) * t132 + m(6) * (-pkin(4) * t132 + t93) + m(7) * (t132 * t179 + t153 * t155 + t93) + 0.2e1 * (t330 + m(7) / 0.2e1) * (t133 * t274 - t134 * t276 + t305) * pkin(11) + (mrSges(6,3) + mrSges(7,2)) * (t305 + t306 + (t133 * t220 - t134 * t217) * qJD(5)); m(7) * (t104 * t179 + t152 * t153) + t179 * t124 + t104 * t183 + t153 * t158 + t152 * t162 - pkin(4) * t125 + t218 * pkin(10) * t163 + t209 - t250 * t221 + ((-Ifges(5,6) + t247) * t218 + (t308 + (-m(6) * pkin(4) + t303) * t221) * pkin(10)) * qJD(4) + (t105 * mrSges(6,3) + t98 * mrSges(7,2) + t249 * t218 + t246 * t277 + (mrSges(7,2) * t137 - mrSges(6,3) * t143 + t218 * t248 + t252) * qJD(5) + (t284 * qJD(5) + m(6) * (-qJD(5) * t143 + t105) + m(7) * (qJD(5) * t137 + t98) + t288) * pkin(11) - t255) * t220 + (-t106 * mrSges(6,3) + t101 * mrSges(7,2) + t251 * t218 + t248 * t277 + (-t136 * mrSges(7,2) - t144 * mrSges(6,3) - t218 * t246 + t253) * qJD(5) + (-t285 * qJD(5) + m(6) * (-qJD(5) * t144 - t106) + m(7) * (-qJD(5) * t136 + t101) + t289) * pkin(11) + t254) * t217; -0.2e1 * pkin(4) * t163 + 0.2e1 * t162 * t179 + (-t165 + t168) * t220 + (t170 + t171) * t217 + 0.2e1 * (m(7) * t179 + t183) * t153 + ((t191 + t192) * t220 + (t186 - t189) * t217) * qJD(5); -pkin(5) * t26 + m(7) * (-pkin(5) * t2 + qJ(6) * t1 + qJD(6) * t6) + qJD(6) * t56 + qJ(6) * t27 + t1 * mrSges(7,3) + t4 * mrSges(6,1) - t2 * mrSges(7,1) - t3 * mrSges(6,2) + t338; m(7) * qJD(6) * t134 + (-m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1)) * t96 + (-mrSges(6,2) + t336) * t95; -pkin(5) * t139 + m(7) * (-pkin(5) * t101 + qJ(6) * t98 + qJD(6) * t136) + qJD(6) * t177 + qJ(6) * t141 + t98 * mrSges(7,3) - Ifges(6,6) * t257 - t101 * mrSges(7,1) - t105 * mrSges(6,2) + t106 * mrSges(6,1) + (-t312 - t343) * t275 + t245 + t283; -t335 * mrSges(7,2) + (m(7) * t272 + (-m(7) * t236 + t183 + t184) * qJD(5)) * pkin(11) + t166 + t167; 0.2e1 * t336 * qJD(6); m(7) * t2 + t26; m(7) * t96; m(7) * t101 + t139; (m(7) * pkin(11) + mrSges(7,2)) * t274; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t15(1) t15(2) t15(4) t15(7) t15(11) t15(16); t15(2) t15(3) t15(5) t15(8) t15(12) t15(17); t15(4) t15(5) t15(6) t15(9) t15(13) t15(18); t15(7) t15(8) t15(9) t15(10) t15(14) t15(19); t15(11) t15(12) t15(13) t15(14) t15(15) t15(20); t15(16) t15(17) t15(18) t15(19) t15(20) t15(21);];
Mq  = res;
