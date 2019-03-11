% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRPPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
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
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPPRR9_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR9_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR9_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR9_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR9_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:28:25
% EndTime: 2019-03-09 09:28:46
% DurationCPUTime: 9.50s
% Computational Cost: add. (6978->665), mult. (18187->883), div. (0->0), fcn. (12464->8), ass. (0->304)
t212 = sin(qJ(6));
t215 = cos(qJ(6));
t217 = cos(qJ(2));
t208 = sin(pkin(6));
t289 = qJD(1) * t208;
t274 = t217 * t289;
t179 = qJD(5) + t274;
t213 = sin(qJ(5));
t216 = cos(qJ(5));
t214 = sin(qJ(2));
t209 = cos(pkin(6));
t288 = qJD(1) * t209;
t280 = pkin(1) * t288;
t160 = pkin(8) * t274 + t214 * t280;
t131 = pkin(3) * t274 + t160;
t107 = -pkin(4) * t274 - t131;
t197 = qJD(2) + t288;
t183 = t197 * qJ(3);
t355 = qJD(4) + t183;
t69 = -t197 * pkin(9) - t107 + t355;
t210 = qJ(3) - pkin(9);
t211 = pkin(2) + qJ(4);
t260 = -t211 * t217 - pkin(1);
t108 = (-t210 * t214 + t260) * t208;
t84 = qJD(1) * t108;
t34 = t213 * t69 + t216 * t84;
t29 = pkin(10) * t179 + t34;
t275 = t214 * t289;
t143 = t197 * t213 - t216 * t275;
t144 = t197 * t216 + t213 * t275;
t195 = t217 * t280;
t337 = pkin(3) + pkin(8);
t282 = pkin(4) + t337;
t270 = t208 * t282;
t255 = t214 * t270;
t106 = -qJD(1) * t255 + t195;
t356 = -t197 * pkin(2) + qJD(3);
t352 = t197 * qJ(4) - t356;
t68 = t106 + t352;
t36 = pkin(5) * t143 - pkin(10) * t144 + t68;
t10 = -t212 * t29 + t215 * t36;
t11 = t212 * t36 + t215 * t29;
t245 = t10 * t215 + t11 * t212;
t252 = mrSges(7,1) * t212 + mrSges(7,2) * t215;
t33 = -t213 * t84 + t216 * t69;
t28 = -pkin(5) * t179 - t33;
t142 = qJD(6) + t143;
t94 = t144 * t215 + t179 * t212;
t328 = Ifges(7,4) * t94;
t93 = -t144 * t212 + t179 * t215;
t31 = Ifges(7,2) * t93 + Ifges(7,6) * t142 + t328;
t85 = Ifges(7,4) * t93;
t32 = Ifges(7,1) * t94 + Ifges(7,5) * t142 + t85;
t329 = t215 / 0.2e1;
t330 = -t212 / 0.2e1;
t369 = -t245 * mrSges(7,3) + t28 * t252 + t31 * t330 + t32 * t329;
t286 = qJD(2) * t217;
t272 = t213 * t286;
t283 = qJD(5) * t216;
t284 = qJD(5) * t213;
t101 = -t197 * t284 + (t214 * t283 + t272) * t289;
t368 = t101 / 0.2e1;
t281 = qJD(1) * qJD(2);
t271 = t208 * t281;
t265 = t217 * t271;
t102 = qJD(5) * t144 - t216 * t265;
t367 = -t102 / 0.2e1;
t366 = mrSges(4,1) + mrSges(3,3);
t365 = mrSges(4,2) - mrSges(3,1);
t308 = t179 * Ifges(6,5);
t140 = Ifges(6,4) * t143;
t310 = t144 * Ifges(6,1);
t65 = -t140 + t308 + t310;
t229 = t33 * mrSges(6,3) - t65 / 0.2e1 - t308 / 0.2e1 - t68 * mrSges(6,2);
t247 = Ifges(7,5) * t215 - Ifges(7,6) * t212;
t317 = Ifges(7,4) * t215;
t249 = -Ifges(7,2) * t212 + t317;
t318 = Ifges(7,4) * t212;
t251 = Ifges(7,1) * t215 - t318;
t333 = -t142 / 0.2e1;
t339 = -t94 / 0.2e1;
t341 = -t93 / 0.2e1;
t364 = t247 * t333 + t249 * t341 + t251 * t339 + t229 - t369;
t347 = -Ifges(6,5) / 0.2e1;
t363 = Ifges(6,6) / 0.2e1;
t334 = t102 / 0.2e1;
t362 = -t197 / 0.2e1;
t361 = t197 / 0.2e1;
t360 = -t289 / 0.2e1;
t359 = t289 / 0.2e1;
t358 = -t101 * Ifges(6,4) / 0.2e1;
t357 = -t310 / 0.2e1 + t140 / 0.2e1;
t180 = t197 * qJD(4);
t327 = pkin(1) * t209;
t201 = t214 * t327;
t225 = (-t217 * t270 - t201) * qJD(2);
t27 = t102 * pkin(5) - t101 * pkin(10) + qJD(1) * t225 + t180;
t285 = qJD(3) * t214;
t224 = (-t285 + (-qJD(2) * t210 - qJD(4)) * t217) * t208;
t266 = t214 * t271;
t177 = pkin(2) * t266;
t293 = qJ(4) * t266 + t177;
t61 = qJD(1) * t224 + t293;
t235 = qJD(2) * t255;
t178 = qJD(2) * t195;
t181 = t197 * qJD(3);
t292 = t178 + t181;
t72 = -qJD(1) * t235 + t292;
t12 = t213 * t72 + t216 * t61 + t69 * t283 - t284 * t84;
t5 = -pkin(10) * t266 + t12;
t1 = qJD(6) * t10 + t212 * t27 + t215 * t5;
t2 = -qJD(6) * t11 - t212 * t5 + t215 * t27;
t257 = t1 * t215 - t2 * t212;
t13 = -qJD(5) * t34 - t213 * t61 + t216 * t72;
t354 = t13 * mrSges(6,1) - t12 * mrSges(6,2) + Ifges(6,5) * t101 - Ifges(6,6) * t102;
t353 = Ifges(6,1) * t368 + Ifges(6,4) * t367;
t100 = t131 + t355;
t302 = qJ(3) * t214;
t120 = (t260 - t302) * t208;
t109 = qJD(1) * t120;
t129 = -t183 - t160;
t148 = (-pkin(2) * t217 - pkin(1) - t302) * t208;
t133 = qJD(1) * t148;
t187 = Ifges(5,6) * t275;
t279 = Ifges(5,4) / 0.2e1 + Ifges(4,5) / 0.2e1;
t350 = -(Ifges(3,6) / 0.2e1 - t279) * t197 + t109 * mrSges(5,3) + t129 * mrSges(4,1) + t34 * mrSges(6,2) + Ifges(3,6) * t362 + (Ifges(3,4) * t214 + Ifges(3,2) * t217) * t360 + (-Ifges(4,6) * t214 - Ifges(4,3) * t217) * t359 - Ifges(5,2) * t274 / 0.2e1 + t187 / 0.2e1 - t179 * Ifges(6,3) - t144 * Ifges(6,5) + t143 * Ifges(6,6) - t100 * mrSges(5,1) - t133 * mrSges(4,2) - t160 * mrSges(3,3) - t33 * mrSges(6,1) + (Ifges(4,5) + Ifges(5,4)) * t361;
t159 = pkin(8) * t275 - t195;
t119 = t159 + t356;
t188 = Ifges(3,4) * t274;
t278 = Ifges(3,5) / 0.2e1 + Ifges(5,5) / 0.2e1;
t79 = t275 * t337 - t195 - t352;
t349 = (-Ifges(4,4) / 0.2e1 + t278) * t197 + t119 * mrSges(4,1) + t159 * mrSges(3,3) + t79 * mrSges(5,1) + Ifges(3,1) * t275 / 0.2e1 + t188 / 0.2e1 + (-Ifges(5,6) * t217 + Ifges(5,3) * t214) * t359 + Ifges(4,4) * t362 + (-Ifges(4,2) * t214 - Ifges(4,6) * t217) * t360 - t109 * mrSges(5,2) - t133 * mrSges(4,3) + (Ifges(3,5) + Ifges(5,5)) * t361;
t299 = t208 * t217;
t166 = pkin(8) * t299 + t201;
t147 = -t209 * qJ(3) - t166;
t121 = pkin(3) * t299 - t147;
t92 = pkin(4) * t299 - t209 * pkin(9) + t121;
t320 = t216 * t108 + t213 * t92;
t287 = qJD(2) * t214;
t273 = t208 * t287;
t192 = pkin(2) * t273;
t291 = qJ(4) * t273 + t192;
t71 = t224 + t291;
t202 = t217 * t327;
t196 = qJD(2) * t202;
t204 = t209 * qJD(3);
t290 = t196 + t204;
t86 = -t235 + t290;
t18 = -qJD(5) * t320 - t213 * t71 + t216 * t86;
t39 = qJD(6) * t93 + t101 * t215 - t212 * t266;
t40 = -qJD(6) * t94 - t101 * t212 - t215 * t266;
t9 = t39 * Ifges(7,1) + t40 * Ifges(7,4) + t102 * Ifges(7,5);
t348 = t9 / 0.2e1;
t346 = -t31 / 0.2e1;
t345 = t39 / 0.2e1;
t344 = t40 / 0.2e1;
t343 = Ifges(6,2) * t334 + t266 * t363 + t358;
t342 = t266 * t347 + t353;
t340 = t93 / 0.2e1;
t338 = t94 / 0.2e1;
t336 = pkin(1) * mrSges(3,1);
t335 = pkin(1) * mrSges(3,2);
t332 = t142 / 0.2e1;
t300 = t208 * t216;
t232 = -t209 * t213 + t214 * t300;
t331 = -t232 / 0.2e1;
t38 = Ifges(7,5) * t39;
t37 = Ifges(7,6) * t40;
t324 = t93 * Ifges(7,6);
t323 = t94 * Ifges(7,5);
t322 = -mrSges(5,1) - mrSges(4,1);
t16 = -mrSges(7,1) * t40 + mrSges(7,2) * t39;
t75 = -mrSges(6,1) * t266 - mrSges(6,3) * t101;
t321 = -t16 + t75;
t184 = qJ(4) * t275;
t190 = pkin(2) * t275;
t103 = -t210 * t274 + t184 + t190;
t51 = t216 * t103 + t213 * t106;
t319 = Ifges(6,4) * t144;
t313 = t142 * Ifges(7,3);
t312 = t143 * Ifges(6,2);
t307 = t179 * Ifges(6,6);
t305 = t213 * t28;
t152 = mrSges(5,1) * t275 - mrSges(5,3) * t197;
t77 = mrSges(6,1) * t143 + mrSges(6,2) * t144;
t304 = t152 - t77;
t105 = mrSges(6,1) * t179 - mrSges(6,3) * t144;
t47 = -mrSges(7,1) * t93 + mrSges(7,2) * t94;
t303 = -t47 + t105;
t301 = t208 * t214;
t298 = t210 * t213;
t297 = t212 * t217;
t296 = t215 * t217;
t295 = -t197 * t365 - t275 * t366;
t153 = -mrSges(4,1) * t274 - t197 * mrSges(4,3);
t154 = mrSges(5,1) * t274 + t197 * mrSges(5,2);
t294 = t153 - t154;
t165 = -pkin(8) * t301 + t202;
t7 = Ifges(7,3) * t102 + t37 + t38;
t149 = -t209 * pkin(2) - t165;
t269 = t337 * t301;
t268 = t209 * qJ(4) - t149;
t267 = t216 * t274;
t262 = 0.3e1 / 0.2e1 * Ifges(5,6) - 0.3e1 / 0.2e1 * Ifges(4,6) - 0.3e1 / 0.2e1 * Ifges(3,4);
t261 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t258 = pkin(5) * t216 + pkin(10) * t213;
t256 = -t1 * t212 - t2 * t215;
t254 = qJD(5) * t258 - qJD(6) * t298 + qJD(4);
t52 = t102 * mrSges(6,1) + t101 * mrSges(6,2);
t253 = mrSges(7,1) * t215 - mrSges(7,2) * t212;
t250 = Ifges(7,1) * t212 + t317;
t248 = Ifges(7,2) * t215 + t318;
t246 = Ifges(7,5) * t212 + Ifges(7,6) * t215;
t244 = t10 * t212 - t11 * t215;
t25 = mrSges(7,1) * t102 - mrSges(7,3) * t39;
t26 = -mrSges(7,2) * t102 + mrSges(7,3) * t40;
t243 = -t212 * t25 + t215 * t26;
t42 = pkin(10) * t299 + t320;
t164 = t209 * t216 + t213 * t301;
t91 = (-pkin(3) - pkin(4)) * t301 + t268;
t55 = -pkin(5) * t232 - pkin(10) * t164 + t91;
t20 = t212 * t55 + t215 * t42;
t19 = -t212 * t42 + t215 * t55;
t56 = -mrSges(7,2) * t142 + mrSges(7,3) * t93;
t57 = mrSges(7,1) * t142 - mrSges(7,3) * t94;
t242 = t212 * t57 - t215 * t56;
t241 = -t212 * t56 - t215 * t57;
t240 = t213 * t33 - t216 * t34;
t48 = -t213 * t108 + t216 * t92;
t236 = qJD(2) * t269;
t50 = -t103 * t213 + t106 * t216;
t161 = -pkin(8) * t273 + t196;
t157 = -qJ(3) * t274 + t190;
t104 = -mrSges(6,2) * t179 - mrSges(6,3) * t143;
t234 = -t104 + t242;
t117 = -t164 * t212 + t208 * t296;
t118 = t164 * t215 + t208 * t297;
t17 = -t108 * t284 + t213 * t86 + t216 * t71 + t92 * t283;
t145 = -pkin(8) * t266 + t178;
t230 = (-qJ(3) * t286 - t285) * t208;
t162 = t166 * qJD(2);
t169 = pkin(5) * t213 - pkin(10) * t216 + t211;
t228 = qJD(3) * t213 + qJD(6) * t169 + t210 * t283;
t227 = (t299 * t337 + t201) * qJD(2);
t226 = (-t285 + (-qJ(3) * qJD(2) - qJD(4)) * t217) * t208;
t203 = t209 * qJD(4);
t87 = t203 + t225;
t114 = -t145 - t181;
t146 = qJD(1) * t162;
t88 = -qJD(1) * t236 + t292;
t89 = qJD(1) * t227 - t180;
t223 = -t145 * mrSges(3,2) + t88 * mrSges(5,2) - t114 * mrSges(4,3) - t89 * mrSges(5,3) + t146 * t365;
t30 = t313 + t323 + t324;
t64 = t307 - t312 + t319;
t221 = -t324 / 0.2e1 - t323 / 0.2e1 - t313 / 0.2e1 + t307 / 0.2e1 - t10 * mrSges(7,1) + t11 * mrSges(7,2) - t30 / 0.2e1 + t64 / 0.2e1 + t34 * mrSges(6,3) - t68 * mrSges(6,1) + t319 / 0.2e1;
t220 = (t312 / 0.2e1 - t221) * t216;
t176 = Ifges(5,4) * t266;
t175 = Ifges(3,5) * t265;
t174 = Ifges(4,5) * t266;
t173 = Ifges(5,5) * t265;
t158 = (-mrSges(5,2) * t214 - mrSges(5,3) * t217) * t289;
t156 = (mrSges(4,2) * t217 - mrSges(4,3) * t214) * t289;
t151 = -t197 * mrSges(3,2) + mrSges(3,3) * t274;
t141 = -t161 - t204;
t139 = t169 * t212 + t215 * t298;
t138 = t169 * t215 - t212 * t298;
t137 = (-t212 * t214 + t213 * t296) * t289;
t136 = (-t213 * t297 - t214 * t215) * t289;
t135 = t212 * t197 - t215 * t267;
t134 = t215 * t197 + t212 * t267;
t132 = t192 + t230;
t130 = -qJD(1) * t269 + t195;
t122 = t157 + t184;
t116 = qJD(5) * t164 - t286 * t300;
t115 = qJD(5) * t232 + t208 * t272;
t113 = qJD(1) * t230 + t177;
t112 = -t203 + t227;
t111 = -t236 + t290;
t110 = pkin(3) * t301 - t268;
t82 = t226 + t291;
t78 = pkin(5) * t144 + pkin(10) * t143;
t76 = mrSges(6,2) * t266 - mrSges(6,3) * t102;
t74 = qJD(1) * t226 + t293;
t73 = -t180 + (t282 * t299 + t201) * t281;
t66 = (-pkin(4) - t258) * t274 - t131;
t60 = -t212 * t228 + t215 * t254;
t59 = t212 * t254 + t215 * t228;
t54 = qJD(6) * t117 + t115 * t215 - t212 * t273;
t53 = -qJD(6) * t118 - t115 * t212 - t215 * t273;
t44 = -pkin(10) * t275 + t51;
t43 = pkin(5) * t275 - t50;
t41 = -pkin(5) * t299 - t48;
t35 = t116 * pkin(5) - t115 * pkin(10) + t87;
t24 = t212 * t66 + t215 * t44;
t23 = -t212 * t44 + t215 * t66;
t22 = t212 * t78 + t215 * t33;
t21 = -t212 * t33 + t215 * t78;
t15 = pkin(5) * t273 - t18;
t14 = -pkin(10) * t273 + t17;
t8 = t39 * Ifges(7,4) + t40 * Ifges(7,2) + t102 * Ifges(7,6);
t6 = pkin(5) * t266 - t13;
t4 = -qJD(6) * t20 - t14 * t212 + t215 * t35;
t3 = qJD(6) * t19 + t14 * t215 + t212 * t35;
t45 = [-t295 * t162 - t143 * (Ifges(6,4) * t115 - Ifges(6,2) * t116) / 0.2e1 + t118 * t348 + (Ifges(7,1) * t54 + Ifges(7,4) * t53 + Ifges(7,5) * t116) * t338 + (Ifges(7,4) * t54 + Ifges(7,2) * t53 + Ifges(7,6) * t116) * t340 + t164 * t342 + m(6) * (t12 * t320 + t13 * t48 + t17 * t34 + t18 * t33 + t68 * t87 - t73 * t91) + t320 * t76 + t7 * t331 + (Ifges(7,5) * t54 + Ifges(7,6) * t53 + Ifges(7,3) * t116) * t332 + (t174 / 0.2e1 + t175 / 0.2e1 + t176 / 0.2e1 + t173 / 0.2e1 + t223) * t209 + t179 * (Ifges(6,5) * t115 - Ifges(6,6) * t116) / 0.2e1 + t144 * (Ifges(6,1) * t115 - Ifges(6,4) * t116) / 0.2e1 + ((t149 * mrSges(4,1) + t110 * mrSges(5,1) - t120 * mrSges(5,2) - t165 * mrSges(3,3) - t148 * mrSges(4,3) + (-t217 * t262 - 0.2e1 * t335) * t208 + (-Ifges(4,4) + t278) * t209) * t217 + (-t148 * mrSges(4,2) + t120 * mrSges(5,3) + t164 * t347 + Ifges(6,6) * t331 + t147 * mrSges(4,1) - t166 * mrSges(3,3) - t121 * mrSges(5,1) + (t214 * t262 - 0.2e1 * t336) * t208 + (-Ifges(3,6) + t279) * t209 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(5,2) + 0.3e1 / 0.2e1 * Ifges(5,3) + 0.3e1 / 0.2e1 * Ifges(4,2) - 0.3e1 / 0.2e1 * Ifges(4,3) - 0.3e1 / 0.2e1 * Ifges(3,2) - Ifges(6,3)) * t299) * t214) * t271 + (Ifges(6,4) * t164 + Ifges(6,2) * t232) * t367 + (Ifges(6,1) * t164 + Ifges(6,4) * t232) * t368 + ((t89 * mrSges(5,1) - t74 * mrSges(5,2) - t113 * mrSges(4,3) + t146 * t366) * t214 + (-t114 * mrSges(4,1) + t88 * mrSges(5,1) + t113 * mrSges(4,2) + t145 * mrSges(3,3) - t74 * mrSges(5,3) + t354) * t217 + (t214 * t350 + t217 * t349) * qJD(2)) * t208 + (Ifges(7,5) * t118 + Ifges(7,6) * t117 - Ifges(7,3) * t232) * t334 + (Ifges(7,4) * t118 + Ifges(7,2) * t117 - Ifges(7,6) * t232) * t344 + (Ifges(7,1) * t118 + Ifges(7,4) * t117 - Ifges(7,5) * t232) * t345 + (-t115 * t33 - t116 * t34 + t12 * t232 - t13 * t164) * mrSges(6,3) + t2 * (-mrSges(7,1) * t232 - mrSges(7,3) * t118) + t1 * (mrSges(7,2) * t232 + mrSges(7,3) * t117) - t73 * (-mrSges(6,1) * t232 + mrSges(6,2) * t164) - t232 * t343 + t161 * t151 + t132 * t156 + t82 * t158 + t112 * t152 + t141 * t153 + t111 * t154 + t6 * (-mrSges(7,1) * t117 + mrSges(7,2) * t118) + t10 * (mrSges(7,1) * t116 - mrSges(7,3) * t54) - t116 * t64 / 0.2e1 + t68 * (mrSges(6,1) * t116 + mrSges(6,2) * t115) + t117 * t8 / 0.2e1 + m(7) * (t1 * t20 + t10 * t4 + t11 * t3 + t15 * t28 + t19 * t2 + t41 * t6) + m(5) * (t100 * t111 + t109 * t82 + t110 * t89 + t112 * t79 + t120 * t74 + t121 * t88) + m(4) * (t113 * t148 + t114 * t147 + t119 * t162 + t129 * t141 + t132 * t133 + t146 * t149) + m(3) * (t145 * t166 - t146 * t165 + t159 * t162 + t160 * t161) + t19 * t25 + t20 * t26 + t41 * t16 + t15 * t47 + t53 * t31 / 0.2e1 + t54 * t32 / 0.2e1 + t28 * (-mrSges(7,1) * t53 + mrSges(7,2) * t54) + t3 * t56 + t4 * t57 + t48 * t75 + t87 * t77 + t91 * t52 + t17 * t104 + t18 * t105 + t115 * t65 / 0.2e1 + t116 * t30 / 0.2e1 + t11 * (-mrSges(7,2) * t116 + mrSges(7,3) * t53); t223 - t294 * qJD(3) + t295 * t160 + (Ifges(7,1) * t137 + Ifges(7,4) * t136) * t339 + (Ifges(7,4) * t137 + Ifges(7,2) * t136) * t341 + t136 * t346 + (m(6) * (qJD(3) * t34 + t12 * t210) + t38 / 0.2e1 + t37 / 0.2e1 - t73 * mrSges(6,1) + t358 + t7 / 0.2e1 + t343 - t12 * mrSges(6,3) + qJD(3) * t104 + t210 * t76 + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t102 + t261) * t213 + (((t322 * qJ(3) + t213 * t363 + t216 * t347 - Ifges(3,6)) * qJD(2) - t187 / 0.2e1 + (t336 + (Ifges(4,6) / 0.2e1 + Ifges(3,4) / 0.2e1) * t214) * t289 - t350) * t214 + (-t188 / 0.2e1 + (t335 + (Ifges(5,6) / 0.2e1 - Ifges(4,6) / 0.2e1) * t217) * t289 + (-pkin(2) * mrSges(4,1) - t211 * mrSges(5,1) - Ifges(4,4)) * qJD(2) + (Ifges(3,2) / 0.2e1 - Ifges(5,3) / 0.2e1 + Ifges(5,2) / 0.2e1 - Ifges(4,2) / 0.2e1 + Ifges(4,3) / 0.2e1 - Ifges(3,1) / 0.2e1) * t275 + (t229 + t357) * t213 + t220 - t349) * t217) * t289 + (Ifges(7,5) * t137 + Ifges(7,6) * t136) * t333 + (t10 * t137 - t11 * t136) * mrSges(7,3) + (qJ(3) * t88 - t109 * t122 - t211 * t89 + (-qJD(4) - t131) * t79 + (qJD(3) - t130) * t100) * m(5) - m(7) * (t10 * t23 + t11 * t24 + t28 * t43) - m(6) * (t107 * t68 + t33 * t50 + t34 * t51) - t304 * qJD(4) + (-pkin(2) * t146 - qJ(3) * t114 - t119 * t160 - t133 * t157 + (-qJD(3) - t159) * t129) * m(4) + t176 + t173 + t174 + t175 + (-t23 + t60) * t57 + (t342 + t251 * t345 + t249 * t344 + t247 * t334 + t6 * t252 + t9 * t329 + t8 * t330 - t13 * mrSges(6,3) - t73 * mrSges(6,2) + t321 * t210 + t303 * qJD(3) + t256 * mrSges(7,3) + m(7) * (-qJD(3) * t28 - t210 * t6) + m(6) * (qJD(3) * t33 + t13 * t210) + (mrSges(7,3) * t244 + t215 * t346 + t246 * t333 + t248 * t341 + t250 * t339 + t253 * t28 + t32 * t330) * qJD(6) + t353) * t216 + m(7) * (t1 * t139 + t10 * t60 + t11 * t59 + t138 * t2) + (-t24 + t59) * t56 + (t220 + (t357 + t364) * t213 + (-m(6) * t240 + m(7) * t305 + t216 * t104 - t303 * t213) * t210) * qJD(5) + t211 * t52 - t157 * t156 - t122 * t158 - t131 * t152 - t130 * t154 - t28 * (-mrSges(7,1) * t136 + mrSges(7,2) * t137) - t137 * t32 / 0.2e1 + t138 * t25 + t139 * t26 + (-t153 + t151) * t159 + m(6) * (qJD(4) * t68 - t211 * t73) - t43 * t47 - t51 * t104 - t50 * t105 - t107 * t77; -t212 * t26 - t215 * t25 + t294 * t197 - t303 * t144 + t242 * qJD(6) + t234 * t143 + ((t156 + t158) * t214 - t322 * t286) * t289 - t52 + (t142 * t244 + t144 * t28 + t256) * m(7) + (-t143 * t34 - t144 * t33 + t73) * m(6) + (-t100 * t197 + t109 * t275 + t89) * m(5) + (t129 * t197 + t133 * t275 + t146) * m(4); -t134 * t57 - t135 * t56 + t304 * t197 + (-mrSges(5,1) * t287 + t158 * t217) * t289 + (-qJD(5) * t234 + t104 * t274 + t321) * t216 + (t241 * qJD(6) - t179 * t303 + t243 + t76) * t213 + (-t10 * t134 - t11 * t135 + t274 * t305 + (-qJD(5) * t244 - t6) * t216 + (qJD(5) * t28 - qJD(6) * t245 + t257) * t213) * m(7) + (t12 * t213 + t13 * t216 - t179 * t240 - t68 * t197) * m(6) + (t109 * t274 + t79 * t197 + t88) * m(5); ((-m(7) * t245 + t241) * qJD(6) + m(7) * t257 + t243) * pkin(10) + (-t140 / 0.2e1 + (Ifges(6,1) / 0.2e1 - Ifges(6,2) / 0.2e1) * t144 - t364) * t143 + t246 * t334 + t248 * t344 + t250 * t345 + t212 * t348 + (-pkin(5) * t6 - t10 * t21 - t11 * t22 - t28 * t34) * m(7) + t221 * t144 + t8 * t329 - t6 * t253 + t303 * t34 - Ifges(6,3) * t266 + (t247 * t332 + t249 * t340 + t251 * t338 + t369) * qJD(6) + t257 * mrSges(7,3) + t354 - pkin(5) * t16 - t22 * t56 - t21 * t57 - t33 * t104; -t28 * (mrSges(7,1) * t94 + mrSges(7,2) * t93) + (Ifges(7,1) * t93 - t328) * t339 + t31 * t338 + (Ifges(7,5) * t93 - Ifges(7,6) * t94) * t333 - t10 * t56 + t11 * t57 + (t10 * t93 + t11 * t94) * mrSges(7,3) + t261 + t7 + (-Ifges(7,2) * t94 + t32 + t85) * t341;];
tauc  = t45(:);
