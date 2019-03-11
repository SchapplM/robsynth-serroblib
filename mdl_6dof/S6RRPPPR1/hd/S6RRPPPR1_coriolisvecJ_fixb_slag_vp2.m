% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta4]';
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
% Datum: 2019-03-09 08:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPPPR1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR1_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR1_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR1_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR1_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:06:10
% EndTime: 2019-03-09 08:06:32
% DurationCPUTime: 13.39s
% Computational Cost: add. (6506->587), mult. (17269->779), div. (0->0), fcn. (12401->8), ass. (0->269)
t371 = Ifges(5,1) + Ifges(6,1);
t210 = sin(pkin(9));
t213 = sin(qJ(2));
t290 = cos(pkin(9));
t309 = cos(qJ(2));
t192 = t210 * t309 + t213 * t290;
t175 = t192 * qJD(1);
t209 = sin(pkin(10));
t211 = cos(pkin(10));
t140 = -t211 * qJD(2) + t175 * t209;
t325 = t140 / 0.2e1;
t242 = t290 * t309;
t275 = qJD(1) * t213;
t173 = -qJD(1) * t242 + t210 * t275;
t317 = t173 / 0.2e1;
t227 = qJD(2) * t209 + t211 * t175;
t323 = t227 / 0.2e1;
t363 = Ifges(5,5) + Ifges(6,4);
t362 = Ifges(6,2) + Ifges(5,3);
t370 = Ifges(5,6) - Ifges(6,6);
t299 = Ifges(6,5) * t209;
t301 = Ifges(5,4) * t209;
t345 = t211 * t371 + t299 - t301;
t174 = t192 * qJD(2);
t161 = qJD(1) * t174;
t322 = -t161 / 0.2e1;
t321 = t161 / 0.2e1;
t217 = -t210 * t213 + t242;
t176 = t217 * qJD(2);
t162 = qJD(1) * t176;
t369 = -t162 / 0.2e1;
t368 = t371 * t323 + t363 * t317 + (-Ifges(5,4) + Ifges(6,5)) * t325;
t212 = sin(qJ(6));
t214 = cos(qJ(6));
t228 = t212 * t209 + t211 * t214;
t108 = t228 * t173;
t177 = t228 * qJD(6);
t348 = t177 - t108;
t193 = t209 * t214 - t212 * t211;
t107 = t193 * t173;
t178 = t193 * qJD(6);
t347 = t178 - t107;
t327 = pkin(4) + pkin(5);
t263 = t327 * t209;
t289 = qJ(5) * t211;
t367 = t263 - t289;
t366 = t212 * t140 + t214 * t227;
t82 = t140 * t214 - t212 * t227;
t170 = qJD(6) - t173;
t308 = Ifges(7,4) * t366;
t24 = t82 * Ifges(7,2) + t170 * Ifges(7,6) + t308;
t335 = t24 / 0.2e1;
t81 = Ifges(7,4) * t82;
t25 = Ifges(7,1) * t366 + t170 * Ifges(7,5) + t81;
t334 = t25 / 0.2e1;
t42 = qJD(6) * t82 + t162 * t228;
t333 = t42 / 0.2e1;
t43 = -qJD(6) * t366 + t162 * t193;
t332 = t43 / 0.2e1;
t307 = pkin(2) * t210;
t203 = qJ(4) + t307;
t306 = -pkin(8) + t203;
t185 = t306 * t209;
t186 = t306 * t211;
t126 = t212 * t185 + t186 * t214;
t266 = pkin(2) * t275;
t115 = pkin(3) * t175 + qJ(4) * t173 + t266;
t267 = t309 * pkin(7);
t201 = qJ(3) * t309 + t267;
t196 = t201 * qJD(1);
t179 = t210 * t196;
t200 = (-qJ(3) - pkin(7)) * t213;
t195 = qJD(1) * t200;
t134 = t195 * t290 - t179;
t128 = t209 * t134;
t26 = t128 + (pkin(8) * t173 - t115) * t211 - t327 * t175;
t286 = t173 * t209;
t60 = t209 * t115 + t211 * t134;
t49 = t175 * qJ(5) + t60;
t36 = -pkin(8) * t286 + t49;
t365 = qJD(4) * t193 - qJD(6) * t126 + t212 * t36 - t214 * t26;
t125 = t185 * t214 - t212 * t186;
t364 = qJD(4) * t228 + qJD(6) * t125 - t212 * t26 - t214 * t36;
t361 = -t140 * t370 + t173 * t362 + t227 * t363;
t359 = t161 * t363 + t162 * t345;
t37 = -mrSges(7,1) * t82 + mrSges(7,2) * t366;
t91 = mrSges(6,1) * t140 - mrSges(6,3) * t227;
t358 = -t91 + t37;
t294 = t175 * Ifges(4,4);
t353 = Ifges(4,6) * qJD(2);
t356 = Ifges(7,5) * t366 - t173 * Ifges(4,2) + t82 * Ifges(7,6) + t170 * Ifges(7,3) + t294 + t353;
t295 = t175 * mrSges(4,3);
t355 = qJD(2) * mrSges(4,1) - mrSges(5,1) * t140 - mrSges(5,2) * t227 - t295;
t354 = Ifges(4,5) * qJD(2);
t288 = t162 * t209;
t109 = -mrSges(5,2) * t161 - mrSges(5,3) * t288;
t112 = -mrSges(6,2) * t288 + mrSges(6,3) * t161;
t350 = t109 + t112;
t287 = t162 * t211;
t110 = mrSges(5,1) * t161 - mrSges(5,3) * t287;
t111 = -t161 * mrSges(6,1) + mrSges(6,2) * t287;
t349 = t110 - t111;
t346 = -t209 * t370 + t211 * t363;
t253 = qJD(1) * t309;
t344 = t213 * pkin(7) * (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t253) + (qJD(2) * mrSges(3,1) - mrSges(3,3) * t275) * t267;
t270 = qJD(1) * qJD(2);
t251 = t213 * t270;
t246 = pkin(2) * t251;
t78 = pkin(3) * t161 - qJ(4) * t162 - qJD(4) * t175 + t246;
t171 = qJD(2) * t200 + qJD(3) * t309;
t154 = t171 * qJD(1);
t172 = -qJD(2) * t201 - t213 * qJD(3);
t215 = qJD(1) * t172;
t95 = t290 * t154 + t210 * t215;
t93 = qJD(2) * qJD(4) + t95;
t86 = t209 * t93;
t30 = t211 * t78 - t86;
t31 = t209 * t78 + t211 * t93;
t343 = -t209 * t30 + t211 * t31;
t298 = Ifges(6,5) * t211;
t232 = t209 * Ifges(6,3) + t298;
t300 = Ifges(5,4) * t211;
t235 = -t209 * Ifges(5,2) + t300;
t342 = Ifges(5,6) * t322 + t235 * t369 + Ifges(6,6) * t321 + t162 * t232 / 0.2e1;
t14 = t86 + (-pkin(8) * t162 - t78) * t211 - t327 * t161;
t16 = t161 * qJ(5) + t173 * qJD(5) + t31;
t15 = pkin(8) * t288 + t16;
t188 = qJD(2) * pkin(2) + t195;
t249 = t290 * t196;
t132 = t210 * t188 + t249;
t124 = qJD(2) * qJ(4) + t132;
t206 = -pkin(2) * t309 - pkin(1);
t276 = qJD(1) * t206;
t197 = qJD(3) + t276;
t99 = t173 * pkin(3) - t175 * qJ(4) + t197;
t53 = -t209 * t124 + t211 * t99;
t240 = qJD(5) - t53;
t19 = -pkin(8) * t227 - t173 * t327 + t240;
t54 = t211 * t124 + t209 * t99;
t44 = t173 * qJ(5) + t54;
t27 = pkin(8) * t140 + t44;
t5 = t19 * t214 - t212 * t27;
t1 = qJD(6) * t5 + t212 * t14 + t15 * t214;
t6 = t212 * t19 + t214 * t27;
t2 = -qJD(6) * t6 + t14 * t214 - t212 * t15;
t341 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t131 = t290 * t188 - t179;
t226 = qJD(2) * pkin(3) - qJD(4) + t131;
t238 = mrSges(6,1) * t209 - mrSges(6,3) * t211;
t239 = mrSges(5,1) * t209 + mrSges(5,2) * t211;
t310 = t209 / 0.2e1;
t216 = qJ(5) * t227 + t226;
t52 = pkin(4) * t140 - t216;
t340 = -t209 * (Ifges(5,4) * t227 - t140 * Ifges(5,2) + Ifges(5,6) * t173) / 0.2e1 + (Ifges(6,5) * t227 + Ifges(6,6) * t173 + t140 * Ifges(6,3)) * t310 + t197 * mrSges(4,2) + t52 * t238 + t354 / 0.2e1 - t226 * t239;
t41 = -pkin(4) * t173 + t240;
t339 = t41 * mrSges(6,1) + t5 * mrSges(7,1) + t54 * mrSges(5,2) + t353 / 0.2e1 - t197 * mrSges(4,1) - t44 * mrSges(6,3) - t53 * mrSges(5,1) - t6 * mrSges(7,2);
t337 = Ifges(7,4) * t333 + Ifges(7,2) * t332 + Ifges(7,6) * t322;
t336 = Ifges(7,1) * t333 + Ifges(7,4) * t332 + Ifges(7,5) * t322;
t331 = -t82 / 0.2e1;
t330 = t82 / 0.2e1;
t329 = -t366 / 0.2e1;
t328 = t366 / 0.2e1;
t326 = -t140 / 0.2e1;
t324 = -t227 / 0.2e1;
t320 = -t170 / 0.2e1;
t319 = t170 / 0.2e1;
t318 = -t173 / 0.2e1;
t314 = -t175 / 0.2e1;
t313 = t175 / 0.2e1;
t304 = mrSges(4,3) * t161;
t303 = mrSges(4,3) * t162;
t302 = Ifges(3,4) * t213;
t138 = -t290 * t200 + t201 * t210;
t94 = t154 * t210 - t290 * t215;
t297 = t138 * t94;
t296 = t173 * mrSges(4,3);
t117 = t171 * t290 + t210 * t172;
t274 = qJD(2) * t213;
t265 = pkin(2) * t274;
t88 = pkin(3) * t174 - qJ(4) * t176 - qJD(4) * t192 + t265;
t46 = t211 * t117 + t209 * t88;
t285 = t173 * t211;
t284 = t176 * t209;
t283 = t176 * t211;
t282 = t192 * t209;
t281 = t192 * t211;
t100 = -mrSges(6,2) * t140 + mrSges(6,3) * t173;
t101 = -mrSges(5,2) * t173 - mrSges(5,3) * t140;
t278 = t100 + t101;
t102 = mrSges(5,1) * t173 - mrSges(5,3) * t227;
t103 = -mrSges(6,1) * t173 + mrSges(6,2) * t227;
t277 = t102 - t103;
t130 = -pkin(3) * t217 - t192 * qJ(4) + t206;
t139 = t210 * t200 + t201 * t290;
t72 = t209 * t130 + t211 * t139;
t98 = mrSges(5,1) * t288 + mrSges(5,2) * t287;
t273 = qJD(4) * t209;
t272 = qJD(4) * t211;
t271 = qJD(5) * t209;
t268 = Ifges(7,5) * t42 + Ifges(7,6) * t43 - Ifges(7,3) * t161;
t264 = Ifges(3,4) * t309;
t57 = -qJ(5) * t217 + t72;
t260 = t290 * pkin(2);
t259 = -t288 / 0.2e1;
t258 = t288 / 0.2e1;
t257 = t287 / 0.2e1;
t252 = qJD(2) * t309;
t248 = t161 * mrSges(4,1) + t162 * mrSges(4,2);
t105 = t209 * t117;
t45 = t211 * t88 - t105;
t59 = t115 * t211 - t128;
t135 = t209 * t139;
t71 = t130 * t211 - t135;
t116 = t171 * t210 - t290 * t172;
t133 = t195 * t210 + t249;
t22 = t174 * qJ(5) - qJD(5) * t217 + t46;
t245 = qJD(1) * t252;
t205 = -t260 - pkin(3);
t97 = mrSges(6,1) * t288 - mrSges(6,3) * t287;
t13 = -t43 * mrSges(7,1) + t42 * mrSges(7,2);
t231 = pkin(4) * t209 - t289;
t230 = t209 * t53 - t211 * t54;
t38 = t135 + (-pkin(8) * t192 - t130) * t211 + t327 * t217;
t47 = pkin(8) * t282 + t57;
t11 = -t212 * t47 + t214 * t38;
t12 = t212 * t38 + t214 * t47;
t62 = -mrSges(7,2) * t170 + mrSges(7,3) * t82;
t63 = mrSges(7,1) * t170 - mrSges(7,3) * t366;
t229 = -t212 * t63 + t214 * t62;
t225 = -qJ(5) * t287 - qJD(5) * t227 + t94;
t224 = Ifges(3,5) * t309 - Ifges(3,6) * t213;
t119 = t193 * t192;
t222 = -t209 * qJ(5) + t205;
t221 = pkin(1) * (mrSges(3,1) * t213 + mrSges(3,2) * t309);
t220 = -qJD(5) * t281 + t116;
t219 = t213 * (Ifges(3,1) * t309 - t302);
t218 = (Ifges(3,2) * t309 + t302) * qJD(1);
t207 = Ifges(3,4) * t253;
t184 = Ifges(3,1) * t275 + Ifges(3,5) * qJD(2) + t207;
t183 = Ifges(3,6) * qJD(2) + t218;
t182 = -t211 * pkin(4) + t222;
t169 = Ifges(4,4) * t173;
t159 = t211 * t327 - t222;
t152 = -qJD(2) * mrSges(4,2) - t296;
t127 = mrSges(4,1) * t173 + mrSges(4,2) * t175;
t122 = t175 * Ifges(4,1) - t169 + t354;
t120 = t228 * t192;
t79 = t192 * t231 + t138;
t64 = -t173 * t231 + t133;
t61 = -t192 * t367 - t138;
t58 = pkin(4) * t217 - t71;
t56 = t176 * t193 - t177 * t192;
t55 = qJD(6) * t119 + t176 * t228;
t51 = t173 * t367 - t133;
t50 = -pkin(4) * t175 - t59;
t48 = t176 * t231 + t220;
t35 = -t140 * t327 + t216;
t34 = -t176 * t367 - t220;
t33 = -pkin(4) * t174 - t45;
t32 = pkin(4) * t288 + t225;
t29 = mrSges(7,2) * t161 + mrSges(7,3) * t43;
t28 = -mrSges(7,1) * t161 - mrSges(7,3) * t42;
t21 = t162 * t263 + t225;
t20 = -pkin(4) * t161 - t30;
t18 = pkin(8) * t284 + t22;
t17 = t105 + (-pkin(8) * t176 - t88) * t211 - t327 * t174;
t4 = -qJD(6) * t12 + t17 * t214 - t212 * t18;
t3 = qJD(6) * t11 + t212 * t17 + t18 * t214;
t7 = [(t235 * t326 + t232 * t325 + Ifges(4,1) * t313 + Ifges(4,4) * t318 + t122 / 0.2e1 - t131 * mrSges(4,3) + t345 * t323 + t346 * t317 + t340) * t176 + (Ifges(7,1) * t55 + Ifges(7,4) * t56) * t328 + (Ifges(7,1) * t120 + Ifges(7,4) * t119) * t333 + (t224 * qJD(2) / 0.2e1 - t344) * qJD(2) + t55 * t334 + t56 * t335 + t120 * t336 + t119 * t337 + m(7) * (t1 * t12 + t11 * t2 - t21 * t61 + t3 * t6 + t34 * t35 + t4 * t5) + m(6) * (t16 * t57 + t20 * t58 + t22 * t44 + t32 * t79 + t33 * t41 + t48 * t52) + (-t356 / 0.2e1 + t361 / 0.2e1 - t132 * mrSges(4,3) - Ifges(4,4) * t313 - Ifges(7,5) * t328 - Ifges(4,2) * t318 + Ifges(5,6) * t326 + Ifges(6,6) * t325 - Ifges(7,6) * t330 - Ifges(7,3) * t319 + t362 * t317 + t363 * t323 - t339) * t174 + (Ifges(7,4) * t55 + Ifges(7,2) * t56) * t330 + (Ifges(7,4) * t120 + Ifges(7,2) * t119) * t332 + t342 * t282 + t127 * t265 + ((t239 + mrSges(4,3)) * t94 + mrSges(4,2) * t246 + Ifges(4,1) * t162 + 0.2e1 * Ifges(4,4) * t322 + t232 * t258 + t235 * t259 + t32 * t238 + t345 * t257 + t346 * t321) * t192 + (Ifges(7,5) * t55 + Ifges(7,6) * t56) * t319 + (Ifges(7,5) * t120 + Ifges(7,6) * t119) * t322 + (-m(4) * t131 - m(5) * t226 - t355) * t116 + (t303 + t98) * t138 - t139 * t304 + (-Ifges(3,2) * t213 + t264) * t245 + (t1 * t119 - t120 * t2 - t5 * t55 + t56 * t6) * mrSges(7,3) + t206 * t248 + m(5) * (t30 * t71 + t31 * t72 + t45 * t53 + t46 * t54 + t297) + m(4) * (t117 * t132 + t297 + t139 * t95 + (t197 + t276) * t265) + t359 * t281 / 0.2e1 + (-t16 * t282 + t20 * t281 + t283 * t41 - t284 * t44) * mrSges(6,2) + t283 * t368 + (t268 / 0.2e1 + Ifges(4,4) * t162 + t346 * t369 + Ifges(7,6) * t332 + Ifges(7,5) * t333 - Ifges(6,6) * t258 - Ifges(5,6) * t259 - t30 * mrSges(5,1) + t20 * mrSges(6,1) + t31 * mrSges(5,2) - t16 * mrSges(6,3) + t95 * mrSges(4,3) - mrSges(4,1) * t246 - (-Ifges(7,3) - Ifges(4,2)) * t322 - t362 * t321 - t363 * t257 + t341 + (-Ifges(4,2) / 0.2e1 - t362 / 0.2e1) * t161) * t217 + (-0.2e1 * t221 + t219) * t270 + t11 * t28 + t12 * t29 + t34 * t37 + t35 * (-mrSges(7,1) * t56 + mrSges(7,2) * t55) + (-t281 * t30 - t282 * t31 - t283 * t53 - t284 * t54) * mrSges(5,3) + t61 * t13 + t3 * t62 + (t184 + qJD(1) * (Ifges(3,1) * t213 + t264)) * t252 / 0.2e1 - (t218 + t183) * t274 / 0.2e1 + t4 * t63 + t48 * t91 + t79 * t97 + t22 * t100 + t46 * t101 + t45 * t102 + t33 * t103 + t72 * t109 + t71 * t110 + t58 * t111 + t57 * t112 - t21 * (-mrSges(7,1) * t119 + mrSges(7,2) * t120) + t117 * t152; (-Ifges(4,1) * t314 - t232 * t326 - t235 * t325 - t318 * t346 - t324 * t345 + t340) * t173 + (mrSges(7,1) * t347 - mrSges(7,2) * t348) * t35 + (t182 * t32 - t41 * t50 - t44 * t49 - t52 * t64) * m(6) + (t272 - t49) * t100 + (-t285 * t53 - t286 * t54 + t343) * mrSges(5,3) + t193 * t336 - t347 * t335 - t348 * t334 + (-t273 - t59) * t102 + (-t298 + t300) * t257 + (t285 * t41 - t286 * t44) * mrSges(6,2) + (t122 - t169) * t317 + (t344 + (-t219 / 0.2e1 + t221) * qJD(1)) * qJD(1) + (t350 * t203 + t16 * mrSges(6,2) + Ifges(5,2) * t259 - Ifges(6,3) * t258 - t32 * mrSges(6,1) - t94 * mrSges(5,1) + t370 * t321 - t342 + (qJD(4) * t44 + t16 * t203) * m(6)) * t211 + (-t349 * t203 + t20 * mrSges(6,2) + (qJD(4) * t41 - qJD(5) * t52 + t20 * t203) * m(6) + t94 * mrSges(5,2) - t32 * mrSges(6,3) + t363 * t321 + t371 * t257) * t209 + (t272 - t60) * t101 + (-mrSges(3,1) * t245 + mrSges(3,2) * t251) * pkin(7) + (-t230 * qJD(4) + t133 * t226 + t203 * t343 + t205 * t94 - t53 * t59 - t54 * t60) * m(5) + (-t1 * t228 - t193 * t2 - t347 * t6 + t348 * t5) * mrSges(7,3) - t228 * t337 + (Ifges(7,1) * t193 - Ifges(7,4) * t228) * t333 + (Ifges(7,4) * t193 - Ifges(7,2) * t228) * t332 + (Ifges(7,5) * t193 - Ifges(7,6) * t228) * t322 - t21 * (mrSges(7,1) * t228 + mrSges(7,2) * t193) + t301 * t259 + t299 * t258 - t304 * t307 - t260 * t303 - t131 * t296 + t183 * t275 / 0.2e1 - t224 * t270 / 0.2e1 - t127 * t266 + t205 * t98 + t182 * t97 - Ifges(3,6) * t251 + Ifges(4,5) * t162 + t159 * t13 - Ifges(4,6) * t161 + t132 * t295 + t356 * t313 + t358 * t271 + t359 * t310 + (-t294 + t361) * t314 + (-Ifges(7,5) * t329 - Ifges(4,2) * t317 + Ifges(5,6) * t325 + Ifges(6,6) * t326 - Ifges(7,6) * t331 - Ifges(7,3) * t320 + t318 * t362 + t324 * t363 + t339) * t175 + t355 * t133 + Ifges(3,5) * t245 + t364 * t62 + t365 * t63 + (t1 * t126 + t125 * t2 - t159 * t21 + t364 * t6 + t365 * t5 + (t271 - t51) * t35) * m(7) + (t131 * t133 - t132 * t134 - t197 * t266 + (t210 * t95 - t290 * t94) * pkin(2)) * m(4) + t285 * t368 + (t273 - t50) * t103 - t51 * t37 - (-Ifges(3,2) * t275 + t184 + t207) * t253 / 0.2e1 + (-Ifges(7,4) * t108 - Ifges(7,2) * t107) * t331 + (-Ifges(7,5) * t108 - Ifges(7,6) * t107) * t320 + (-Ifges(7,1) * t108 - Ifges(7,4) * t107) * t329 + (-Ifges(7,1) * t177 - Ifges(7,4) * t178) * t328 + (-Ifges(7,4) * t177 - Ifges(7,2) * t178) * t330 + (-Ifges(7,5) * t177 - Ifges(7,6) * t178) * t319 - t64 * t91 - t94 * mrSges(4,1) - t95 * mrSges(4,2) + t125 * t28 + t126 * t29 - t134 * t152; -t228 * t28 + t193 * t29 - t347 * t63 - t348 * t62 + t349 * t211 + t350 * t209 + (t358 + t355) * t175 + (-t209 * t277 + t211 * t278 + t152) * t173 + t248 + (t1 * t193 + t175 * t35 - t2 * t228 - t347 * t5 - t348 * t6) * m(7) + (t16 * t209 - t20 * t211 - t175 * t52 - (-t209 * t41 - t211 * t44) * t173) * m(6) + (-t173 * t230 + t175 * t226 + t209 * t31 + t211 * t30) * m(5) + (t131 * t175 + t132 * t173 + t246) * m(4); t277 * t227 + t278 * t140 + t82 * t62 - t366 * t63 - t13 + t97 + t98 + (-t366 * t5 + t6 * t82 + t21) * m(7) + (t140 * t44 - t227 * t41 + t32) * m(6) + (t140 * t54 + t227 * t53 + t94) * m(5); t212 * t29 + t214 * t28 - t358 * t227 + t229 * qJD(6) + (-t100 - t229) * t173 + t111 + (t1 * t212 - t35 * t227 + t2 * t214 + t170 * (-t212 * t5 + t214 * t6)) * m(7) + (-t173 * t44 + t227 * t52 + t20) * m(6); -t35 * (mrSges(7,1) * t366 + mrSges(7,2) * t82) + (Ifges(7,1) * t82 - t308) * t329 + t24 * t328 + (Ifges(7,5) * t82 - Ifges(7,6) * t366) * t320 - t5 * t62 + t6 * t63 + (t366 * t6 + t5 * t82) * mrSges(7,3) + t268 + (-Ifges(7,2) * t366 + t25 + t81) * t331 + t341;];
tauc  = t7(:);
