% Calculate vector of inverse dynamics joint torques for
% S5RRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPP4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP4_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP4_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP4_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP4_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP4_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP4_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP4_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:54:53
% EndTime: 2019-12-31 20:55:11
% DurationCPUTime: 8.81s
% Computational Cost: add. (5016->467), mult. (11802->589), div. (0->0), fcn. (8122->12), ass. (0->218)
t352 = -mrSges(5,1) - mrSges(6,1);
t340 = Ifges(5,1) + Ifges(6,1);
t339 = -Ifges(5,4) + Ifges(6,5);
t338 = Ifges(6,4) + Ifges(5,5);
t351 = mrSges(6,2) + mrSges(5,3);
t215 = sin(qJ(2));
t218 = cos(qJ(2));
t178 = -mrSges(3,1) * t218 + mrSges(3,2) * t215;
t212 = qJ(2) + qJ(3);
t205 = pkin(8) + t212;
t192 = sin(t205);
t193 = cos(t205);
t206 = sin(t212);
t207 = cos(t212);
t323 = -t207 * mrSges(4,1) + t206 * mrSges(4,2) + t352 * t193 + (mrSges(5,2) - mrSges(6,3)) * t192;
t350 = t178 + t323;
t214 = sin(qJ(3));
t217 = cos(qJ(3));
t159 = -t214 * t215 + t217 * t218;
t146 = t159 * qJD(1);
t160 = t214 * t218 + t215 * t217;
t147 = t160 * qJD(1);
t213 = sin(pkin(8));
t281 = cos(pkin(8));
t102 = -t281 * t146 + t147 * t213;
t316 = t102 / 0.2e1;
t349 = -Ifges(5,6) + Ifges(6,6);
t211 = qJD(2) + qJD(3);
t232 = t213 * t146 + t147 * t281;
t348 = t339 * t102 + t338 * t211 + t340 * t232;
t220 = -pkin(7) - pkin(6);
t182 = t220 * t218;
t166 = qJD(1) * t182;
t151 = t217 * t166;
t181 = t220 * t215;
t165 = qJD(1) * t181;
t115 = -t165 * t214 + t151;
t280 = qJ(4) * t146;
t235 = t115 - t280;
t250 = t281 * t214;
t148 = t214 * t166;
t116 = t217 * t165 + t148;
t140 = t147 * qJ(4);
t89 = -t140 + t116;
t332 = t213 * t89 - t235 * t281 - (t213 * t217 + t250) * qJD(3) * pkin(2);
t347 = mrSges(4,1) * t206 + mrSges(5,1) * t192 + mrSges(4,2) * t207 + mrSges(5,2) * t193;
t267 = qJD(1) * qJD(2);
t169 = qJDD(1) * t218 - t215 * t267;
t216 = sin(qJ(1));
t219 = cos(qJ(1));
t327 = g(1) * t219 + g(2) * t216;
t346 = t193 * (m(6) * qJ(5) + mrSges(6,3));
t208 = t218 * pkin(2);
t197 = t208 + pkin(1);
t180 = t197 * qJD(1);
t121 = -pkin(3) * t146 + qJD(4) - t180;
t155 = qJD(2) * pkin(2) + t165;
t108 = t217 * t155 + t148;
t84 = t108 - t140;
t75 = pkin(3) * t211 + t84;
t109 = t155 * t214 - t151;
t85 = t109 + t280;
t78 = t281 * t85;
t28 = t213 * t75 + t78;
t26 = qJ(5) * t211 + t28;
t36 = pkin(4) * t102 - qJ(5) * t232 + t121;
t345 = t121 * mrSges(5,1) + t36 * mrSges(6,1) - t26 * mrSges(6,2) - t28 * mrSges(5,3);
t282 = t213 * t85;
t27 = t281 * t75 - t282;
t25 = -t211 * pkin(4) + qJD(5) - t27;
t344 = mrSges(5,2) * t121 + t25 * mrSges(6,2) - t27 * mrSges(5,3) - mrSges(6,3) * t36;
t170 = qJDD(1) * t215 + t218 * t267;
t230 = t159 * qJD(3);
t90 = qJD(1) * t230 + t169 * t214 + t170 * t217;
t231 = t160 * qJD(3);
t91 = -qJD(1) * t231 + t169 * t217 - t170 * t214;
t48 = t213 * t91 + t281 * t90;
t319 = t48 / 0.2e1;
t343 = m(5) + m(6);
t342 = t169 / 0.2e1;
t209 = qJDD(2) + qJDD(3);
t309 = t209 / 0.2e1;
t306 = t218 / 0.2e1;
t313 = t232 / 0.2e1;
t92 = -mrSges(5,2) * t211 - mrSges(5,3) * t102;
t95 = -mrSges(6,2) * t102 + mrSges(6,3) * t211;
t337 = -t92 - t95;
t336 = t352 * t211 + t351 * t232;
t335 = Ifges(5,4) * t232;
t334 = Ifges(6,5) * t232;
t333 = t218 * Ifges(3,2);
t124 = t214 * t181 - t217 * t182;
t237 = t193 * pkin(4) + t192 * qJ(5);
t271 = qJD(1) * t218;
t272 = qJD(1) * t215;
t297 = pkin(6) * t218;
t298 = pkin(6) * t215;
t329 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t272) * t297 + (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t271) * t298;
t157 = t169 * pkin(6);
t158 = t170 * pkin(6);
t328 = t157 * t218 + t158 * t215;
t325 = 0.2e1 * t309;
t322 = -m(3) * pkin(6) + m(4) * t220 + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - t351;
t321 = m(3) * pkin(1) + m(4) * t197 + mrSges(2,1) - t350;
t317 = -t102 / 0.2e1;
t310 = t147 / 0.2e1;
t307 = t211 / 0.2e1;
t305 = pkin(2) * t214;
t304 = pkin(2) * t215;
t303 = pkin(2) * t217;
t302 = pkin(3) * t147;
t301 = pkin(3) * t206;
t195 = pkin(3) * t207;
t300 = pkin(3) * t213;
t299 = pkin(4) * t192;
t120 = qJDD(2) * pkin(2) - pkin(7) * t170 - t158;
t122 = pkin(7) * t169 + t157;
t52 = -qJD(3) * t109 + t217 * t120 - t122 * t214;
t12 = pkin(3) * t209 - qJ(4) * t90 - qJD(4) * t147 + t52;
t268 = qJD(3) * t217;
t269 = qJD(3) * t214;
t51 = t214 * t120 + t217 * t122 + t155 * t268 + t166 * t269;
t14 = qJ(4) * t91 + qJD(4) * t146 + t51;
t6 = t213 * t12 + t281 * t14;
t290 = Ifges(3,4) * t215;
t289 = Ifges(3,4) * t218;
t288 = t108 * mrSges(4,3);
t287 = t109 * mrSges(4,3);
t286 = t147 * Ifges(4,4);
t285 = t192 * mrSges(6,1);
t279 = qJDD(1) * pkin(1);
t196 = pkin(3) + t303;
t139 = pkin(2) * t250 + t213 * t196;
t273 = t195 + t208;
t270 = qJD(2) * t215;
t263 = pkin(2) * t269;
t262 = pkin(2) * t268;
t201 = pkin(2) * t270;
t200 = pkin(2) * t272;
t260 = t195 + t237;
t259 = qJD(2) * t220;
t258 = t281 * pkin(3);
t47 = t213 * t90 - t281 * t91;
t256 = t47 * mrSges(5,1) + t48 * mrSges(5,2);
t255 = t47 * mrSges(6,1) - t48 * mrSges(6,3);
t252 = t216 * t346;
t251 = t219 * t346;
t31 = -t209 * mrSges(6,1) + t48 * mrSges(6,2);
t118 = -qJD(2) * t160 - t231;
t105 = -pkin(3) * t118 + t201;
t249 = t267 / 0.2e1;
t123 = t217 * t181 + t182 * t214;
t128 = -pkin(3) * t159 - t197;
t245 = -g(1) * t216 + g(2) * t219;
t244 = mrSges(3,1) * t215 + mrSges(3,2) * t218;
t239 = t290 + t333;
t238 = Ifges(3,5) * t218 - Ifges(3,6) * t215;
t141 = -pkin(2) * t169 - t279;
t137 = -t213 * t263 + t281 * t262;
t236 = pkin(1) * t244;
t234 = -qJ(4) * t160 + t123;
t5 = t12 * t281 - t213 * t14;
t233 = t215 * (Ifges(3,1) * t218 - t290);
t167 = t215 * t259;
t168 = t218 * t259;
t71 = t217 * t167 + t214 * t168 + t181 * t268 + t182 * t269;
t53 = pkin(4) * t232 + qJ(5) * t102 + t302;
t138 = t196 * t281 - t213 * t305;
t171 = -t301 - t304;
t226 = m(6) * (t171 - t299) - t285;
t67 = -pkin(3) * t91 + qJDD(4) + t141;
t223 = m(6) * (-t299 - t301) - t285;
t72 = -qJD(3) * t124 - t167 * t214 + t217 * t168;
t117 = qJD(2) * t159 + t230;
t222 = -qJ(4) * t117 - qJD(4) * t160 + t72;
t142 = Ifges(4,4) * t146;
t2 = qJ(5) * t209 + qJD(5) * t211 + t6;
t3 = -t209 * pkin(4) + qJDD(5) - t5;
t59 = t211 * Ifges(6,6) + t102 * Ifges(6,3) + t334;
t60 = -t102 * Ifges(5,2) + t211 * Ifges(5,6) + t335;
t98 = t146 * Ifges(4,2) + t211 * Ifges(4,6) + t286;
t99 = t147 * Ifges(4,1) + t211 * Ifges(4,5) + t142;
t221 = t147 * t287 + t98 * t310 + t146 * t288 - t51 * mrSges(4,2) + t52 * mrSges(4,1) + Ifges(4,5) * t90 + Ifges(4,6) * t91 + t60 * t313 - t147 * (Ifges(4,1) * t146 - t286) / 0.2e1 + t180 * (mrSges(4,1) * t147 + mrSges(4,2) * t146) + t2 * mrSges(6,3) - t3 * mrSges(6,1) + t5 * mrSges(5,1) - t6 * mrSges(5,2) + t338 * t48 + t349 * t47 + t348 * t316 - (-Ifges(4,2) * t147 + t142 + t99) * t146 / 0.2e1 + (-Ifges(5,2) * t316 + Ifges(6,3) * t317 - t345) * t232 + (-Ifges(5,4) * t316 - Ifges(6,5) * t317 + t344) * t102 - (-t340 * t102 + t334 - t335 + t59) * t232 / 0.2e1 - (Ifges(4,5) * t146 - Ifges(4,6) * t147 - t102 * t338 + t232 * t349) * t211 / 0.2e1 + (Ifges(6,2) + Ifges(4,3) + Ifges(5,3)) * t209;
t210 = -qJ(4) + t220;
t199 = Ifges(3,4) * t271;
t194 = -t258 - pkin(4);
t190 = qJ(5) + t300;
t164 = pkin(1) + t273;
t145 = Ifges(3,1) * t272 + Ifges(3,5) * qJD(2) + t199;
t144 = Ifges(3,6) * qJD(2) + qJD(1) * t239;
t134 = -pkin(4) - t138;
t133 = qJ(5) + t139;
t130 = qJD(5) + t137;
t127 = mrSges(4,1) * t211 - mrSges(4,3) * t147;
t126 = -mrSges(4,2) * t211 + mrSges(4,3) * t146;
t125 = t200 + t302;
t112 = t213 * t159 + t160 * t281;
t111 = -t159 * t281 + t160 * t213;
t107 = -mrSges(4,1) * t146 + mrSges(4,2) * t147;
t100 = qJ(4) * t159 + t124;
t77 = -mrSges(4,2) * t209 + mrSges(4,3) * t91;
t76 = mrSges(4,1) * t209 - mrSges(4,3) * t90;
t69 = t117 * t281 + t213 * t118;
t68 = t117 * t213 - t118 * t281;
t66 = mrSges(5,1) * t102 + mrSges(5,2) * t232;
t65 = mrSges(6,1) * t102 - mrSges(6,3) * t232;
t54 = pkin(4) * t111 - qJ(5) * t112 + t128;
t49 = t200 + t53;
t46 = t213 * t235 + t281 * t89;
t44 = qJ(4) * t118 + qJD(4) * t159 + t71;
t34 = t281 * t84 - t282;
t33 = t213 * t84 + t78;
t32 = -mrSges(6,2) * t47 + mrSges(6,3) * t209;
t30 = mrSges(5,1) * t209 - mrSges(5,3) * t48;
t29 = -mrSges(5,2) * t209 - mrSges(5,3) * t47;
t10 = pkin(4) * t68 - qJ(5) * t69 - qJD(5) * t112 + t105;
t7 = pkin(4) * t47 - qJ(5) * t48 - qJD(5) * t232 + t67;
t1 = [t170 * t289 / 0.2e1 + (t218 * t289 + t233) * t249 + (Ifges(5,4) * t317 + Ifges(6,5) * t316 + t307 * t338 + t313 * t340 + t344 + t348 / 0.2e1) * t69 + t118 * t287 + (Ifges(4,1) * t117 + Ifges(4,4) * t118) * t310 + (Ifges(6,3) * t316 - Ifges(5,2) * t317 + t59 / 0.2e1 - t60 / 0.2e1 + t339 * t313 + t349 * t307 + t345) * t68 + (mrSges(5,1) * t67 + mrSges(6,1) * t7 - mrSges(6,2) * t2 - mrSges(5,3) * t6 + (Ifges(6,3) + Ifges(5,2)) * t47 + 0.2e1 * t339 * t319 + t349 * t325) * t111 + t107 * t201 + (-mrSges(3,1) * t298 - mrSges(3,2) * t297 + 0.2e1 * Ifges(3,6) * t306) * qJDD(2) + (Ifges(3,1) * t170 + Ifges(3,4) * t342 + Ifges(3,5) * qJDD(2) - t249 * t333) * t215 + ((t210 * t343 + t322) * t219 + (m(5) * t164 - m(6) * (-t164 - t237) + t321) * t216) * g(1) + (-t343 * (t219 * t164 - t210 * t216) + (-m(6) * t237 - t321) * t219 + t322 * t216) * g(2) + (Ifges(4,5) * t117 + Ifges(4,6) * t118) * t307 + (Ifges(3,4) * t170 + Ifges(3,2) * t169) * t306 + (t209 * t338 + t339 * t47 + t340 * t48) * t112 / 0.2e1 + (-m(5) * t27 + m(6) * t25 + t336) * (t213 * t44 - t222 * t281) + (m(5) * t28 + m(6) * t26 - t337) * (t213 * t222 + t281 * t44) + (mrSges(5,2) * t67 + mrSges(6,2) * t3 - mrSges(5,3) * t5 - mrSges(6,3) * t7 + (Ifges(6,5) / 0.2e1 - Ifges(5,4) / 0.2e1) * t47 + t340 * t319 + t338 * t309) * t112 + m(6) * (t10 * t36 + t54 * t7) + m(5) * (t105 * t121 + t128 * t67) + (t169 * t297 + t170 * t298 + t328) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(6) * t328) + (t145 * t306 + t238 * qJD(2) / 0.2e1 - t329) * qJD(2) + (-t141 * mrSges(4,1) + t51 * mrSges(4,3) + Ifges(4,4) * t90 + Ifges(4,2) * t91 + Ifges(4,6) * t325) * t159 + (t141 * mrSges(4,2) - t52 * mrSges(4,3) + Ifges(4,1) * t90 + Ifges(4,4) * t91 + Ifges(4,5) * t325) * t160 + m(4) * (t108 * t72 + t109 * t71 + t123 * t52 + t124 * t51 - t141 * t197 - t180 * t201) + t10 * t65 + t239 * t342 + t105 * t66 + t117 * t99 / 0.2e1 + t118 * t98 / 0.2e1 + t54 * t255 + t128 * t256 + t123 * t76 + t124 * t77 + t71 * t126 + t72 * t127 + t146 * (Ifges(4,4) * t117 + Ifges(4,2) * t118) / 0.2e1 - pkin(1) * (-mrSges(3,1) * t169 + mrSges(3,2) * t170) - t180 * (-mrSges(4,1) * t118 + mrSges(4,2) * t117) - t236 * t267 - t144 * t270 / 0.2e1 - t197 * (-mrSges(4,1) * t91 + mrSges(4,2) * t90) - t178 * t279 + Ifges(2,3) * qJDD(1) - t117 * t288 + (-m(5) * t5 + m(6) * t3 - t30 + t31) * (t100 * t213 - t234 * t281) + (m(5) * t6 + m(6) * t2 + t29 + t32) * (t100 * t281 + t213 * t234); -t332 * t336 + t221 + t327 * (m(4) * t304 - m(5) * t171 + t244 + t347) + (t262 - t116) * t126 + (-t121 * t125 + t138 * t5 + t139 * t6 + (-t46 + t137) * t28 + t332 * t27) * m(5) + (t133 * t2 + t134 * t3 - t36 * t49 + (-t46 + t130) * t26 - t332 * t25) * m(6) + (-t115 - t263) * t127 + t76 * t303 + t77 * t305 + t337 * t46 + (-m(5) * t273 - m(6) * (t208 + t260) - m(4) * t208 + t350) * g(3) + ((t214 * t51 + t217 * t52 + (-t108 * t214 + t109 * t217) * qJD(3)) * pkin(2) - t108 * t115 - t109 * t116 + t180 * t200) * m(4) + Ifges(3,3) * qJDD(2) - t49 * t65 - g(1) * (t219 * t226 + t251) - g(2) * (t216 * t226 + t252) - t125 * t66 + t130 * t95 + t133 * t32 + t134 * t31 + t137 * t92 + t138 * t30 + t139 * t29 - t107 * t200 - t157 * mrSges(3,2) - t158 * mrSges(3,1) + Ifges(3,6) * t169 + Ifges(3,5) * t170 - t238 * t267 / 0.2e1 + t144 * t272 / 0.2e1 - (-Ifges(3,2) * t272 + t145 + t199) * t271 / 0.2e1 + (t329 + (-t233 / 0.2e1 + t236) * qJD(1)) * qJD(1); t30 * t258 + t221 + t29 * t300 - t53 * t65 + qJD(5) * t95 - g(1) * (t219 * t223 + t251) - g(2) * (t216 * t223 + t252) - t108 * t126 + t109 * t127 + t190 * t32 + t194 * t31 - t66 * t302 + t337 * t34 - t336 * t33 + (t190 * t2 + t194 * t3 - t25 * t33 - t36 * t53 + (-t34 + qJD(5)) * t26) * m(6) + ((t213 * t6 + t281 * t5) * pkin(3) - t121 * t302 + t27 * t33 - t28 * t34) * m(5) + (-m(5) * t195 - m(6) * t260 + t323) * g(3) + (m(5) * t301 + t347) * t327; -t337 * t102 - t336 * t232 + t255 + t256 + (t102 * t26 - t232 * t25 + t245 + t7) * m(6) + (t102 * t28 + t232 * t27 + t245 + t67) * m(5); t232 * t65 - t211 * t95 + (g(3) * t193 - t192 * t327 - t26 * t211 + t36 * t232 + t3) * m(6) + t31;];
tau = t1;
