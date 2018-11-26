% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6]';
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
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 17:34
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRRPPR3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR3_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR3_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPPR3_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR3_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:34:02
% EndTime: 2018-11-23 17:34:11
% DurationCPUTime: 8.79s
% Computational Cost: add. (6260->523), mult. (15258->649), div. (0->0), fcn. (10016->6), ass. (0->239)
t373 = Ifges(4,1) + Ifges(5,1);
t359 = Ifges(6,1) + Ifges(5,3);
t371 = -Ifges(4,5) - Ifges(5,4);
t372 = -mrSges(5,1) - mrSges(4,1);
t358 = -Ifges(4,4) + Ifges(5,5);
t370 = -Ifges(4,6) + Ifges(5,6);
t216 = qJD(2) + qJD(3);
t322 = cos(qJ(3));
t323 = cos(qJ(2));
t252 = t322 * t323;
t237 = qJD(1) * t252;
t221 = sin(qJ(3));
t222 = sin(qJ(2));
t283 = t221 * t222;
t248 = t216 * t283;
t124 = qJD(1) * t248 - t216 * t237;
t187 = -t252 + t283;
t327 = t187 / 0.2e1;
t369 = t124 * t327;
t188 = t221 * t323 + t222 * t322;
t174 = t188 * qJD(1);
t165 = Ifges(5,5) * t174;
t281 = qJD(1) * t222;
t173 = t221 * t281 - t237;
t302 = t174 * Ifges(6,4);
t368 = t165 - t302 + (-Ifges(6,5) + Ifges(5,6)) * t216 + t359 * t173;
t220 = sin(qJ(6));
t223 = cos(qJ(6));
t263 = qJD(1) * t323;
t195 = -qJD(1) * pkin(1) - pkin(2) * t263;
t95 = t173 * pkin(3) - t174 * qJ(4) + t195;
t239 = qJD(5) - t95;
t340 = -pkin(4) - pkin(9);
t34 = pkin(5) * t174 + t173 * t340 + t239;
t217 = -pkin(3) + t340;
t273 = t323 * pkin(7);
t197 = pkin(8) * t323 + t273;
t191 = t197 * qJD(1);
t175 = t221 * t191;
t196 = (-pkin(8) - pkin(7)) * t222;
t190 = qJD(1) * t196;
t181 = qJD(2) * pkin(2) + t190;
t282 = -t322 * t181 + t175;
t287 = t174 * qJ(5);
t86 = t282 - t287;
t236 = qJD(4) + t86;
t61 = t216 * t217 + t236;
t14 = -t220 * t61 + t223 * t34;
t15 = t220 * t34 + t223 * t61;
t243 = t14 * t223 + t15 * t220;
t143 = t216 * t188;
t125 = t143 * qJD(1);
t145 = -t173 * t220 - t216 * t223;
t62 = qJD(6) * t145 + t125 * t223;
t29 = -mrSges(7,1) * t124 - mrSges(7,3) * t62;
t265 = t322 * t191;
t135 = t221 * t181 + t265;
t229 = qJD(2) * t197;
t180 = t322 * t229;
t192 = qJD(2) * t196;
t238 = qJD(1) * t192;
t58 = qJD(1) * t180 + qJD(3) * t135 + t221 * t238;
t25 = t124 * qJ(5) - t174 * qJD(5) + t58;
t276 = qJD(1) * qJD(2);
t260 = t222 * t276;
t253 = pkin(2) * t260;
t33 = t125 * pkin(3) + t124 * qJ(4) - t174 * qJD(4) + t253;
t6 = -pkin(5) * t124 + t125 * t340 - t33;
t3 = -qJD(6) * t15 - t220 * t25 + t223 * t6;
t146 = t173 * t223 - t216 * t220;
t63 = -qJD(6) * t146 - t125 * t220;
t30 = mrSges(7,2) * t124 + mrSges(7,3) * t63;
t2 = qJD(6) * t14 + t220 * t6 + t223 * t25;
t319 = t2 * t223;
t367 = m(7) * (-qJD(6) * t243 - t3 * t220 + t319) - t220 * t29 + t223 * t30;
t366 = qJD(4) + t282;
t167 = Ifges(4,4) * t173;
t307 = t173 * Ifges(5,5);
t170 = qJD(6) + t174;
t356 = t170 * Ifges(7,3);
t357 = t145 * Ifges(7,6);
t365 = t146 * Ifges(7,5) + t373 * t174 - t371 * t216 - t167 + t307 + t356 + t357;
t214 = t216 * qJ(4);
t112 = t214 + t135;
t65 = -pkin(4) * t173 + t239;
t364 = t195 * mrSges(4,1) + t95 * mrSges(5,1) - t112 * mrSges(5,2) + t65 * mrSges(6,2);
t244 = Ifges(7,5) * t223 - Ifges(7,6) * t220;
t312 = Ifges(7,4) * t223;
t245 = -Ifges(7,2) * t220 + t312;
t313 = Ifges(7,4) * t220;
t246 = Ifges(7,1) * t223 - t313;
t247 = mrSges(7,1) * t220 + mrSges(7,2) * t223;
t332 = -t170 / 0.2e1;
t334 = -t146 / 0.2e1;
t335 = -t145 / 0.2e1;
t314 = Ifges(7,4) * t146;
t52 = Ifges(7,2) * t145 + Ifges(7,6) * t170 + t314;
t293 = qJ(5) * t173;
t87 = t135 + t293;
t76 = -t214 - t87;
t66 = pkin(5) * t216 - t76;
t363 = t244 * t332 + t245 * t335 + t246 * t334 + t220 * t52 / 0.2e1 - t247 * t66;
t362 = t65 * mrSges(6,1) + t14 * mrSges(7,1) + t195 * mrSges(4,2) - t15 * mrSges(7,2) - t95 * mrSges(5,3);
t361 = -t125 / 0.2e1;
t295 = -mrSges(6,1) * t216 + mrSges(7,1) * t145 - mrSges(7,2) * t146 - t173 * mrSges(6,3);
t318 = mrSges(4,3) * t173;
t150 = -mrSges(4,2) * t216 - t318;
t154 = -t173 * mrSges(5,2) + mrSges(5,3) * t216;
t354 = -t150 - t154;
t304 = t174 * mrSges(4,3);
t353 = -mrSges(5,2) * t174 - t216 * t372 - t304;
t317 = mrSges(6,3) * t174;
t153 = mrSges(6,2) * t216 - t317;
t92 = -mrSges(7,2) * t170 + mrSges(7,3) * t145;
t93 = mrSges(7,1) * t170 - mrSges(7,3) * t146;
t241 = -t220 * t93 + t223 * t92;
t352 = t153 + t241;
t240 = -t220 * t92 - t223 * t93;
t351 = qJD(6) * t240 + t367;
t350 = t222 * pkin(7) * (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t263) + (qJD(2) * mrSges(3,1) - mrSges(3,3) * t281) * t273;
t147 = -t322 * t196 + t197 * t221;
t148 = t221 * t196 + t197 * t322;
t348 = -t124 * t147 - t125 * t148 + t188 * t58;
t347 = m(5) * t112 + t154;
t96 = -pkin(3) * t216 + t366;
t346 = m(5) * t96 - t353;
t345 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t344 = m(6) * t76 - m(7) * t66 + t295;
t342 = t62 / 0.2e1;
t341 = t63 / 0.2e1;
t224 = -pkin(3) - pkin(4);
t339 = -t124 / 0.2e1;
t338 = t125 / 0.2e1;
t333 = t146 / 0.2e1;
t331 = -t173 / 0.2e1;
t330 = t173 / 0.2e1;
t329 = -t174 / 0.2e1;
t328 = t174 / 0.2e1;
t325 = -t216 / 0.2e1;
t324 = t216 / 0.2e1;
t321 = pkin(2) * t221;
t320 = pkin(4) * t174;
t210 = -t323 * pkin(2) - pkin(1);
t316 = mrSges(7,3) * t220;
t315 = Ifges(3,4) * t222;
t111 = t187 * qJ(5) + t148;
t176 = t221 * t229;
t261 = qJD(3) * t322;
t279 = qJD(3) * t221;
t57 = -qJD(1) * t176 + t181 * t261 - t191 * t279 + t322 * t238;
t48 = t216 * qJD(4) + t57;
t23 = -t125 * qJ(5) - t173 * qJD(5) - t48;
t311 = t111 * t23;
t309 = t124 * mrSges(6,3);
t308 = t147 * t58;
t64 = t216 * t224 + t236;
t306 = t173 * t64;
t305 = t173 * t96;
t303 = t174 * Ifges(4,4);
t301 = t187 * t23;
t294 = qJ(4) * t125;
t207 = qJ(4) + t321;
t290 = t125 * t207;
t289 = t143 * t220;
t288 = t143 * t223;
t286 = t174 * t223;
t285 = t187 * t220;
t284 = t187 * t223;
t127 = t174 * pkin(3) + t173 * qJ(4);
t280 = qJD(2) * t222;
t278 = qJD(6) * t220;
t277 = qJD(6) * t223;
t274 = Ifges(7,5) * t62 + Ifges(7,6) * t63 - Ifges(7,3) * t124;
t272 = t322 * pkin(2);
t271 = pkin(2) * t281;
t270 = pkin(2) * t280;
t268 = -t154 + t295;
t267 = Ifges(3,4) * t323;
t262 = qJD(2) * t323;
t258 = -t277 / 0.2e1;
t257 = -t124 * mrSges(6,1) + t125 * mrSges(6,2);
t139 = t190 * t221 + t265;
t256 = -t187 * pkin(3) + t188 * qJ(4) - t210;
t254 = pkin(2) * t261;
t209 = -t272 - pkin(3);
t251 = qJD(1) * t262;
t206 = -pkin(4) + t209;
t242 = t14 * t220 - t15 * t223;
t110 = -qJ(5) * t188 + t147;
t49 = pkin(5) * t188 + t187 * t340 + t256;
t27 = t110 * t223 + t220 * t49;
t26 = -t110 * t220 + t223 * t49;
t97 = t271 + t127;
t235 = Ifges(3,5) * t323 - Ifges(3,6) * t222;
t140 = t190 * t322 - t175;
t234 = t187 * t277 + t289;
t233 = t187 * t278 - t288;
t142 = -t216 * t252 + t248;
t42 = t143 * pkin(3) + t142 * qJ(4) - t188 * qJD(4) + t270;
t67 = t322 * t192 + t196 * t261 - t197 * t279 - t176;
t232 = pkin(1) * (mrSges(3,1) * t222 + mrSges(3,2) * t323);
t231 = t222 * (Ifges(3,1) * t323 - t315);
t230 = (Ifges(3,2) * t323 + t315) * qJD(1);
t39 = -pkin(5) * t173 + t174 * t340 - t127;
t68 = qJD(3) * t148 + t221 * t192 + t180;
t10 = t62 * Ifges(7,4) + t63 * Ifges(7,2) - t124 * Ifges(7,6);
t166 = Ifges(6,4) * t173;
t104 = -t174 * Ifges(6,2) - t216 * Ifges(6,6) + t166;
t105 = -t173 * Ifges(4,2) + t216 * Ifges(4,6) + t303;
t11 = t62 * Ifges(7,1) + t63 * Ifges(7,4) - Ifges(7,5) * t124;
t144 = Ifges(7,4) * t145;
t53 = t146 * Ifges(7,1) + t170 * Ifges(7,5) + t144;
t225 = -t223 * t10 / 0.2e1 - t220 * t11 / 0.2e1 - t57 * mrSges(4,2) + t48 * mrSges(5,3) + t25 * mrSges(6,2) + t3 * t316 + t76 * t317 + t282 * t318 + (-Ifges(7,5) * t220 - Ifges(7,6) * t223) * t339 + (-Ifges(7,2) * t223 - t313) * t341 + (-Ifges(7,1) * t220 - t312) * t342 + t372 * t58 + (t258 - t286 / 0.2e1) * t53 + (t302 + t105) * t328 + (-mrSges(7,1) * t223 + mrSges(7,2) * t220 - mrSges(6,1)) * t23 + (-t307 + t166 + t104) * t331 + (-Ifges(7,5) * t334 + Ifges(6,2) * t328 + Ifges(6,6) * t324 - Ifges(7,6) * t335 - Ifges(7,3) * t332 + t325 * t371 + t362) * t173 + (Ifges(6,5) * t324 + t15 * t316 + t325 * t370 + t331 * t359 + t363 - t364) * t174 + (-Ifges(6,5) + t370) * t125 + (-Ifges(6,6) + t371) * t124 + (t15 * t278 - t319 + (t277 + t286) * t14) * mrSges(7,3) + (-Ifges(4,2) * t174 - t167 + t365) * t330 + (-t173 * t373 + t165 - t303 + t368) * t329 + t363 * qJD(6);
t219 = qJ(4) + pkin(5);
t211 = Ifges(3,4) * t263;
t205 = pkin(5) + t207;
t172 = Ifges(3,1) * t281 + Ifges(3,5) * qJD(2) + t211;
t171 = Ifges(3,6) * qJD(2) + t230;
t132 = mrSges(4,1) * t173 + mrSges(4,2) * t174;
t130 = mrSges(5,1) * t173 - mrSges(5,3) * t174;
t126 = mrSges(6,1) * t174 + mrSges(6,2) * t173;
t91 = t140 + t287;
t90 = t139 + t293;
t89 = -pkin(4) * t187 + t256;
t85 = -t127 - t320;
t75 = -t97 - t320;
t38 = t39 - t271;
t32 = t142 * qJ(5) - t188 * qJD(5) + t68;
t28 = -pkin(4) * t143 - t42;
t22 = -pkin(4) * t125 - t33;
t21 = -mrSges(7,1) * t63 + mrSges(7,2) * t62;
t19 = t220 * t39 + t223 * t87;
t18 = -t220 * t87 + t223 * t39;
t17 = t220 * t38 + t223 * t90;
t16 = -t220 * t90 + t223 * t38;
t9 = -pkin(5) * t142 + t143 * t340 - t42;
t5 = -qJD(6) * t27 - t220 * t32 + t223 * t9;
t4 = qJD(6) * t26 + t220 * t9 + t223 * t32;
t1 = [t145 * (-Ifges(7,4) * t233 - Ifges(7,2) * t234) / 0.2e1 + (t235 * qJD(2) / 0.2e1 - t350) * qJD(2) + ((t361 - t338) * t188 + 0.2e1 * t369) * Ifges(6,4) + t344 * (-t187 * qJD(5) - t67) + t66 * (mrSges(7,1) * t234 - mrSges(7,2) * t233) + (t172 + qJD(1) * (Ifges(3,1) * t222 + t267)) * t262 / 0.2e1 - (t230 + t171) * t280 / 0.2e1 - (qJD(6) * t53 + t10) * t285 / 0.2e1 + m(5) * (t148 * t48 - t256 * t33 + t42 * t95 + t308) - t256 * (mrSges(5,1) * t125 + mrSges(5,3) * t124) + (m(4) * t282 + t346) * t68 + m(7) * (t14 * t5 + t15 * t4 + t2 * t27 + t26 * t3 - t311) + m(6) * (t110 * t25 + t22 * t89 + t28 * t65 + t32 * t64 - t311) + (-t187 * t48 + t348) * mrSges(5,2) + (m(4) * t135 + t150 + t347) * t67 + (t22 * mrSges(6,1) + mrSges(4,2) * t253 - t33 * mrSges(5,3) + Ifges(5,5) * t338 + Ifges(7,5) * t342 - Ifges(6,2) * t124 + Ifges(7,6) * t341 + t345) * t188 + (t231 - 0.2e1 * t232) * t276 + t132 * t270 + m(4) * (t308 + t148 * t57 + (qJD(1) * t210 + t195) * t270) + (-t357 / 0.2e1 - t356 / 0.2e1 + t104 / 0.2e1 + Ifges(6,6) * t325 + Ifges(6,2) * t329 - Ifges(7,5) * t333 + (-Ifges(5,5) + Ifges(6,4)) * t330 - t373 * t328 + t371 * t324 - t362 - Ifges(4,4) * t331 + t64 * mrSges(6,3) - t365 / 0.2e1 - t282 * mrSges(4,3) - t96 * mrSges(5,2)) * t142 + ((Ifges(7,3) + t373) * t188 + (t244 + t358) * t187) * t339 + (-t124 * t373 + t125 * t358 + t274) * t188 / 0.2e1 + (-Ifges(5,5) * t124 + t125 * t359) * t327 + (mrSges(4,1) * t253 + t33 * mrSges(5,1) + t22 * mrSges(6,2) + Ifges(4,2) * t125 + t245 * t341 + t246 * t342 + t52 * t258 + t338 * t359) * t187 + t170 * (-Ifges(7,5) * t233 - Ifges(7,6) * t234) / 0.2e1 + (-Ifges(7,1) * t233 - Ifges(7,4) * t234) * t333 - t247 * t301 - t52 * t289 / 0.2e1 + t53 * t288 / 0.2e1 + t11 * t284 / 0.2e1 + t210 * (mrSges(4,1) * t125 - mrSges(4,2) * t124) + (-Ifges(3,2) * t222 + t267) * t251 + t32 * t153 + (t188 * t361 + t369) * Ifges(4,4) + (t111 * t125 - t188 * t25 - t301) * mrSges(6,3) + t28 * t126 + t42 * t130 + t89 * t257 + t111 * t21 + t5 * t93 + t4 * t92 + t26 * t29 + t27 * t30 + (-t187 * t57 + t348) * mrSges(4,3) + (t14 * t233 - t15 * t234 - t2 * t285 - t284 * t3) * mrSges(7,3) + t110 * t309 + (-t105 / 0.2e1 - Ifges(4,2) * t331 + Ifges(6,5) * t325 + Ifges(6,4) * t329 + t359 * t330 + t370 * t324 + t364 - t76 * mrSges(6,3) + t368 / 0.2e1 - t135 * mrSges(4,3) - t344 * qJ(5) + t358 * t328) * t143; (t206 * t25 - t207 * t23 - t64 * t90 - t65 * t75 + t76 * t91) * m(6) + t354 * t140 + t353 * t139 + (-t112 * t140 - t139 * t96 + t207 * t48 + t209 * t58 - t95 * t97) * m(5) + t225 + t295 * t91 - (-Ifges(3,2) * t281 + t172 + t211) * t263 / 0.2e1 + (-t277 * t93 - t278 * t92 + t367) * (-pkin(9) + t206) + (-t14 * t16 - t15 * t17 - t205 * t23 - t66 * t91) * m(7) + (-t124 * t209 - t290 + t305) * mrSges(5,2) + (-t282 * t139 - t135 * t140 - t195 * t271 + (-t322 * t58 + t221 * t57 + (t135 * t322 + t221 * t282) * qJD(3)) * pkin(2)) * m(4) + (t124 * t272 - t125 * t321) * mrSges(4,3) + (-t344 + t347) * (t254 + qJD(4)) + (-mrSges(3,1) * t251 + mrSges(3,2) * t260) * pkin(7) + (m(6) * t64 - m(7) * t242 + t346 + t352) * pkin(2) * t279 + t150 * t254 + Ifges(3,5) * t251 + (-t306 + t290) * mrSges(6,3) + (t350 + (t232 - t231 / 0.2e1) * qJD(1)) * qJD(1) + t171 * t281 / 0.2e1 - t235 * t276 / 0.2e1 + t205 * t21 - t132 * t271 - t90 * t153 - t75 * t126 - t97 * t130 - Ifges(3,6) * t260 - t16 * t93 - t17 * t92 + t135 * t304 + t206 * t309; (t124 * t224 + t294 - t306) * mrSges(6,3) + (pkin(3) * t124 - t294 + t305) * mrSges(5,2) + t225 - t295 * t86 + (t304 + t353) * t135 - t354 * t282 - t268 * qJD(4) + t219 * t21 - t87 * t153 - t85 * t126 - t127 * t130 - t19 * t92 - t18 * t93 + t351 * t217 + (-t14 * t18 - t15 * t19 - t219 * t23 + t236 * t66) * m(7) + (-t23 * qJ(4) + t224 * t25 - t236 * t76 - t64 * t87 - t65 * t85) * m(6) + (-pkin(3) * t58 + qJ(4) * t48 + t112 * t366 - t127 * t95 - t135 * t96) * m(5); (-mrSges(5,2) + mrSges(6,3)) * t124 + t268 * t216 + (-t126 + t130 + t240) * t174 - m(7) * (t174 * t243 + t216 * t66) + (-t174 * t65 + t216 * t76 + t25) * m(6) + (-t112 * t216 + t174 * t95 + t58) * m(5) + t351; t220 * t30 + t223 * t29 + t295 * t173 + t241 * qJD(6) + t352 * t174 + t257 + (-t170 * t242 - t173 * t66 + t2 * t220 + t223 * t3) * m(7) + (t173 * t76 + t174 * t64 + t22) * m(6); -t66 * (mrSges(7,1) * t146 + mrSges(7,2) * t145) + (Ifges(7,1) * t145 - t314) * t334 + t52 * t333 + (Ifges(7,5) * t145 - Ifges(7,6) * t146) * t332 - t14 * t92 + t15 * t93 + (t14 * t145 + t146 * t15) * mrSges(7,3) + t274 + (-Ifges(7,2) * t146 + t144 + t53) * t335 + t345;];
tauc  = t1(:);
