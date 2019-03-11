% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
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
% Datum: 2019-03-09 06:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRRP10_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP10_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP10_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP10_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP10_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP10_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP10_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:30:19
% EndTime: 2019-03-09 06:30:45
% DurationCPUTime: 14.89s
% Computational Cost: add. (7325->618), mult. (16137->844), div. (0->0), fcn. (10025->6), ass. (0->264)
t360 = Ifges(6,1) + Ifges(7,1);
t359 = Ifges(7,4) + Ifges(6,5);
t358 = Ifges(7,5) - Ifges(6,4);
t202 = sin(qJ(5));
t203 = sin(qJ(4));
t205 = cos(qJ(5));
t206 = cos(qJ(4));
t166 = t202 * t203 - t205 * t206;
t264 = qJD(4) + qJD(5);
t111 = t264 * t166;
t204 = sin(qJ(3));
t351 = t166 * t204;
t133 = qJD(1) * t351;
t372 = t111 + t133;
t167 = t202 * t206 + t203 * t205;
t149 = t167 * qJD(1);
t132 = t204 * t149;
t366 = t264 * t167;
t279 = t366 + t132;
t370 = Ifges(6,6) - Ifges(7,6);
t356 = Ifges(6,3) + Ifges(7,2);
t371 = -Ifges(4,6) * qJD(3) / 0.2e1;
t207 = cos(qJ(3));
t266 = qJD(1) * qJD(3);
t246 = t207 * t266;
t276 = qJD(3) * t203;
t277 = qJD(1) * t207;
t164 = t206 * t277 + t276;
t274 = qJD(3) * t206;
t218 = t203 * t277 - t274;
t107 = t202 * t164 + t205 * t218;
t269 = qJD(4) * t207;
t216 = -t203 * t269 - t204 * t274;
t265 = qJD(3) * qJD(4);
t123 = qJD(1) * t216 + t206 * t265;
t275 = qJD(3) * t204;
t257 = t203 * t275;
t365 = t206 * t269 - t257;
t124 = -qJD(1) * t365 - t203 * t265;
t41 = -qJD(5) * t107 + t205 * t123 + t202 * t124;
t213 = t205 * t164 - t202 * t218;
t42 = qJD(5) * t213 + t202 * t123 - t205 * t124;
t369 = t359 * t246 + t358 * t42 + t360 * t41;
t208 = -pkin(1) - pkin(7);
t188 = qJD(1) * t208 + qJD(2);
t172 = t204 * t188;
t201 = t204 * qJD(1);
t258 = t203 * t201;
t134 = -pkin(4) * t258 + t172;
t271 = qJD(4) * t203;
t368 = pkin(4) * t271 + t279 * pkin(5) + t372 * qJ(6) - qJD(6) * t167 - t134;
t101 = Ifges(6,4) * t107;
t186 = t201 + t264;
t298 = Ifges(7,5) * t107;
t355 = t359 * t186 + t360 * t213 - t101 + t298;
t337 = t42 / 0.2e1;
t367 = Ifges(7,3) * t337;
t281 = t207 * t188;
t156 = -qJD(3) * pkin(3) - t281;
t117 = pkin(4) * t218 + t156;
t173 = pkin(3) * t204 - pkin(8) * t207 + qJ(2);
t151 = t173 * qJD(1);
t155 = qJD(3) * pkin(8) + t172;
t98 = t203 * t151 + t206 * t155;
t83 = -pkin(9) * t218 + t98;
t293 = t202 * t83;
t195 = t201 + qJD(4);
t97 = t206 * t151 - t155 * t203;
t82 = -pkin(9) * t164 + t97;
t72 = pkin(4) * t195 + t82;
t21 = t205 * t72 - t293;
t350 = qJD(6) - t21;
t19 = -pkin(5) * t186 + t350;
t33 = t107 * pkin(5) - qJ(6) * t213 + t117;
t364 = t117 * mrSges(6,2) + t19 * mrSges(7,2) - t21 * mrSges(6,3) - t33 * mrSges(7,3) + t355 / 0.2e1;
t62 = pkin(5) * t213 + qJ(6) * t107;
t339 = t41 / 0.2e1;
t327 = t123 / 0.2e1;
t326 = t124 / 0.2e1;
t362 = t246 / 0.2e1;
t141 = t166 * t207;
t354 = -qJD(3) * t141 - t204 * t366 - t149;
t239 = t264 * t206;
t268 = qJD(5) * t202;
t249 = t203 * t268;
t252 = t204 * t271;
t273 = qJD(3) * t207;
t254 = t203 * t273;
t255 = t206 * t273;
t353 = t204 * t249 - (t204 * t239 + t254) * t205 - (t255 - t252) * t202 + t166 * qJD(1);
t352 = qJ(2) * (m(3) + m(4));
t284 = t204 * t208;
t185 = t206 * t284;
t129 = t203 * t173 + t185;
t349 = t356 * t246 + t359 * t41 - t370 * t42;
t235 = pkin(3) * t207 + pkin(8) * t204;
t162 = qJD(3) * t235 + qJD(2);
t137 = t162 * qJD(1);
t270 = qJD(4) * t206;
t53 = t203 * t137 + t151 * t270 - t155 * t271 + t188 * t255;
t214 = -qJD(4) * t98 + t206 * t137;
t54 = -t188 * t254 + t214;
t225 = -t203 * t54 + t206 * t53;
t126 = -t195 * mrSges(5,2) - mrSges(5,3) * t218;
t127 = mrSges(5,1) * t195 - mrSges(5,3) * t164;
t220 = -t203 * t126 - t206 * t127;
t224 = t203 * t98 + t206 * t97;
t348 = -m(5) * t224 + t220;
t291 = t205 * t83;
t22 = t202 * t72 + t291;
t26 = -pkin(9) * t123 + (pkin(4) * qJD(1) - t188 * t203) * t273 + t214;
t28 = pkin(9) * t124 + t53;
t6 = -qJD(5) * t22 - t202 * t28 + t205 * t26;
t161 = t206 * t173;
t241 = -t203 * t208 + pkin(4);
t282 = t206 * t207;
t105 = -pkin(9) * t282 + t204 * t241 + t161;
t285 = t203 * t207;
t116 = -pkin(9) * t285 + t129;
t289 = t202 * t105 + t205 * t116;
t143 = t206 * t162;
t263 = pkin(9) * t204 * t206;
t58 = t143 + (-t185 + (pkin(9) * t207 - t173) * t203) * qJD(4) + (t207 * t241 + t263) * qJD(3);
t272 = qJD(3) * t208;
t253 = t207 * t272;
t80 = t203 * t162 + t173 * t270 + t206 * t253 - t208 * t252;
t65 = -pkin(9) * t365 + t80;
t11 = -qJD(5) * t289 - t202 * t65 + t205 * t58;
t297 = Ifges(5,6) * t203;
t228 = Ifges(5,5) * t206 - t297;
t301 = Ifges(5,4) * t203;
t232 = Ifges(5,1) * t206 - t301;
t233 = mrSges(5,1) * t203 + mrSges(5,2) * t206;
t316 = t206 / 0.2e1;
t318 = -t203 / 0.2e1;
t319 = t195 / 0.2e1;
t323 = t164 / 0.2e1;
t302 = Ifges(5,4) * t164;
t94 = -Ifges(5,2) * t218 + Ifges(5,6) * t195 + t302;
t158 = Ifges(5,4) * t218;
t95 = t164 * Ifges(5,1) + t195 * Ifges(5,5) - t158;
t347 = -mrSges(5,3) * t224 + t156 * t233 + t228 * t319 + t232 * t323 + t316 * t95 + t318 * t94;
t20 = qJ(6) * t186 + t22;
t100 = Ifges(7,5) * t213;
t47 = Ifges(7,6) * t186 + Ifges(7,3) * t107 + t100;
t299 = Ifges(6,4) * t213;
t50 = -Ifges(6,2) * t107 + Ifges(6,6) * t186 + t299;
t345 = t117 * mrSges(6,1) + t33 * mrSges(7,1) + t47 / 0.2e1 - t50 / 0.2e1 - t20 * mrSges(7,2) - t22 * mrSges(6,3);
t344 = t97 * mrSges(5,1) + t20 * mrSges(7,3) + t21 * mrSges(6,1) + t371 - (Ifges(4,4) * t207 - t204 * Ifges(4,2)) * qJD(1) / 0.2e1 - t98 * mrSges(5,2) - t19 * mrSges(7,1) - t22 * mrSges(6,2);
t343 = m(5) / 0.2e1;
t342 = m(7) / 0.2e1;
t341 = Ifges(7,5) * t339 + Ifges(7,6) * t362 + t367;
t340 = -t41 * Ifges(6,4) / 0.2e1 + Ifges(6,2) * t337 - Ifges(6,6) * t246 / 0.2e1;
t338 = -t42 / 0.2e1;
t335 = Ifges(5,1) * t327 + Ifges(5,4) * t326 + Ifges(5,5) * t362;
t333 = -pkin(9) - pkin(8);
t332 = -t107 / 0.2e1;
t331 = t107 / 0.2e1;
t329 = -t213 / 0.2e1;
t328 = t213 / 0.2e1;
t324 = -t164 / 0.2e1;
t321 = -t186 / 0.2e1;
t320 = t186 / 0.2e1;
t315 = m(6) * t117;
t314 = t53 * mrSges(5,2);
t313 = t54 * mrSges(5,1);
t30 = mrSges(6,1) * t246 - mrSges(6,3) * t41;
t31 = -mrSges(7,1) * t246 + t41 * mrSges(7,2);
t309 = t31 - t30;
t29 = -mrSges(7,2) * t42 + mrSges(7,3) * t246;
t32 = -mrSges(6,2) * t246 - mrSges(6,3) * t42;
t308 = t32 + t29;
t85 = -mrSges(7,2) * t107 + mrSges(7,3) * t186;
t305 = mrSges(6,3) * t107;
t86 = -mrSges(6,2) * t186 - t305;
t307 = t85 + t86;
t304 = mrSges(6,3) * t213;
t87 = mrSges(6,1) * t186 - t304;
t88 = -mrSges(7,1) * t186 + mrSges(7,2) * t213;
t306 = -t88 + t87;
t170 = t235 * qJD(1);
t113 = t206 * t170 - t203 * t281;
t89 = (pkin(4) * t207 + t263) * qJD(1) + t113;
t114 = t203 * t170 + t206 * t281;
t96 = pkin(9) * t258 + t114;
t44 = t202 * t89 + t205 * t96;
t303 = Ifges(4,4) * t204;
t300 = Ifges(5,4) * t206;
t296 = Ifges(5,3) * t195;
t295 = qJ(2) * mrSges(4,1);
t294 = qJ(2) * mrSges(4,2);
t287 = qJD(3) * mrSges(4,2);
t278 = qJD(3) * mrSges(4,1) - mrSges(5,1) * t218 - t164 * mrSges(5,2) - mrSges(4,3) * t277;
t267 = qJD(5) * t205;
t260 = Ifges(5,5) * t123 + Ifges(5,6) * t124 + Ifges(5,3) * t246;
t200 = -pkin(4) * t206 - pkin(3);
t259 = qJD(4) * t333;
t165 = t188 * t275;
t194 = t204 * t272;
t245 = t276 / 0.2e1;
t244 = -t275 / 0.2e1;
t242 = -t269 / 0.2e1;
t91 = -pkin(4) * t124 + t165;
t240 = t156 - t281;
t163 = pkin(4) * t285 - t207 * t208;
t238 = t206 * t259;
t234 = mrSges(4,1) * t204 + mrSges(4,2) * t207;
t231 = Ifges(5,1) * t203 + t300;
t230 = -Ifges(5,2) * t203 + t300;
t229 = Ifges(5,2) * t206 + t301;
t43 = -t202 * t96 + t205 * t89;
t223 = t203 * t97 - t206 * t98;
t103 = mrSges(5,1) * t246 - mrSges(5,3) * t123;
t104 = -mrSges(5,2) * t246 + mrSges(5,3) * t124;
t222 = -t203 * t103 + t206 * t104;
t60 = t105 * t205 - t116 * t202;
t183 = t333 * t203;
t184 = t333 * t206;
t219 = t205 * t183 + t184 * t202;
t122 = t183 * t202 - t184 * t205;
t5 = t202 * t26 + t205 * t28 + t72 * t267 - t268 * t83;
t10 = t105 * t267 - t116 * t268 + t202 * t58 + t205 * t65;
t125 = pkin(4) * t365 + t194;
t2 = qJ(6) * t246 + qJD(6) * t186 + t5;
t3 = -pkin(5) * t246 - t6;
t212 = t6 * mrSges(6,1) - t3 * mrSges(7,1) - t5 * mrSges(6,2) + t2 * mrSges(7,3) + t349;
t199 = -pkin(4) * t205 - pkin(5);
t197 = pkin(4) * t202 + qJ(6);
t189 = pkin(4) * t267 + qJD(6);
t181 = -mrSges(4,3) * t201 - t287;
t171 = t203 * t259;
t169 = qJD(1) * t234;
t147 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t207 - t303) * qJD(1);
t139 = t167 * t207;
t138 = t167 * t204;
t128 = -t203 * t284 + t161;
t102 = pkin(5) * t166 - qJ(6) * t167 + t200;
t93 = Ifges(5,5) * t164 - Ifges(5,6) * t218 + t296;
t81 = -qJD(4) * t129 - t203 * t253 + t143;
t79 = pkin(5) * t139 + qJ(6) * t141 + t163;
t78 = -mrSges(5,1) * t124 + mrSges(5,2) * t123;
t76 = qJD(5) * t122 + t171 * t202 - t205 * t238;
t75 = qJD(5) * t219 + t205 * t171 + t202 * t238;
t73 = t123 * Ifges(5,4) + t124 * Ifges(5,2) + Ifges(5,6) * t246;
t71 = -t207 * t249 + (t207 * t239 - t257) * t205 + t216 * t202;
t69 = qJD(3) * t351 - t207 * t366;
t64 = mrSges(6,1) * t107 + mrSges(6,2) * t213;
t63 = mrSges(7,1) * t107 - mrSges(7,3) * t213;
t57 = -pkin(5) * t204 - t60;
t56 = qJ(6) * t204 + t289;
t49 = Ifges(7,4) * t213 + t186 * Ifges(7,2) + t107 * Ifges(7,6);
t48 = Ifges(6,5) * t213 - t107 * Ifges(6,6) + t186 * Ifges(6,3);
t46 = pkin(4) * t164 + t62;
t35 = -pkin(5) * t277 - t43;
t34 = qJ(6) * t277 + t44;
t25 = t205 * t82 - t293;
t24 = t202 * t82 + t291;
t18 = pkin(5) * t71 - qJ(6) * t69 + qJD(6) * t141 + t125;
t17 = mrSges(6,1) * t42 + mrSges(6,2) * t41;
t16 = mrSges(7,1) * t42 - mrSges(7,3) * t41;
t9 = -pkin(5) * t273 - t11;
t8 = qJ(6) * t273 + qJD(6) * t204 + t10;
t7 = pkin(5) * t42 - qJ(6) * t41 - qJD(6) * t213 + t91;
t1 = [(Ifges(6,4) * t332 + Ifges(7,5) * t331 + t320 * t359 + t328 * t360 + t364) * t69 + (t141 * t7 + t2 * t204) * mrSges(7,3) + (-t141 * t91 - t204 * t5) * mrSges(6,2) + (-Ifges(6,4) * t141 + Ifges(6,6) * t204) * t338 + (-Ifges(7,5) * t141 + Ifges(7,6) * t204) * t337 + t3 * (-mrSges(7,1) * t204 - mrSges(7,2) * t141) + t6 * (mrSges(6,1) * t204 + mrSges(6,3) * t141) + m(6) * (t10 * t22 + t11 * t21 + t117 * t125 + t163 * t91 + t289 * t5 + t6 * t60) + t289 * t32 + (Ifges(6,6) * t332 + Ifges(7,6) * t331 + t320 * t356 + t328 * t359 + t344) * t273 + (-t141 * t360 + t204 * t359) * t339 + (t169 + ((2 * mrSges(3,3)) + t234 + 0.2e1 * t352) * qJD(1)) * qJD(2) + (t260 + t349) * t204 / 0.2e1 + (-0.2e1 * t204 * t294 + (-0.3e1 / 0.2e1 * t207 ^ 2 + 0.3e1 / 0.2e1 * t204 ^ 2) * Ifges(4,4)) * t266 + ((-0.3e1 / 0.2e1 * Ifges(4,1) + 0.3e1 / 0.2e1 * Ifges(4,2)) * t204 + 0.2e1 * t295) * t246 + t156 * (mrSges(5,1) * t365 + mrSges(5,2) * t216) + ((t207 * t228 - t359 * t141 - t370 * t139 + (Ifges(5,3) + t356) * t204) * qJD(1) + t93 + t49 + t48) * t273 / 0.2e1 + (-Ifges(6,2) * t332 + Ifges(7,3) * t331 - t370 * t320 + t328 * t358 + t345) * t71 - t369 * t141 / 0.2e1 + (t91 * mrSges(6,1) + t7 * mrSges(7,1) - t2 * mrSges(7,2) - t5 * mrSges(6,3) - Ifges(6,2) * t338 + t339 * t358 + t340 + t341 + t367) * t139 - t73 * t285 / 0.2e1 + (-t207 * t78 + (t181 * t207 - t278 * t204) * qJD(3)) * t208 - t204 * t314 + m(7) * (t18 * t33 + t19 * t9 + t2 * t56 + t20 * t8 + t3 * t57 + t7 * t79) + (t224 * t275 + (qJD(4) * t223 - t53 * t203 - t54 * t206) * t207) * mrSges(5,3) - t218 * (-t229 * t269 + (Ifges(5,6) * t207 - t204 * t230) * qJD(3)) / 0.2e1 + (t204 * t245 + t206 * t242) * t94 + (t203 * t242 + t206 * t244) * t95 + qJD(3) ^ 2 * (-Ifges(4,5) * t204 - Ifges(4,6) * t207) / 0.2e1 + t163 * t17 + t125 * t64 + t80 * t126 + t81 * t127 + t128 * t103 + t129 * t104 + t8 * t85 + t10 * t86 + t11 * t87 + t9 * t88 + t79 * t16 + t18 * t63 + m(5) * (t128 * t54 + t129 * t53 + t194 * t240 + t80 * t98 + t81 * t97) + t56 * t29 + t57 * t31 + t60 * t30 + t147 * t244 + t204 * t313 + ((-Ifges(5,5) * t203 - Ifges(5,6) * t206) * t269 + (Ifges(5,3) * t207 - t204 * t228) * qJD(3)) * t319 + (-t231 * t269 + (Ifges(5,5) * t207 - t204 * t232) * qJD(3)) * t323 + (Ifges(5,6) * t204 + t207 * t230) * t326 + (Ifges(5,5) * t204 + t207 * t232) * t327 + t282 * t335 + t233 * t207 * t165; -t308 * t351 + t309 * t138 + (-t16 - t17 - t78 + (t126 * t206 - t127 * t203 + t181) * qJD(3)) * t207 + (t220 * qJD(4) + (t63 + t64 - t278) * qJD(3) + t222) * t204 - m(5) * t223 * t273 + 0.2e1 * ((-qJD(4) * t224 + t225) * t343 + (t33 * t342 + t315 / 0.2e1 + t240 * t343) * qJD(3)) * t204 + t354 * t307 + t353 * t306 + (t138 * t3 - t19 * t353 - t2 * t351 + t20 * t354 - t207 * t7) * m(7) + (-t138 * t6 - t207 * t91 + t21 * t353 + t22 * t354 - t351 * t5) * m(6) + (-t169 + (-mrSges(3,3) - t352) * qJD(1) + t348) * qJD(1); (mrSges(7,1) * t279 + mrSges(7,3) * t372) * t33 + (-t166 * t5 - t167 * t6 + t21 * t372 - t22 * t279) * mrSges(6,3) + (t47 - t50) * (t132 / 0.2e1 + t366 / 0.2e1) + (-Ifges(6,4) * t111 + Ifges(7,5) * t133 - Ifges(6,2) * t366 - Ifges(7,3) * t132) * t332 + (Ifges(6,4) * t133 - Ifges(7,5) * t111 + Ifges(6,2) * t132 + Ifges(7,3) * t366) * t331 + (-t111 * t359 - t366 * t370) * t320 + (-t111 * t360 + t358 * t366) * t328 + (-t117 * t134 + t219 * t6 + t122 * t5 + t200 * t91 + (-t44 + t75) * t22 + (-t43 - t76) * t21) * m(6) - t309 * t219 + (t230 * t274 / 0.2e1 + (t64 + t315) * t203 * pkin(4) + t348 * pkin(8) + t347) * qJD(4) + ((t147 / 0.2e1 + (t230 * t316 - Ifges(4,5) / 0.2e1) * qJD(3) + (t294 - t303 / 0.2e1) * qJD(1) + t347) * t204 + ((-t166 * t370 + t359 * t167) * qJD(3) / 0.2e1 + (-Ifges(6,5) / 0.2e1 - Ifges(7,4) / 0.2e1) * t213 + (Ifges(6,6) / 0.2e1 - Ifges(7,6) / 0.2e1) * t107 + t371 + (t230 * t318 - Ifges(4,2) / 0.2e1 + Ifges(4,1) / 0.2e1) * t201 + ((t297 / 0.2e1 + Ifges(4,4) / 0.2e1) * t207 - t295) * qJD(1) - t296 / 0.2e1 - t230 * t271 / 0.2e1 - t344 + (t245 + t324) * Ifges(5,5) - t93 / 0.2e1 - t49 / 0.2e1 - t48 / 0.2e1 + (-Ifges(6,3) / 0.2e1 - Ifges(7,2) / 0.2e1) * t186) * t207) * qJD(1) + (t166 * t358 + t167 * t360) * t339 + (-t132 * t358 + t133 * t360) * t329 + t355 * (-t133 / 0.2e1 - t111 / 0.2e1) + (t102 * t7 - t219 * t3 + t122 * t2 + t368 * t33 + (-t34 + t75) * t20 + (-t35 + t76) * t19) * m(7) + t368 * t63 + t369 * t167 / 0.2e1 + (t132 * t370 + t133 * t359) * t321 + ((-t181 - t287) * t207 + ((-mrSges(5,1) * t206 + mrSges(5,2) * t203 - mrSges(4,1)) * qJD(3) + t278) * t204) * t188 - t306 * t76 + t307 * t75 + t308 * t122 + (mrSges(6,1) * t279 - mrSges(6,2) * t372) * t117 + (-t166 * t2 + t167 * t3 - t19 * t372 - t20 * t279) * mrSges(7,2) + (-pkin(3) * t165 + pkin(8) * t225 - t113 * t97 - t114 * t98 - t156 * t172) * m(5) + t200 * t17 + t7 * (mrSges(7,1) * t166 - mrSges(7,3) * t167) + t91 * (mrSges(6,1) * t166 + mrSges(6,2) * t167) - t134 * t64 - t114 * t126 - t113 * t127 + t102 * t16 - t34 * t85 - t44 * t86 - t43 * t87 - t35 * t88 - pkin(3) * t78 + t225 * mrSges(5,3) + t222 * pkin(8) + t73 * t316 + t229 * t326 + t231 * t327 + t203 * t335 + (Ifges(7,5) * t167 + Ifges(7,3) * t166) * t337 + (Ifges(6,4) * t167 - Ifges(6,2) * t166) * t338 + t166 * t340 + t166 * t341; (-Ifges(6,2) * t331 + Ifges(7,3) * t332 - t370 * t321 + t329 * t358 - t345) * t213 + t260 - t307 * t25 + t313 - t314 + t306 * t24 + (-t19 * t24 + t197 * t2 + t199 * t3 - t33 * t46 + (t189 - t25) * t20) * m(7) + t212 - m(6) * (-t21 * t24 + t22 * t25) + t197 * t29 + t199 * t31 + t189 * t85 - t97 * t126 + t98 * t127 - t46 * t63 + (-t164 * t64 + t202 * t32 + t205 * t30 + (-t202 * t306 + t205 * t86) * qJD(5) + 0.2e1 * t315 * t324 + 0.2e1 * t19 * t268 * t342 + m(6) * (t202 * t5 + t205 * t6 - t21 * t268 + t22 * t267)) * pkin(4) + t94 * t323 + (-Ifges(5,1) * t218 - t302) * t324 + (t98 * t164 - t218 * t97) * mrSges(5,3) - t195 * (-Ifges(5,5) * t218 - Ifges(5,6) * t164) / 0.2e1 - t156 * (t164 * mrSges(5,1) - mrSges(5,2) * t218) + (-Ifges(5,2) * t164 - t158 + t95) * t218 / 0.2e1 + (-Ifges(6,4) * t331 - Ifges(7,5) * t332 - t359 * t321 - t360 * t329 + t364) * t107; (t304 + t306) * t22 + (-t305 - t307) * t21 + t212 + (t107 * t19 + t20 * t213) * mrSges(7,2) - t117 * (mrSges(6,1) * t213 - mrSges(6,2) * t107) + t50 * t328 + (Ifges(7,3) * t213 - t298) * t332 - t33 * (mrSges(7,1) * t213 + mrSges(7,3) * t107) + qJD(6) * t85 - t62 * t63 - pkin(5) * t31 + qJ(6) * t29 + (-t107 * t359 - t213 * t370) * t321 + (-pkin(5) * t3 + qJ(6) * t2 - t19 * t22 + t20 * t350 - t33 * t62) * m(7) + (-Ifges(6,2) * t213 - t101 + t355) * t331 + (-t107 * t360 + t100 - t299 + t47) * t329; t213 * t63 - t186 * t85 + 0.2e1 * (t3 / 0.2e1 + t33 * t328 + t20 * t321) * m(7) + t31;];
tauc  = t1(:);
