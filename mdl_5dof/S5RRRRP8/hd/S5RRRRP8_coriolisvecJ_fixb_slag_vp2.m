% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 22:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRRP8_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP8_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP8_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP8_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP8_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP8_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP8_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:59:40
% EndTime: 2019-12-31 22:00:07
% DurationCPUTime: 11.70s
% Computational Cost: add. (5616->496), mult. (14230->682), div. (0->0), fcn. (9345->6), ass. (0->225)
t207 = sin(qJ(3));
t208 = sin(qJ(2));
t259 = qJD(1) * t208;
t242 = t207 * t259;
t210 = cos(qJ(3));
t256 = qJD(2) * t210;
t166 = -t242 + t256;
t241 = t210 * t259;
t167 = qJD(2) * t207 + t241;
t206 = sin(qJ(4));
t209 = cos(qJ(4));
t114 = t166 * t206 + t167 * t209;
t182 = -qJD(2) * pkin(2) + pkin(6) * t259;
t135 = -pkin(3) * t166 + t182;
t211 = cos(qJ(2));
t258 = qJD(1) * t211;
t196 = qJD(3) - t258;
t178 = -pkin(2) * t211 - t208 * pkin(7) - pkin(1);
t158 = t178 * qJD(1);
t203 = pkin(6) * t258;
t183 = qJD(2) * pkin(7) + t203;
t117 = t210 * t158 - t183 * t207;
t86 = -pkin(8) * t167 + t117;
t78 = pkin(3) * t196 + t86;
t118 = t158 * t207 + t183 * t210;
t87 = pkin(8) * t166 + t118;
t81 = t206 * t87;
t23 = t209 * t78 - t81;
t347 = qJ(5) * t114;
t18 = t23 - t347;
t186 = qJD(4) + t196;
t17 = pkin(4) * t186 + t18;
t83 = t209 * t87;
t24 = t206 * t78 + t83;
t233 = t209 * t166 - t167 * t206;
t325 = qJ(5) * t233;
t19 = t24 + t325;
t293 = -t186 / 0.2e1;
t303 = t114 / 0.2e1;
t304 = -t114 / 0.2e1;
t257 = qJD(2) * t208;
t234 = qJD(1) * t257;
t250 = qJD(2) * qJD(3);
t254 = qJD(3) * t207;
t255 = qJD(2) * t211;
t129 = t210 * t250 + (-t208 * t254 + t210 * t255) * qJD(1);
t253 = qJD(3) * t210;
t345 = t207 * t255 + t208 * t253;
t130 = -qJD(1) * t345 - t207 * t250;
t44 = qJD(4) * t233 + t129 * t209 + t130 * t206;
t230 = pkin(2) * t208 - pkin(7) * t211;
t176 = t230 * qJD(2);
t159 = qJD(1) * t176;
t232 = pkin(6) * t234;
t68 = -qJD(3) * t118 + t210 * t159 + t207 * t232;
t37 = pkin(3) * t234 - pkin(8) * t129 + t68;
t67 = t158 * t253 + t207 * t159 - t183 * t254 - t210 * t232;
t47 = pkin(8) * t130 + t67;
t6 = -qJD(4) * t24 - t206 * t47 + t209 * t37;
t2 = pkin(4) * t234 - qJ(5) * t44 - qJD(5) * t114 + t6;
t45 = -qJD(4) * t114 - t129 * t206 + t130 * t209;
t251 = qJD(4) * t209;
t252 = qJD(4) * t206;
t5 = t206 * t37 + t209 * t47 + t78 * t251 - t252 * t87;
t3 = qJ(5) * t45 + qJD(5) * t233 + t5;
t319 = t6 * mrSges(5,1) + t2 * mrSges(6,1) - t5 * mrSges(5,2) - t3 * mrSges(6,2);
t335 = -Ifges(6,6) - Ifges(5,6);
t337 = -Ifges(6,5) - Ifges(5,5);
t359 = Ifges(6,3) + Ifges(5,3);
t323 = t234 * t359 - t335 * t45 - t337 * t44;
t336 = Ifges(5,2) + Ifges(6,2);
t339 = Ifges(5,1) + Ifges(6,1);
t340 = -t233 / 0.2e1;
t338 = Ifges(5,4) + Ifges(6,4);
t358 = t338 * t233;
t350 = t339 * t114 - t337 * t186 + t358;
t357 = t114 * t338;
t351 = -t335 * t186 + t336 * t233 + t357;
t73 = -pkin(4) * t233 + qJD(5) + t135;
t365 = t319 + t323 + (t114 * t335 - t337 * t233) * t293 + (t114 * t19 + t17 * t233) * mrSges(6,3) + (t114 * t24 + t23 * t233) * mrSges(5,3) - t135 * (mrSges(5,1) * t114 + mrSges(5,2) * t233) - t73 * (mrSges(6,1) * t114 + mrSges(6,2) * t233) + t351 * t303 + (-t114 * t336 + t350 + t358) * t340 + (t339 * t233 - t357) * t304;
t173 = t230 * qJD(1);
t131 = pkin(6) * t242 + t210 * t173;
t263 = t210 * t211;
t219 = pkin(3) * t208 - pkin(8) * t263;
t313 = -pkin(8) - pkin(7);
t243 = qJD(3) * t313;
t364 = -qJD(1) * t219 + t210 * t243 - t131;
t152 = t207 * t173;
t264 = t208 * t210;
t265 = t207 * t211;
t363 = t152 + (-pkin(6) * t264 - pkin(8) * t265) * qJD(1) - t207 * t243;
t220 = t206 * t207 - t209 * t210;
t320 = qJD(3) + qJD(4);
t121 = t320 * t220;
t217 = t220 * t211;
t140 = qJD(1) * t217;
t362 = t121 - t140;
t169 = t206 * t210 + t207 * t209;
t122 = t320 * t169;
t218 = t211 * t169;
t139 = qJD(1) * t218;
t361 = t122 - t139;
t238 = Ifges(3,5) * qJD(2) / 0.2e1;
t360 = -t336 * t45 / 0.2e1 - t338 * t44 / 0.2e1 + t335 * t234 / 0.2e1;
t184 = t313 * t207;
t185 = t313 * t210;
t128 = t206 * t184 - t209 * t185;
t349 = -qJD(4) * t128 + t206 * t363 + t209 * t364;
t348 = t184 * t251 + t185 * t252 + t206 * t364 - t209 * t363;
t354 = -t337 * t234 + t338 * t45 + t339 * t44;
t353 = -pkin(4) * t259 + t362 * qJ(5) - qJD(5) * t169 + t349;
t352 = -t361 * qJ(5) - qJD(5) * t220 + t348;
t286 = pkin(3) * t207;
t161 = t258 * t286 + t203;
t346 = pkin(3) * t254 + t361 * pkin(4) - t161;
t201 = Ifges(3,4) * t258;
t221 = t117 * t210 + t118 * t207;
t279 = Ifges(4,4) * t210;
t225 = -Ifges(4,2) * t207 + t279;
t280 = Ifges(4,4) * t207;
t227 = Ifges(4,1) * t210 - t280;
t228 = mrSges(4,1) * t207 + mrSges(4,2) * t210;
t275 = Ifges(4,6) * t207;
t276 = Ifges(4,5) * t210;
t289 = t210 / 0.2e1;
t290 = -t207 / 0.2e1;
t296 = t167 / 0.2e1;
t273 = t167 * Ifges(4,4);
t97 = t166 * Ifges(4,2) + t196 * Ifges(4,6) + t273;
t162 = Ifges(4,4) * t166;
t98 = Ifges(4,1) * t167 + Ifges(4,5) * t196 + t162;
t212 = -t221 * mrSges(4,3) + t182 * t228 + t166 * t225 / 0.2e1 + t227 * t296 + t196 * (-t275 + t276) / 0.2e1 + t97 * t290 + t98 * t289;
t343 = t212 + Ifges(3,1) * t259 / 0.2e1 + t201 / 0.2e1 + t238;
t292 = t186 / 0.2e1;
t306 = t233 / 0.2e1;
t237 = -Ifges(3,6) * qJD(2) / 0.2e1;
t61 = -mrSges(5,1) * t233 + mrSges(5,2) * t114;
t326 = m(5) * t135 + t61;
t147 = t220 * t208;
t165 = t210 * t178;
t285 = pkin(6) * t207;
t116 = -pkin(8) * t264 + t165 + (-pkin(3) - t285) * t211;
t198 = pkin(6) * t263;
t138 = t207 * t178 + t198;
t266 = t207 * t208;
t124 = -pkin(8) * t266 + t138;
t64 = t206 * t116 + t209 * t124;
t324 = Ifges(4,5) * t129 + Ifges(4,6) * t130;
t322 = -t68 * mrSges(4,1) + t67 * mrSges(4,2);
t321 = pkin(1) * mrSges(3,2) * qJD(1);
t245 = Ifges(5,3) / 0.2e1 + Ifges(6,3) / 0.2e1;
t246 = Ifges(6,6) / 0.2e1 + Ifges(5,6) / 0.2e1;
t247 = Ifges(6,5) / 0.2e1 + Ifges(5,5) / 0.2e1;
t281 = Ifges(3,4) * t208;
t318 = t247 * t114 + t245 * t186 + t246 * t233 + t117 * mrSges(4,1) + t17 * mrSges(6,1) + t23 * mrSges(5,1) + t237 - (t211 * Ifges(3,2) + t281) * qJD(1) / 0.2e1 + t196 * Ifges(4,3) + t167 * Ifges(4,5) + t166 * Ifges(4,6) - t118 * mrSges(4,2) - t19 * mrSges(6,2) - t24 * mrSges(5,2) - t335 * t306 - t337 * t303 + t359 * t292;
t317 = t44 / 0.2e1;
t316 = t45 / 0.2e1;
t312 = pkin(1) * mrSges(3,1);
t302 = t129 / 0.2e1;
t301 = t130 / 0.2e1;
t298 = -t166 / 0.2e1;
t297 = -t167 / 0.2e1;
t291 = -t196 / 0.2e1;
t34 = t209 * t86 - t81;
t268 = qJD(2) * mrSges(3,2);
t260 = t210 * t176 + t257 * t285;
t177 = pkin(3) * t266 + t208 * pkin(6);
t244 = Ifges(4,3) * t234 + t324;
t136 = pkin(3) * t345 + pkin(6) * t255;
t200 = -pkin(3) * t210 - pkin(2);
t15 = -t45 * mrSges(6,1) + t44 * mrSges(6,2);
t110 = -pkin(3) * t130 + qJD(2) * t203;
t33 = -t206 * t86 - t83;
t63 = t209 * t116 - t206 * t124;
t127 = t209 * t184 + t185 * t206;
t231 = m(4) * t182 - qJD(2) * mrSges(3,1) - mrSges(4,1) * t166 + mrSges(4,2) * t167 + mrSges(3,3) * t259;
t229 = mrSges(4,1) * t210 - mrSges(4,2) * t207;
t226 = Ifges(4,1) * t207 + t279;
t224 = Ifges(4,2) * t210 + t280;
t223 = Ifges(4,5) * t207 + Ifges(4,6) * t210;
t222 = -t207 * t68 + t210 * t67;
t62 = t219 * qJD(2) + (-t198 + (pkin(8) * t208 - t178) * t207) * qJD(3) + t260;
t84 = t207 * t176 + t178 * t253 + (-t208 * t256 - t211 * t254) * pkin(6);
t66 = -pkin(8) * t345 + t84;
t9 = t116 * t251 - t124 * t252 + t206 * t62 + t209 * t66;
t10 = -qJD(4) * t64 - t206 * t66 + t209 * t62;
t199 = pkin(3) * t209 + pkin(4);
t180 = mrSges(3,3) * t258 - t268;
t146 = t169 * t208;
t141 = pkin(4) * t220 + t200;
t137 = -pkin(6) * t265 + t165;
t134 = mrSges(4,1) * t196 - mrSges(4,3) * t167;
t133 = -mrSges(4,2) * t196 + mrSges(4,3) * t166;
t132 = -pkin(6) * t241 + t152;
t119 = pkin(4) * t146 + t177;
t109 = -mrSges(4,2) * t234 + mrSges(4,3) * t130;
t108 = mrSges(4,1) * t234 - mrSges(4,3) * t129;
t95 = -qJ(5) * t220 + t128;
t94 = -qJ(5) * t169 + t127;
t91 = mrSges(5,1) * t186 - mrSges(5,3) * t114;
t90 = mrSges(6,1) * t186 - mrSges(6,3) * t114;
t89 = -mrSges(5,2) * t186 + mrSges(5,3) * t233;
t88 = -mrSges(6,2) * t186 + mrSges(6,3) * t233;
t85 = -qJD(3) * t138 + t260;
t80 = pkin(3) * t167 + pkin(4) * t114;
t79 = -mrSges(4,1) * t130 + mrSges(4,2) * t129;
t72 = t129 * Ifges(4,1) + t130 * Ifges(4,4) + Ifges(4,5) * t234;
t71 = t129 * Ifges(4,4) + t130 * Ifges(4,2) + Ifges(4,6) * t234;
t70 = -qJD(2) * t218 + t147 * t320;
t69 = -qJD(2) * t217 - t122 * t208;
t60 = -mrSges(6,1) * t233 + mrSges(6,2) * t114;
t49 = -pkin(4) * t70 + t136;
t48 = -qJ(5) * t146 + t64;
t46 = -pkin(4) * t211 + t147 * qJ(5) + t63;
t31 = -mrSges(5,2) * t234 + mrSges(5,3) * t45;
t30 = -mrSges(6,2) * t234 + mrSges(6,3) * t45;
t29 = mrSges(5,1) * t234 - mrSges(5,3) * t44;
t28 = mrSges(6,1) * t234 - mrSges(6,3) * t44;
t22 = -pkin(4) * t45 + t110;
t21 = t34 - t347;
t20 = t33 - t325;
t16 = -mrSges(5,1) * t45 + mrSges(5,2) * t44;
t8 = qJ(5) * t70 - qJD(5) * t146 + t9;
t7 = pkin(4) * t257 - qJ(5) * t69 + qJD(5) * t147 + t10;
t1 = [(t337 * t317 + t335 * t316 + (0.3e1 / 0.2e1 * Ifges(3,4) * t255 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2) - Ifges(4,3) / 0.2e1 + (m(4) * pkin(6) + t228) * pkin(6) - t245) * t257) * qJD(1) - t319 + t322) * t211 + (-t335 * t70 - t337 * t69) * t292 + (-t146 * t336 - t147 * t338) * t316 - (t244 + t323 + t324) * t211 / 0.2e1 + (pkin(6) * t79 + t71 * t290 + t72 * t289 + t227 * t302 + t225 * t301 + (-t207 * t67 - t210 * t68) * mrSges(4,3) + (t223 * t291 + t182 * t229 + t224 * t298 + t226 * t297 - t210 * t97 / 0.2e1 + t98 * t290 + (t117 * t207 - t118 * t210) * mrSges(4,3)) * qJD(3) + (((-0.3e1 / 0.2e1 * Ifges(3,4) + t276 / 0.2e1 - t275 / 0.2e1) * t208 - 0.2e1 * t312 - t247 * t147 - t246 * t146) * qJD(1) - pkin(6) * t180 + t237 + t318) * qJD(2)) * t208 + m(6) * (t119 * t22 + t17 * t7 + t19 * t8 + t2 * t46 + t3 * t48 + t49 * t73) + m(5) * (t10 * t23 + t110 * t177 + t135 * t136 + t24 * t9 + t5 * t64 + t6 * t63) + (-t146 * t3 + t147 * t2 - t17 * t69 + t19 * t70) * mrSges(6,3) + (-t146 * t5 + t147 * t6 - t23 * t69 + t24 * t70) * mrSges(5,3) + t22 * (mrSges(6,1) * t146 - mrSges(6,2) * t147) + t110 * (mrSges(5,1) * t146 - mrSges(5,2) * t147) + t350 * t69 / 0.2e1 + t351 * t70 / 0.2e1 - t354 * t147 / 0.2e1 + t146 * t360 + m(4) * (t117 * t85 + t118 * t84 + t68 * t137 + t67 * t138) + (t336 * t70 + t338 * t69) * t306 + (-t146 * t338 - t147 * t339) * t317 + (t338 * t70 + t339 * t69) * t303 + (t231 * pkin(6) + t238 - 0.2e1 * t321 + t343) * t255 + t177 * t16 + t46 * t28 + t48 * t30 + t49 * t60 + t63 * t29 + t64 * t31 + t73 * (-mrSges(6,1) * t70 + mrSges(6,2) * t69) + t8 * t88 + t9 * t89 + t7 * t90 + t10 * t91 + t119 * t15 + t84 * t133 + t85 * t134 + t135 * (-mrSges(5,1) * t70 + mrSges(5,2) * t69) + t136 * t61 + t137 * t108 + t138 * t109; -m(4) * (t117 * t131 + t118 * t132) + (t326 * t286 + t212) * qJD(3) + t348 * t89 + (t110 * t200 + t127 * t6 + t128 * t5 - t135 * t161 + t23 * t349 + t24 * t348) * m(5) + t349 * t91 + t350 * (-t121 / 0.2e1 + t140 / 0.2e1) + t351 * (-t122 / 0.2e1 + t139 / 0.2e1) + t352 * t88 + (t141 * t22 + t17 * t353 + t19 * t352 + t2 * t94 + t3 * t95 + t346 * t73) * m(6) + t353 * t90 + t354 * t169 / 0.2e1 + t220 * t360 + t346 * t60 + (t121 * t337 + t122 * t335) * t292 + (t139 * t335 + t140 * t337) * t293 + (-t121 * t338 - t122 * t336) * t306 + (-t139 * t336 - t140 * t338) * t340 + (t169 * t338 - t220 * t336) * t316 + (-t121 * t339 - t122 * t338) * t303 + (-t139 * t338 - t140 * t339) * t304 + (t169 * t339 - t220 * t338) * t317 + t222 * mrSges(4,3) + ((t238 - t201 / 0.2e1 + t321 + ((-m(4) * pkin(2) - mrSges(3,1) - t229) * qJD(2) - t231) * pkin(6) - t343) * t211 + ((t312 + t281 / 0.2e1 + (Ifges(3,2) / 0.2e1 - Ifges(3,1) / 0.2e1) * t211) * qJD(1) + (t180 + t268) * pkin(6) + t237 - t318) * t208 + (-t169 * t337 + t220 * t335 + t223) * t257 / 0.2e1) * qJD(1) + (-t169 * t2 + t17 * t362 - t19 * t361 - t220 * t3) * mrSges(6,3) + (mrSges(5,1) * t361 - mrSges(5,2) * t362) * t135 + (mrSges(6,1) * t361 - mrSges(6,2) * t362) * t73 + (-t169 * t6 - t220 * t5 + t23 * t362 - t24 * t361) * mrSges(5,3) + t207 * t72 / 0.2e1 + t200 * t16 + t71 * t289 + t224 * t301 + t226 * t302 + t22 * (mrSges(6,1) * t220 + mrSges(6,2) * t169) + t110 * (mrSges(5,1) * t220 + mrSges(5,2) * t169) + (-t108 * t207 + t109 * t210 + m(4) * t222 + (-m(4) * t221 - t207 * t133 - t210 * t134) * qJD(3)) * pkin(7) - pkin(2) * t79 + t94 * t28 + t95 * t30 + t127 * t29 + t128 * t31 - t132 * t133 - t131 * t134 + t141 * t15 - t161 * t61; -t322 + (-Ifges(4,2) * t167 + t162 + t98) * t298 + (t117 * t166 + t118 * t167) * mrSges(4,3) + t244 - m(5) * (t23 * t33 + t24 * t34) + t199 * t28 - t182 * (mrSges(4,1) * t167 + mrSges(4,2) * t166) + (Ifges(4,5) * t166 - Ifges(4,6) * t167) * t291 + t97 * t296 + (Ifges(4,1) * t166 - t273) * t297 + (t209 * t29 + (t30 + t31) * t206 + ((t88 + t89) * t209 + (-t90 - t91) * t206) * qJD(4) + m(5) * (t206 * t5 + t209 * t6 - t23 * t252 + t24 * t251) - t326 * t167) * pkin(3) + (-t17 * t20 - t19 * t21 - t73 * t80 + t199 * t2 + (-t17 * t252 + t19 * t251 + t206 * t3) * pkin(3)) * m(6) - t80 * t60 - t21 * t88 - t34 * t89 - t20 * t90 - t33 * t91 - t117 * t133 + t118 * t134 + t365; (-(-t17 + t18) * t19 + (-t114 * t73 + t2) * pkin(4)) * m(6) + (-t114 * t60 + t28) * pkin(4) - t18 * t88 - t23 * t89 + t19 * t90 + t24 * t91 + t365; -t233 * t88 + t114 * t90 + 0.2e1 * (t22 / 0.2e1 + t19 * t340 + t17 * t303) * m(6) + t15;];
tauc = t1(:);
