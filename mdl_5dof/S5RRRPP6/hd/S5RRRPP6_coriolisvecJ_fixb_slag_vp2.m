% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRPP6_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP6_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP6_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP6_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP6_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP6_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP6_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:59:59
% EndTime: 2019-12-31 21:00:20
% DurationCPUTime: 9.57s
% Computational Cost: add. (3987->499), mult. (10459->675), div. (0->0), fcn. (6675->6), ass. (0->232)
t192 = cos(qJ(3));
t235 = qJD(2) * qJD(3);
t190 = sin(qJ(3));
t191 = sin(qJ(2));
t239 = qJD(3) * t191;
t193 = cos(qJ(2));
t241 = qJD(2) * t193;
t314 = t190 * t239 - t192 * t241;
t114 = -qJD(1) * t314 + t192 * t235;
t238 = qJD(3) * t192;
t327 = t190 * t241 + t191 * t238;
t115 = -qJD(1) * t327 - t190 * t235;
t189 = sin(pkin(8));
t257 = cos(pkin(8));
t67 = t114 * t189 - t115 * t257;
t303 = -t67 / 0.2e1;
t68 = t114 * t257 + t189 * t115;
t301 = t68 / 0.2e1;
t322 = Ifges(5,1) + Ifges(6,1);
t321 = Ifges(5,4) - Ifges(6,5);
t320 = Ifges(6,4) + Ifges(5,5);
t216 = t257 * t190;
t149 = t189 * t192 + t216;
t236 = t193 * qJD(1);
t121 = t149 * t236;
t135 = t149 * qJD(3);
t248 = t121 - t135;
t215 = t257 * t192;
t253 = t189 * t190;
t199 = t215 - t253;
t197 = t193 * t199;
t122 = qJD(1) * t197;
t136 = t199 * qJD(3);
t247 = t122 - t136;
t222 = Ifges(3,5) * qJD(2) / 0.2e1;
t243 = qJD(2) * t191;
t218 = qJD(1) * t243;
t335 = t322 * t301 + t321 * t303 + t320 * t218 / 0.2e1;
t319 = Ifges(5,6) - Ifges(6,6);
t211 = pkin(2) * t191 - pkin(7) * t193;
t155 = t211 * qJD(1);
t139 = t190 * t155;
t250 = t191 * t192;
t251 = t190 * t193;
t101 = t139 + (-pkin(6) * t250 - qJ(4) * t251) * qJD(1);
t244 = qJD(1) * t191;
t228 = t190 * t244;
t116 = pkin(6) * t228 + t192 * t155;
t249 = t192 * t193;
t201 = pkin(3) * t191 - qJ(4) * t249;
t87 = qJD(1) * t201 + t116;
t38 = -t189 * t101 + t257 * t87;
t269 = -qJ(4) - pkin(7);
t217 = qJD(3) * t269;
t237 = qJD(4) * t192;
t132 = t190 * t217 + t237;
t198 = -qJD(4) * t190 + t192 * t217;
t84 = t132 * t189 - t198 * t257;
t334 = -t84 - t38;
t39 = t257 * t101 + t189 * t87;
t85 = t132 * t257 + t189 * t198;
t333 = t85 - t39;
t332 = Ifges(5,3) + Ifges(6,2) + Ifges(4,3);
t177 = qJD(3) - t236;
t242 = qJD(2) * t192;
t153 = -t228 + t242;
t227 = t192 * t244;
t154 = qJD(2) * t190 + t227;
t200 = t189 * t153 + t154 * t257;
t98 = -t257 * t153 + t154 * t189;
t318 = t320 * t177 + t322 * t200 - t321 * t98;
t186 = pkin(6) * t236;
t274 = pkin(3) * t190;
t146 = t236 * t274 + t186;
t240 = qJD(3) * t190;
t330 = pkin(3) * t240 - t248 * pkin(4) + t247 * qJ(5) - qJD(5) * t149 - t146;
t329 = pkin(4) * t244 - t334;
t328 = -qJ(5) * t244 + t333;
t184 = Ifges(3,4) * t236;
t166 = -qJD(2) * pkin(2) + pkin(6) * t244;
t161 = -pkin(2) * t193 - t191 * pkin(7) - pkin(1);
t143 = t161 * qJD(1);
t167 = qJD(2) * pkin(7) + t186;
t104 = t192 * t143 - t167 * t190;
t105 = t143 * t190 + t167 * t192;
t202 = t104 * t192 + t105 * t190;
t263 = Ifges(4,4) * t192;
t206 = -Ifges(4,2) * t190 + t263;
t264 = Ifges(4,4) * t190;
t208 = Ifges(4,1) * t192 - t264;
t209 = mrSges(4,1) * t190 + mrSges(4,2) * t192;
t261 = Ifges(4,6) * t190;
t262 = Ifges(4,5) * t192;
t278 = t192 / 0.2e1;
t279 = -t190 / 0.2e1;
t280 = t177 / 0.2e1;
t282 = t154 / 0.2e1;
t265 = Ifges(4,4) * t154;
t89 = Ifges(4,2) * t153 + Ifges(4,6) * t177 + t265;
t147 = Ifges(4,4) * t153;
t90 = t154 * Ifges(4,1) + t177 * Ifges(4,5) + t147;
t194 = -t202 * mrSges(4,3) + t89 * t279 + t90 * t278 + t153 * t206 / 0.2e1 + t208 * t282 + t166 * t209 + (-t261 + t262) * t280;
t326 = t194 + Ifges(3,1) * t244 / 0.2e1 + t184 / 0.2e1 + t222;
t113 = -pkin(3) * t153 + qJD(4) + t166;
t77 = qJ(4) * t153 + t105;
t258 = t189 * t77;
t76 = -qJ(4) * t154 + t104;
t69 = pkin(3) * t177 + t76;
t21 = t257 * t69 - t258;
t15 = -t177 * pkin(4) + qJD(5) - t21;
t27 = pkin(4) * t98 - qJ(5) * t200 + t113;
t325 = t113 * mrSges(5,2) + t15 * mrSges(6,2) - t21 * mrSges(5,3) - t27 * mrSges(6,3) + t318 / 0.2e1;
t324 = -Ifges(5,6) / 0.2e1;
t323 = Ifges(6,6) / 0.2e1;
t302 = t67 / 0.2e1;
t299 = -t98 / 0.2e1;
t298 = t98 / 0.2e1;
t289 = t200 / 0.2e1;
t221 = -Ifges(3,6) * qJD(2) / 0.2e1;
t78 = -mrSges(5,2) * t177 - mrSges(5,3) * t98;
t81 = -mrSges(6,2) * t98 + mrSges(6,3) * t177;
t317 = -t78 - t81;
t79 = mrSges(5,1) * t177 - mrSges(5,3) * t200;
t80 = -mrSges(6,1) * t177 + mrSges(6,2) * t200;
t316 = t79 - t80;
t315 = Ifges(4,5) * t114 + Ifges(4,6) * t115;
t179 = pkin(6) * t249;
t124 = t190 * t161 + t179;
t313 = qJD(1) * pkin(1) * mrSges(3,2);
t312 = t218 * t332 - t319 * t67 + t320 * t68 + t315;
t156 = t211 * qJD(2);
t144 = qJD(1) * t156;
t214 = pkin(6) * t218;
t54 = -qJD(3) * t105 + t192 * t144 + t190 * t214;
t20 = pkin(3) * t218 - qJ(4) * t114 - qJD(4) * t154 + t54;
t53 = t143 * t238 + t190 * t144 - t167 * t240 - t192 * t214;
t24 = qJ(4) * t115 + qJD(4) * t153 + t53;
t4 = t189 * t20 + t257 * t24;
t1 = qJ(5) * t218 + qJD(5) * t177 + t4;
t3 = -t189 * t24 + t20 * t257;
t2 = -pkin(4) * t218 - t3;
t311 = -t54 * mrSges(4,1) - t3 * mrSges(5,1) + t2 * mrSges(6,1) + t53 * mrSges(4,2) + t4 * mrSges(5,2) - t1 * mrSges(6,3);
t72 = t257 * t77;
t22 = t189 * t69 + t72;
t16 = qJ(5) * t177 + t22;
t213 = Ifges(4,3) / 0.2e1 + Ifges(5,3) / 0.2e1 + Ifges(6,2) / 0.2e1;
t231 = t323 + t324;
t232 = Ifges(5,5) / 0.2e1 + Ifges(6,4) / 0.2e1;
t266 = Ifges(3,4) * t191;
t310 = t213 * t177 + t232 * t200 + t231 * t98 + t104 * mrSges(4,1) + t16 * mrSges(6,3) + t21 * mrSges(5,1) + t221 - (t193 * Ifges(3,2) + t266) * qJD(1) / 0.2e1 + Ifges(5,6) * t299 + Ifges(6,6) * t298 + t154 * Ifges(4,5) + t153 * Ifges(4,6) - t105 * mrSges(4,2) - t15 * mrSges(6,1) - t22 * mrSges(5,2) + t320 * t289 + t332 * t280;
t32 = Ifges(6,5) * t200 + Ifges(6,6) * t177 + Ifges(6,3) * t98;
t35 = Ifges(5,4) * t200 - Ifges(5,2) * t98 + Ifges(5,6) * t177;
t308 = -t113 * mrSges(5,1) - t27 * mrSges(6,1) + t16 * mrSges(6,2) + t22 * mrSges(5,3) - t32 / 0.2e1 + t35 / 0.2e1;
t307 = Ifges(6,5) * t301 + Ifges(6,3) * t302 + t218 * t323;
t306 = -t68 * Ifges(5,4) / 0.2e1 + Ifges(5,2) * t302 + t218 * t324;
t296 = pkin(1) * mrSges(3,1);
t290 = -t200 / 0.2e1;
t288 = t114 / 0.2e1;
t287 = t115 / 0.2e1;
t284 = -t153 / 0.2e1;
t283 = -t154 / 0.2e1;
t281 = -t177 / 0.2e1;
t276 = pkin(3) * t154;
t275 = pkin(3) * t189;
t273 = pkin(6) * t190;
t245 = t192 * t156 + t243 * t273;
t40 = -t191 * t237 + t201 * qJD(2) + (-t179 + (qJ(4) * t191 - t161) * t190) * qJD(3) + t245;
t246 = t190 * t156 + t161 * t238;
t44 = (-qJD(2) * pkin(6) - qJ(4) * qJD(3)) * t250 + (-qJD(4) * t191 + (-pkin(6) * qJD(3) - qJ(4) * qJD(2)) * t193) * t190 + t246;
t9 = t189 * t40 + t257 * t44;
t268 = mrSges(4,3) * t153;
t267 = mrSges(4,3) * t154;
t151 = t192 * t161;
t102 = -qJ(4) * t250 + t151 + (-pkin(3) - t273) * t193;
t252 = t190 * t191;
t107 = -qJ(4) * t252 + t124;
t48 = t189 * t102 + t257 * t107;
t254 = qJD(2) * mrSges(3,2);
t157 = pkin(3) * t252 + t191 * pkin(6);
t120 = pkin(3) * t327 + pkin(6) * t241;
t183 = -pkin(3) * t192 - pkin(2);
t229 = t257 * pkin(3);
t18 = t67 * mrSges(5,1) + t68 * mrSges(5,2);
t17 = t67 * mrSges(6,1) - t68 * mrSges(6,3);
t96 = -pkin(3) * t115 + qJD(2) * t186;
t212 = m(4) * t166 - qJD(2) * mrSges(3,1) - mrSges(4,1) * t153 + mrSges(4,2) * t154 + mrSges(3,3) * t244;
t210 = mrSges(4,1) * t192 - mrSges(4,2) * t190;
t207 = Ifges(4,1) * t190 + t263;
t205 = Ifges(4,2) * t192 + t264;
t204 = Ifges(4,5) * t190 + Ifges(4,6) * t192;
t203 = -t190 * t54 + t192 * t53;
t52 = -mrSges(6,1) * t218 + t68 * mrSges(6,2);
t8 = -t189 * t44 + t257 * t40;
t47 = t102 * t257 - t189 * t107;
t182 = -t229 - pkin(4);
t180 = qJ(5) + t275;
t165 = t269 * t192;
t163 = mrSges(3,3) * t236 - t254;
t128 = t199 * t191;
t127 = t149 * t191;
t123 = -pkin(6) * t251 + t151;
t119 = mrSges(4,1) * t177 - t267;
t118 = -mrSges(4,2) * t177 + t268;
t117 = -pkin(6) * t227 + t139;
t110 = -t165 * t257 + t253 * t269;
t109 = -t165 * t189 - t216 * t269;
t95 = -mrSges(4,2) * t218 + mrSges(4,3) * t115;
t94 = mrSges(4,1) * t218 - mrSges(4,3) * t114;
t92 = -pkin(4) * t199 - qJ(5) * t149 + t183;
t83 = qJD(2) * t197 - t149 * t239;
t82 = t189 * t314 - t215 * t239 - t216 * t241;
t75 = -qJD(3) * t124 + t245;
t74 = (-t191 * t242 - t193 * t240) * pkin(6) + t246;
t71 = pkin(4) * t127 - qJ(5) * t128 + t157;
t70 = -mrSges(4,1) * t115 + mrSges(4,2) * t114;
t57 = t114 * Ifges(4,1) + t115 * Ifges(4,4) + Ifges(4,5) * t218;
t56 = t114 * Ifges(4,4) + t115 * Ifges(4,2) + Ifges(4,6) * t218;
t51 = mrSges(5,1) * t218 - mrSges(5,3) * t68;
t50 = -mrSges(5,2) * t218 - mrSges(5,3) * t67;
t49 = -mrSges(6,2) * t67 + mrSges(6,3) * t218;
t46 = mrSges(5,1) * t98 + mrSges(5,2) * t200;
t45 = mrSges(6,1) * t98 - mrSges(6,3) * t200;
t43 = t193 * pkin(4) - t47;
t42 = -qJ(5) * t193 + t48;
t28 = pkin(4) * t200 + qJ(5) * t98 + t276;
t26 = t257 * t76 - t258;
t25 = t189 * t76 + t72;
t10 = -pkin(4) * t82 - qJ(5) * t83 - qJD(5) * t128 + t120;
t7 = -pkin(4) * t243 - t8;
t6 = qJ(5) * t243 - qJD(5) * t193 + t9;
t5 = pkin(4) * t67 - qJ(5) * t68 - qJD(5) * t200 + t96;
t11 = [t9 * t78 + t8 * t79 + t7 * t80 + t6 * t81 + t71 * t17 + t10 * t45 + t42 * t49 + t48 * t50 + t47 * t51 + t43 * t52 + (-t1 * t127 + t128 * t2) * mrSges(6,2) - (t312 + t315) * t193 / 0.2e1 + (Ifges(5,2) * t299 - Ifges(6,3) * t298 + t280 * t319 + t289 * t321 + t308) * t82 + (-t127 * t321 + t128 * t322) * t301 + m(4) * (t104 * t75 + t105 * t74 + t54 * t123 + t53 * t124) + (-t127 * t4 - t128 * t3) * mrSges(5,3) + (Ifges(6,5) * t128 + Ifges(6,3) * t127) * t302 + (-Ifges(5,6) * t303 - Ifges(6,6) * t302 - t320 * t301 + (0.3e1 / 0.2e1 * Ifges(3,4) * t241 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2) + (m(4) * pkin(6) + t209) * pkin(6) - t213) * t243) * qJD(1) + t311) * t193 + t128 * t335 + m(5) * (t113 * t120 + t157 * t96 + t21 * t8 + t22 * t9 + t3 * t47 + t4 * t48) + m(6) * (t1 * t42 + t10 * t27 + t15 * t7 + t16 * t6 + t2 * t43 + t5 * t71) + (Ifges(5,4) * t299 + Ifges(6,5) * t298 + t320 * t280 + t322 * t289 + t325) * t83 + (t212 * pkin(6) + t222 - 0.2e1 * t313 + t326) * t241 + t74 * t118 + t75 * t119 + t120 * t46 + t123 * t94 + t124 * t95 + t5 * (mrSges(6,1) * t127 - mrSges(6,3) * t128) + t96 * (mrSges(5,1) * t127 + mrSges(5,2) * t128) + t157 * t18 + t127 * t306 + t127 * t307 + (t208 * t288 + t206 * t287 + pkin(6) * t70 + t56 * t279 + t57 * t278 + (-t190 * t53 - t192 * t54) * mrSges(4,3) + (t205 * t284 + t207 * t283 + t166 * t210 + t204 * t281 + t90 * t279 - t192 * t89 / 0.2e1 + (t104 * t190 - t105 * t192) * mrSges(4,3)) * qJD(3) + (-pkin(6) * t163 + ((-0.3e1 / 0.2e1 * Ifges(3,4) + t262 / 0.2e1 - t261 / 0.2e1) * t191 - 0.2e1 * t296 + t232 * t128 + t231 * t127) * qJD(1) + t221 + t310) * qJD(2)) * t191 + (Ifges(5,4) * t128 - Ifges(5,2) * t127) * t303; t92 * t17 - pkin(2) * t70 + t333 * t78 + t334 * t79 + (-t190 * t94 + t192 * t95 + m(4) * t203 + (-m(4) * t202 - t190 * t118 - t192 * t119) * qJD(3)) * pkin(7) + (t49 + t50) * t110 + (Ifges(6,5) * t149 - Ifges(6,3) * t199) * t302 + (Ifges(5,4) * t149 + Ifges(5,2) * t199) * t303 + t5 * (-mrSges(6,1) * t199 - mrSges(6,3) * t149) + t96 * (-mrSges(5,1) * t199 + mrSges(5,2) * t149) + (-t149 * t3 + t199 * t4 + t21 * t247 + t22 * t248) * mrSges(5,3) + (t1 * t199 + t149 * t2 - t15 * t247 + t16 * t248) * mrSges(6,2) - t199 * t306 - t199 * t307 + (Ifges(5,4) * t122 + Ifges(6,5) * t136 - Ifges(5,2) * t121 + Ifges(6,3) * t135) * t298 + (Ifges(5,4) * t136 + Ifges(6,5) * t122 - Ifges(5,2) * t135 + Ifges(6,3) * t121) * t299 + ((m(5) * t113 + t46) * t274 + t194) * qJD(3) - m(5) * (t113 * t146 + t21 * t38 + t22 * t39) + t205 * t287 + t207 * t288 + t56 * t278 + t149 * t335 + t203 * mrSges(4,3) - m(4) * (t104 * t116 + t105 * t117) + (t52 - t51) * t109 + m(5) * (-t109 * t3 + t110 * t4 + t183 * t96 - t21 * t84 + t22 * t85) + ((t222 + t313 - t184 / 0.2e1 + ((-m(4) * pkin(2) - mrSges(3,1) - t210) * qJD(2) - t212) * pkin(6) - t326) * t193 + ((t163 + t254) * pkin(6) + (t296 + t266 / 0.2e1 + (Ifges(3,2) / 0.2e1 - Ifges(3,1) / 0.2e1) * t193) * qJD(1) + t221 - t310) * t191 + (t149 * t320 + t199 * t319 + t204) * t243 / 0.2e1) * qJD(1) - t117 * t118 - t116 * t119 + t328 * t81 + t329 * t80 + t330 * t45 + (t1 * t110 + t109 * t2 + t15 * t329 + t16 * t328 + t27 * t330 + t5 * t92) * m(6) + t318 * (-t122 / 0.2e1 + t136 / 0.2e1) + (-t135 * t319 + t136 * t320) * t280 + (-t121 * t319 + t122 * t320) * t281 + (t149 * t322 + t321 * t199) * t301 + (-t321 * t135 + t136 * t322) * t289 + (-t321 * t121 + t122 * t322) * t290 - t146 * t46 + t183 * t18 + t190 * t57 / 0.2e1 + (t35 - t32) * (t121 / 0.2e1 - t135 / 0.2e1) + (-mrSges(6,1) * t248 + mrSges(6,3) * t247) * t27 + (-mrSges(5,1) * t248 - mrSges(5,2) * t247) * t113; qJD(5) * t81 - t28 * t45 + (t268 - t118) * t104 + t312 + (-Ifges(5,2) * t298 + Ifges(6,3) * t299 - t281 * t319 - t290 * t321 + t308) * t200 + (-Ifges(4,2) * t154 + t147 + t90) * t284 - t311 + ((t189 * t4 + t257 * t3) * pkin(3) - t113 * t276 + t21 * t25 - t22 * t26) * m(5) + (Ifges(4,5) * t153 - Ifges(4,6) * t154) * t281 + t89 * t282 + (Ifges(4,1) * t153 - t265) * t283 + t50 * t275 + t51 * t229 + (t267 + t119) * t105 + (t1 * t180 - t15 * t25 + t182 * t2 - t27 * t28 + (-t26 + qJD(5)) * t16) * m(6) + t316 * t25 + t317 * t26 - t166 * (mrSges(4,1) * t154 + mrSges(4,2) * t153) + t180 * t49 + t182 * t52 - t46 * t276 + (-Ifges(5,4) * t298 - Ifges(6,5) * t299 - t320 * t281 - t322 * t290 + t325) * t98; -t317 * t98 + t316 * t200 + t17 + t18 + (-t15 * t200 + t16 * t98 + t5) * m(6) + (t200 * t21 + t22 * t98 + t96) * m(5); t200 * t45 - t177 * t81 + 0.2e1 * (t2 / 0.2e1 + t27 * t289 + t16 * t281) * m(6) + t52;];
tauc = t11(:);
