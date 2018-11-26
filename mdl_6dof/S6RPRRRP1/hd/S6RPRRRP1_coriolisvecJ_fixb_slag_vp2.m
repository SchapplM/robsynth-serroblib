% Calculate vector of centrifugal and coriolis load on the joints for
% S6RPRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2018-11-23 16:24
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RPRRRP1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP1_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP1_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP1_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP1_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:24:04
% EndTime: 2018-11-23 16:24:10
% DurationCPUTime: 6.24s
% Computational Cost: add. (7664->497), mult. (18013->643), div. (0->0), fcn. (11861->8), ass. (0->229)
t357 = Ifges(6,1) + Ifges(7,1);
t351 = Ifges(7,4) + Ifges(6,5);
t356 = Ifges(6,6) - Ifges(7,6);
t210 = sin(qJ(5));
t213 = cos(qJ(5));
t346 = -t356 * t210 + t351 * t213;
t298 = Ifges(7,5) * t210;
t301 = Ifges(6,4) * t210;
t344 = t357 * t213 + t298 - t301;
t271 = qJD(3) + qJD(4);
t355 = Ifges(5,6) * t271 / 0.2e1;
t215 = cos(qJ(3));
t263 = -cos(pkin(10)) * pkin(1) - pkin(2);
t193 = -pkin(3) * t215 + t263;
t187 = qJD(1) * t193;
t255 = Ifges(5,5) * t271;
t354 = t187 * mrSges(5,2) + t255 / 0.2e1;
t211 = sin(qJ(4));
t212 = sin(qJ(3));
t214 = cos(qJ(4));
t191 = t211 * t212 - t214 * t215;
t185 = t191 * qJD(1);
t192 = t211 * t215 + t214 * t212;
t186 = t192 * qJD(1);
t303 = Ifges(5,4) * t186;
t353 = t355 - Ifges(5,2) * t185 / 0.2e1 + t303 / 0.2e1;
t156 = t271 * t192;
t148 = t156 * qJD(1);
t155 = t271 * t191;
t147 = t155 * qJD(1);
t163 = t210 * t186 - t213 * t271;
t96 = -qJD(5) * t163 - t213 * t147;
t164 = t213 * t186 + t210 * t271;
t97 = qJD(5) * t164 - t210 * t147;
t350 = (-Ifges(6,4) + Ifges(7,5)) * t97 + t357 * t96 + t351 * t148;
t159 = Ifges(6,4) * t163;
t182 = qJD(5) + t185;
t299 = Ifges(7,5) * t163;
t339 = t164 * t357 + t351 * t182 - t159 + t299;
t248 = mrSges(7,1) * t210 - mrSges(7,3) * t213;
t250 = mrSges(6,1) * t210 + mrSges(6,2) * t213;
t200 = sin(pkin(10)) * pkin(1) + pkin(7);
t195 = t200 * qJD(1);
t256 = pkin(8) * qJD(1) + t195;
t272 = t212 * qJD(2);
t166 = t215 * t256 + t272;
t160 = t211 * t166;
t207 = t215 * qJD(2);
t230 = t256 * t212;
t165 = t207 - t230;
t162 = qJD(3) * pkin(3) + t165;
t102 = t214 * t162 - t160;
t98 = -pkin(4) * t271 - t102;
t38 = t163 * pkin(5) - t164 * qJ(6) + t98;
t348 = t38 * t248 + t98 * t250;
t347 = t210 * t351 + t213 * t356;
t297 = Ifges(7,5) * t213;
t300 = Ifges(6,4) * t213;
t345 = t210 * t357 - t297 + t300;
t114 = t185 * pkin(4) - t186 * pkin(9) + t187;
t273 = qJD(5) * t213;
t274 = qJD(5) * t210;
t201 = qJD(3) * t207;
t157 = -qJD(3) * t230 + t201;
t218 = qJD(3) * t166;
t32 = qJD(4) * t102 + t214 * t157 - t211 * t218;
t275 = qJD(3) * t212;
t270 = pkin(3) * t275;
t65 = pkin(4) * t148 + pkin(9) * t147 + qJD(1) * t270;
t161 = t214 * t166;
t103 = t211 * t162 + t161;
t99 = pkin(9) * t271 + t103;
t6 = t114 * t273 + t210 * t65 + t213 * t32 - t274 * t99;
t37 = t114 * t210 + t213 * t99;
t7 = -qJD(5) * t37 - t210 * t32 + t213 * t65;
t343 = -t210 * t7 + t213 * t6;
t2 = qJ(6) * t148 + qJD(6) * t182 + t6;
t36 = t114 * t213 - t210 * t99;
t336 = qJD(6) - t36;
t26 = -pkin(5) * t182 + t336;
t4 = -pkin(5) * t148 - t7;
t342 = t2 * t213 + t210 * t4 + t26 * t273;
t237 = Ifges(7,3) * t210 + t297;
t243 = -Ifges(6,2) * t210 + t300;
t319 = t210 / 0.2e1;
t320 = -t210 / 0.2e1;
t325 = t164 / 0.2e1;
t327 = t163 / 0.2e1;
t328 = -t163 / 0.2e1;
t158 = Ifges(7,5) * t164;
t82 = t182 * Ifges(7,6) + t163 * Ifges(7,3) + t158;
t302 = Ifges(6,4) * t164;
t85 = -t163 * Ifges(6,2) + t182 * Ifges(6,6) + t302;
t341 = (-t210 * t37 - t213 * t36) * mrSges(6,3) + t237 * t327 + t243 * t328 + t85 * t320 + t82 * t319 + t348 + t344 * t325 + t346 * t182 / 0.2e1;
t27 = qJ(6) * t182 + t37;
t340 = -t187 * mrSges(5,1) - t36 * mrSges(6,1) + t26 * mrSges(7,1) + t37 * mrSges(6,2) - t27 * mrSges(7,3) + t353;
t181 = pkin(5) * t274 - qJ(6) * t273 - qJD(6) * t210;
t285 = t185 * t213;
t286 = t185 * t210;
t252 = -pkin(5) * t286 + qJ(6) * t285;
t338 = t181 - t103 - t252;
t25 = mrSges(6,1) * t97 + mrSges(6,2) * t96;
t33 = qJD(4) * t103 + t157 * t211 + t214 * t218;
t337 = m(6) * t33 + t25;
t135 = t191 * pkin(4) - t192 * pkin(9) + t193;
t308 = pkin(8) + t200;
t189 = t308 * t212;
t190 = t308 * t215;
t146 = -t189 * t211 + t190 * t214;
t335 = t210 * t135 + t213 * t146;
t334 = -t214 * t189 - t190 * t211;
t257 = qJD(3) * t308;
t178 = t212 * t257;
t179 = t215 * t257;
t72 = qJD(4) * t334 - t178 * t214 - t179 * t211;
t95 = pkin(4) * t156 + pkin(9) * t155 + t270;
t13 = -qJD(5) * t335 - t210 * t72 + t213 * t95;
t333 = t96 / 0.2e1;
t332 = -t97 / 0.2e1;
t331 = t97 / 0.2e1;
t329 = t148 / 0.2e1;
t326 = -t164 / 0.2e1;
t324 = -t182 / 0.2e1;
t322 = t185 / 0.2e1;
t318 = -t213 / 0.2e1;
t317 = t213 / 0.2e1;
t316 = m(5) * t187;
t315 = pkin(3) * t211;
t314 = pkin(3) * t214;
t313 = pkin(5) * t186;
t307 = mrSges(6,3) * t163;
t306 = mrSges(6,3) * t164;
t305 = Ifges(5,1) * t186;
t304 = Ifges(4,4) * t212;
t180 = Ifges(5,4) * t185;
t295 = t334 * t33;
t294 = t185 * mrSges(5,3);
t293 = t186 * mrSges(5,3);
t292 = t191 * t33;
t197 = t263 * qJD(1);
t291 = t197 * mrSges(4,2);
t290 = Ifges(4,5) * qJD(3);
t289 = Ifges(4,6) * qJD(3);
t288 = t155 * t210;
t287 = t155 * t213;
t283 = t210 * t214;
t282 = t213 * t214;
t151 = pkin(4) * t186 + pkin(9) * t185;
t48 = t213 * t102 + t210 * t151;
t110 = t165 * t214 - t160;
t277 = qJD(1) * t212;
t129 = pkin(3) * t277 + t151;
t46 = t213 * t110 + t210 * t129;
t115 = -t163 * mrSges(7,2) + mrSges(7,3) * t182;
t116 = -mrSges(6,2) * t182 - t307;
t280 = t115 + t116;
t117 = mrSges(6,1) * t182 - t306;
t118 = -mrSges(7,1) * t182 + mrSges(7,2) * t164;
t279 = -t117 + t118;
t278 = mrSges(5,1) * t271 - mrSges(6,1) * t163 - mrSges(6,2) * t164 - t293;
t276 = qJD(1) * t215;
t269 = Ifges(7,4) / 0.2e1 + Ifges(6,5) / 0.2e1;
t268 = Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1;
t267 = -Ifges(6,6) / 0.2e1 + Ifges(7,6) / 0.2e1;
t205 = Ifges(4,4) * t276;
t262 = t290 / 0.2e1;
t261 = -t289 / 0.2e1;
t259 = m(4) * t200 + mrSges(4,3);
t43 = -t148 * mrSges(7,1) + t96 * mrSges(7,2);
t109 = t165 * t211 + t161;
t253 = t289 / 0.2e1 + (Ifges(4,2) * t215 + t304) * qJD(1) / 0.2e1 - t197 * mrSges(4,1);
t251 = mrSges(6,1) * t213 - mrSges(6,2) * t210;
t249 = mrSges(7,1) * t213 + mrSges(7,3) * t210;
t242 = Ifges(6,2) * t213 + t301;
t236 = -Ifges(7,3) * t213 + t298;
t235 = pkin(5) * t213 + qJ(6) * t210;
t234 = pkin(5) * t210 - qJ(6) * t213;
t47 = -t102 * t210 + t151 * t213;
t45 = -t110 * t210 + t129 * t213;
t66 = t135 * t213 - t146 * t210;
t175 = t195 * t215 + t272;
t194 = -pkin(4) - t235;
t12 = t135 * t273 - t146 * t274 + t210 * t95 + t213 * t72;
t220 = t7 * mrSges(6,1) - t4 * mrSges(7,1) - t6 * mrSges(6,2) + t2 * mrSges(7,3);
t167 = -mrSges(5,2) * t271 - t294;
t219 = t210 * t279 + t213 * t280 + t167;
t73 = qJD(4) * t146 - t178 * t211 + t214 * t179;
t41 = -mrSges(7,2) * t97 + mrSges(7,3) * t148;
t42 = mrSges(6,1) * t148 - mrSges(6,3) * t96;
t44 = -mrSges(6,2) * t148 - mrSges(6,3) * t97;
t217 = (t41 + t44) * t213 + (-t42 + t43) * t210 + (-t210 * t280 + t213 * t279) * qJD(5) + m(7) * (-t27 * t274 + t342) + m(6) * (-t273 * t36 - t274 * t37 + t343);
t131 = -t180 + t255 + t305;
t19 = t96 * Ifges(7,5) + t148 * Ifges(7,6) + t97 * Ifges(7,3);
t20 = t96 * Ifges(6,4) - t97 * Ifges(6,2) + t148 * Ifges(6,6);
t83 = t164 * Ifges(6,5) - t163 * Ifges(6,6) + t182 * Ifges(6,3);
t84 = t164 * Ifges(7,4) + t182 * Ifges(7,2) + t163 * Ifges(7,6);
t9 = pkin(5) * t97 - qJ(6) * t96 - qJD(6) * t164 + t33;
t216 = -t32 * mrSges(5,2) - Ifges(5,5) * t147 - Ifges(5,6) * t148 - t102 * t294 + t19 * t318 + t20 * t317 + t236 * t331 + t242 * t332 - t9 * t249 + t345 * t333 + (-t251 - mrSges(5,1)) * t33 + t347 * t329 + (t131 - t180) * t322 + t350 * t319 + (-t85 / 0.2e1 + t82 / 0.2e1) * t286 + (-t237 * t328 - t243 * t327 - t346 * t324 - t344 * t326 + t348 + t354) * t185 + (Ifges(6,6) * t327 + Ifges(7,6) * t328 + t355 - Ifges(5,2) * t322 + t351 * t326 + (Ifges(6,3) + Ifges(7,2)) * t324 + t340) * t186 + (-t285 * t36 - t286 * t37 + t343) * mrSges(6,3) - (-Ifges(5,1) * t185 - t303 + t83 + t84) * t186 / 0.2e1 + (t26 * t285 + (-t274 - t286) * t27 + t342) * mrSges(7,2) + t341 * qJD(5) + (t285 / 0.2e1 + t273 / 0.2e1) * t339;
t198 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t276;
t196 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t277;
t188 = t194 - t314;
t184 = Ifges(4,1) * t277 + t205 + t290;
t177 = t186 * qJ(6);
t174 = -t195 * t212 + t207;
t171 = t175 * qJD(3);
t170 = -t195 * t275 + t201;
t169 = qJD(4) * t315 + t181;
t150 = mrSges(5,1) * t185 + mrSges(5,2) * t186;
t143 = Ifges(7,2) * t148;
t141 = Ifges(6,3) * t148;
t107 = mrSges(7,1) * t163 - mrSges(7,3) * t164;
t106 = pkin(5) * t164 + qJ(6) * t163;
t94 = Ifges(7,4) * t96;
t93 = Ifges(6,5) * t96;
t92 = Ifges(6,6) * t97;
t91 = Ifges(7,6) * t97;
t89 = t192 * t234 - t334;
t56 = t109 + t252;
t51 = -pkin(5) * t191 - t66;
t50 = qJ(6) * t191 + t335;
t40 = -t47 - t313;
t39 = t177 + t48;
t35 = -t45 - t313;
t34 = t177 + t46;
t24 = mrSges(7,1) * t97 - mrSges(7,3) * t96;
t14 = -t234 * t155 + (qJD(5) * t235 - qJD(6) * t213) * t192 + t73;
t11 = -pkin(5) * t156 - t13;
t10 = qJ(6) * t156 + qJD(6) * t191 + t12;
t1 = [t193 * (mrSges(5,1) * t148 - mrSges(5,2) * t147) + t72 * t167 - t334 * t25 + t10 * t115 + t12 * t116 + t13 * t117 + t11 * t118 + t14 * t107 + t89 * t24 + t66 * t42 + t335 * t44 + t50 * t41 + t51 * t43 - t278 * t73 + (-t146 * t148 + t147 * t334) * mrSges(5,3) + m(5) * (-t102 * t73 + t103 * t72 + t146 * t32 - t295) + m(7) * (t10 * t27 + t11 * t26 + t14 * t38 + t2 * t50 + t4 * t51 + t89 * t9) + m(6) * (t12 * t37 + t13 * t36 + t335 * t6 + t66 * t7 + t73 * t98 - t295) + (t83 / 0.2e1 + t84 / 0.2e1 - t103 * mrSges(5,3) + t268 * t182 + t269 * t164 + t267 * t163 - t340 - t353) * t156 - (t339 * t317 - t180 / 0.2e1 - t102 * mrSges(5,3) + (-t210 * t27 + t213 * t26) * mrSges(7,2) + t131 / 0.2e1 + t305 / 0.2e1 + t341 + t354) * t155 + (t259 * t170 + (-t200 * t196 + t184 / 0.2e1 + 0.3e1 / 0.2e1 * t205 + t262 - t259 * t174 + 0.2e1 * t291) * qJD(3)) * t215 + (t93 / 0.2e1 - t92 / 0.2e1 + t141 / 0.2e1 + t94 / 0.2e1 + t143 / 0.2e1 + t91 / 0.2e1 + Ifges(5,4) * t147 - t32 * mrSges(5,3) + t267 * t97 + t269 * t96 + (Ifges(5,2) + t268) * t148 + t220) * t191 + (-Ifges(5,1) * t147 - Ifges(5,4) * t148 + t237 * t331 + t243 * t332 + t19 * t319 + t20 * t320 + t9 * t248 + (mrSges(5,3) + t250) * t33 + (-t210 * t6 - t213 * t7) * mrSges(6,3) + (-t2 * t210 + t213 * t4) * mrSges(7,2) + (t85 * t318 + t242 * t327 + t236 * t328 + t98 * t251 + t38 * t249 + (t210 * t36 - t213 * t37) * mrSges(6,3) + (-t210 * t26 - t213 * t27) * mrSges(7,2) + t345 * t326 + t347 * t324 + t339 * t320) * qJD(5) + t344 * t333 + t346 * t329 + (qJD(5) * t82 + t350) * t317) * t192 + (t259 * t171 + (-t200 * t198 + t261 - t259 * t175 + (t263 * mrSges(4,1) - 0.3e1 / 0.2e1 * t304 + (-0.3e1 / 0.2e1 * Ifges(4,2) + 0.3e1 / 0.2e1 * Ifges(4,1)) * t215) * qJD(1) + (t150 + qJD(1) * (mrSges(5,1) * t191 + mrSges(5,2) * t192) + 0.2e1 * t316) * pkin(3) - t253) * qJD(3)) * t212; (-t147 * mrSges(5,3) + t24 + t25) * t191 + (t107 - t278) * t156 + (-t212 * t196 + t215 * t198 + (-t212 ^ 2 - t215 ^ 2) * qJD(1) * mrSges(4,3)) * qJD(3) - t219 * t155 + m(5) * (-t102 * t156 - t103 * t155 + t292) + m(6) * (t156 * t98 - t287 * t37 + t288 * t36 + t292) + m(7) * (t156 * t38 + t191 * t9 - t26 * t288 - t27 * t287) + m(4) * (t170 * t212 - t171 * t215 + (-t174 * t212 + t175 * t215) * qJD(3)) + (m(5) * t32 - t148 * mrSges(5,3) + t217) * t192; t278 * t109 + (t169 - t56) * t107 + t216 - m(7) * (t26 * t35 + t27 * t34 + t38 * t56) - m(6) * (t109 * t98 + t36 * t45 + t37 * t46) + t217 * (pkin(9) + t315) + (m(5) * (t211 * t32 - t214 * t33) + (t214 * t147 - t211 * t148) * mrSges(5,3) + (-t278 * t211 + t219 * t214 + m(5) * (-t102 * t211 + t103 * t214) + m(7) * (t26 * t283 + t27 * t282) + m(6) * (t211 * t98 + t282 * t37 - t283 * t36)) * qJD(4)) * pkin(3) + ((t174 * mrSges(4,3) + t262 - t291 - t184 / 0.2e1 - t205 / 0.2e1) * t215 + (t175 * mrSges(4,3) + t261 + (t304 / 0.2e1 + (-Ifges(4,1) / 0.2e1 + Ifges(4,2) / 0.2e1) * t215) * qJD(1) + (-t150 - t316) * pkin(3) + t253) * t212) * qJD(1) - m(5) * (-t102 * t109 + t103 * t110) + m(7) * (t169 * t38 + t188 * t9) + t175 * t196 - t174 * t198 + t188 * t24 - t110 * t167 - t170 * mrSges(4,2) - t171 * mrSges(4,1) - t34 * t115 - t46 * t116 - t45 * t117 - t35 * t118 + t103 * t293 + t337 * (-pkin(4) - t314); (t278 + t293) * t103 + t216 + t217 * pkin(9) + t338 * t107 - m(6) * (t103 * t98 + t36 * t47 + t37 * t48) + t194 * t24 - t102 * t167 - t39 * t115 - t48 * t116 - t47 * t117 - t40 * t118 - t337 * pkin(4) + (t194 * t9 - t26 * t40 - t27 * t39 + t338 * t38) * m(7); t143 + t141 + t220 + t94 + t93 - t92 + t91 - t38 * (mrSges(7,1) * t164 + mrSges(7,3) * t163) - t98 * (mrSges(6,1) * t164 - mrSges(6,2) * t163) + qJD(6) * t115 - t106 * t107 - pkin(5) * t43 + qJ(6) * t41 + t85 * t325 + (Ifges(7,3) * t164 - t299) * t328 + (-t279 + t306) * t37 + (-t280 - t307) * t36 + (t163 * t26 + t164 * t27) * mrSges(7,2) + (-t351 * t163 - t164 * t356) * t324 + (-pkin(5) * t4 + qJ(6) * t2 - t106 * t38 - t26 * t37 + t27 * t336) * m(7) + (-Ifges(6,2) * t164 - t159 + t339) * t327 + (-t163 * t357 + t158 - t302 + t82) * t326; t164 * t107 - t182 * t115 + 0.2e1 * (t4 / 0.2e1 + t38 * t325 + t27 * t324) * m(7) + t43;];
tauc  = t1(:);
