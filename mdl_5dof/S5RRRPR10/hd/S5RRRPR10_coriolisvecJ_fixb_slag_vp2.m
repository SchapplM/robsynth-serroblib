% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-31 21:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRPR10_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR10_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR10_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR10_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR10_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR10_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR10_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:26:39
% EndTime: 2019-12-31 21:27:16
% DurationCPUTime: 16.41s
% Computational Cost: add. (9697->652), mult. (26299->938), div. (0->0), fcn. (20240->10), ass. (0->294)
t229 = sin(pkin(10));
t231 = cos(pkin(10));
t234 = sin(qJ(3));
t237 = cos(qJ(3));
t209 = t229 * t237 + t231 * t234;
t238 = cos(qJ(2));
t230 = sin(pkin(5));
t285 = qJD(1) * t230;
t273 = t238 * t285;
t157 = t209 * t273;
t202 = t209 * qJD(3);
t287 = t157 - t202;
t208 = t229 * t234 - t231 * t237;
t158 = t208 * t273;
t203 = t208 * qJD(3);
t286 = -t158 + t203;
t318 = -qJ(4) - pkin(8);
t265 = qJD(3) * t318;
t201 = qJD(4) * t237 + t234 * t265;
t243 = -qJD(4) * t234 + t237 * t265;
t137 = t231 * t201 + t229 * t243;
t235 = sin(qJ(2));
t274 = t235 * t285;
t232 = cos(pkin(5));
t284 = qJD(1) * t232;
t277 = pkin(1) * t284;
t195 = -pkin(7) * t274 + t238 * t277;
t245 = (pkin(2) * t235 - pkin(8) * t238) * t230;
t196 = qJD(1) * t245;
t134 = -t234 * t195 + t237 * t196;
t111 = (-qJ(4) * t237 * t238 + pkin(3) * t235) * t285 + t134;
t135 = t237 * t195 + t234 * t196;
t263 = t234 * t273;
t119 = -qJ(4) * t263 + t135;
t67 = t229 * t111 + t231 * t119;
t378 = -t67 + t137;
t221 = qJD(2) + t284;
t282 = qJD(2) * t238;
t271 = t237 * t282;
t280 = qJD(3) * t237;
t281 = qJD(3) * t234;
t147 = t221 * t280 + (-t235 * t281 + t271) * t285;
t178 = t221 * t234 + t237 * t274;
t283 = qJD(2) * t230;
t267 = qJD(1) * t283;
t262 = t235 * t267;
t198 = pkin(7) * t273 + t235 * t277;
t165 = pkin(8) * t221 + t198;
t191 = (-pkin(2) * t238 - pkin(8) * t235 - pkin(1)) * t230;
t172 = qJD(1) * t191;
t118 = t165 * t237 + t172 * t234;
t197 = qJD(2) * t245;
t185 = qJD(1) * t197;
t290 = t230 * t235;
t222 = pkin(7) * t290;
t326 = pkin(1) * t238;
t206 = t232 * t326 - t222;
t199 = t206 * qJD(2);
t186 = qJD(1) * t199;
t75 = -qJD(3) * t118 + t237 * t185 - t186 * t234;
t37 = pkin(3) * t262 - qJ(4) * t147 - qJD(4) * t178 + t75;
t272 = t234 * t282;
t148 = -t221 * t281 + (-t235 * t280 - t272) * t285;
t177 = t221 * t237 - t234 * t274;
t74 = -t165 * t281 + t172 * t280 + t234 * t185 + t237 * t186;
t43 = qJ(4) * t148 + qJD(4) * t177 + t74;
t11 = t229 * t37 + t231 * t43;
t289 = t230 * t238;
t207 = t232 * t235 * pkin(1) + pkin(7) * t289;
t200 = t207 * qJD(2);
t187 = qJD(1) * t200;
t115 = -t148 * pkin(3) + t187;
t233 = sin(qJ(5));
t236 = cos(qJ(5));
t215 = qJD(3) - t273;
t117 = -t165 * t234 + t237 * t172;
t92 = -qJ(4) * t178 + t117;
t80 = pkin(3) * t215 + t92;
t93 = qJ(4) * t177 + t118;
t86 = t231 * t93;
t36 = t229 * t80 + t86;
t32 = pkin(9) * t215 + t36;
t164 = -t221 * pkin(2) - t195;
t123 = -t177 * pkin(3) + qJD(4) + t164;
t246 = t177 * t229 + t231 * t178;
t264 = t231 * t177 - t178 * t229;
t51 = -pkin(4) * t264 - pkin(9) * t246 + t123;
t12 = -t233 * t32 + t236 * t51;
t94 = t147 * t229 - t231 * t148;
t95 = t147 * t231 + t148 * t229;
t27 = t94 * pkin(4) - t95 * pkin(9) + t115;
t9 = pkin(9) * t262 + t11;
t1 = qJD(5) * t12 + t233 * t27 + t236 * t9;
t13 = t233 * t51 + t236 * t32;
t2 = -qJD(5) * t13 - t233 * t9 + t236 * t27;
t260 = t2 * mrSges(6,1) - t1 * mrSges(6,2);
t320 = t95 * Ifges(5,4);
t103 = t215 * t233 + t236 * t246;
t47 = -qJD(5) * t103 - t233 * t95 + t236 * t262;
t44 = Ifges(6,6) * t47;
t102 = t215 * t236 - t233 * t246;
t46 = qJD(5) * t102 + t233 * t262 + t236 * t95;
t45 = Ifges(6,5) * t46;
t5 = Ifges(6,3) * t94 + t44 + t45;
t377 = t260 + t115 * mrSges(5,1) - t11 * mrSges(5,3) + t5 / 0.2e1 - t320 / 0.2e1;
t375 = t94 * Ifges(5,2);
t374 = -pkin(9) * t274 + t378;
t155 = pkin(3) * t263 + t198;
t373 = pkin(3) * t281 - t287 * pkin(4) + t286 * pkin(9) - t155;
t122 = qJD(5) - t264;
t308 = t122 * Ifges(6,3);
t309 = t103 * Ifges(6,5);
t310 = t102 * Ifges(6,6);
t28 = t308 + t309 + t310;
t296 = t215 * Ifges(5,6);
t303 = t246 * Ifges(5,4);
t306 = t264 * Ifges(5,2);
t72 = t296 + t303 + t306;
t372 = -t72 / 0.2e1 + t28 / 0.2e1;
t370 = -t123 * mrSges(5,1) - t12 * mrSges(6,1) + t13 * mrSges(6,2) + t36 * mrSges(5,3);
t330 = t215 / 0.2e1;
t369 = t246 / 0.2e1;
t368 = t264 / 0.2e1;
t367 = Ifges(4,3) + Ifges(5,3);
t136 = t201 * t229 - t231 * t243;
t66 = t111 * t231 - t119 * t229;
t366 = -t66 - t136;
t228 = -pkin(3) * t237 - pkin(2);
t149 = pkin(4) * t208 - pkin(9) * t209 + t228;
t216 = t318 * t237;
t268 = t318 * t234;
t163 = -t231 * t216 + t229 * t268;
t101 = t149 * t233 + t163 * t236;
t365 = -qJD(5) * t101 - t233 * t374 + t236 * t373;
t100 = t149 * t236 - t163 * t233;
t364 = qJD(5) * t100 + t233 * t373 + t236 * t374;
t363 = pkin(4) * t274 - t366;
t350 = Ifges(5,1) * t369 + Ifges(5,4) * t368 + Ifges(5,5) * t330;
t250 = t12 * t236 + t13 * t233;
t362 = t250 * mrSges(6,3);
t190 = pkin(8) * t232 + t207;
t133 = t237 * t190 + t234 * t191;
t361 = Ifges(4,5) * t147 + Ifges(5,5) * t95 + Ifges(4,6) * t148 - Ifges(5,6) * t94 + t262 * t367;
t360 = -t234 * t75 + t237 * t74;
t359 = t1 * t236 - t2 * t233;
t300 = t178 * Ifges(4,4);
t113 = t177 * Ifges(4,2) + t215 * Ifges(4,6) + t300;
t173 = Ifges(4,4) * t177;
t114 = t178 * Ifges(4,1) + t215 * Ifges(4,5) + t173;
t247 = t117 * t237 + t118 * t234;
t315 = Ifges(4,4) * t237;
t316 = Ifges(4,4) * t234;
t327 = t237 / 0.2e1;
t334 = t178 / 0.2e1;
t335 = t177 / 0.2e1;
t358 = -t247 * mrSges(4,3) + (Ifges(4,5) * t237 - Ifges(4,6) * t234) * t330 + t164 * (mrSges(4,1) * t234 + mrSges(4,2) * t237) + (-Ifges(4,2) * t234 + t315) * t335 + (Ifges(4,1) * t237 - t316) * t334 - t234 * t113 / 0.2e1 + t114 * t327;
t7 = t46 * Ifges(6,1) + t47 * Ifges(6,4) + t94 * Ifges(6,5);
t356 = t7 / 0.2e1;
t314 = Ifges(6,4) * t103;
t29 = Ifges(6,2) * t102 + Ifges(6,6) * t122 + t314;
t354 = -t29 / 0.2e1;
t353 = t46 / 0.2e1;
t352 = t47 / 0.2e1;
t349 = t94 / 0.2e1;
t348 = pkin(1) * mrSges(3,1);
t347 = pkin(1) * mrSges(3,2);
t346 = -t102 / 0.2e1;
t345 = t102 / 0.2e1;
t344 = -t103 / 0.2e1;
t343 = t103 / 0.2e1;
t342 = -t122 / 0.2e1;
t341 = t122 / 0.2e1;
t204 = t232 * t237 - t234 * t290;
t205 = t232 * t234 + t237 * t290;
t141 = -t231 * t204 + t205 * t229;
t340 = -t141 / 0.2e1;
t142 = t204 * t229 + t205 * t231;
t339 = t142 / 0.2e1;
t338 = t147 / 0.2e1;
t337 = t148 / 0.2e1;
t333 = t204 / 0.2e1;
t332 = t205 / 0.2e1;
t331 = -t215 / 0.2e1;
t329 = -t233 / 0.2e1;
t328 = t236 / 0.2e1;
t325 = pkin(3) * t178;
t322 = t94 * Ifges(5,4);
t321 = t95 * Ifges(5,1);
t154 = qJD(3) * t204 + t230 * t271;
t270 = t235 * t283;
t82 = -qJD(3) * t133 + t237 * t197 - t199 * t234;
t55 = pkin(3) * t270 - qJ(4) * t154 - qJD(4) * t205 + t82;
t153 = -qJD(3) * t205 - t230 * t272;
t81 = -t190 * t281 + t191 * t280 + t234 * t197 + t237 * t199;
t59 = qJ(4) * t153 + qJD(4) * t204 + t81;
t20 = t229 * t55 + t231 * t59;
t317 = Ifges(3,4) * t235;
t313 = Ifges(6,4) * t233;
t312 = Ifges(6,4) * t236;
t311 = Ifges(3,5) * t238;
t305 = t264 * Ifges(5,6);
t302 = t246 * Ifges(5,5);
t301 = t177 * Ifges(4,6);
t299 = t178 * Ifges(4,5);
t298 = t186 * mrSges(3,2);
t295 = t221 * Ifges(3,5);
t294 = t229 * t93;
t110 = qJ(4) * t204 + t133;
t132 = -t234 * t190 + t237 * t191;
t98 = -pkin(3) * t289 - t205 * qJ(4) + t132;
t58 = t231 * t110 + t229 * t98;
t109 = mrSges(5,1) * t215 - mrSges(5,3) * t246;
t60 = -mrSges(6,1) * t102 + mrSges(6,2) * t103;
t291 = t109 - t60;
t288 = -mrSges(3,1) * t221 - mrSges(4,1) * t177 + mrSges(4,2) * t178 + mrSges(3,3) * t274;
t276 = -Ifges(5,3) / 0.2e1 - Ifges(4,3) / 0.2e1;
t48 = t94 * mrSges(5,1) + t95 * mrSges(5,2);
t259 = -t1 * t233 - t2 * t236;
t258 = mrSges(6,1) * t236 - mrSges(6,2) * t233;
t257 = mrSges(6,1) * t233 + mrSges(6,2) * t236;
t256 = Ifges(6,1) * t236 - t313;
t255 = Ifges(6,1) * t233 + t312;
t254 = -Ifges(6,2) * t233 + t312;
t253 = Ifges(6,2) * t236 + t313;
t252 = Ifges(6,5) * t236 - Ifges(6,6) * t233;
t251 = Ifges(6,5) * t233 + Ifges(6,6) * t236;
t249 = t12 * t233 - t13 * t236;
t10 = -t229 * t43 + t231 * t37;
t19 = -t229 * t59 + t231 * t55;
t35 = t231 * t80 - t294;
t54 = -pkin(9) * t289 + t58;
t189 = t222 + (-pkin(2) - t326) * t232;
t146 = -t204 * pkin(3) + t189;
t70 = t141 * pkin(4) - t142 * pkin(9) + t146;
t22 = t233 * t70 + t236 * t54;
t21 = -t233 * t54 + t236 * t70;
t68 = -mrSges(6,2) * t122 + mrSges(6,3) * t102;
t69 = mrSges(6,1) * t122 - mrSges(6,3) * t103;
t248 = -t233 * t69 + t236 * t68;
t57 = -t229 * t110 + t231 * t98;
t120 = -t233 * t142 - t236 * t289;
t244 = -t236 * t142 + t233 * t289;
t130 = -t153 * pkin(3) + t200;
t242 = -t296 / 0.2e1 - t303 / 0.2e1 - t306 / 0.2e1 + t310 / 0.2e1 + t309 / 0.2e1 + t308 / 0.2e1 + t372;
t99 = Ifges(6,4) * t102;
t30 = Ifges(6,1) * t103 + Ifges(6,5) * t122 + t99;
t31 = -pkin(4) * t215 - t35;
t241 = t252 * t341 + t254 * t345 + t256 * t343 + t257 * t31 + t29 * t329 + t30 * t328;
t239 = 0.2e1 * t350 + t241;
t227 = -pkin(3) * t231 - pkin(4);
t217 = Ifges(3,4) * t273;
t214 = t267 * t311;
t194 = -t221 * mrSges(3,2) + mrSges(3,3) * t273;
t162 = -t216 * t229 - t231 * t268;
t160 = Ifges(3,1) * t274 + t217 + t295;
t159 = Ifges(3,6) * t221 + (Ifges(3,2) * t238 + t317) * t285;
t151 = mrSges(4,1) * t215 - mrSges(4,3) * t178;
t150 = -mrSges(4,2) * t215 + mrSges(4,3) * t177;
t140 = -t158 * t236 + t233 * t274;
t139 = t158 * t233 + t236 * t274;
t129 = -mrSges(4,2) * t262 + mrSges(4,3) * t148;
t128 = mrSges(4,1) * t262 - mrSges(4,3) * t147;
t112 = t215 * Ifges(4,3) + t299 + t301;
t108 = -mrSges(5,2) * t215 + mrSges(5,3) * t264;
t107 = t153 * t229 + t154 * t231;
t106 = -t231 * t153 + t154 * t229;
t96 = -mrSges(4,1) * t148 + mrSges(4,2) * t147;
t85 = t147 * Ifges(4,1) + t148 * Ifges(4,4) + Ifges(4,5) * t262;
t84 = t147 * Ifges(4,4) + t148 * Ifges(4,2) + Ifges(4,6) * t262;
t79 = mrSges(5,1) * t262 - mrSges(5,3) * t95;
t78 = -mrSges(5,2) * t262 - mrSges(5,3) * t94;
t76 = -mrSges(5,1) * t264 + mrSges(5,2) * t246;
t71 = t215 * Ifges(5,3) + t302 + t305;
t65 = pkin(4) * t246 - pkin(9) * t264 + t325;
t62 = qJD(5) * t244 - t233 * t107 + t236 * t270;
t61 = qJD(5) * t120 + t236 * t107 + t233 * t270;
t53 = pkin(4) * t289 - t57;
t42 = t231 * t92 - t294;
t41 = t229 * t92 + t86;
t39 = Ifges(5,5) * t262 + t321 - t322;
t38 = Ifges(5,6) * t262 + t320 - t375;
t33 = t106 * pkin(4) - t107 * pkin(9) + t130;
t24 = -mrSges(6,2) * t94 + mrSges(6,3) * t47;
t23 = mrSges(6,1) * t94 - mrSges(6,3) * t46;
t18 = pkin(9) * t270 + t20;
t17 = -pkin(4) * t270 - t19;
t16 = t233 * t65 + t236 * t42;
t15 = -t233 * t42 + t236 * t65;
t14 = -mrSges(6,1) * t47 + mrSges(6,2) * t46;
t8 = -pkin(4) * t262 - t10;
t6 = t46 * Ifges(6,4) + t47 * Ifges(6,2) + t94 * Ifges(6,6);
t4 = -qJD(5) * t22 - t18 * t233 + t236 * t33;
t3 = qJD(5) * t21 + t18 * t236 + t233 * t33;
t25 = [(Ifges(4,5) * t154 + Ifges(5,5) * t107 + Ifges(4,6) * t153 + t270 * t367) * t330 + ((Ifges(3,5) * t232 / 0.2e1 - t206 * mrSges(3,3) + (-0.2e1 * t347 + 0.3e1 / 0.2e1 * Ifges(3,4) * t238) * t230) * t238 + (-Ifges(3,6) * t232 + Ifges(4,5) * t332 + Ifges(4,6) * t333 + Ifges(5,5) * t339 + Ifges(5,6) * t340 - t207 * mrSges(3,3) + (-0.2e1 * t348 - 0.3e1 / 0.2e1 * t317) * t230 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2) + t276) * t289) * t235) * t267 - t94 * (Ifges(5,4) * t142 - Ifges(5,6) * t289) / 0.2e1 + (Ifges(6,5) * t61 + Ifges(6,6) * t62) * t341 + (t1 * t120 - t12 * t61 + t13 * t62 + t2 * t244) * mrSges(6,3) + (t186 * t238 + t187 * t235 + (-t195 * t238 - t198 * t235) * qJD(2)) * t230 * mrSges(3,3) + (-Ifges(6,4) * t244 + Ifges(6,2) * t120) * t352 + m(6) * (t1 * t22 + t12 * t4 + t13 * t3 + t17 * t31 + t2 * t21 + t53 * t8) + m(5) * (t10 * t57 + t11 * t58 + t115 * t146 + t123 * t130 + t19 * t35 + t20 * t36) + m(4) * (t117 * t82 + t118 * t81 + t132 * t75 + t133 * t74 + t164 * t200 + t187 * t189) + m(3) * (t186 * t207 - t187 * t206 - t195 * t200 + t198 * t199) + (t107 * t123 + t11 * t289 + t115 * t142 - t270 * t36) * mrSges(5,2) + (Ifges(6,4) * t61 + Ifges(6,2) * t62) * t345 + t95 * (Ifges(5,1) * t142 - Ifges(5,5) * t289) / 0.2e1 - t361 * t289 / 0.2e1 + (-Ifges(6,5) * t244 + Ifges(6,6) * t120) * t349 + (t214 / 0.2e1 - t298 - t187 * mrSges(3,1)) * t232 + t10 * (-mrSges(5,1) * t289 - t142 * mrSges(5,3)) + t75 * (-mrSges(4,1) * t289 - t205 * mrSges(4,3)) + t74 * (mrSges(4,2) * t289 + t204 * mrSges(4,3)) + (Ifges(5,4) * t107 + Ifges(5,6) * t270) * t368 + t288 * t200 + (t238 * t160 + t221 * (-Ifges(3,6) * t235 + t311) + (t112 + t71) * t235) * t283 / 0.2e1 + t8 * (-mrSges(6,1) * t120 - mrSges(6,2) * t244) - t244 * t356 - t159 * t270 / 0.2e1 + t118 * (-mrSges(4,2) * t270 + mrSges(4,3) * t153) + t35 * (mrSges(5,1) * t270 - mrSges(5,3) * t107) + t117 * (mrSges(4,1) * t270 - mrSges(4,3) * t154) + (Ifges(6,1) * t61 + Ifges(6,4) * t62) * t343 + t187 * (-mrSges(4,1) * t204 + mrSges(4,2) * t205) + (Ifges(5,1) * t107 + Ifges(5,5) * t270) * t369 + (-Ifges(5,4) * t369 + Ifges(6,5) * t343 - Ifges(5,2) * t368 - Ifges(5,6) * t330 + Ifges(6,6) * t345 + Ifges(6,3) * t341 - t370 + t372) * t106 + (t375 / 0.2e1 + Ifges(6,3) * t349 + Ifges(6,6) * t352 + Ifges(6,5) * t353 + t377) * t141 + t21 * t23 + t22 * t24 + t53 * t14 + t17 * t60 + t61 * t30 / 0.2e1 + t62 * t29 / 0.2e1 + t31 * (-mrSges(6,1) * t62 + mrSges(6,2) * t61) + t3 * t68 + t4 * t69 + t58 * t78 + t57 * t79 + t20 * t108 + t19 * t109 + t120 * t6 / 0.2e1 + t130 * t76 + t132 * t128 + t133 * t129 + t146 * t48 + t81 * t150 + t82 * t151 + t153 * t113 / 0.2e1 + t154 * t114 / 0.2e1 + t164 * (-mrSges(4,1) * t153 + mrSges(4,2) * t154) + t189 * t96 + (-Ifges(6,1) * t244 + Ifges(6,4) * t120) * t353 + t199 * t194 + t85 * t332 + t84 * t333 + (Ifges(4,1) * t154 + Ifges(4,4) * t153 + Ifges(4,5) * t270) * t334 + (Ifges(4,4) * t154 + Ifges(4,2) * t153 + Ifges(4,6) * t270) * t335 + (Ifges(4,4) * t205 + Ifges(4,2) * t204 - Ifges(4,6) * t289) * t337 + (Ifges(4,1) * t205 + Ifges(4,4) * t204 - Ifges(4,5) * t289) * t338 + t39 * t339 + t38 * t340 + t107 * t350; t242 * t202 + (-mrSges(4,1) * t237 + mrSges(4,2) * t234 - mrSges(3,1)) * t187 + ((m(5) * t123 + t76) * t234 * pkin(3) + (-m(4) * t247 - t234 * t150 - t237 * t151) * pkin(8) + t358) * qJD(3) + ((t195 * mrSges(3,3) + t285 * t347 - t160 / 0.2e1 - t295 / 0.2e1 - t217 / 0.2e1 - t358) * t238 + (t36 * mrSges(5,2) - t35 * mrSges(5,1) - t305 / 0.2e1 - t302 / 0.2e1 + t198 * mrSges(3,3) - t301 / 0.2e1 - t299 / 0.2e1 - t117 * mrSges(4,1) + t118 * mrSges(4,2) + t159 / 0.2e1 - t71 / 0.2e1 - t112 / 0.2e1 + t276 * t215 + (t348 + t317 / 0.2e1 + (Ifges(3,2) / 0.2e1 - Ifges(3,1) / 0.2e1) * t238) * t285 + (-qJD(2) + t221 / 0.2e1) * Ifges(3,6) + (Ifges(4,5) * t234 + Ifges(5,5) * t209 + Ifges(4,6) * t237 - Ifges(5,6) * t208) * qJD(2) / 0.2e1) * t235) * t285 + (t256 * t353 + t254 * t352 + t252 * t349 + t8 * t257 - t322 / 0.2e1 + t321 / 0.2e1 + t115 * mrSges(5,2) + t39 / 0.2e1 - t10 * mrSges(5,3) + t6 * t329 + t7 * t328 + t259 * mrSges(6,3) + (mrSges(6,3) * t249 + t236 * t354 + t251 * t342 + t253 * t346 + t255 * t344 + t258 * t31 + t30 * t329) * qJD(5)) * t209 - m(5) * (t123 * t155 + t35 * t66 + t36 * t67) - m(4) * (t117 * t134 + t118 * t135 + t164 * t198) - t298 + t214 - t264 * (-Ifges(5,4) * t158 - Ifges(5,2) * t157) / 0.2e1 - t246 * (-Ifges(5,1) * t158 - Ifges(5,4) * t157) / 0.2e1 + (-Ifges(5,5) * t158 - Ifges(5,6) * t157) * t331 + ((t233 * t203 - t139) * mrSges(6,3) + t287 * mrSges(6,2)) * t13 + ((t236 * t203 + t140) * mrSges(6,3) - t287 * mrSges(6,1)) * t12 - t239 * t203 + m(4) * (-pkin(2) * t187 + pkin(8) * t360) + t360 * mrSges(4,3) - t288 * t198 + (t286 * t35 + t287 * t36) * mrSges(5,3) + (-mrSges(5,1) * t287 - mrSges(5,2) * t286) * t123 + t364 * t68 + t365 * t69 + (t1 * t101 + t100 * t2 + t12 * t365 + t364 * t13 + t162 * t8 + t363 * t31) * m(6) + t363 * t60 + m(5) * (-t10 * t162 + t11 * t163 + t115 * t228 - t136 * t35 + t137 * t36) + (t14 - t79) * t162 + (-t128 * t234 + t129 * t237) * pkin(8) + t234 * t85 / 0.2e1 + t228 * t48 + t366 * t109 + t378 * t108 + (t45 / 0.2e1 + t44 / 0.2e1 - t38 / 0.2e1 + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t94 + t377) * t208 - pkin(2) * t96 + t100 * t23 + t101 * t24 - t140 * t30 / 0.2e1 - t31 * (-mrSges(6,1) * t139 + mrSges(6,2) * t140) - t135 * t150 - t134 * t151 - t155 * t76 - t157 * t28 / 0.2e1 + t157 * t72 / 0.2e1 + t163 * t78 - t195 * t194 + t84 * t327 + (Ifges(4,2) * t237 + t316) * t337 + (Ifges(4,1) * t234 + t315) * t338 + (Ifges(6,5) * t140 + Ifges(6,6) * t139 + Ifges(6,3) * t157) * t342 + (Ifges(6,1) * t140 + Ifges(6,4) * t139 + Ifges(6,5) * t157) * t344 + (Ifges(6,4) * t140 + Ifges(6,2) * t139 + Ifges(6,6) * t157) * t346 + t139 * t354 + t158 * t350; (-t178 * t76 + t229 * t78 + t231 * t79) * pkin(3) - t8 * t258 - (-Ifges(4,2) * t178 + t114 + t173) * t177 / 0.2e1 + (t117 * t177 + t118 * t178) * mrSges(4,3) + t359 * mrSges(6,3) + (m(6) * t359 + (-m(6) * t250 - t233 * t68 - t236 * t69) * qJD(5) - t23 * t233 + t236 * t24) * (pkin(3) * t229 + pkin(9)) + (-t123 * t325 + t35 * t41 - t36 * t42 + (t10 * t231 + t11 * t229) * pkin(3)) * m(5) + t361 - t178 * (Ifges(4,1) * t177 - t300) / 0.2e1 + t291 * t41 + (t241 - t362) * qJD(5) + (-t123 * mrSges(5,2) + t35 * mrSges(5,3) - t239 + t362) * t264 + (-t12 * t15 - t13 * t16 + t227 * t8 - t31 * t41) * m(6) + t227 * t14 + t10 * mrSges(5,1) - t11 * mrSges(5,2) + (-t242 + t370) * t246 - t16 * t68 - t15 * t69 - t74 * mrSges(4,2) + t75 * mrSges(4,1) - t42 * t108 - t117 * t150 + t118 * t151 - t164 * (mrSges(4,1) * t178 + mrSges(4,2) * t177) + t6 * t328 + (Ifges(4,5) * t177 - Ifges(4,6) * t178) * t331 + t113 * t334 + t251 * t349 + t253 * t352 + t255 * t353 + t233 * t356; t236 * t23 + t233 * t24 + t291 * t246 + t248 * qJD(5) + (-t108 - t248) * t264 + t48 + (-t122 * t249 - t246 * t31 - t259) * m(6) + (t246 * t35 - t264 * t36 + t115) * m(5); -t31 * (mrSges(6,1) * t103 + mrSges(6,2) * t102) + (Ifges(6,1) * t102 - t314) * t344 + t29 * t343 + (Ifges(6,5) * t102 - Ifges(6,6) * t103) * t342 - t12 * t68 + t13 * t69 + (t102 * t12 + t103 * t13) * mrSges(6,3) + t260 + t5 + (-Ifges(6,2) * t103 + t30 + t99) * t346;];
tauc = t25(:);
