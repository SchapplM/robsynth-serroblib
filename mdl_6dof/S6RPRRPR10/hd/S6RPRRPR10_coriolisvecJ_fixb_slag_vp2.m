% Calculate vector of centrifugal and coriolis load on the joints for
% S6RPRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6]';
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
% Datum: 2018-11-23 16:21
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RPRRPR10_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR10_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR10_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPR10_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR10_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR10_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR10_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:21:27
% EndTime: 2018-11-23 16:21:35
% DurationCPUTime: 7.66s
% Computational Cost: add. (5558->571), mult. (11879->751), div. (0->0), fcn. (6889->6), ass. (0->274)
t383 = Ifges(5,1) + Ifges(6,1);
t374 = Ifges(6,4) + Ifges(5,5);
t387 = Ifges(5,6) - Ifges(6,6);
t192 = sin(qJ(4));
t266 = qJD(4) * t192;
t193 = sin(qJ(3));
t270 = qJD(1) * t193;
t326 = pkin(8) - pkin(9);
t196 = cos(qJ(3));
t269 = qJD(1) * t196;
t245 = pkin(3) * t196 + pkin(8) * t193;
t160 = t245 * qJD(1);
t195 = cos(qJ(4));
t198 = -pkin(1) - pkin(7);
t181 = qJD(1) * t198 + qJD(2);
t275 = t181 * t196;
t92 = t192 * t160 + t195 * t275;
t77 = qJ(5) * t269 + t92;
t386 = -pkin(9) * t192 * t270 + t326 * t266 + t77;
t174 = t326 * t195;
t327 = pkin(4) + pkin(5);
t209 = pkin(9) * t193 * t195 - t196 * t327;
t274 = t192 * t196;
t91 = t160 * t195 - t181 * t274;
t385 = -qJD(1) * t209 + qJD(4) * t174 + t91;
t166 = pkin(3) * t193 - pkin(8) * t196 + qJ(2);
t142 = t166 * qJD(1);
t163 = t193 * t181;
t144 = qJD(3) * pkin(8) + t163;
t74 = t195 * t142 - t192 * t144;
t347 = qJD(5) - t74;
t145 = -qJD(3) * pkin(3) - t275;
t75 = t192 * t142 + t195 * t144;
t220 = t192 * t75 + t195 * t74;
t185 = qJD(4) + t270;
t60 = -pkin(4) * t185 + t347;
t177 = t185 * qJ(5);
t61 = t177 + t75;
t222 = t192 * t61 - t195 * t60;
t294 = Ifges(6,5) * t195;
t229 = Ifges(6,3) * t192 + t294;
t299 = Ifges(5,4) * t195;
t235 = -Ifges(5,2) * t192 + t299;
t240 = mrSges(6,1) * t192 - mrSges(6,3) * t195;
t242 = mrSges(5,1) * t192 + mrSges(5,2) * t195;
t292 = Ifges(6,6) * t192;
t293 = Ifges(5,6) * t192;
t297 = Ifges(5,5) * t195;
t298 = Ifges(6,4) * t195;
t311 = t195 / 0.2e1;
t313 = t192 / 0.2e1;
t314 = -t192 / 0.2e1;
t315 = t185 / 0.2e1;
t156 = qJD(3) * t192 + t195 * t269;
t319 = t156 / 0.2e1;
t261 = t195 * qJD(3);
t155 = t192 * t269 - t261;
t321 = t155 / 0.2e1;
t322 = -t155 / 0.2e1;
t295 = Ifges(6,5) * t192;
t300 = Ifges(5,4) * t192;
t345 = t195 * t383 + t295 - t300;
t150 = Ifges(5,4) * t155;
t296 = Ifges(6,5) * t155;
t354 = t156 * t383 + t374 * t185 - t150 + t296;
t208 = qJ(5) * t156 - t145;
t65 = pkin(4) * t155 - t208;
t149 = Ifges(6,5) * t156;
t67 = t185 * Ifges(6,6) + t155 * Ifges(6,3) + t149;
t301 = Ifges(5,4) * t156;
t70 = -t155 * Ifges(5,2) + t185 * Ifges(5,6) + t301;
t384 = t311 * t354 + t313 * t67 + t314 * t70 + t145 * t242 + t229 * t321 + t235 * t322 + t65 * t240 + (t292 + t298 - t293 + t297) * t315 + t345 * t319 - t222 * mrSges(6,2) - t220 * mrSges(5,3);
t191 = sin(qJ(6));
t194 = cos(qJ(6));
t218 = t155 * t191 + t156 * t194;
t86 = t155 * t194 - t156 * t191;
t79 = Ifges(7,4) * t86;
t382 = Ifges(7,2) * t218 - t79;
t381 = -pkin(9) * t156 + t347;
t380 = qJD(5) * t192 + t163;
t277 = qJ(5) * t195;
t214 = t192 * t327 - t277;
t379 = t198 - t214;
t281 = Ifges(4,5) * qJD(3);
t303 = Ifges(4,4) * t193;
t376 = qJD(1) / 0.2e1;
t378 = t281 / 0.2e1 + (t196 * Ifges(4,1) - t303) * t376 + t384;
t267 = qJD(3) * t196;
t248 = qJD(1) * t267;
t180 = qJD(6) - t185;
t318 = -t180 / 0.2e1;
t50 = -t155 * t327 + t208;
t35 = -t185 * t327 + t381;
t57 = pkin(9) * t155 + t75;
t39 = t177 + t57;
t8 = -t191 * t39 + t194 * t35;
t9 = t191 * t35 + t194 * t39;
t377 = (t218 * t9 + t8 * t86) * mrSges(7,3) - Ifges(7,3) * t248 + (Ifges(7,5) * t86 - Ifges(7,6) * t218) * t318 - t50 * (mrSges(7,1) * t218 + mrSges(7,2) * t86);
t375 = -qJD(3) / 0.2e1;
t173 = t326 * t192;
t106 = t173 * t191 + t174 * t194;
t373 = -qJD(6) * t106 + t191 * t386 + t194 * t385;
t105 = t173 * t194 - t174 * t191;
t372 = qJD(6) * t105 + t191 * t385 - t194 * t386;
t264 = qJD(4) * t196;
t252 = t192 * t264;
t107 = qJD(4) * t261 + (-t193 * t261 - t252) * qJD(1);
t268 = qJD(3) * t193;
t254 = t192 * t268;
t108 = -qJD(1) * t254 + qJD(4) * t156;
t371 = t374 * t248 + (-Ifges(5,4) + Ifges(6,5)) * t108 + t383 * t107;
t370 = -t185 * t214 + t380;
t226 = pkin(4) * t192 - t277;
t369 = t185 * t226 - t380;
t217 = t191 * t195 - t192 * t194;
t344 = qJD(4) - qJD(6);
t368 = t344 * t217;
t154 = qJD(3) * t245 + qJD(2);
t273 = t193 * t198;
t367 = -qJD(4) * t273 + t154;
t366 = t192 * t374 + t195 * t387;
t365 = t192 * t383 - t294 + t299;
t129 = t154 * qJD(1);
t255 = t181 * t267;
t265 = qJD(4) * t195;
t33 = t192 * t129 + t142 * t265 - t144 * t266 + t195 * t255;
t19 = qJ(5) * t248 + t185 * qJD(5) + t33;
t10 = pkin(9) * t108 + t19;
t34 = t129 * t195 - t142 * t266 - t144 * t265 - t192 * t255;
t11 = -pkin(9) * t107 - t248 * t327 - t34;
t1 = qJD(6) * t8 + t10 * t194 + t11 * t191;
t2 = -qJD(6) * t9 - t10 * t191 + t11 * t194;
t22 = qJD(6) * t86 + t107 * t194 + t108 * t191;
t23 = -qJD(6) * t218 - t107 * t191 + t108 * t194;
t364 = t2 * mrSges(7,1) - t1 * mrSges(7,2) + Ifges(7,5) * t22 + Ifges(7,6) * t23;
t363 = qJD(1) * (qJ(2) * (m(3) + m(4)) + mrSges(3,3));
t310 = Ifges(7,4) * t218;
t360 = Ifges(7,1) * t86 - t310;
t334 = t22 / 0.2e1;
t333 = t23 / 0.2e1;
t28 = Ifges(7,1) * t218 + Ifges(7,5) * t180 + t79;
t358 = t28 / 0.2e1;
t330 = -t218 / 0.2e1;
t357 = -t248 / 0.2e1;
t249 = Ifges(4,6) * t375;
t165 = qJ(5) * t194 - t191 * t327;
t350 = -qJD(6) * t165 - t191 * t381 - t194 * t57;
t164 = -qJ(5) * t191 - t194 * t327;
t349 = qJD(6) * t164 - t191 * t57 + t194 * t381;
t348 = t217 * t193;
t223 = -t192 * t34 + t195 * t33;
t24 = -pkin(4) * t248 - t34;
t225 = t19 * t195 + t192 * t24;
t343 = t34 * mrSges(5,1) - t24 * mrSges(6,1) - t33 * mrSges(5,2) + t19 * mrSges(6,3) - t364;
t282 = -qJD(3) * mrSges(4,1) + mrSges(5,1) * t155 + mrSges(5,2) * t156 + mrSges(4,3) * t269;
t342 = m(5) * t145 + t282;
t278 = qJ(5) * t192;
t340 = -t195 * t327 - t278;
t304 = mrSges(5,3) * t156;
t110 = mrSges(5,1) * t185 - t304;
t111 = -mrSges(6,1) * t185 + mrSges(6,2) * t156;
t271 = -t110 + t111;
t305 = mrSges(5,3) * t155;
t109 = -mrSges(5,2) * t185 - t305;
t112 = -mrSges(6,2) * t155 + mrSges(6,3) * t185;
t272 = t109 + t112;
t206 = -t192 * t272 + t195 * t271;
t339 = -m(5) * t220 - m(6) * t222 + t206;
t257 = Ifges(6,6) / 0.2e1 - Ifges(5,6) / 0.2e1;
t258 = Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1;
t259 = Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1;
t302 = Ifges(4,4) * t196;
t338 = -t257 * t155 - t259 * t156 - t258 * t185 - t61 * mrSges(6,3) - t74 * mrSges(5,1) - t9 * mrSges(7,2) - t249 + (-Ifges(4,2) * t193 + t302) * t376 + t180 * Ifges(7,3) + t218 * Ifges(7,5) + t86 * Ifges(7,6) - Ifges(5,6) * t322 - Ifges(6,6) * t321 + t60 * mrSges(6,1) + t75 * mrSges(5,2) + t8 * mrSges(7,1) - t374 * t319 - (Ifges(5,3) + Ifges(6,2)) * t315;
t336 = Ifges(7,4) * t334 + Ifges(7,2) * t333 + Ifges(7,6) * t357;
t335 = Ifges(7,1) * t334 + Ifges(7,4) * t333 + Ifges(7,5) * t357;
t332 = -t86 / 0.2e1;
t331 = t86 / 0.2e1;
t329 = t218 / 0.2e1;
t325 = t107 / 0.2e1;
t324 = -t108 / 0.2e1;
t323 = t108 / 0.2e1;
t320 = -t156 / 0.2e1;
t317 = t180 / 0.2e1;
t316 = -t185 / 0.2e1;
t312 = -t195 / 0.2e1;
t36 = -mrSges(7,1) * t86 + mrSges(7,2) * t218;
t94 = mrSges(6,1) * t155 - mrSges(6,3) * t156;
t306 = -t36 + t94;
t291 = qJ(2) * mrSges(4,1);
t290 = qJ(2) * mrSges(4,2);
t117 = qJD(1) * t348;
t284 = t117 + t368;
t216 = t191 * t192 + t194 * t195;
t138 = t216 * qJD(1);
t118 = t193 * t138;
t89 = t344 * t216;
t283 = -t118 - t89;
t279 = qJ(5) * t155;
t276 = qJD(3) * mrSges(4,2);
t114 = t192 * t166 + t195 * t273;
t262 = qJD(5) * t195;
t260 = qJD(1) * qJD(2);
t98 = t193 * qJ(5) + t114;
t256 = t181 * t268;
t253 = t198 * t267;
t250 = -t281 / 0.2e1;
t178 = t192 * t273;
t113 = t166 * t195 - t178;
t243 = mrSges(5,1) * t195 - mrSges(5,2) * t192;
t241 = mrSges(6,1) * t195 + mrSges(6,3) * t192;
t234 = Ifges(5,2) * t195 + t300;
t228 = -Ifges(6,3) * t195 + t295;
t227 = pkin(4) * t195 + t278;
t62 = -mrSges(7,2) * t180 + mrSges(7,3) * t86;
t63 = mrSges(7,1) * t180 - mrSges(7,3) * t218;
t224 = -t191 * t63 + t194 * t62;
t73 = t178 + (-pkin(9) * t196 - t166) * t195 - t327 * t193;
t80 = pkin(9) * t274 + t98;
t30 = -t191 * t80 + t194 * t73;
t31 = t191 * t73 + t194 * t80;
t221 = t192 * t60 + t195 * t61;
t219 = t192 * t74 - t195 * t75;
t55 = -t166 * t266 - t192 * t253 + t195 * t367;
t83 = -mrSges(6,1) * t248 + t107 * mrSges(6,2);
t215 = -t198 + t226;
t132 = t217 * t196;
t133 = t216 * t196;
t81 = -mrSges(6,2) * t108 + mrSges(6,3) * t248;
t82 = mrSges(5,1) * t248 - mrSges(5,3) * t107;
t84 = -mrSges(5,2) * t248 - mrSges(5,3) * t108;
t207 = (t81 + t84) * t195 + (-t82 + t83) * t192;
t54 = t166 * t265 + t192 * t367 + t195 * t253;
t205 = qJ(5) * t107 + qJD(5) * t156 - t256;
t38 = qJ(5) * t267 + t193 * qJD(5) + t54;
t184 = Ifges(6,2) * t248;
t183 = Ifges(5,3) * t248;
t171 = -mrSges(4,3) * t270 - t276;
t167 = -pkin(3) - t227;
t159 = (mrSges(4,1) * t193 + mrSges(4,2) * t196) * qJD(1);
t148 = pkin(3) - t340;
t139 = t217 * qJD(1);
t131 = t216 * t193;
t119 = t215 * t196;
t104 = Ifges(6,4) * t107;
t103 = Ifges(5,5) * t107;
t102 = Ifges(5,6) * t108;
t101 = Ifges(6,6) * t108;
t99 = -pkin(4) * t193 - t113;
t97 = t379 * t196;
t93 = pkin(4) * t156 + t279;
t78 = -pkin(4) * t269 - t91;
t66 = -t156 * t327 - t279;
t59 = (qJD(4) * t227 - t262) * t196 - t215 * t268;
t53 = mrSges(5,1) * t108 + mrSges(5,2) * t107;
t52 = mrSges(6,1) * t108 - mrSges(6,3) * t107;
t51 = -pkin(4) * t267 - t55;
t45 = t107 * Ifges(5,4) - t108 * Ifges(5,2) + Ifges(5,6) * t248;
t44 = t107 * Ifges(6,5) + Ifges(6,6) * t248 + t108 * Ifges(6,3);
t43 = t196 * t368 - t216 * t268;
t42 = qJD(3) * t348 + t196 * t89;
t41 = -qJD(3) * t132 + t193 * t89;
t40 = qJD(3) * t133 + t193 * t368;
t37 = (qJD(4) * t340 + t262) * t196 - t379 * t268;
t32 = (t195 * t264 - t254) * pkin(9) + t38;
t29 = pkin(9) * t252 + qJD(3) * t209 - t55;
t27 = Ifges(7,2) * t86 + Ifges(7,6) * t180 + t310;
t25 = pkin(4) * t108 - t205;
t18 = mrSges(7,2) * t248 + mrSges(7,3) * t23;
t17 = -mrSges(7,1) * t248 - mrSges(7,3) * t22;
t14 = -t108 * t327 + t205;
t7 = -mrSges(7,1) * t23 + mrSges(7,2) * t22;
t4 = -qJD(6) * t31 - t191 * t32 + t194 * t29;
t3 = qJD(6) * t30 + t191 * t29 + t194 * t32;
t5 = [m(6) * (t119 * t25 + t19 * t98 + t24 * t99 + t38 * t61 + t51 * t60 + t59 * t65) + m(7) * (t1 * t31 + t14 * t97 + t2 * t30 + t3 * t9 + t37 * t50 + t4 * t8) + (Ifges(7,4) * t133 - Ifges(7,2) * t132) * t333 + (Ifges(7,1) * t133 - Ifges(7,4) * t132) * t334 + (-t1 * t132 - t133 * t2 + t42 * t9 - t43 * t8) * mrSges(7,3) + t14 * (mrSges(7,1) * t132 + mrSges(7,2) * t133) + (Ifges(7,1) * t43 + Ifges(7,4) * t42) * t329 + (Ifges(7,4) * t43 + Ifges(7,2) * t42) * t331 + t133 * t335 - t132 * t336 + (Ifges(7,5) * t43 + Ifges(7,6) * t42) * t317 + m(5) * (t113 * t34 + t114 * t33 + t54 * t75 + t55 * t74) + (mrSges(4,1) * t260 + t259 * t107 + ((0.3e1 / 0.2e1 * t303 - 0.2e1 * t290) * qJD(1) + t342 * t198 + t250 - t378) * qJD(3) + t257 * t108 + t183 / 0.2e1 + t184 / 0.2e1 + t103 / 0.2e1 + t104 / 0.2e1 + t101 / 0.2e1 - t102 / 0.2e1 + t343) * t193 + t43 * t358 + (t159 + 0.2e1 * t363) * qJD(2) + t30 * t17 + t31 * t18 + t37 * t36 + t42 * t27 / 0.2e1 + t50 * (-mrSges(7,1) * t42 + mrSges(7,2) * t43) + t3 * t62 + t4 * t63 + (t25 * t240 + t229 * t323 + t235 * t324 + t44 * t313 - t198 * t53 + mrSges(4,2) * t260 + (-t192 * t33 - t195 * t34) * mrSges(5,3) + (-t19 * t192 + t195 * t24) * mrSges(6,2) + (-mrSges(6,2) * t221 + mrSges(5,3) * t219 + t145 * t243 + t228 * t322 + t234 * t321 + t241 * t65 + t312 * t70 + t316 * t366 + t320 * t365) * qJD(4) + (t198 * t171 + (-m(5) * t198 + t242) * t163 + (0.2e1 * t291 + (t298 / 0.2e1 + t292 / 0.2e1 + t297 / 0.2e1 - t293 / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(4,4)) * t196 - Ifges(7,5) * t133 / 0.2e1 + Ifges(7,6) * t132 / 0.2e1 + (-0.3e1 / 0.2e1 * Ifges(4,1) + 0.3e1 / 0.2e1 * Ifges(4,2) + Ifges(7,3) + t258) * t193) * qJD(1) + t249 - t338) * qJD(3) + t345 * t325 + (qJD(4) * t354 + t45) * t314 + (qJD(4) * t67 + t371) * t311) * t196 + t59 * t94 + t97 * t7 + t98 * t81 + t99 * t83 + t54 * t109 + t55 * t110 + t51 * t111 + t38 * t112 + t113 * t82 + t114 * t84 + t119 * t52; -t348 * t17 + t131 * t18 + (t41 + t138) * t63 + (t40 + t139) * t62 + m(7) * (t1 * t131 - t2 * t348 + t40 * t9 + t41 * t8) - m(7) * (-t138 * t8 - t139 * t9) + (t7 - t52 - t53 - m(6) * t25 + m(7) * t14 + (-m(5) * t219 + m(6) * t221 + t192 * t271 + t195 * t272 + t171) * qJD(3)) * t196 + (t206 * qJD(4) + m(6) * (t265 * t60 - t266 * t61 + t225) + m(5) * (-t265 * t74 - t266 * t75 + t223 - t255) + t207 + (m(6) * t65 - m(7) * t50 + t306 + t342) * qJD(3)) * t193 + (-t159 + t339 - t363) * qJD(1); (-t368 / 0.2e1 - t117 / 0.2e1) * t27 + (Ifges(7,1) * t89 - Ifges(7,4) * t368) * t329 + t228 * t323 + t234 * t324 - t216 * t336 - t217 * t335 + (-Ifges(7,4) * t217 - Ifges(7,2) * t216) * t333 + (-Ifges(7,1) * t217 - Ifges(7,4) * t216) * t334 + (-t1 * t216 + t2 * t217 + t283 * t8 - t284 * t9) * mrSges(7,3) + t14 * (mrSges(7,1) * t216 - mrSges(7,2) * t217) + t207 * pkin(8) + (Ifges(7,4) * t89 - Ifges(7,2) * t368) * t331 + (Ifges(7,5) * t89 - Ifges(7,6) * t368) * t317 + (-Ifges(7,1) * t118 + Ifges(7,4) * t117) * t330 + (-Ifges(7,4) * t118 + Ifges(7,2) * t117) * t332 + (-Ifges(7,5) * t118 + Ifges(7,6) * t117) * t318 + (t89 / 0.2e1 + t118 / 0.2e1) * t28 + (((-Ifges(7,5) * t217 - Ifges(7,6) * t216) * t375 + t249 + (t302 / 0.2e1 - t291) * qJD(1) + t366 * qJD(3) / 0.2e1 + t338) * t196 + ((t290 - t303 / 0.2e1 + (-Ifges(4,2) / 0.2e1 + Ifges(4,1) / 0.2e1) * t196) * qJD(1) + t250 + t378) * t193) * qJD(1) + (mrSges(7,1) * t284 - mrSges(7,2) * t283) * t50 + ((-t171 - t276) * t196 + ((-mrSges(4,1) - t243) * qJD(3) - t282) * t193) * t181 + t45 * t311 + t44 * t312 + t365 * t325 + (pkin(8) * t339 + t384) * qJD(4) + t223 * mrSges(5,3) + t225 * mrSges(6,2) + t369 * t94 + (pkin(8) * t225 + t167 * t25 + t369 * t65 - t60 * t78 - t61 * t77) * m(6) + t370 * t36 - pkin(3) * t53 - t25 * t241 + t371 * t313 + t372 * t62 + t373 * t63 + (t1 * t106 + t105 * t2 + t14 * t148 + t370 * t50 + t372 * t9 + t373 * t8) * m(7) + t105 * t17 + t106 * t18 - t92 * t109 - t91 * t110 - t78 * t111 - t77 * t112 + t148 * t7 + t167 * t52 + (-pkin(3) * t256 + pkin(8) * t223 - t145 * t163 - t74 * t91 - t75 * t92) * m(5); (t155 * t60 + t156 * t61) * mrSges(6,2) + t70 * t319 + (Ifges(6,3) * t156 - t296) * t322 + (t27 - t360) * t330 + t86 * t358 + t343 + t183 + t184 - t377 + t103 + t104 + t101 - t102 + (-Ifges(5,2) * t156 - t150 + t354) * t321 + t349 * t62 + t350 * t63 + (t1 * t165 + t164 * t2 + t349 * t9 + t350 * t8 - t50 * t66) * m(7) + (-pkin(4) * t24 + qJ(5) * t19 + t347 * t61 - t60 * t75 - t65 * t93) * m(6) + (-t155 * t383 + t149 - t301 + t67) * t320 + t382 * t332 + (-t271 + t304) * t75 + (-t272 - t305) * t74 - t66 * t36 + (-t374 * t155 - t156 * t387) * t316 + qJ(5) * t81 - pkin(4) * t83 - t93 * t94 + qJD(5) * t112 - t65 * (mrSges(6,1) * t156 + mrSges(6,3) * t155) - t145 * (mrSges(5,1) * t156 - mrSges(5,2) * t155) + t164 * t17 + t165 * t18; t194 * t17 + t191 * t18 + t306 * t156 + t224 * qJD(6) + (-t112 - t224) * t185 + t83 + (t1 * t191 - t156 * t50 + t194 * t2 + t180 * (-t191 * t8 + t194 * t9)) * m(7) + (t156 * t65 - t185 * t61 + t24) * m(6); t360 * t330 + t27 * t329 - t8 * t62 + t9 * t63 + (t28 - t382) * t332 + t364 + t377;];
tauc  = t5(:);
