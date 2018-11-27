% Calculate joint inertia matrix for
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
% m_mdh [8x1]
%   mass of all robot links (including the base)
% mrSges [8x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [8x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [7x7]
%   inertia matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-26 21:21
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S7RRRRRRR1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(4,1),zeros(8,1),zeros(8,3),zeros(8,6)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_inertiaJ_slag_vp2: qJ has to be [7x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_inertiaJ_slag_vp2: pkin has to be [4x1] (double)');
assert( isreal(m) && all(size(m) == [8 1]), ...
  'S7RRRRRRR1_inertiaJ_slag_vp2: m has to be [8x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [8,3]), ...
  'S7RRRRRRR1_inertiaJ_slag_vp2: mrSges has to be [8x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [8 6]), ...
  'S7RRRRRRR1_inertiaJ_slag_vp2: Ifges has to be [8x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-26 19:14:29
% EndTime: 2018-11-26 19:14:51
% DurationCPUTime: 7.68s
% Computational Cost: add. (4602->748), mult. (12544->1037), div. (0->0), fcn. (13896->12), ass. (0->300)
t375 = 2 * pkin(2);
t224 = sin(qJ(5));
t230 = cos(qJ(5));
t232 = cos(qJ(3));
t226 = sin(qJ(3));
t231 = cos(qJ(4));
t289 = t226 * t231;
t146 = t224 * t289 - t230 * t232;
t151 = -t224 * t232 - t230 * t289;
t225 = sin(qJ(4));
t292 = t225 * t226;
t380 = (-mrSges(5,1) * t232 + mrSges(6,1) * t146 - mrSges(6,2) * t151 + mrSges(5,3) * t289) * t225 - (mrSges(5,2) * t232 + mrSges(5,3) * t292) * t231;
t223 = sin(qJ(6));
t297 = t223 * t224;
t158 = mrSges(7,2) * t230 - mrSges(7,3) * t297;
t229 = cos(qJ(6));
t222 = sin(qJ(7));
t228 = cos(qJ(7));
t293 = t224 * t229;
t144 = -t222 * t293 + t228 * t230;
t149 = t222 * t230 + t228 * t293;
t306 = -mrSges(7,1) * t230 - mrSges(8,1) * t144 + mrSges(8,2) * t149 - mrSges(7,3) * t293;
t379 = t158 * t229 - t306 * t223;
t291 = t225 * t230;
t145 = t223 * t291 + t229 * t231;
t294 = t224 * t225;
t109 = mrSges(7,2) * t294 + mrSges(7,3) * t145;
t159 = -mrSges(6,2) * t231 + mrSges(6,3) * t294;
t287 = t229 * t230;
t150 = t223 * t231 - t225 * t287;
t101 = -t150 * t222 + t228 * t294;
t103 = t150 * t228 + t222 * t294;
t307 = -mrSges(7,1) * t294 - mrSges(8,1) * t101 + mrSges(8,2) * t103 - mrSges(7,3) * t150;
t378 = t109 * t229 - t307 * t223 + t159;
t233 = cos(qJ(2));
t227 = sin(qJ(2));
t288 = t227 * t232;
t147 = t225 * t288 + t231 * t233;
t290 = t226 * t227;
t376 = mrSges(4,2) * t233 + mrSges(4,3) * t290 - (mrSges(5,2) * t290 - mrSges(5,3) * t147) * t231;
t260 = -pkin(2) * t231 - pkin(3);
t132 = t260 * t290;
t148 = -t225 * t233 + t231 * t288;
t240 = -pkin(2) * t288 - pkin(3) * t148;
t58 = t132 * t224 - t230 * t240;
t374 = t58 ^ 2;
t170 = t260 * t232;
t248 = (pkin(3) * t231 + pkin(2)) * t226;
t89 = t170 * t224 - t230 * t248;
t373 = t89 ^ 2;
t105 = t148 * t224 + t230 * t290;
t106 = t148 * t230 - t224 * t290;
t56 = t106 * t229 + t147 * t223;
t24 = -t105 * t228 - t222 * t56;
t372 = t24 / 0.2e1;
t25 = -t105 * t222 + t228 * t56;
t371 = t25 / 0.2e1;
t102 = -t151 * t223 - t229 * t292;
t104 = t151 * t229 - t223 * t292;
t37 = Ifges(7,1) * t104 + Ifges(7,4) * t102 - Ifges(7,5) * t146;
t370 = t37 / 0.2e1;
t53 = -t104 * t222 + t146 * t228;
t369 = t53 / 0.2e1;
t54 = t104 * t228 + t146 * t222;
t368 = t54 / 0.2e1;
t75 = Ifges(8,4) * t149 + Ifges(8,2) * t144 - Ifges(8,6) * t297;
t367 = t75 / 0.2e1;
t79 = Ifges(8,1) * t149 + Ifges(8,4) * t144 - Ifges(8,5) * t297;
t366 = t79 / 0.2e1;
t80 = Ifges(7,1) * t150 + Ifges(7,4) * t145 - Ifges(7,5) * t294;
t365 = t80 / 0.2e1;
t81 = Ifges(6,1) * t151 + Ifges(6,4) * t146 - Ifges(6,5) * t292;
t364 = t81 / 0.2e1;
t82 = Ifges(5,1) * t148 - Ifges(5,4) * t147 - Ifges(5,5) * t290;
t363 = -t82 / 0.2e1;
t362 = t101 / 0.2e1;
t361 = t103 / 0.2e1;
t360 = t104 / 0.2e1;
t320 = Ifges(8,4) * t228;
t121 = Ifges(8,6) * t229 + (-Ifges(8,2) * t222 + t320) * t223;
t359 = t121 / 0.2e1;
t321 = Ifges(8,4) * t222;
t126 = Ifges(8,5) * t229 + (Ifges(8,1) * t228 - t321) * t223;
t358 = t126 / 0.2e1;
t323 = Ifges(7,4) * t223;
t127 = -Ifges(7,5) * t230 + (Ifges(7,1) * t229 - t323) * t224;
t357 = t127 / 0.2e1;
t325 = Ifges(6,4) * t224;
t128 = Ifges(6,5) * t231 + (-Ifges(6,1) * t230 + t325) * t225;
t356 = t128 / 0.2e1;
t355 = t144 / 0.2e1;
t354 = t148 / 0.2e1;
t353 = t149 / 0.2e1;
t352 = t150 / 0.2e1;
t351 = t151 / 0.2e1;
t176 = -Ifges(8,5) * t222 - Ifges(8,6) * t228;
t350 = t176 / 0.2e1;
t177 = Ifges(7,5) * t223 + Ifges(7,6) * t229;
t349 = -t177 / 0.2e1;
t178 = Ifges(6,5) * t224 + Ifges(6,6) * t230;
t348 = t178 / 0.2e1;
t180 = -Ifges(8,2) * t228 - t321;
t347 = t180 / 0.2e1;
t185 = -Ifges(8,1) * t222 - t320;
t346 = t185 / 0.2e1;
t322 = Ifges(7,4) * t229;
t186 = Ifges(7,1) * t223 + t322;
t345 = t186 / 0.2e1;
t324 = Ifges(6,4) * t230;
t187 = Ifges(6,1) * t224 + t324;
t344 = t187 / 0.2e1;
t343 = -t222 / 0.2e1;
t342 = -t225 / 0.2e1;
t341 = -t228 / 0.2e1;
t340 = t228 / 0.2e1;
t339 = t229 / 0.2e1;
t338 = -t230 / 0.2e1;
t337 = mrSges(7,1) * t105 - mrSges(8,1) * t24 + mrSges(8,2) * t25 - mrSges(7,3) * t56;
t336 = pkin(2) * t232;
t335 = pkin(3) * t230;
t334 = t58 * t89;
t333 = -mrSges(7,1) * t146 - mrSges(8,1) * t53 + mrSges(8,2) * t54 - mrSges(7,3) * t104;
t55 = -t106 * t223 + t147 * t229;
t332 = -mrSges(6,1) * t147 - mrSges(7,1) * t55 + mrSges(7,2) * t56 + mrSges(6,3) * t106;
t331 = mrSges(5,2) * t231;
t330 = mrSges(7,2) * t229;
t329 = Ifges(4,4) * t226;
t328 = Ifges(4,4) * t232;
t327 = Ifges(5,4) * t225;
t326 = Ifges(5,4) * t231;
t319 = Ifges(5,5) * t225;
t318 = Ifges(4,6) * t233;
t317 = t224 * t58;
t316 = t224 * t89;
t315 = t230 * t58;
t314 = t230 * t89;
t313 = t233 * Ifges(4,5);
t312 = mrSges(5,1) * t231 - mrSges(5,2) * t225 + mrSges(4,1);
t311 = mrSges(7,1) * t229 - mrSges(7,2) * t223 + mrSges(6,1);
t310 = -mrSges(6,1) * t292 + mrSges(7,1) * t102 - mrSges(7,2) * t104 - mrSges(6,3) * t151;
t309 = -mrSges(4,1) * t233 - mrSges(5,1) * t147 - mrSges(5,2) * t148 + mrSges(4,3) * t288;
t172 = mrSges(8,1) * t228 - mrSges(8,2) * t222;
t308 = t172 + mrSges(7,1);
t305 = -mrSges(6,1) * t231 - mrSges(7,1) * t145 + mrSges(7,2) * t150 - mrSges(6,3) * t291;
t213 = t224 ^ 2;
t235 = pkin(3) ^ 2;
t209 = t213 * t235;
t215 = t226 ^ 2;
t236 = pkin(2) ^ 2;
t210 = t215 * t236;
t219 = t230 ^ 2;
t300 = t219 * t235;
t221 = t232 ^ 2;
t299 = t221 * t236;
t298 = t222 * t223;
t296 = t223 * t228;
t295 = t223 * t230;
t286 = t230 * t235;
t285 = Ifges(7,5) * t150 + Ifges(7,6) * t145;
t284 = Ifges(5,5) * t148 - Ifges(5,6) * t147;
t283 = Ifges(6,6) * t294 + Ifges(6,3) * t231;
t282 = Ifges(8,5) * t296 + Ifges(8,3) * t229;
t281 = Ifges(4,5) * t288 + Ifges(4,3) * t233;
t280 = t222 ^ 2 + t228 ^ 2;
t279 = 2 * pkin(4);
t278 = -Ifges(4,6) + Ifges(5,6) * t231 / 0.2e1 + t319 / 0.2e1;
t1 = Ifges(8,5) * t25 + Ifges(8,6) * t24 + Ifges(8,3) * t55;
t9 = Ifges(8,5) * t54 + Ifges(8,6) * t53 + Ifges(8,3) * t102;
t10 = Ifges(7,5) * t56 + Ifges(7,6) * t55 + Ifges(7,3) * t105;
t277 = t225 * t336;
t276 = pkin(2) * t290;
t12 = Ifges(7,4) * t56 + Ifges(7,2) * t55 + Ifges(7,6) * t105;
t275 = t1 / 0.2e1 + t12 / 0.2e1;
t34 = Ifges(7,4) * t104 + Ifges(7,2) * t102 - Ifges(7,6) * t146;
t274 = t9 / 0.2e1 + t34 / 0.2e1;
t30 = Ifges(8,5) * t103 + Ifges(8,6) * t101 + Ifges(8,3) * t145;
t31 = Ifges(7,5) * t104 + Ifges(7,6) * t102 - Ifges(7,3) * t146;
t251 = t225 * t276;
t61 = t230 * t132 + t224 * t240;
t40 = -t223 * t61 - t229 * t251;
t273 = t40 * t297;
t92 = t230 * t170 + t224 * t248;
t69 = -t223 * t92 - t229 * t277;
t272 = t69 * t297;
t271 = t40 * t295;
t270 = t69 * t295;
t35 = Ifges(6,4) * t106 - Ifges(6,2) * t105 + Ifges(6,6) * t147;
t269 = -t10 / 0.2e1 + t35 / 0.2e1;
t76 = Ifges(7,4) * t150 + Ifges(7,2) * t145 - Ifges(7,6) * t294;
t268 = t30 / 0.2e1 + t76 / 0.2e1;
t77 = Ifges(6,4) * t151 + Ifges(6,2) * t146 - Ifges(6,6) * t292;
t267 = t31 / 0.2e1 - t77 / 0.2e1;
t32 = Ifges(6,5) * t106 - Ifges(6,6) * t105 + Ifges(6,3) * t147;
t78 = Ifges(5,4) * t148 - Ifges(5,2) * t147 - Ifges(5,6) * t290;
t266 = t32 / 0.2e1 - t78 / 0.2e1;
t214 = t225 ^ 2;
t265 = t214 * t209;
t216 = t227 ^ 2;
t264 = t216 * t210;
t122 = -Ifges(7,6) * t230 + (-Ifges(7,2) * t223 + t322) * t224;
t71 = Ifges(8,5) * t149 + Ifges(8,6) * t144 - Ifges(8,3) * t297;
t263 = t71 / 0.2e1 + t122 / 0.2e1;
t123 = Ifges(6,6) * t231 + (Ifges(6,2) * t224 - t324) * t225;
t72 = -Ifges(7,3) * t294 + t285;
t262 = t72 / 0.2e1 - t123 / 0.2e1;
t124 = -Ifges(5,6) * t232 + (Ifges(5,2) * t225 - t326) * t226;
t73 = Ifges(6,5) * t151 + Ifges(6,6) * t146 - Ifges(6,3) * t292;
t261 = t73 / 0.2e1 - t124 / 0.2e1;
t259 = pkin(3) * t229 + pkin(4);
t117 = -Ifges(8,6) * t298 + t282;
t181 = Ifges(7,2) * t229 + t323;
t258 = t117 / 0.2e1 + t181 / 0.2e1;
t195 = Ifges(7,5) * t293;
t118 = -Ifges(7,6) * t297 - Ifges(7,3) * t230 + t195;
t182 = Ifges(6,2) * t230 + t325;
t257 = t118 / 0.2e1 - t182 / 0.2e1;
t119 = -Ifges(6,5) * t291 + t283;
t183 = -Ifges(5,2) * t231 - t327;
t256 = t119 / 0.2e1 - t183 / 0.2e1;
t255 = (-mrSges(5,1) * t290 - mrSges(6,1) * t105 - mrSges(6,2) * t106 - mrSges(5,3) * t148) * t225;
t250 = t226 * t236 * t288;
t249 = m(8) * t280 * (pkin(4) ^ 2);
t20 = pkin(4) * t56 + t58;
t41 = -t223 * t251 + t229 * t61;
t21 = -pkin(4) * t105 + t41;
t4 = -t20 * t228 - t21 * t222;
t5 = -t20 * t222 + t21 * t228;
t247 = t222 * t4 - t228 * t5;
t220 = t231 ^ 2;
t246 = mrSges(4,2) + (t214 + t220) * mrSges(5,3);
t212 = t223 ^ 2;
t218 = t229 ^ 2;
t245 = -mrSges(6,2) + (t212 + t218) * mrSges(7,3);
t244 = mrSges(8,1) * t222 + mrSges(8,2) * t228;
t70 = -t223 * t277 + t229 * t92;
t42 = pkin(4) * t146 + t70;
t43 = pkin(4) * t104 + t89;
t15 = -t222 * t42 - t228 * t43;
t16 = -t222 * t43 + t228 * t42;
t243 = t15 * t222 - t16 * t228;
t116 = -pkin(3) * t291 + pkin(4) * t150;
t131 = t259 * t294;
t57 = -t116 * t228 - t131 * t222;
t60 = -t116 * t222 + t131 * t228;
t242 = t222 * t57 - t228 * t60;
t169 = t259 * t230;
t171 = (pkin(4) * t229 + pkin(3)) * t224;
t88 = -t169 * t222 - t171 * t228;
t91 = t169 * t228 - t171 * t222;
t241 = t222 * t88 - t228 * t91;
t28 = -mrSges(7,2) * t105 + mrSges(7,3) * t55;
t66 = -mrSges(6,2) * t147 - mrSges(6,3) * t105;
t239 = -t223 * t337 + t229 * t28 + t66;
t110 = mrSges(6,2) * t292 + mrSges(6,3) * t146;
t63 = mrSges(7,2) * t146 + mrSges(7,3) * t102;
t238 = -t223 * t333 + t229 * t63 + t110;
t237 = -Ifges(5,3) + (t331 + (t230 * mrSges(6,1) - mrSges(6,2) * t224 + mrSges(5,1)) * t225) * pkin(2);
t200 = t216 * t299;
t199 = t214 * t300;
t198 = t214 * t299;
t197 = t212 * t300;
t193 = Ifges(5,6) * t292;
t191 = t214 * t264;
t190 = t212 * t265;
t189 = -Ifges(4,1) * t226 - t328;
t188 = -Ifges(5,1) * t225 - t326;
t184 = -Ifges(4,2) * t232 - t329;
t168 = t214 * t250;
t167 = t212 * t286 * t294;
t162 = mrSges(8,1) * t229 - mrSges(8,3) * t296;
t157 = -mrSges(8,2) * t229 - mrSges(8,3) * t298;
t156 = (-mrSges(5,1) * t225 - t331) * t226;
t155 = (-mrSges(6,1) * t224 - mrSges(6,2) * t230) * t225;
t154 = (mrSges(7,1) * t223 + t330) * t224;
t153 = t244 * t223;
t130 = t313 + (Ifges(4,1) * t232 - t329) * t227;
t129 = -Ifges(5,5) * t232 + (-Ifges(5,1) * t231 + t327) * t226;
t125 = t318 + (-Ifges(4,2) * t226 + t328) * t227;
t120 = -Ifges(5,5) * t289 - Ifges(5,3) * t232 + t193;
t113 = -mrSges(8,1) * t297 - mrSges(8,3) * t149;
t108 = mrSges(8,2) * t297 + mrSges(8,3) * t144;
t74 = -Ifges(5,3) * t290 + t284;
t68 = t69 ^ 2;
t64 = mrSges(8,1) * t145 - mrSges(8,3) * t103;
t62 = -mrSges(8,2) * t145 + mrSges(8,3) * t101;
t39 = t40 ^ 2;
t38 = Ifges(6,1) * t106 - Ifges(6,4) * t105 + Ifges(6,5) * t147;
t36 = Ifges(8,1) * t103 + Ifges(8,4) * t101 + Ifges(8,5) * t145;
t33 = Ifges(8,4) * t103 + Ifges(8,2) * t101 + Ifges(8,6) * t145;
t27 = mrSges(8,1) * t102 - mrSges(8,3) * t54;
t26 = -mrSges(8,2) * t102 + mrSges(8,3) * t53;
t17 = t69 * t40;
t14 = Ifges(7,1) * t56 + Ifges(7,4) * t55 + Ifges(7,5) * t105;
t13 = Ifges(8,1) * t54 + Ifges(8,4) * t53 + Ifges(8,5) * t102;
t11 = Ifges(8,4) * t54 + Ifges(8,2) * t53 + Ifges(8,6) * t102;
t8 = mrSges(8,1) * t55 - mrSges(8,3) * t25;
t7 = -mrSges(8,2) * t55 + mrSges(8,3) * t24;
t3 = Ifges(8,1) * t25 + Ifges(8,4) * t24 + Ifges(8,5) * t55;
t2 = Ifges(8,4) * t25 + Ifges(8,2) * t24 + Ifges(8,6) * t55;
t6 = [t106 * t38 + t56 * t14 + t148 * t82 + t24 * t2 + t25 * t3 + 0.2e1 * t41 * t28 + 0.2e1 * t4 * t8 + 0.2e1 * t5 * t7 + 0.2e1 * t61 * t66 + Ifges(2,3) + (t1 + t12) * t55 + (Ifges(3,2) * t233 + t281) * t233 + (t32 - t78) * t147 + (t10 - t35) * t105 + m(8) * (t4 ^ 2 + t5 ^ 2 + t39) + m(7) * (t41 ^ 2 + t374 + t39) + m(6) * (t61 ^ 2 + t191 + t374) + m(4) * (t200 + t264) + m(5) * (t220 * t264 + t191 + t200) + 0.2e1 * t332 * t58 + 0.2e1 * t337 * t40 + (0.2e1 * Ifges(3,4) * t233 + t232 * t130 + (t375 * t376 - t125 - t318 - t74) * t226 + (t226 * t255 + t232 * t309) * t375 + Ifges(3,1) * t227) * t227; m(8) * (t15 * t4 + t16 * t5 + t17) + m(5) * (t168 + (t220 - 0.1e1) * t250) + ((-t156 * t227 + t255 + t376) * t232 + (t227 * t380 - t309) * t226) * pkin(2) + Ifges(3,6) * t233 + Ifges(3,5) * t227 + t261 * t147 + t38 * t351 + t129 * t354 + t14 * t360 + t106 * t364 + t3 * t368 + t2 * t369 + t56 * t370 + t13 * t371 + t11 * t372 + t267 * t105 + t269 * t146 + t274 * t55 + t275 * t102 - t310 * t58 + (-t74 / 0.2e1 - t125 / 0.2e1 - t318 / 0.2e1 + t227 * t189 / 0.2e1) * t232 + t15 * t8 + t16 * t7 + t5 * t26 + t4 * t27 + t332 * t89 + t333 * t40 + m(7) * (t41 * t70 + t17 + t334) + m(6) * (t61 * t92 + t168 + t334) + t337 * t69 + t41 * t63 + t70 * t28 + t92 * t66 + (-t313 / 0.2e1 - t130 / 0.2e1 + t231 * t363 + (-t120 / 0.2e1 - t184 / 0.2e1) * t227 - t266 * t225) * t226 + t61 * t110; t104 * t37 + t53 * t11 + 0.2e1 * t92 * t110 + t54 * t13 + 0.2e1 * t15 * t27 + t151 * t81 + 0.2e1 * t16 * t26 + 0.2e1 * t70 * t63 + Ifges(3,3) - 0.2e1 * t310 * t89 + 0.2e1 * t333 * t69 + (-t120 - t184) * t232 + (t77 - t31) * t146 + (t9 + t34) * t102 + (-t231 * t129 - t189 + (t124 - t73) * t225) * t226 + (t226 * t156 + (t215 + t221) * mrSges(4,3) + t380 * t232) * t375 + m(8) * (t15 ^ 2 + t16 ^ 2 + t68) + m(7) * (t70 ^ 2 + t373 + t68) + m(6) * (t92 ^ 2 + t198 + t373) + m(4) * (t210 + t299) + m(5) * (t220 * t299 + t198 + t210); m(8) * (t4 * t57 + t5 * t60) + t281 + t256 * t147 + t262 * t105 + t14 * t352 + t188 * t354 + t106 * t356 + t3 * t361 + t2 * t362 + t56 * t365 + t36 * t371 + t33 * t372 + t266 * t231 + t268 * t55 + t275 * t145 + t305 * t58 + t307 * t40 + (t278 * t226 + (t226 * t246 - t232 * t312) * pkin(2)) * t227 + t57 * t8 + t60 * t7 + t5 * t62 + t4 * t64 + (-t155 * t276 + t363 + t38 * t338 + t269 * t224 + (-t332 * t230 + t239 * t224 - m(8) * t273 + m(6) * (t224 * t61 - t315) + m(7) * (t293 * t41 - t273 - t315)) * pkin(3)) * t225 + t41 * t109 + t61 * t159; -Ifges(4,5) * t226 + m(8) * (t15 * t57 + t16 * t60) + t33 * t369 + t36 * t368 + t57 * t27 + t60 * t26 + t16 * t62 + t15 * t64 + t11 * t362 + t13 * t361 + t80 * t360 + t70 * t109 + t37 * t352 + t128 * t351 + t92 * t159 + t305 * t89 + t307 * t69 + t278 * t232 - t262 * t146 + t274 * t145 + t268 * t102 + (-t226 * t188 / 0.2e1 + t261) * t231 + (t226 * t312 + t232 * t246) * pkin(2) + (-t129 / 0.2e1 - t155 * t336 + t81 * t338 - t256 * t226 - t267 * t224 + (t310 * t230 + t238 * t224 - m(8) * t272 + m(6) * (t224 * t92 - t314) + m(7) * (t293 * t70 - t272 - t314)) * pkin(3)) * t225; t101 * t33 + t103 * t36 + t150 * t80 + 0.2e1 * t57 * t64 + 0.2e1 * t60 * t62 + Ifges(4,3) + (t119 - t183) * t231 + (t76 + t30) * t145 + m(7) * (t218 * t265 + t190 + t199) + m(6) * (t199 + t265) + m(8) * (t57 ^ 2 + t60 ^ 2 + t190) + (-t188 + (-0.2e1 * pkin(3) * t305 - t128) * t230 + (0.2e1 * pkin(3) * t378 + t123 - t72) * t224) * t225; t237 * t290 + m(8) * (t4 * t88 + t5 * t91) + t306 * t40 + t284 + t257 * t105 + t263 * t55 + (t61 * mrSges(6,3) + t269) * t230 + (t332 * t224 + t239 * t230 - m(8) * t271 + m(6) * (t230 * t61 + t317) + m(7) * (t287 * t41 - t271 + t317)) * pkin(3) + (t38 / 0.2e1 + t58 * mrSges(6,3) + t14 * t339 - t275 * t223) * t224 + t24 * t367 + t25 * t366 + t88 * t8 + t91 * t7 + t5 * t108 + t4 * t113 + t56 * t357 + t2 * t355 + t3 * t353 + t58 * t154 + t41 * t158 + t147 * t348 + t106 * t344; t193 + m(8) * (t15 * t88 + t16 * t91) + t53 * t367 + t54 * t366 + t88 * t27 + t91 * t26 + t16 * t108 + t15 * t113 + t104 * t357 + t11 * t355 + t13 * t353 + t89 * t154 + t70 * t158 + t151 * t344 + t306 * t69 + (-Ifges(5,5) * t231 + t178 * t342) * t226 - t257 * t146 + t263 * t102 + (t92 * mrSges(6,3) - t267) * t230 + t237 * t232 + (t89 * mrSges(6,3) - t223 * t274 + t339 * t37 + t364) * t224 + (-t310 * t224 + t238 * t230 - m(8) * t270 + m(6) * (t230 * t92 + t316) + m(7) * (t287 * t70 - t270 + t316)) * pkin(3); -t319 + t88 * t64 + t91 * t62 + t75 * t362 + t79 * t361 + t60 * t108 + t57 * t113 + t33 * t355 + t36 * t353 + t127 * t352 + (t348 - Ifges(5,6)) * t231 + t263 * t145 + m(8) * (t57 * t88 + t60 * t91 + t167) + m(7) * t167 + (t187 * t342 - t262) * t230 + (t356 + t80 * t339 - t268 * t223 + (m(7) * (t218 - 0.1e1) * t286 - t257) * t225) * t224 + ((-t154 * t225 + t378) * t230 + (t225 * t379 + t305) * t224) * pkin(3); 0.2e1 * t91 * t108 + 0.2e1 * t88 * t113 + t144 * t75 + t149 * t79 + Ifges(5,3) + (t182 - t118) * t230 + m(7) * (t218 * t300 + t197 + t209) + m(6) * (t209 + t300) + m(8) * (t88 ^ 2 + t91 ^ 2 + t197) + (t127 * t229 + t187 + (-t122 - t71) * t223) * t224 + 0.2e1 * (t154 * t224 + (t213 + t219) * mrSges(6,3) + t379 * t230) * pkin(3); -t61 * mrSges(6,2) + t24 * t359 + t25 * t358 + t40 * t153 + t5 * t157 + t4 * t162 + t105 * t177 / 0.2e1 + t56 * t345 - t311 * t58 + t258 * t55 + (t41 * mrSges(7,3) + t275) * t229 + (t14 / 0.2e1 - t40 * mrSges(7,3) + t2 * t343 + t3 * t340 + (-t222 * t7 - t228 * t8 + m(8) * (-t222 * t5 - t228 * t4)) * pkin(4)) * t223 + t32; -t92 * mrSges(6,2) + t53 * t359 + t54 * t358 + t69 * t153 + t16 * t157 + t15 * t162 + t146 * t349 + t104 * t345 - t311 * t89 + t258 * t102 + (t70 * mrSges(7,3) + t274) * t229 + (t370 - t69 * mrSges(7,3) + t11 * t343 + t13 * t340 + (-t222 * t26 - t228 * t27 + m(8) * (-t15 * t228 - t16 * t222)) * pkin(4)) * t223 + t73; t101 * t359 + t103 * t358 + t60 * t157 + t57 * t162 + t150 * t345 + t268 * t229 + t258 * t145 + (t224 * t349 - Ifges(6,5) * t230 + (t224 * t245 + t230 * t311) * pkin(3)) * t225 + (-pkin(3) * t153 * t294 + t365 + t33 * t343 + t36 * t340 + (-t222 * t62 - t228 * t64 + m(8) * (-t222 * t60 - t228 * t57)) * pkin(4)) * t223 + t283; t88 * t162 + t126 * t353 + t121 * t355 + t91 * t157 + t177 * t338 + (t224 * t345 + t263) * t229 + (-t224 * t311 + t230 * t245) * pkin(3) + (t357 - t153 * t335 + t75 * t343 + t79 * t340 - t258 * t224 + (m(8) * (-t222 * t91 - t228 * t88) - t228 * t113 - t222 * t108) * pkin(4)) * t223 + t178; Ifges(6,3) + (t117 + t181) * t229 + t212 * t249 + (-t222 * t121 + t228 * t126 + t186 + (-t157 * t222 - t162 * t228) * t279) * t223; t3 * t343 + t2 * t341 + t25 * t346 + t24 * t347 + t55 * t350 - t41 * mrSges(7,2) + t308 * t40 + t247 * mrSges(8,3) + (m(8) * t247 + t222 * t8 - t228 * t7) * pkin(4) + t10; t54 * t346 + t13 * t343 + t11 * t341 + t53 * t347 + t102 * t350 - t70 * mrSges(7,2) + t308 * t69 + t243 * mrSges(8,3) + (m(8) * t243 + t222 * t27 - t228 * t26) * pkin(4) + t31; t33 * t341 + t36 * t343 + t145 * t350 + t103 * t346 + t101 * t347 + t242 * mrSges(8,3) + (-Ifges(7,3) + (-t223 * t308 - t330) * pkin(3)) * t294 + (m(8) * t242 + t222 * t64 - t228 * t62) * pkin(4) + t285; t195 + t149 * t346 + t144 * t347 + t75 * t341 + t79 * t343 + (-pkin(3) * t330 - Ifges(7,3)) * t230 + t241 * mrSges(8,3) + ((-t176 / 0.2e1 - Ifges(7,6)) * t224 - t308 * t335) * t223 + (m(8) * t241 - t228 * t108 + t222 * t113) * pkin(4); t176 * t339 + (t223 * t346 - pkin(4) * t157 - t121 / 0.2e1) * t228 + (-t223 * t180 / 0.2e1 - t126 / 0.2e1 + pkin(4) * t162) * t222 + t177; mrSges(8,3) * t279 * t280 - t228 * t180 - t222 * t185 + Ifges(7,3) + t249; mrSges(8,1) * t4 - mrSges(8,2) * t5 + t1; mrSges(8,1) * t15 - mrSges(8,2) * t16 + t9; mrSges(8,1) * t57 - mrSges(8,2) * t60 + t30; mrSges(8,1) * t88 - mrSges(8,2) * t91 + t71; (-Ifges(8,6) * t222 - pkin(4) * t172) * t223 + t282; pkin(4) * t244 + t176; Ifges(8,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_7_matlab.m
res = [t6(1) t6(2) t6(4) t6(7) t6(11) t6(16) t6(22); t6(2) t6(3) t6(5) t6(8) t6(12) t6(17) t6(23); t6(4) t6(5) t6(6) t6(9) t6(13) t6(18) t6(24); t6(7) t6(8) t6(9) t6(10) t6(14) t6(19) t6(25); t6(11) t6(12) t6(13) t6(14) t6(15) t6(20) t6(26); t6(16) t6(17) t6(18) t6(19) t6(20) t6(21) t6(27); t6(22) t6(23) t6(24) t6(25) t6(26) t6(27) t6(28);];
Mq  = res;
