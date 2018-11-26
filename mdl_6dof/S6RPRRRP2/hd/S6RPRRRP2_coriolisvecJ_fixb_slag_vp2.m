% Calculate vector of centrifugal and coriolis load on the joints for
% S6RPRRRP2
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

function tauc = S6RPRRRP2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP2_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP2_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP2_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP2_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:24:32
% EndTime: 2018-11-23 16:24:44
% DurationCPUTime: 11.79s
% Computational Cost: add. (7570->553), mult. (18081->736), div. (0->0), fcn. (11728->8), ass. (0->239)
t208 = sin(pkin(10)) * pkin(1) + pkin(7);
t192 = t208 * qJD(1);
t217 = sin(qJ(3));
t220 = cos(qJ(3));
t153 = qJD(2) * t220 - t217 * t192;
t144 = -qJD(3) * pkin(3) - t153;
t216 = sin(qJ(4));
t219 = cos(qJ(4));
t267 = qJD(3) * t219;
t271 = qJD(1) * t217;
t179 = -t216 * t271 + t267;
t116 = -t179 * pkin(4) + t144;
t180 = qJD(3) * t216 + t219 * t271;
t215 = sin(qJ(5));
t218 = cos(qJ(5));
t124 = t179 * t215 + t180 * t218;
t270 = qJD(1) * t220;
t206 = qJD(4) - t270;
t154 = t217 * qJD(2) + t220 * t192;
t145 = qJD(3) * pkin(8) + t154;
t255 = -cos(pkin(10)) * pkin(1) - pkin(2);
t173 = -pkin(3) * t220 - t217 * pkin(8) + t255;
t148 = t173 * qJD(1);
t90 = -t145 * t216 + t219 * t148;
t76 = -pkin(9) * t180 + t90;
t69 = pkin(4) * t206 + t76;
t91 = t145 * t219 + t148 * t216;
t77 = pkin(9) * t179 + t91;
t73 = t215 * t77;
t22 = t218 * t69 - t73;
t351 = qJ(6) * t124;
t18 = t22 - t351;
t199 = qJD(5) + t206;
t15 = pkin(5) * t199 + t18;
t75 = t218 * t77;
t23 = t215 * t69 + t75;
t244 = t218 * t179 - t180 * t215;
t335 = qJ(6) * t244;
t19 = t23 + t335;
t268 = qJD(3) * t217;
t245 = qJD(1) * t268;
t203 = Ifges(7,3) * t245;
t204 = Ifges(6,3) * t245;
t309 = -t199 / 0.2e1;
t319 = t124 / 0.2e1;
t320 = -t124 / 0.2e1;
t261 = qJD(3) * qJD(4);
t265 = qJD(4) * t216;
t266 = qJD(3) * t220;
t138 = t219 * t261 + (-t217 * t265 + t219 * t266) * qJD(1);
t264 = qJD(4) * t219;
t349 = t216 * t266 + t217 * t264;
t139 = -qJD(1) * t349 - t216 * t261;
t51 = qJD(5) * t244 + t138 * t218 + t139 * t215;
t146 = t153 * qJD(3);
t242 = pkin(3) * t217 - pkin(8) * t220;
t190 = t242 * qJD(3);
t172 = qJD(1) * t190;
t45 = -qJD(4) * t91 - t146 * t216 + t219 * t172;
t30 = pkin(4) * t245 - pkin(9) * t138 + t45;
t44 = -t145 * t265 + t219 * t146 + t148 * t264 + t216 * t172;
t33 = pkin(9) * t139 + t44;
t6 = -qJD(5) * t23 - t215 * t33 + t218 * t30;
t2 = pkin(5) * t245 - qJ(6) * t51 - qJD(6) * t124 + t6;
t262 = qJD(5) * t218;
t263 = qJD(5) * t215;
t5 = t215 * t30 + t218 * t33 + t69 * t262 - t263 * t77;
t52 = -qJD(5) * t124 - t138 * t215 + t139 * t218;
t3 = qJ(6) * t52 + qJD(6) * t244 + t5;
t333 = t6 * mrSges(6,1) + t2 * mrSges(7,1) - t5 * mrSges(6,2) - t3 * mrSges(7,2);
t345 = -t244 / 0.2e1;
t363 = Ifges(6,5) + Ifges(7,5);
t365 = Ifges(6,1) + Ifges(7,1);
t364 = Ifges(6,4) + Ifges(7,4);
t368 = t364 * t244;
t355 = t365 * t124 + t363 * t199 + t368;
t361 = Ifges(6,6) + Ifges(7,6);
t362 = Ifges(6,2) + Ifges(7,2);
t367 = t124 * t364;
t356 = t361 * t199 + t362 * t244 + t367;
t47 = Ifges(7,6) * t52;
t48 = Ifges(6,6) * t52;
t49 = Ifges(7,5) * t51;
t50 = Ifges(6,5) * t51;
t70 = -pkin(5) * t244 + qJD(6) + t116;
t376 = t203 + t204 + t47 + t48 + t49 + t50 + t333 + (-t124 * t361 + t363 * t244) * t309 + (t124 * t19 + t15 * t244) * mrSges(7,3) + (t124 * t23 + t22 * t244) * mrSges(6,3) - t116 * (mrSges(6,1) * t124 + mrSges(6,2) * t244) - t70 * (mrSges(7,1) * t124 + mrSges(7,2) * t244) + t356 * t319 + (-t124 * t362 + t355 + t368) * t345 + (t365 * t244 - t367) * t320;
t187 = t242 * qJD(1);
t113 = t219 * t153 + t216 * t187;
t253 = t216 * t270;
t325 = -pkin(9) - pkin(8);
t254 = qJD(4) * t325;
t375 = pkin(9) * t253 + t216 * t254 - t113;
t112 = -t216 * t153 + t219 * t187;
t277 = t219 * t220;
t228 = pkin(4) * t217 - pkin(9) * t277;
t374 = -qJD(1) * t228 + t219 * t254 - t112;
t229 = t215 * t216 - t218 * t219;
t334 = qJD(4) + qJD(5);
t127 = t334 * t229;
t226 = t229 * t220;
t150 = qJD(1) * t226;
t373 = t127 - t150;
t182 = t215 * t219 + t216 * t218;
t128 = t334 * t182;
t227 = t182 * t220;
t149 = qJD(1) * t227;
t372 = t128 - t149;
t371 = -Ifges(4,1) / 0.2e1;
t211 = Ifges(4,4) * t270;
t370 = -t211 / 0.2e1;
t369 = -t362 * t52 / 0.2e1 - t364 * t51 / 0.2e1 - t361 * t245 / 0.2e1;
t197 = t325 * t216;
t198 = t325 * t219;
t354 = t197 * t262 + t198 * t263 + t215 * t374 + t218 * t375;
t135 = t215 * t197 - t218 * t198;
t353 = -qJD(5) * t135 - t215 * t375 + t218 * t374;
t250 = Ifges(4,5) * qJD(3) / 0.2e1;
t359 = t363 * t245 + t364 * t52 + t365 * t51;
t358 = -t372 * qJ(6) - qJD(6) * t229 + t354;
t357 = -pkin(5) * t271 + t373 * qJ(6) - qJD(6) * t182 + t353;
t131 = pkin(4) * t253 + t154;
t352 = pkin(4) * t265 + t372 * pkin(5) - t131;
t296 = Ifges(5,4) * t180;
t109 = Ifges(5,2) * t179 + Ifges(5,6) * t206 + t296;
t174 = Ifges(5,4) * t179;
t110 = t180 * Ifges(5,1) + t206 * Ifges(5,5) + t174;
t233 = t216 * t91 + t219 * t90;
t294 = Ifges(5,4) * t219;
t237 = -Ifges(5,2) * t216 + t294;
t295 = Ifges(5,4) * t216;
t239 = Ifges(5,1) * t219 - t295;
t240 = mrSges(5,1) * t216 + mrSges(5,2) * t219;
t290 = Ifges(5,6) * t216;
t291 = Ifges(5,5) * t219;
t305 = t219 / 0.2e1;
t306 = -t216 / 0.2e1;
t312 = t180 / 0.2e1;
t221 = -t233 * mrSges(5,3) + t179 * t237 / 0.2e1 + t239 * t312 + t144 * t240 + t206 * (-t290 + t291) / 0.2e1 + t109 * t306 + t110 * t305;
t350 = t153 * mrSges(4,3) + t271 * t371 - t221 - t250 + t370;
t308 = t199 / 0.2e1;
t322 = t244 / 0.2e1;
t249 = -Ifges(4,6) * qJD(3) / 0.2e1;
t72 = -mrSges(6,1) * t244 + mrSges(6,2) * t124;
t336 = m(6) * t116 + t72;
t158 = t229 * t217;
t160 = t219 * t173;
t279 = t208 * t216;
t303 = pkin(9) * t217;
t104 = -t219 * t303 + t160 + (-pkin(4) - t279) * t220;
t186 = t208 * t277;
t126 = t216 * t173 + t186;
t278 = t216 * t217;
t114 = -pkin(9) * t278 + t126;
t57 = t215 * t104 + t218 * t114;
t234 = -t216 * t45 + t219 * t44;
t331 = -t45 * mrSges(5,1) + t44 * mrSges(5,2) - Ifges(5,5) * t138 - Ifges(5,6) * t139;
t194 = t255 * qJD(1);
t258 = Ifges(7,3) / 0.2e1 + Ifges(6,3) / 0.2e1;
t259 = Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1;
t260 = Ifges(7,5) / 0.2e1 + Ifges(6,5) / 0.2e1;
t297 = Ifges(4,4) * t217;
t330 = t260 * t124 + t258 * t199 + t259 * t244 + t15 * mrSges(7,1) + t194 * mrSges(4,1) + t22 * mrSges(6,1) + t90 * mrSges(5,1) + t206 * Ifges(5,3) + t180 * Ifges(5,5) + t179 * Ifges(5,6) + t249 - (t220 * Ifges(4,2) + t297) * qJD(1) / 0.2e1 - t19 * mrSges(7,2) - t23 * mrSges(6,2) - t91 * mrSges(5,2) + t361 * t322 + t363 * t319 + (Ifges(7,3) + Ifges(6,3)) * t308;
t329 = t51 / 0.2e1;
t328 = t52 / 0.2e1;
t318 = t138 / 0.2e1;
t317 = t139 / 0.2e1;
t314 = -t179 / 0.2e1;
t313 = -t180 / 0.2e1;
t307 = -t206 / 0.2e1;
t42 = -mrSges(7,2) * t245 + mrSges(7,3) * t52;
t43 = -mrSges(6,2) * t245 + mrSges(6,3) * t52;
t299 = t42 + t43;
t26 = t218 * t76 - t73;
t98 = -mrSges(7,2) * t199 + mrSges(7,3) * t244;
t99 = -mrSges(6,2) * t199 + mrSges(6,3) * t244;
t298 = t98 + t99;
t286 = t194 * mrSges(4,2);
t100 = mrSges(7,1) * t199 - mrSges(7,3) * t124;
t101 = mrSges(6,1) * t199 - mrSges(6,3) * t124;
t276 = t100 + t101;
t273 = t219 * t190 + t268 * t279;
t257 = mrSges(4,3) * t271;
t272 = qJD(3) * mrSges(4,1) + mrSges(5,1) * t179 - mrSges(5,2) * t180 - t257;
t161 = pkin(4) * t278 + t217 * t208;
t256 = mrSges(4,3) * t270;
t130 = pkin(4) * t349 + t208 * t266;
t210 = -pkin(4) * t219 - pkin(3);
t246 = m(4) * t208 + mrSges(4,3);
t16 = -t52 * mrSges(7,1) + t51 * mrSges(7,2);
t25 = -t215 * t76 - t75;
t56 = t218 * t104 - t215 * t114;
t134 = t218 * t197 + t198 * t215;
t241 = mrSges(5,1) * t219 - mrSges(5,2) * t216;
t238 = Ifges(5,1) * t216 + t294;
t236 = Ifges(5,2) * t219 + t295;
t235 = Ifges(5,5) * t216 + Ifges(5,6) * t219;
t232 = t216 * t90 - t219 * t91;
t119 = mrSges(5,1) * t245 - mrSges(5,3) * t138;
t120 = -mrSges(5,2) * t245 + mrSges(5,3) * t139;
t231 = -t216 * t119 + t219 * t120;
t141 = -mrSges(5,2) * t206 + mrSges(5,3) * t179;
t142 = mrSges(5,1) * t206 - mrSges(5,3) * t180;
t230 = -t216 * t141 - t219 * t142;
t58 = t228 * qJD(3) + (-t186 + (-t173 + t303) * t216) * qJD(4) + t273;
t78 = t173 * t264 + t216 * t190 + (-t217 * t267 - t220 * t265) * t208;
t66 = -pkin(9) * t349 + t78;
t9 = t104 * t262 - t114 * t263 + t215 * t58 + t218 * t66;
t147 = t154 * qJD(3);
t95 = -t139 * pkin(4) + t147;
t10 = -qJD(5) * t57 - t215 * t66 + t218 * t58;
t209 = pkin(4) * t218 + pkin(5);
t205 = Ifges(5,3) * t245;
t195 = -qJD(3) * mrSges(4,2) + t256;
t157 = t182 * t217;
t151 = pkin(5) * t229 + t210;
t125 = -t220 * t279 + t160;
t115 = pkin(5) * t157 + t161;
t107 = -qJ(6) * t229 + t135;
t106 = -qJ(6) * t182 + t134;
t92 = pkin(4) * t180 + pkin(5) * t124;
t88 = -mrSges(5,1) * t139 + mrSges(5,2) * t138;
t83 = t138 * Ifges(5,1) + t139 * Ifges(5,4) + Ifges(5,5) * t245;
t82 = t138 * Ifges(5,4) + t139 * Ifges(5,2) + Ifges(5,6) * t245;
t81 = -qJD(3) * t227 + t158 * t334;
t80 = -qJD(3) * t226 - t128 * t217;
t79 = -qJD(4) * t126 + t273;
t71 = -mrSges(7,1) * t244 + mrSges(7,2) * t124;
t53 = -pkin(5) * t81 + t130;
t41 = mrSges(6,1) * t245 - mrSges(6,3) * t51;
t40 = mrSges(7,1) * t245 - mrSges(7,3) * t51;
t35 = -qJ(6) * t157 + t57;
t34 = -pkin(5) * t220 + t158 * qJ(6) + t56;
t24 = -t52 * pkin(5) + t95;
t21 = t26 - t351;
t20 = t25 - t335;
t17 = -mrSges(6,1) * t52 + mrSges(6,2) * t51;
t8 = qJ(6) * t81 - qJD(6) * t157 + t9;
t7 = pkin(5) * t268 - qJ(6) * t80 + qJD(6) * t158 + t10;
t1 = [t356 * t81 / 0.2e1 - t359 * t158 / 0.2e1 + (t361 * t81 + t363 * t80) * t308 + (-t203 / 0.2e1 - t204 / 0.2e1 - t205 / 0.2e1 - t50 / 0.2e1 - t49 / 0.2e1 - t48 / 0.2e1 - t47 / 0.2e1 - t259 * t52 - t260 * t51 + t246 * t146 + (0.3e1 / 0.2e1 * t211 + t250 + 0.2e1 * t286 + (-m(4) * t153 + m(5) * t144 - t272) * t208 - t350) * qJD(3) + t331 - t333) * t220 + t355 * t80 / 0.2e1 + m(5) * (t45 * t125 + t44 * t126 + t91 * t78 + t90 * t79) + m(7) * (t115 * t24 + t15 * t7 + t19 * t8 + t2 * t34 + t3 * t35 + t53 * t70) + m(6) * (t10 * t22 + t116 * t130 + t161 * t95 + t23 * t9 + t5 * t57 + t56 * t6) + (t83 * t305 + t82 * t306 + t239 * t318 + t237 * t317 + (-t216 * t44 - t219 * t45) * mrSges(5,3) + (mrSges(4,3) + t240) * t147 + (t110 * t306 - t219 * t109 / 0.2e1 + t235 * t307 + t236 * t314 + t238 * t313 + t144 * t241 + t232 * mrSges(5,3)) * qJD(4) + (-t246 * t154 + ((-0.3e1 / 0.2e1 * Ifges(4,4) + t291 / 0.2e1 - t290 / 0.2e1) * t217 + t255 * mrSges(4,1) - t260 * t158 - t259 * t157 + (0.3e1 / 0.2e1 * Ifges(4,1) - Ifges(5,3) / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(4,2) - t258) * t220) * qJD(1) + t249 + t330) * qJD(3) + (t88 + (m(5) + m(4)) * t147 - t195 * qJD(3)) * t208) * t217 + (t362 * t81 + t364 * t80) * t322 + (-t157 * t362 - t158 * t364) * t328 + (t364 * t81 + t365 * t80) * t319 + (-t157 * t364 - t158 * t365) * t329 + t161 * t17 + t79 * t142 + t78 * t141 + t130 * t72 + t125 * t119 + t126 * t120 + t115 * t16 + t116 * (-mrSges(6,1) * t81 + mrSges(6,2) * t80) + t10 * t101 + t8 * t98 + t9 * t99 + t7 * t100 + t70 * (-mrSges(7,1) * t81 + mrSges(7,2) * t80) + t53 * t71 + t56 * t41 + t57 * t43 + t35 * t42 + t34 * t40 + t157 * t369 + (-t15 * t80 - t157 * t3 + t158 * t2 + t19 * t81) * mrSges(7,3) + (-t157 * t5 + t158 * t6 - t22 * t80 + t23 * t81) * mrSges(6,3) + t95 * (mrSges(6,1) * t157 - mrSges(6,2) * t158) + t24 * (mrSges(7,1) * t157 - mrSges(7,2) * t158); t276 * t81 + t298 * t80 - t299 * t158 - (t40 + t41) * t157 + (-t16 - t17 - t88) * t220 + (t230 * qJD(4) + t231) * t217 + ((t141 * t219 - t142 * t216 + t195 - t256) * t220 + (t71 + t72 - t257 - t272) * t217) * qJD(3) + m(4) * (t146 * t217 - t147 * t220 + (-t153 * t217 + t154 * t220) * qJD(3)) + m(7) * (t15 * t81 - t2 * t157 - t3 * t158 + t19 * t80 - t220 * t24 + t268 * t70) + m(6) * (t116 * t268 - t6 * t157 - t5 * t158 + t22 * t81 - t220 * t95 + t23 * t80) + m(5) * ((-qJD(3) * t232 - t147) * t220 + (qJD(3) * t144 - qJD(4) * t233 + t234) * t217); t352 * t71 + t234 * mrSges(5,3) + t356 * (t149 / 0.2e1 - t128 / 0.2e1) + t357 * t100 + (t106 * t2 + t107 * t3 + t15 * t357 + t151 * t24 + t19 * t358 + t352 * t70) * m(7) + t358 * t98 + t359 * t182 / 0.2e1 + (-t127 * t363 - t128 * t361) * t308 + (-t149 * t361 - t150 * t363) * t309 + t353 * t101 + (-t116 * t131 + t134 * t6 + t135 * t5 + t210 * t95 + t22 * t353 + t23 * t354) * m(6) + t354 * t99 + t355 * (t150 / 0.2e1 - t127 / 0.2e1) + ((t250 + t370 - t286 + t350) * t220 + ((t297 / 0.2e1 + (t371 + Ifges(4,2) / 0.2e1) * t220) * qJD(1) + t154 * mrSges(4,3) + t249 - t330) * t217 + (t182 * t363 - t229 * t361 + t235) * t268 / 0.2e1) * qJD(1) + (-t182 * t6 + t22 * t373 - t229 * t5 - t23 * t372) * mrSges(6,3) + (mrSges(7,1) * t372 - mrSges(7,2) * t373) * t70 + (mrSges(6,1) * t372 - mrSges(6,2) * t373) * t116 + (t15 * t373 - t182 * t2 - t19 * t372 - t229 * t3) * mrSges(7,3) + t24 * (mrSges(7,1) * t229 + mrSges(7,2) * t182) + t95 * (mrSges(6,1) * t229 + mrSges(6,2) * t182) + t216 * t83 / 0.2e1 + t231 * pkin(8) + (-t127 * t364 - t128 * t362) * t322 + (-t149 * t362 - t150 * t364) * t345 + (t182 * t364 - t229 * t362) * t328 + (-t127 * t365 - t128 * t364) * t319 + (-t149 * t364 - t150 * t365) * t320 + (t182 * t365 - t229 * t364) * t329 + t210 * t17 - m(5) * (t112 * t90 + t113 * t91 + t144 * t154) - t153 * t195 + t151 * t16 - t112 * t142 - t146 * mrSges(4,2) - t113 * t141 - t131 * t72 + t134 * t41 + t135 * t43 + t106 * t40 + t107 * t42 - pkin(3) * t88 + t229 * t369 + t236 * t317 + t238 * t318 + m(5) * (-pkin(3) * t147 + pkin(8) * t234) + (t336 * t216 * pkin(4) + (-m(5) * t233 + t230) * pkin(8) + t221) * qJD(4) + t272 * t154 + (-mrSges(4,1) - t241) * t147 + t82 * t305; t205 - m(6) * (t22 * t25 + t23 * t26) + (t218 * t41 + t299 * t215 + (-t215 * t276 + t218 * t298) * qJD(5) + m(6) * (t215 * t5 + t218 * t6 - t22 * t263 + t23 * t262) - t336 * t180) * pkin(4) + ((-t15 * t263 + t19 * t262 + t215 * t3) * pkin(4) + t2 * t209 - t15 * t20 - t19 * t21 - t70 * t92) * m(7) - t331 + (t179 * t90 + t180 * t91) * mrSges(5,3) + (-Ifges(5,2) * t180 + t110 + t174) * t314 + t209 * t40 - t144 * (mrSges(5,1) * t180 + mrSges(5,2) * t179) + t91 * t142 - t90 * t141 - t21 * t98 - t26 * t99 - t20 * t100 - t25 * t101 - t92 * t71 + (Ifges(5,5) * t179 - Ifges(5,6) * t180) * t307 + t109 * t312 + (Ifges(5,1) * t179 - t296) * t313 + t376; (-t124 * t71 + t40) * pkin(5) + (-(-t15 + t18) * t19 + (-t124 * t70 + t2) * pkin(5)) * m(7) - t18 * t98 - t22 * t99 + t19 * t100 + t23 * t101 + t376; t124 * t100 - t244 * t98 + 0.2e1 * (t24 / 0.2e1 + t19 * t345 + t15 * t319) * m(7) + t16;];
tauc  = t1(:);
