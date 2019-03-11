% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
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
% Datum: 2019-03-09 04:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRPP4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP4_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP4_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP4_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP4_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:38:47
% EndTime: 2019-03-09 04:39:08
% DurationCPUTime: 10.50s
% Computational Cost: add. (9052->540), mult. (24180->700), div. (0->0), fcn. (18089->8), ass. (0->240)
t358 = Ifges(7,2) + Ifges(6,3);
t363 = Ifges(5,3) + t358;
t362 = Ifges(4,2) / 0.2e1;
t214 = sin(pkin(9));
t215 = cos(pkin(9));
t217 = sin(qJ(3));
t219 = cos(qJ(3));
t198 = t214 * t219 + t215 * t217;
t191 = t198 * qJD(1);
t216 = sin(qJ(4));
t218 = cos(qJ(4));
t164 = qJD(3) * t218 - t191 * t216;
t196 = t214 * t217 - t219 * t215;
t192 = t196 * qJD(3);
t178 = qJD(1) * t192;
t119 = qJD(4) * t164 - t178 * t218;
t165 = qJD(3) * t216 + t191 * t218;
t120 = -qJD(4) * t165 + t178 * t216;
t213 = sin(pkin(10));
t279 = cos(pkin(10));
t76 = t119 * t213 - t120 * t279;
t329 = -t76 / 0.2e1;
t77 = t119 * t279 + t213 * t120;
t327 = t77 / 0.2e1;
t193 = t198 * qJD(3);
t179 = qJD(1) * t193;
t310 = t179 / 0.2e1;
t353 = Ifges(6,1) + Ifges(7,1);
t352 = Ifges(6,4) - Ifges(7,5);
t351 = Ifges(7,4) + Ifges(6,5);
t350 = Ifges(6,6) - Ifges(7,6);
t189 = t196 * qJD(1);
t248 = t279 * t216;
t197 = t213 * t218 + t248;
t128 = t197 * t189;
t188 = t197 * qJD(4);
t361 = t128 + t188;
t247 = t279 * t218;
t269 = t213 * t216;
t227 = t247 - t269;
t129 = t227 * t189;
t190 = t227 * qJD(4);
t265 = -t129 - t190;
t298 = pkin(7) + qJ(2);
t203 = t298 * t214;
t200 = qJD(1) * t203;
t204 = t298 * t215;
t201 = qJD(1) * t204;
t154 = -t200 * t219 - t217 * t201;
t148 = -qJD(3) * pkin(3) - t154;
t254 = -pkin(2) * t215 - pkin(1);
t202 = qJD(1) * t254 + qJD(2);
t125 = pkin(3) * t189 - pkin(8) * t191 + t202;
t155 = -t217 * t200 + t219 * t201;
t149 = qJD(3) * pkin(8) + t155;
t85 = t218 * t125 - t149 * t216;
t86 = t125 * t216 + t149 * t218;
t233 = t86 * t216 + t85 * t218;
t236 = Ifges(5,5) * t218 - Ifges(5,6) * t216;
t289 = Ifges(5,4) * t218;
t238 = -Ifges(5,2) * t216 + t289;
t290 = Ifges(5,4) * t216;
t240 = Ifges(5,1) * t218 - t290;
t241 = mrSges(5,1) * t216 + mrSges(5,2) * t218;
t304 = t218 / 0.2e1;
t305 = -t216 / 0.2e1;
t182 = qJD(4) + t189;
t308 = t182 / 0.2e1;
t311 = t165 / 0.2e1;
t287 = t165 * Ifges(5,4);
t98 = t164 * Ifges(5,2) + t182 * Ifges(5,6) + t287;
t163 = Ifges(5,4) * t164;
t99 = t165 * Ifges(5,1) + t182 * Ifges(5,5) + t163;
t360 = t98 * t305 + t99 * t304 + t148 * t241 + t164 * t238 / 0.2e1 + t240 * t311 + t236 * t308 - t233 * mrSges(5,3);
t359 = t189 * t362;
t181 = Ifges(4,4) * t189;
t342 = -t181 / 0.2e1 + t191 * Ifges(4,1) / 0.2e1;
t357 = t202 * mrSges(4,2) + Ifges(4,5) * qJD(3) + t342 + t360;
t104 = -pkin(4) * t164 + qJD(5) + t148;
t67 = qJ(5) * t164 + t86;
t282 = t213 * t67;
t66 = -qJ(5) * t165 + t85;
t53 = pkin(4) * t182 + t66;
t17 = t279 * t53 - t282;
t13 = -t182 * pkin(5) + qJD(6) - t17;
t106 = -t279 * t164 + t165 * t213;
t228 = t213 * t164 + t165 * t279;
t348 = -t352 * t106 + t351 * t182 + t353 * t228;
t41 = pkin(5) * t106 - qJ(6) * t228 + t104;
t356 = t104 * mrSges(6,2) + t13 * mrSges(7,2) - t17 * mrSges(6,3) - t41 * mrSges(7,3) + t348 / 0.2e1;
t225 = t196 * qJD(2);
t111 = -qJD(1) * t225 + qJD(3) * t154;
t135 = pkin(3) * t179 + pkin(8) * t178;
t40 = -qJD(4) * t86 - t111 * t216 + t218 * t135;
t12 = pkin(4) * t179 - qJ(5) * t119 - qJD(5) * t165 + t40;
t261 = qJD(4) * t218;
t262 = qJD(4) * t216;
t39 = t218 * t111 + t125 * t261 + t216 * t135 - t149 * t262;
t16 = qJ(5) * t120 + qJD(5) * t164 + t39;
t3 = t12 * t279 - t213 * t16;
t2 = -t179 * pkin(5) - t3;
t328 = t76 / 0.2e1;
t226 = t198 * qJD(2);
t112 = qJD(1) * t226 + qJD(3) * t155;
t79 = -pkin(4) * t120 + t112;
t9 = pkin(5) * t76 - qJ(6) * t77 - qJD(6) * t228 + t79;
t355 = mrSges(6,2) * t79 + mrSges(7,2) * t2 - mrSges(6,3) * t3 - mrSges(7,3) * t9 + Ifges(7,5) * t328 + 0.2e1 * t351 * t310 + 0.2e1 * t327 * t353 + (Ifges(6,4) + t352) * t329;
t354 = -Ifges(6,6) / 0.2e1;
t321 = -t106 / 0.2e1;
t320 = t106 / 0.2e1;
t309 = -t182 / 0.2e1;
t318 = -t228 / 0.2e1;
t91 = -mrSges(7,2) * t106 + mrSges(7,3) * t182;
t92 = -mrSges(6,2) * t182 - mrSges(6,3) * t106;
t294 = t91 + t92;
t93 = mrSges(6,1) * t182 - mrSges(6,3) * t228;
t94 = -mrSges(7,1) * t182 + mrSges(7,2) * t228;
t293 = t94 - t93;
t274 = t189 * t216;
t114 = -pkin(4) * t274 + t155;
t347 = pkin(4) * t262 + t361 * pkin(5) + t265 * qJ(6) - qJD(6) * t197 - t114;
t346 = Ifges(5,5) * t119 + Ifges(5,6) * t120;
t153 = pkin(3) * t196 - pkin(8) * t198 + t254;
t160 = -t203 * t217 + t204 * t219;
t156 = t218 * t160;
t103 = t216 * t153 + t156;
t345 = -t219 * t203 - t204 * t217;
t344 = -t216 * t40 + t218 * t39;
t343 = t363 * t179 - t350 * t76 + t351 * t77 + t346;
t341 = (m(3) * qJ(2) + mrSges(3,3)) * (t214 ^ 2 + t215 ^ 2);
t4 = t213 * t12 + t279 * t16;
t1 = qJ(6) * t179 + qJD(6) * t182 + t4;
t340 = t40 * mrSges(5,1) + t3 * mrSges(6,1) - t2 * mrSges(7,1) - t39 * mrSges(5,2) - t4 * mrSges(6,2) + t1 * mrSges(7,3);
t61 = t279 * t67;
t18 = t213 * t53 + t61;
t14 = qJ(6) * t182 + t18;
t334 = Ifges(5,3) / 0.2e1;
t339 = (Ifges(7,6) / 0.2e1 + t354) * t106 + (t334 + Ifges(6,3) / 0.2e1 + Ifges(7,2) / 0.2e1) * t182 + (Ifges(6,5) / 0.2e1 + Ifges(7,4) / 0.2e1) * t228 - t13 * mrSges(7,1) - t18 * mrSges(6,2) - t86 * mrSges(5,2) - Ifges(4,6) * qJD(3) - t191 * Ifges(4,4) + t359 - Ifges(6,6) * t320 - Ifges(7,6) * t321 + t165 * Ifges(5,5) + t164 * Ifges(5,6) + t14 * mrSges(7,3) + t17 * mrSges(6,1) + t202 * mrSges(4,1) + t85 * mrSges(5,1) - t351 * t318 - t363 * t309;
t47 = Ifges(7,5) * t228 + t182 * Ifges(7,6) + t106 * Ifges(7,3);
t50 = Ifges(6,4) * t228 - t106 * Ifges(6,2) + t182 * Ifges(6,6);
t337 = -t104 * mrSges(6,1) - t41 * mrSges(7,1) + t14 * mrSges(7,2) + t18 * mrSges(6,3) + t50 / 0.2e1 - t47 / 0.2e1;
t335 = mrSges(6,1) * t79 + mrSges(7,1) * t9 - mrSges(7,2) * t1 - mrSges(6,3) * t4 + 0.2e1 * Ifges(7,3) * t328 - t77 * Ifges(6,4) / 0.2e1 + t179 * t354 + Ifges(7,6) * t310 + (-t352 + Ifges(7,5)) * t327 + (-t329 + t328) * Ifges(6,2);
t317 = t228 / 0.2e1;
t316 = t119 / 0.2e1;
t315 = t120 / 0.2e1;
t313 = -t164 / 0.2e1;
t312 = -t165 / 0.2e1;
t303 = pkin(4) * t165;
t302 = pkin(4) * t213;
t301 = pkin(4) * t216;
t297 = -qJ(5) - pkin(8);
t230 = qJ(5) * t192 - qJD(5) * t198;
t126 = qJD(3) * t345 - t225;
t151 = pkin(3) * t193 + pkin(8) * t192;
t245 = -t126 * t216 + t218 * t151;
t28 = pkin(4) * t193 + t230 * t218 + (-t156 + (qJ(5) * t198 - t153) * t216) * qJD(4) + t245;
t252 = t198 * t261;
t255 = t218 * t126 + t216 * t151 + t153 * t261;
t32 = -qJ(5) * t252 + (-qJD(4) * t160 + t230) * t216 + t255;
t8 = t213 * t28 + t279 * t32;
t54 = -mrSges(7,2) * t76 + mrSges(7,3) * t179;
t55 = -mrSges(6,2) * t179 - mrSges(6,3) * t76;
t296 = t54 + t55;
t56 = mrSges(6,1) * t179 - mrSges(6,3) * t77;
t57 = -t179 * mrSges(7,1) + t77 * mrSges(7,2);
t295 = t57 - t56;
t276 = qJ(5) * t218;
t150 = pkin(3) * t191 + pkin(8) * t189;
t95 = t218 * t150 - t154 * t216;
t68 = pkin(4) * t191 + t189 * t276 + t95;
t96 = t216 * t150 + t218 * t154;
t83 = qJ(5) * t274 + t96;
t31 = t213 * t68 + t279 * t83;
t102 = t218 * t153 - t160 * t216;
t81 = pkin(4) * t196 - t198 * t276 + t102;
t272 = t198 * t216;
t87 = -qJ(5) * t272 + t103;
t38 = t213 * t81 + t279 * t87;
t292 = mrSges(5,3) * t164;
t291 = mrSges(5,3) * t165;
t275 = t112 * t345;
t267 = -qJD(3) * mrSges(4,1) - mrSges(5,1) * t164 + mrSges(5,2) * t165 + mrSges(4,3) * t191;
t263 = qJD(4) * t198;
t210 = -pkin(4) * t218 - pkin(3);
t253 = t279 * pkin(4);
t34 = t76 * mrSges(6,1) + t77 * mrSges(6,2);
t33 = t76 * mrSges(7,1) - t77 * mrSges(7,3);
t249 = qJD(4) * t297;
t246 = t179 * mrSges(4,1) - t178 * mrSges(4,2);
t130 = pkin(4) * t272 - t345;
t242 = mrSges(5,1) * t218 - mrSges(5,2) * t216;
t239 = Ifges(5,1) * t216 + t289;
t237 = Ifges(5,2) * t218 + t290;
t235 = Ifges(5,5) * t216 + Ifges(5,6) * t218;
t234 = -t216 * t39 - t218 * t40;
t232 = t216 * t85 - t218 * t86;
t123 = -mrSges(5,2) * t182 + t292;
t124 = mrSges(5,1) * t182 - t291;
t231 = t123 * t218 - t124 * t216;
t7 = -t213 * t32 + t279 * t28;
t30 = -t213 * t83 + t279 * t68;
t37 = -t213 * t87 + t279 * t81;
t224 = -qJD(5) * t216 + t218 * t249;
t127 = qJD(3) * t160 + t226;
t90 = pkin(4) * t252 - t192 * t301 + t127;
t209 = -t253 - pkin(5);
t206 = qJ(6) + t302;
t205 = t297 * t218;
t184 = qJD(5) * t218 + t216 * t249;
t168 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t189;
t162 = -t205 * t279 + t269 * t297;
t161 = -t205 * t213 - t248 * t297;
t152 = -pkin(5) * t227 - qJ(6) * t197 + t210;
t138 = t227 * t198;
t137 = t197 * t198;
t134 = t184 * t279 + t213 * t224;
t133 = t184 * t213 - t224 * t279;
t101 = -mrSges(5,2) * t179 + mrSges(5,3) * t120;
t100 = mrSges(5,1) * t179 - mrSges(5,3) * t119;
t89 = -t192 * t227 - t197 * t263;
t88 = t198 * t213 * t262 + t192 * t197 - t247 * t263;
t80 = -mrSges(5,1) * t120 + mrSges(5,2) * t119;
t64 = mrSges(6,1) * t106 + mrSges(6,2) * t228;
t63 = mrSges(7,1) * t106 - mrSges(7,3) * t228;
t60 = t119 * Ifges(5,1) + t120 * Ifges(5,4) + t179 * Ifges(5,5);
t59 = t119 * Ifges(5,4) + t120 * Ifges(5,2) + t179 * Ifges(5,6);
t58 = pkin(5) * t137 - qJ(6) * t138 + t130;
t45 = pkin(5) * t228 + qJ(6) * t106 + t303;
t43 = -qJD(4) * t103 + t245;
t42 = -t160 * t262 + t255;
t36 = -t196 * pkin(5) - t37;
t35 = qJ(6) * t196 + t38;
t27 = -t191 * pkin(5) - t30;
t26 = qJ(6) * t191 + t31;
t20 = t279 * t66 - t282;
t19 = t213 * t66 + t61;
t10 = -pkin(5) * t88 - qJ(6) * t89 - qJD(6) * t138 + t90;
t6 = -t193 * pkin(5) - t7;
t5 = qJ(6) * t193 + qJD(6) * t196 + t8;
t11 = [(-t310 * t350 + t335) * t137 - (t342 + t357) * t192 + (Ifges(6,4) * t321 + Ifges(7,5) * t320 + t351 * t308 + t353 * t317 + t356) * t89 + 0.2e1 * t341 * qJD(2) * qJD(1) + t355 * t138 + (t59 * t305 + t60 * t304 - Ifges(4,1) * t178 - Ifges(4,4) * t179 + t240 * t316 + t238 * t315 + t236 * t310 + (mrSges(4,3) + t241) * t112 + t234 * mrSges(5,3) + (t99 * t305 - t218 * t98 / 0.2e1 + t148 * t242 + t237 * t313 + t239 * t312 + t235 * t309 + t232 * mrSges(5,3)) * qJD(4)) * t198 + m(5) * (t102 * t40 + t103 * t39 + t127 * t148 + t42 * t86 + t43 * t85 - t275) + m(4) * (t111 * t160 + t126 * t155 - t127 * t154 - t275) + (-mrSges(4,3) * t111 + Ifges(4,4) * t178 + Ifges(6,6) * t329 + Ifges(7,6) * t328 + t351 * t327 + t358 * t310 + t340 + t343 / 0.2e1 + t346 / 0.2e1 + (t334 + Ifges(4,2)) * t179) * t196 + (t154 * t192 - t155 * t193 - t160 * t179) * mrSges(4,3) + m(7) * (t1 * t35 + t10 * t41 + t13 * t6 + t14 * t5 + t2 * t36 + t58 * t9) + m(6) * (t104 * t90 + t130 * t79 + t17 * t7 + t18 * t8 + t3 * t37 + t38 * t4) + (t359 + t339) * t193 - (-mrSges(4,3) * t178 + t80) * t345 + (Ifges(6,2) * t321 - Ifges(7,3) * t320 + t350 * t308 + t352 * t317 + t337) * t88 + t254 * t246 + t267 * t127 + t35 * t54 + t38 * t55 + t37 * t56 + t36 * t57 + t58 * t33 + t10 * t63 + t90 * t64 + t5 * t91 + t8 * t92 + t7 * t93 + t6 * t94 + t102 * t100 + t103 * t101 + t42 * t123 + t43 * t124 + t130 * t34 + t126 * t168; t296 * t197 - t295 * t227 + (-t63 - t64 - t267) * t191 + t246 - m(4) * (-t154 * t191 - t155 * t189) - (-t168 - t231) * t189 + t231 * qJD(4) + t216 * t101 + t218 * t100 - t265 * t294 + t361 * t293 - t341 * qJD(1) ^ 2 + (t1 * t197 + t13 * t361 - t14 * t265 - t191 * t41 - t2 * t227) * m(7) + (-t104 * t191 - t17 * t361 - t18 * t265 + t197 * t4 + t227 * t3) * m(6) + (-t148 * t191 - t182 * t232 - t234) * m(5); (-t188 * t350 + t190 * t351) * t308 + (-t188 * t352 + t190 * t353) * t317 + (t1 * t162 + t152 * t9 + t161 * t2 + t347 * t41 + (t134 - t26) * t14 + (t133 - t27) * t13) * m(7) + t347 * t63 - (t154 * mrSges(4,3) + t181 / 0.2e1 + (t362 - Ifges(4,1) / 0.2e1) * t191 - t357) * t189 + t355 * t197 + (t227 * t350 + t235) * t310 - t335 * t227 + (t50 - t47) * (-t128 / 0.2e1 - t188 / 0.2e1) + t348 * (t129 / 0.2e1 + t190 / 0.2e1) + (t128 * t350 - t129 * t351) * t309 + (t128 * t352 - t129 * t353) * t318 + (Ifges(6,4) * t190 - Ifges(7,5) * t129 - Ifges(6,2) * t188 - Ifges(7,3) * t128) * t321 + (-Ifges(6,4) * t129 + Ifges(7,5) * t190 + Ifges(6,2) * t128 + Ifges(7,3) * t188) * t320 + t237 * t315 + t239 * t316 + t59 * t304 + (-mrSges(4,1) - t242) * t112 + m(6) * (-t133 * t17 + t134 * t18 - t161 * t3 + t162 * t4 + t210 * t79) + ((m(6) * t104 + t64) * t301 + t360) * qJD(4) + (t155 * mrSges(4,3) - t339) * t191 - m(6) * (t104 * t114 + t17 * t30 + t18 * t31) - t267 * t155 + t293 * t133 + t294 * t134 + t295 * t161 + t296 * t162 + (-pkin(3) * t112 - t148 * t155 - t85 * t95 - t86 * t96) * m(5) - pkin(3) * t80 - t26 * t91 - t31 * t92 - t30 * t93 - t27 * t94 - t111 * mrSges(4,2) - t114 * t64 - t96 * t123 - t95 * t124 + t152 * t33 - t154 * t168 - Ifges(4,5) * t178 - Ifges(4,6) * t179 + t210 * t34 + t216 * t60 / 0.2e1 + (-t13 * t265 - t14 * t361) * mrSges(7,2) + (t17 * t265 - t18 * t361) * mrSges(6,3) + (mrSges(7,1) * t361 + mrSges(7,3) * t265) * t41 + (mrSges(6,1) * t361 - mrSges(6,2) * t265) * t104 + (-t100 * t216 + t101 * t218 + m(5) * t344 + (-m(5) * t233 - t216 * t123 - t218 * t124) * qJD(4)) * pkin(8) + t344 * mrSges(5,3); -t294 * t20 - t293 * t19 + (t1 * t206 - t13 * t19 + t2 * t209 - t41 * t45 + (qJD(6) - t20) * t14) * m(7) + (Ifges(5,5) * t164 - Ifges(5,6) * t165) * t309 + t98 * t311 + (Ifges(5,1) * t164 - t287) * t312 + t55 * t302 + t56 * t253 + (t291 + t124) * t86 + (-Ifges(6,2) * t320 + Ifges(7,3) * t321 - t350 * t309 - t352 * t318 + t337) * t228 + (-Ifges(5,2) * t165 + t163 + t99) * t313 + t340 + ((t213 * t4 + t279 * t3) * pkin(4) - t104 * t303 + t17 * t19 - t18 * t20) * m(6) + (t292 - t123) * t85 + t343 - t64 * t303 - t45 * t63 + qJD(6) * t91 - t148 * (mrSges(5,1) * t165 + mrSges(5,2) * t164) + t206 * t54 + t209 * t57 + (-Ifges(6,4) * t320 - Ifges(7,5) * t321 - t351 * t309 - t353 * t318 + t356) * t106; -t293 * t228 + t294 * t106 + t33 + t34 + (t106 * t14 - t13 * t228 + t9) * m(7) + (t106 * t18 + t17 * t228 + t79) * m(6); t228 * t63 - t182 * t91 + 0.2e1 * (t2 / 0.2e1 + t41 * t317 + t14 * t309) * m(7) + t57;];
tauc  = t11(:);
