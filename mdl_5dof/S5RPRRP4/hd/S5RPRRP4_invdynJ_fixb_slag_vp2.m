% Calculate vector of inverse dynamics joint torques for
% S5RPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% m [6x1]
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:33
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRP4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP4_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP4_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP4_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP4_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP4_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:31:51
% EndTime: 2022-01-23 09:32:11
% DurationCPUTime: 11.12s
% Computational Cost: add. (3827->432), mult. (9228->566), div. (0->0), fcn. (6276->10), ass. (0->227)
t334 = Ifges(5,4) + Ifges(6,4);
t335 = Ifges(5,1) + Ifges(6,1);
t333 = Ifges(5,5) + Ifges(6,5);
t332 = Ifges(5,2) + Ifges(6,2);
t331 = Ifges(5,6) + Ifges(6,6);
t195 = sin(pkin(8));
t198 = sin(qJ(3));
t201 = cos(qJ(3));
t248 = qJD(1) * qJD(3);
t137 = (qJDD(1) * t201 - t198 * t248) * t195;
t361 = Ifges(4,5) * t137;
t138 = (-qJDD(1) * t198 - t201 * t248) * t195;
t360 = Ifges(4,6) * t138;
t196 = cos(pkin(8));
t245 = qJDD(1) * t196;
t172 = qJDD(3) - t245;
t359 = Ifges(4,3) * t172;
t197 = sin(qJ(4));
t200 = cos(qJ(4));
t214 = t197 * t198 - t200 * t201;
t210 = t214 * qJD(4);
t255 = qJD(1) * t195;
t52 = -t137 * t197 + t138 * t200 + t210 * t255;
t358 = t331 * t52;
t153 = t197 * t201 + t198 * t200;
t211 = qJD(1) * t153;
t109 = t195 * t211;
t51 = -qJD(4) * t109 + t137 * t200 + t138 * t197;
t357 = t333 * t51;
t164 = qJDD(4) + t172;
t356 = (Ifges(5,3) + Ifges(6,3)) * t164;
t247 = qJDD(1) * qJ(2);
t249 = qJD(1) * qJD(2);
t169 = t247 + t249;
t355 = t169 + t249;
t232 = t201 * t255;
t233 = t198 * t255;
t111 = -t197 * t233 + t200 * t232;
t161 = pkin(2) * t196 + pkin(6) * t195 + pkin(1);
t136 = -qJDD(1) * t161 + qJDD(2);
t270 = t196 * t198;
t139 = -qJD(1) * t161 + qJD(2);
t254 = qJD(1) * t196;
t235 = qJ(2) * t254;
t97 = t139 * t198 + t201 * t235;
t50 = -qJD(3) * t97 + t201 * t136 - t169 * t270;
t21 = pkin(3) * t172 - pkin(7) * t137 + t50;
t174 = qJD(3) - t254;
t123 = t201 * t139;
t237 = qJ(2) * t270;
t272 = t195 * t201;
t244 = pkin(7) * t272;
t209 = -t237 - t244;
t79 = qJD(1) * t209 + t123;
t62 = pkin(3) * t174 + t79;
t80 = -pkin(7) * t233 + t97;
t74 = t200 * t80;
t24 = t197 * t62 + t74;
t222 = qJD(3) * t237;
t252 = qJD(3) * t201;
t268 = t196 * t201;
t49 = -qJD(1) * t222 + t198 * t136 + t139 * t252 + t169 * t268;
t27 = pkin(7) * t138 + t49;
t6 = -qJD(4) * t24 - t197 * t27 + t200 * t21;
t2 = pkin(4) * t164 - qJ(5) * t51 - qJD(5) * t111 + t6;
t250 = qJD(4) * t200;
t251 = qJD(4) * t197;
t5 = t197 * t21 + t200 * t27 + t62 * t250 - t251 * t80;
t3 = qJ(5) * t52 - qJD(5) * t109 + t5;
t354 = t6 * mrSges(5,1) + t2 * mrSges(6,1) - t5 * mrSges(5,2) - t3 * mrSges(6,2);
t353 = -t50 * mrSges(4,1) + t49 * mrSges(4,2);
t338 = -mrSges(6,1) - mrSges(5,1);
t337 = -mrSges(6,2) - mrSges(5,2);
t351 = t334 * t109;
t350 = t334 * t111;
t194 = qJ(3) + qJ(4);
t185 = sin(t194);
t294 = pkin(3) * t198;
t157 = pkin(4) * t185 + t294;
t343 = -m(6) - m(3);
t319 = -m(4) + t343;
t349 = t319 * qJ(2) - m(5) * (qJ(2) + t294) - m(6) * t157 + mrSges(2,2) - mrSges(3,3);
t165 = qJD(4) + t174;
t348 = -t332 * t109 + t331 * t165 + t350;
t347 = t335 * t111 + t333 * t165 - t351;
t192 = t195 ^ 2;
t193 = t196 ^ 2;
t346 = mrSges(3,3) * (t193 + t192);
t345 = t355 * t192;
t189 = t201 * pkin(3);
t220 = mrSges(3,1) * t196 - mrSges(3,2) * t195;
t306 = pkin(7) + pkin(6);
t344 = m(4) * t161 + m(5) * (pkin(1) + (pkin(2) + t189) * t196) + mrSges(2,1) + t220 + (m(5) * t306 + mrSges(4,3) + mrSges(5,3) + mrSges(6,3)) * t195;
t36 = -mrSges(6,2) * t164 + mrSges(6,3) * t52;
t37 = -mrSges(5,2) * t164 + mrSges(5,3) * t52;
t330 = t36 + t37;
t286 = mrSges(6,3) * t109;
t83 = -mrSges(6,2) * t165 - t286;
t287 = mrSges(5,3) * t109;
t84 = -mrSges(5,2) * t165 - t287;
t329 = t84 + t83;
t282 = t111 * mrSges(6,3);
t85 = mrSges(6,1) * t165 - t282;
t283 = t111 * mrSges(5,3);
t86 = mrSges(5,1) * t165 - t283;
t328 = t85 + t86;
t103 = t111 * qJ(5);
t72 = t197 * t80;
t23 = t200 * t62 - t72;
t13 = -t103 + t23;
t312 = m(5) * pkin(3);
t327 = t312 + mrSges(4,1);
t318 = -qJD(3) - qJD(4);
t326 = (t254 + t318) * t214;
t95 = t318 * t153;
t325 = t196 * t211 + t95;
t173 = qJ(2) * t268;
t107 = -t198 * t161 + t173;
t323 = t356 + t357 + t358;
t186 = cos(t194);
t199 = sin(qJ(1));
t202 = cos(qJ(1));
t267 = t196 * t202;
t127 = -t185 * t267 + t186 * t199;
t128 = t185 * t199 + t186 * t267;
t322 = t127 * t338 - t337 * t128;
t269 = t196 * t199;
t125 = t185 * t269 + t186 * t202;
t126 = t185 * t202 - t186 * t269;
t321 = -t125 * t338 + t337 * t126;
t320 = mrSges(5,1) * t185 - t186 * t337;
t317 = t193 * t355;
t284 = Ifges(4,4) * t201;
t285 = Ifges(4,4) * t198;
t316 = (-qJ(2) * (mrSges(4,1) * t201 - mrSges(4,2) * t198) - t201 * (-Ifges(4,1) * t198 - t284) / 0.2e1 + t198 * (-Ifges(4,2) * t201 - t285) / 0.2e1) * t192;
t212 = (-t198 * Ifges(4,2) + t284) * t195;
t213 = (t201 * Ifges(4,1) - t285) * t195;
t315 = t201 * (t174 * Ifges(4,6) + qJD(1) * t212) + t198 * (t174 * Ifges(4,5) + qJD(1) * t213);
t203 = qJD(1) ^ 2;
t311 = m(6) * pkin(4);
t302 = t111 / 0.2e1;
t293 = pkin(3) * t200;
t292 = pkin(4) * t111;
t291 = g(3) * t195;
t29 = t200 * t79 - t72;
t150 = t201 * t161;
t87 = -t244 - t150 + (-qJ(2) * t198 - pkin(3)) * t196;
t274 = t195 * t198;
t98 = -pkin(7) * t274 + t107;
t39 = t197 * t87 + t200 * t98;
t277 = qJ(2) * t203;
t276 = qJ(5) * t109;
t184 = t195 * qJ(2);
t156 = t195 * t169;
t266 = t198 * t199;
t265 = t198 * t202;
t264 = t199 * t201;
t263 = t201 * t202;
t258 = t345 * qJ(2);
t253 = qJD(2) * t196;
t257 = -t161 * t252 + t201 * t253;
t140 = pkin(3) * t233 + qJ(2) * t255;
t158 = pkin(4) * t186 + t189;
t183 = t195 * qJD(2);
t230 = t195 * t252;
t147 = pkin(3) * t230 + t183;
t151 = pkin(3) * t274 + t184;
t246 = qJDD(1) * t195;
t241 = mrSges(4,3) * t274;
t236 = t359 + t360 + t361;
t231 = t198 * t253;
t229 = -t52 * mrSges(6,1) + t51 * mrSges(6,2);
t93 = -pkin(3) * t138 + t156;
t28 = -t197 * t79 - t74;
t38 = -t197 * t98 + t200 * t87;
t225 = pkin(3) * t232;
t224 = mrSges(4,3) * t233;
t223 = mrSges(4,3) * t232;
t221 = -mrSges(3,1) * t245 + mrSges(3,2) * t246;
t219 = mrSges(4,1) * t198 + mrSges(4,2) * t201;
t133 = -mrSges(4,2) * t174 - t224;
t134 = mrSges(4,1) * t174 - t223;
t216 = t133 * t201 - t134 * t198;
t215 = (pkin(2) + t158) * t196 - (-qJ(5) - t306) * t195;
t143 = -t196 * t265 + t264;
t141 = t196 * t266 + t263;
t75 = qJD(3) * t209 + t257;
t76 = -t231 + (-t173 + (pkin(7) * t195 + t161) * t198) * qJD(3);
t9 = t197 * t76 + t200 * t75 + t87 * t250 - t251 * t98;
t207 = t174 * t195 * (-Ifges(4,5) * t198 - Ifges(4,6) * t201);
t10 = -qJD(4) * t39 - t197 * t75 + t200 * t76;
t12 = pkin(4) * t165 + t13;
t14 = t24 - t276;
t77 = pkin(4) * t109 + qJD(5) + t140;
t204 = t14 * t282 + t24 * t283 - t77 * (mrSges(6,1) * t111 - mrSges(6,2) * t109) - t140 * (mrSges(5,1) * t111 - mrSges(5,2) * t109) + t323 - t12 * t286 - t23 * t287 - (-t335 * t109 - t350) * t111 / 0.2e1 + t348 * t302 - (-t109 * t333 - t111 * t331) * t165 / 0.2e1 + (-t332 * t111 + t347 - t351) * t109 / 0.2e1 + t354;
t182 = -qJDD(1) * pkin(1) + qJDD(2);
t181 = pkin(4) + t293;
t176 = t192 * t277;
t145 = t219 * t195;
t144 = t196 * t263 + t266;
t142 = -t196 * t264 + t265;
t135 = t219 * t255;
t132 = t214 * t195;
t131 = t153 * t195;
t106 = -t150 - t237;
t100 = -mrSges(4,2) * t172 + mrSges(4,3) * t138;
t99 = mrSges(4,1) * t172 - mrSges(4,3) * t137;
t96 = -t198 * t235 + t123;
t92 = t225 + t292;
t90 = -qJD(3) * t107 - t231;
t89 = -t222 + t257;
t88 = pkin(4) * t131 + t151;
t71 = (qJD(3) * t214 + t210) * t195;
t70 = t95 * t195;
t66 = mrSges(5,1) * t109 + mrSges(5,2) * t111;
t65 = mrSges(6,1) * t109 + mrSges(6,2) * t111;
t53 = -pkin(4) * t71 + t147;
t35 = mrSges(5,1) * t164 - mrSges(5,3) * t51;
t34 = mrSges(6,1) * t164 - mrSges(6,3) * t51;
t25 = -qJ(5) * t131 + t39;
t22 = -pkin(4) * t196 + qJ(5) * t132 + t38;
t18 = -pkin(4) * t52 + qJDD(5) + t93;
t17 = -t103 + t29;
t16 = t28 + t276;
t8 = -qJ(5) * t70 + qJD(5) * t132 + t10;
t7 = qJ(5) * t71 - qJD(5) * t131 + t9;
t1 = [t138 * t212 / 0.2e1 + t137 * t213 / 0.2e1 + (t317 + t345) * mrSges(3,3) - pkin(1) * t221 - t182 * t220 + (t172 * (Ifges(4,5) * t201 - Ifges(4,6) * t198) / 0.2e1 - t315 * qJD(3) / 0.2e1 + Ifges(3,4) * t245 + Ifges(3,1) * t246) * t195 + (t207 / 0.2e1 + t96 * t241) * qJD(3) - (Ifges(4,4) * t137 + Ifges(4,2) * t138 + Ifges(4,6) * t172) * t274 / 0.2e1 + t145 * t156 + (-t360 / 0.2e1 - t361 / 0.2e1 - t356 / 0.2e1 - t359 / 0.2e1 - t358 / 0.2e1 - t357 / 0.2e1 - t236 / 0.2e1 - t323 / 0.2e1 + Ifges(3,2) * t245 + Ifges(3,4) * t246 + t353 - t354) * t196 + (-t12 * t70 + t14 * t71) * mrSges(6,3) + (-t23 * t70 + t24 * t71) * mrSges(5,3) + (-t97 * t230 - t50 * t272) * mrSges(4,3) + t147 * t66 + t151 * (-mrSges(5,1) * t52 + mrSges(5,2) * t51) + t140 * (-mrSges(5,1) * t71 + mrSges(5,2) * t70) + t89 * t133 + t90 * t134 + t106 * t99 + t107 * t100 + t9 * t84 + t8 * t85 + t10 * t86 + t77 * (-mrSges(6,1) * t71 + mrSges(6,2) * t70) + t7 * t83 + t53 * t65 + t22 * t34 + t25 * t36 + t38 * t35 + t39 * t37 + t347 * t70 / 0.2e1 + t348 * t71 / 0.2e1 + (mrSges(5,1) * t93 + mrSges(6,1) * t18 - mrSges(5,3) * t5 - mrSges(6,3) * t3 - t164 * t331 - t332 * t52 - t334 * t51) * t131 + m(5) * (t10 * t23 + t140 * t147 + t151 * t93 + t24 * t9 + t38 * t6 + t39 * t5) + m(6) * (t12 * t8 + t14 * t7 + t18 * t88 + t2 * t22 + t25 * t3 + t53 * t77) + (-t142 * mrSges(4,1) - t141 * mrSges(4,2) + t337 * t125 + t338 * t126 + t349 * t202 + (m(3) * pkin(1) - m(6) * (-pkin(1) - t215) + t344) * t199) * g(1) + (-t144 * mrSges(4,1) - t143 * mrSges(4,2) + t337 * t127 + t338 * t128 + (-m(6) * t215 + t343 * pkin(1) - t344) * t202 + t349 * t199) * g(2) + (t331 * t71 + t333 * t70) * t165 / 0.2e1 - (t332 * t71 + t334 * t70) * t109 / 0.2e1 + (t334 * t71 + t335 * t70) * t302 + m(3) * (-pkin(1) * t182 + qJ(2) * t317 + t258) + (-mrSges(5,2) * t93 - mrSges(6,2) * t18 + mrSges(5,3) * t6 + mrSges(6,3) * t2 - t164 * t333 - t334 * t52 - t335 * t51) * t132 - t316 * t248 + t88 * t229 - t49 * t241 + (-mrSges(4,1) * t138 + mrSges(4,2) * t137) * t184 + t247 * t346 + (Ifges(4,1) * t137 + Ifges(4,4) * t138 + Ifges(4,5) * t172) * t272 / 0.2e1 + t135 * t183 + m(4) * (t106 * t50 + t107 * t49 + t89 * t97 + t90 * t96 + t258) + Ifges(2,3) * qJDD(1); t216 * qJD(3) + (-t216 * t196 + (-t135 - t65 - t66) * t195) * qJD(1) + t330 * t153 - (t34 + t35) * t214 - t203 * t346 + t198 * t100 + t201 * t99 + t221 + t326 * t329 + t325 * t328 + (-g(1) * t199 + g(2) * t202) * (m(5) - t319) + (t12 * t325 + t14 * t326 + t153 * t3 - t2 * t214 - t77 * t255) * m(6) + (-t140 * t255 + t153 * t5 - t214 * t6 + t23 * t325 + t24 * t326) * m(5) + (t49 * t198 + t50 * t201 - t176 + t174 * (-t96 * t198 + t97 * t201)) * m(4) + (-t193 * t277 - t176 + t182) * m(3); (t197 * t5 + t200 * t6 + (-t197 * t23 + t200 * t24) * qJD(4)) * t312 - t353 - qJD(1) * t207 / 0.2e1 + t236 + t181 * t34 + (t134 + t223) * t97 + g(3) * t145 + (-t224 - t133) * t96 - t29 * t84 - t16 * t85 - t28 * t86 - t92 * t65 - t17 * t83 - t66 * t225 - m(5) * (t140 * t225 + t23 * t28 + t24 * t29) + (-mrSges(4,2) * t142 - m(6) * (-t157 * t269 - t158 * t202) + t327 * t141 + t321) * g(2) + (mrSges(4,2) * t144 - m(6) * (-t157 * t267 + t158 * t199) - t327 * t143 + t322) * g(1) + (m(5) * t294 + mrSges(6,1) * t185 + t320) * t291 + t316 * t203 + t315 * t255 / 0.2e1 + t35 * t293 + t204 + (-t12 * t16 - t14 * t17 + t157 * t291 + t181 * t2 - t77 * t92) * m(6) + ((t197 * t3 + (-t12 * t197 + t14 * t200) * qJD(4)) * m(6) + t330 * t197 - t328 * t251 + t329 * t250) * pkin(3); -m(6) * (t77 * t292 + (-t12 + t13) * t14) - t23 * t84 + t14 * t85 + t24 * t86 - t13 * t83 + pkin(4) * t34 - t65 * t292 + t2 * t311 + t204 + (-(-mrSges(6,1) - t311) * t185 + t320) * t291 + (t125 * t311 + t321) * g(2) + (-t127 * t311 + t322) * g(1); t109 * t83 + t111 * t85 + (g(3) * t196 + t109 * t14 + t12 * t111 + t18 + (-g(1) * t202 - g(2) * t199) * t195) * m(6) + t229;];
tau = t1;
