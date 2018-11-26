% Calculate vector of centrifugal and coriolis load on the joints for
% S6PRRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2018-11-23 15:12
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6PRRPRP2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP2_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP2_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP2_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP2_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRP2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:11:53
% EndTime: 2018-11-23 15:12:02
% DurationCPUTime: 9.30s
% Computational Cost: add. (5381->530), mult. (14187->698), div. (0->0), fcn. (10270->10), ass. (0->260)
t176 = sin(qJ(2));
t172 = sin(pkin(6));
t261 = qJD(1) * t172;
t238 = t176 * t261;
t175 = sin(qJ(3));
t256 = qJD(3) * t175;
t356 = pkin(3) * t256 - t238;
t353 = Ifges(6,1) + Ifges(7,1);
t352 = Ifges(7,4) + Ifges(6,5);
t354 = Ifges(6,6) - Ifges(7,6);
t171 = sin(pkin(11));
t178 = cos(qJ(3));
t271 = cos(pkin(11));
t225 = t271 * t175;
t146 = t171 * t178 + t225;
t137 = t146 * qJD(3);
t190 = -t171 * t175 + t178 * t271;
t138 = t190 * qJD(3);
t355 = pkin(4) * t137 - pkin(9) * t138 + t356;
t295 = -qJ(4) - pkin(8);
t227 = qJD(3) * t295;
t254 = qJD(4) * t178;
t132 = t175 * t227 + t254;
t133 = -qJD(4) * t175 + t178 * t227;
t179 = cos(qJ(2));
t237 = t179 * t261;
t345 = t132 * t271 + t171 * t133 - t190 * t237;
t157 = t295 * t175;
t158 = t295 * t178;
t112 = t171 * t157 - t158 * t271;
t174 = sin(qJ(5));
t177 = cos(qJ(5));
t252 = qJD(5) * t177;
t253 = qJD(5) * t174;
t168 = -pkin(3) * t178 - pkin(2);
t98 = -pkin(4) * t190 - pkin(9) * t146 + t168;
t348 = -t112 * t253 + t355 * t174 + t345 * t177 + t98 * t252;
t336 = t177 * t112 + t174 * t98;
t347 = -qJD(5) * t336 - t345 * t174 + t355 * t177;
t135 = t190 * qJD(2);
t307 = -t135 / 0.2e1;
t351 = Ifges(5,2) * t307;
t330 = -t354 * t174 + t352 * t177;
t277 = Ifges(7,5) * t174;
t281 = Ifges(6,4) * t174;
t329 = t353 * t177 + t277 - t281;
t342 = qJD(5) - t135;
t350 = qJ(6) * t137 - qJD(6) * t190 + t348;
t349 = -pkin(5) * t137 - t347;
t109 = t146 * t237;
t81 = t132 * t171 - t271 * t133;
t346 = t109 - t81;
t344 = t352 * t174 + t354 * t177;
t276 = Ifges(7,5) * t177;
t280 = Ifges(6,4) * t177;
t343 = t353 * t174 - t276 + t280;
t341 = -Ifges(5,6) / 0.2e1;
t126 = qJD(2) * t137;
t127 = qJD(2) * t138;
t257 = qJD(2) * t178;
t136 = -qJD(2) * t225 - t171 * t257;
t194 = t177 * qJD(3) + t136 * t174;
t72 = qJD(5) * t194 + t127 * t177;
t116 = qJD(3) * t174 - t136 * t177;
t73 = qJD(5) * t116 + t127 * t174;
t339 = (-Ifges(6,4) + Ifges(7,5)) * t73 + t353 * t72 + t352 * t126;
t114 = Ifges(6,4) * t194;
t278 = Ifges(7,5) * t194;
t338 = t116 * t353 + t342 * t352 + t114 - t278;
t40 = mrSges(6,1) * t126 - mrSges(6,3) * t72;
t41 = -t126 * mrSges(7,1) + t72 * mrSges(7,2);
t294 = t41 - t40;
t42 = -mrSges(6,2) * t126 - mrSges(6,3) * t73;
t43 = -mrSges(7,2) * t73 + mrSges(7,3) * t126;
t293 = t42 + t43;
t74 = mrSges(7,2) * t194 + mrSges(7,3) * t342;
t286 = mrSges(6,3) * t194;
t75 = -mrSges(6,2) * t342 + t286;
t292 = t74 + t75;
t285 = mrSges(6,3) * t116;
t76 = mrSges(6,1) * t342 - t285;
t77 = -mrSges(7,1) * t342 + mrSges(7,2) * t116;
t291 = t76 - t77;
t128 = qJD(2) * t168 + qJD(4) - t237;
t337 = t128 * mrSges(5,2);
t203 = pkin(5) * t174 - qJ(6) * t177;
t173 = cos(pkin(6));
t262 = t173 * t178;
t163 = qJD(1) * t262;
t150 = qJD(2) * pkin(8) + t238;
t224 = qJ(4) * qJD(2) + t150;
t107 = -t175 * t224 + t163;
t260 = qJD(1) * t175;
t236 = t173 * t260;
t108 = t178 * t224 + t236;
t226 = t271 * t108;
t52 = t107 * t171 + t226;
t335 = -qJD(6) * t174 + t203 * t342 - t52;
t287 = mrSges(5,3) * t136;
t272 = -qJD(3) * mrSges(5,1) - mrSges(6,1) * t194 + mrSges(6,2) * t116 - t287;
t334 = Ifges(5,5) * qJD(3);
t251 = qJD(2) * qJD(3);
t230 = t175 * t251;
t103 = qJD(3) * pkin(3) + t107;
t50 = t171 * t103 + t226;
t45 = qJD(3) * pkin(9) + t50;
t63 = -pkin(4) * t135 + pkin(9) * t136 + t128;
t16 = -t174 * t45 + t177 * t63;
t333 = qJD(6) - t16;
t118 = t150 * t178 + t236;
t255 = qJD(3) * t178;
t332 = -t118 * qJD(3) + (-qJ(4) * t255 + (-qJD(4) - t237) * t175) * qJD(2);
t99 = t171 * t108;
t49 = t103 * t271 - t99;
t44 = -qJD(3) * pkin(4) - t49;
t19 = -pkin(5) * t194 - t116 * qJ(6) + t44;
t217 = mrSges(7,1) * t174 - mrSges(7,3) * t177;
t219 = mrSges(6,1) * t174 + mrSges(6,2) * t177;
t331 = -t19 * t217 - t44 * t219;
t267 = t135 * t177;
t328 = t252 - t267;
t268 = t135 * t174;
t327 = -t253 + t268;
t263 = t172 * t179;
t234 = qJD(2) * t263;
t223 = t178 * t234;
t94 = qJD(1) * t223 + qJD(3) * t163 - t150 * t256;
t71 = (-qJ(4) * t256 + t254) * qJD(2) + t94;
t28 = t171 * t332 + t271 * t71;
t258 = qJD(2) * t176;
t235 = t172 * t258;
t134 = pkin(3) * t230 + qJD(1) * t235;
t58 = pkin(4) * t126 - pkin(9) * t127 + t134;
t3 = t174 * t58 + t177 * t28 + t63 * t252 - t253 * t45;
t17 = t174 * t63 + t177 * t45;
t4 = -qJD(5) * t17 - t174 * t28 + t177 * t58;
t326 = -t174 * t4 + t177 * t3;
t1 = qJ(6) * t126 + qJD(6) * t342 + t3;
t2 = -pkin(5) * t126 - t4;
t325 = t1 * t177 + t174 * t2;
t323 = qJD(3) * t341;
t283 = Ifges(5,4) * t136;
t322 = t351 + t283 / 0.2e1;
t206 = Ifges(7,3) * t174 + t276;
t212 = -Ifges(6,2) * t174 + t280;
t304 = t174 / 0.2e1;
t305 = -t174 / 0.2e1;
t311 = t116 / 0.2e1;
t313 = -t194 / 0.2e1;
t314 = t194 / 0.2e1;
t113 = Ifges(7,5) * t116;
t34 = Ifges(7,6) * t342 - Ifges(7,3) * t194 + t113;
t282 = Ifges(6,4) * t116;
t37 = Ifges(6,2) * t194 + Ifges(6,6) * t342 + t282;
t321 = t206 * t313 + t212 * t314 + t34 * t304 + t37 * t305 - t331 + t329 * t311 + t330 * t342 / 0.2e1;
t10 = -pkin(5) * t342 + t333;
t11 = qJ(6) * t342 + t17;
t320 = t128 * mrSges(5,1) + t16 * mrSges(6,1) - t10 * mrSges(7,1) - t17 * mrSges(6,2) + t11 * mrSges(7,3) + t322 + t323;
t180 = qJD(2) ^ 2;
t318 = t72 / 0.2e1;
t317 = -t73 / 0.2e1;
t316 = t73 / 0.2e1;
t312 = -t116 / 0.2e1;
t310 = t126 / 0.2e1;
t309 = -t342 / 0.2e1;
t306 = t136 / 0.2e1;
t303 = -t177 / 0.2e1;
t302 = t177 / 0.2e1;
t301 = pkin(3) * t171;
t27 = t171 * t71 - t271 * t332;
t264 = t172 * t176;
t139 = -t175 * t264 + t262;
t140 = t173 * t175 + t178 * t264;
t89 = -t139 * t271 + t140 * t171;
t296 = t27 * t89;
t54 = t107 * t271 - t99;
t259 = qJD(2) * t175;
t249 = pkin(3) * t259;
t84 = -pkin(4) * t136 - pkin(9) * t135 + t249;
t25 = t174 * t84 + t177 * t54;
t290 = mrSges(5,3) * t126;
t289 = mrSges(5,3) * t127;
t288 = mrSges(5,3) * t135;
t284 = Ifges(4,4) * t175;
t129 = Ifges(5,4) * t135;
t111 = -t271 * t157 - t158 * t171;
t275 = t111 * t27;
t273 = t136 * Ifges(5,1);
t270 = Ifges(4,5) * qJD(3);
t269 = Ifges(4,6) * qJD(3);
t61 = -mrSges(7,1) * t194 - mrSges(7,3) * t116;
t250 = -t61 - t272;
t247 = Ifges(6,5) / 0.2e1 + Ifges(7,4) / 0.2e1;
t246 = Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1;
t245 = Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1;
t244 = mrSges(4,3) * t259;
t243 = mrSges(4,3) * t257;
t242 = t174 * t263;
t239 = t271 * pkin(3);
t83 = t126 * mrSges(5,1) + t127 * mrSges(5,2);
t167 = -t239 - pkin(4);
t222 = -t1 * t174 + t2 * t177;
t221 = -t3 * t174 - t4 * t177;
t220 = mrSges(6,1) * t177 - mrSges(6,2) * t174;
t218 = mrSges(7,1) * t177 + mrSges(7,3) * t174;
t211 = Ifges(6,2) * t177 + t281;
t205 = -Ifges(7,3) * t177 + t277;
t204 = t177 * pkin(5) + t174 * qJ(6);
t202 = t10 * t177 - t11 * t174;
t201 = t10 * t174 + t11 * t177;
t200 = -t16 * t177 - t17 * t174;
t199 = t16 * t174 - t17 * t177;
t24 = -t174 * t54 + t177 * t84;
t95 = -t150 * t255 + (-qJD(3) * t173 - t234) * t260;
t196 = -t95 * t175 + t94 * t178;
t46 = -t112 * t174 + t177 * t98;
t90 = t171 * t139 + t140 * t271;
t64 = t174 * t90 + t177 * t263;
t191 = (mrSges(4,1) * t175 + mrSges(4,2) * t178) * qJD(2);
t182 = t4 * mrSges(6,1) - t2 * mrSges(7,1) - t3 * mrSges(6,2) + t1 * mrSges(7,3);
t169 = Ifges(4,4) * t257;
t156 = -qJD(3) * mrSges(4,2) + t243;
t155 = qJD(3) * mrSges(4,1) - t244;
t151 = -qJD(2) * pkin(2) - t237;
t144 = -t204 + t167;
t143 = qJD(3) * t191;
t142 = Ifges(4,1) * t259 + t169 + t270;
t141 = t269 + (t178 * Ifges(4,2) + t284) * qJD(2);
t124 = Ifges(7,2) * t126;
t123 = Ifges(6,3) * t126;
t119 = -qJD(3) * mrSges(5,2) + t288;
t117 = -t150 * t175 + t163;
t106 = -qJD(3) * t140 - t175 * t234;
t105 = qJD(3) * t139 + t223;
t96 = -mrSges(5,1) * t135 - mrSges(5,2) * t136;
t92 = t129 - t273 + t334;
t70 = Ifges(7,4) * t72;
t69 = Ifges(6,5) * t72;
t68 = Ifges(6,6) * t73;
t67 = Ifges(7,6) * t73;
t65 = t177 * t90 - t242;
t60 = pkin(5) * t116 - qJ(6) * t194;
t55 = t146 * t203 + t111;
t53 = t105 * t271 + t171 * t106;
t51 = t105 * t171 - t106 * t271;
t36 = t116 * Ifges(7,4) + Ifges(7,2) * t342 - Ifges(7,6) * t194;
t35 = t116 * Ifges(6,5) + Ifges(6,6) * t194 + Ifges(6,3) * t342;
t33 = pkin(5) * t190 - t46;
t32 = -qJ(6) * t190 + t336;
t31 = mrSges(6,1) * t73 + mrSges(6,2) * t72;
t30 = mrSges(7,1) * t73 - mrSges(7,3) * t72;
t21 = t72 * Ifges(6,4) - t73 * Ifges(6,2) + t126 * Ifges(6,6);
t20 = t72 * Ifges(7,5) + t126 * Ifges(7,6) + t73 * Ifges(7,3);
t18 = pkin(5) * t136 - t24;
t15 = -qJ(6) * t136 + t25;
t14 = -qJD(5) * t242 + t174 * t53 - t177 * t235 + t252 * t90;
t13 = -qJD(5) * t64 + t174 * t235 + t177 * t53;
t12 = t203 * t138 + (qJD(5) * t204 - qJD(6) * t177) * t146 + t81;
t5 = pkin(5) * t73 - qJ(6) * t72 - qJD(6) * t116 + t27;
t6 = [-t90 * t290 + t105 * t156 + t106 * t155 + t53 * t119 + t293 * t65 + t294 * t64 - t291 * t14 + t292 * t13 + (-t139 * t178 - t140 * t175) * mrSges(4,3) * t251 + (t30 + t31 + t289) * t89 - t250 * t51 + ((-mrSges(3,2) * t180 - t143 - t83) * t179 + (-mrSges(3,1) * t180 + (qJD(2) * (-mrSges(4,1) * t178 + mrSges(4,2) * t175) + t96) * qJD(2)) * t176) * t172 + m(7) * (t1 * t65 + t10 * t14 + t11 * t13 + t19 * t51 + t2 * t64 + t5 * t89) + m(6) * (t13 * t17 - t14 * t16 + t3 * t65 - t4 * t64 + t44 * t51 + t296) + m(5) * (t296 + t28 * t90 - t49 * t51 + t50 * t53 + (t128 * t258 - t134 * t179) * t172) + m(4) * (t105 * t118 + t106 * t117 + t139 * t95 + t140 * t94 + (t151 - t237) * t235); t345 * t119 + (Ifges(5,1) * t127 - Ifges(5,4) * t126 + t206 * t316 + t5 * t217 + t212 * t317 + t134 * mrSges(5,2) + t20 * t304 + (mrSges(5,3) + t219) * t27 + t221 * mrSges(6,3) + t222 * mrSges(7,2) + (-mrSges(7,2) * t201 + mrSges(6,3) * t199 + t19 * t218 + t205 * t314 + t211 * t313 + t220 * t44 + t303 * t37 + t309 * t344 + t312 * t343) * qJD(5) + t329 * t318 + t330 * t310 + (qJD(5) * t338 + t21) * t305 + (qJD(5) * t34 + t339) * t302) * t146 + t250 * t109 + (t112 * t28 + t128 * t356 + t134 * t168 + t345 * t50 + t346 * t49 + t275) * m(5) + t347 * t76 + (t16 * t347 + t17 * t348 + t3 * t336 - t346 * t44 + t4 * t46 + t275) * m(6) + t348 * t75 + t349 * t77 + (t1 * t32 + t2 * t33 + t5 * t55 + (-t109 + t12) * t19 + t350 * t11 + t349 * t10) * m(7) + t350 * t74 + (t111 * t127 - t112 * t126) * mrSges(5,3) + t196 * mrSges(4,3) + t272 * t81 + (((-t117 * t178 - t118 * t175) * qJD(3) + t196) * pkin(8) - (pkin(2) * t258 + t151 * t176 + (-t117 * t175 + t118 * t178) * t179) * t261) * m(4) + (-t50 * mrSges(5,3) + t35 / 0.2e1 + t36 / 0.2e1 + t246 * t342 + t247 * t116 - t245 * t194 + t320 + t322) * t137 + (-t49 * mrSges(5,3) + t337 + t321 + t202 * mrSges(7,2) + t200 * mrSges(6,3) + t338 * t302 + t129 / 0.2e1 + t92 / 0.2e1 - t273 / 0.2e1) * t138 + (Ifges(5,5) * t138 / 0.2e1 + t137 * t341 + (t151 * mrSges(4,2) - pkin(8) * t155 - t117 * mrSges(4,3) + t142 / 0.2e1 + t270 / 0.2e1) * t178 + (t151 * mrSges(4,1) - pkin(8) * t156 - t118 * mrSges(4,3) + pkin(3) * t96 - t141 / 0.2e1 - t269 / 0.2e1 + (0.3e1 / 0.2e1 * Ifges(4,1) - 0.3e1 / 0.2e1 * Ifges(4,2)) * t257) * t175) * qJD(3) + t336 * t42 - (-Ifges(5,4) * t127 + t134 * mrSges(5,1) + t69 / 0.2e1 - t68 / 0.2e1 + t123 / 0.2e1 + t70 / 0.2e1 + t124 / 0.2e1 + t67 / 0.2e1 - t28 * mrSges(5,3) + t245 * t73 + t247 * t72 + (Ifges(5,2) + t246) * t126 + t182) * t190 + t168 * t83 - pkin(2) * t143 + (0.3e1 / 0.2e1 * t178 ^ 2 - 0.3e1 / 0.2e1 * t175 ^ 2) * Ifges(4,4) * t251 + t33 * t41 + t32 * t43 + t46 * t40 + t55 * t30 + t12 * t61 + ((t155 * t175 - t156 * t178) * t179 - t176 * t96) * t261 + t111 * t31; t343 * t318 + t344 * t310 + t321 * qJD(5) + (-t16 * t24 + t167 * t27 - t17 * t25 - t44 * t52) * m(6) + (-t10 * t18 - t11 * t15 + t144 * t5 + t19 * t335) * m(7) - t96 * t249 + t141 * t259 / 0.2e1 + (t283 + t36 + t35) * t306 - (-Ifges(4,2) * t259 + t142 + t169) * t257 / 0.2e1 - t5 * t218 - t151 * t191 + (t37 / 0.2e1 - t34 / 0.2e1) * t268 - Ifges(4,6) * t230 / 0.2e1 + (t206 * t314 + t212 * t313 + Ifges(5,1) * t306 - t334 / 0.2e1 - t337 + t329 * t312 + t330 * t309 + t331) * t135 + t251 * Ifges(4,5) * t178 / 0.2e1 + (t244 + t155) * t118 + (t243 - t156) * t117 + (t252 / 0.2e1 - t267 / 0.2e1) * t338 + (t129 + t92) * t307 + (t10 * t328 + t11 * t327 + t325) * mrSges(7,2) + (-t16 * t328 + t17 * t327 + t326) * mrSges(6,3) + t339 * t304 - t272 * t52 + t335 * t61 + ((m(6) * t200 + m(7) * t202) * qJD(5) - t292 * t253 + t293 * t177 + t294 * t174 - t291 * t252 + t326 * m(6) + t325 * m(7)) * (pkin(9) + t301) - t50 * t287 - t239 * t289 + (-Ifges(7,6) * t314 - Ifges(6,6) * t313 + t351 + t323 - t352 * t312 + (-Ifges(6,3) - Ifges(7,2)) * t309 + t320) * t136 + t205 * t316 + t211 * t317 - t290 * t301 - t175 * t180 * (Ifges(4,1) * t178 - t284) / 0.2e1 + t167 * t31 + t144 * t30 + (-t128 * t249 + t49 * t52 - t50 * t54 + (t171 * t28 - t27 * t271) * pkin(3)) * m(5) - t28 * mrSges(5,2) + t49 * t288 + t21 * t302 + t20 * t303 - t15 * t74 - t25 * t75 - t24 * t76 - t18 * t77 - t94 * mrSges(4,2) + t95 * mrSges(4,1) + (-t220 - mrSges(5,1)) * t27 - t54 * t119 - Ifges(5,6) * t126 + Ifges(5,5) * t127; -t135 * t119 - t250 * t136 + (t292 * t342 - t294) * t177 + (-t291 * t342 + t293) * t174 + t83 + (t136 * t19 + t201 * t342 - t222) * m(7) + (t136 * t44 - t199 * t342 - t221) * m(6) + (-t135 * t50 - t136 * t49 + t134) * m(5); t182 + (-t10 * t194 + t11 * t116) * mrSges(7,2) + t124 + t123 + t70 - t68 + t69 + t67 + (Ifges(7,3) * t116 + t278) * t314 + (t285 + t291) * t17 + (t286 - t292) * t16 + t37 * t311 - pkin(5) * t41 + qJ(6) * t43 - t60 * t61 + qJD(6) * t74 - t19 * (mrSges(7,1) * t116 - mrSges(7,3) * t194) - t44 * (mrSges(6,1) * t116 + mrSges(6,2) * t194) + (-t354 * t116 + t194 * t352) * t309 + (-pkin(5) * t2 + qJ(6) * t1 - t10 * t17 + t11 * t333 - t19 * t60) * m(7) + (-Ifges(6,2) * t116 + t114 + t338) * t313 + (t194 * t353 + t113 - t282 + t34) * t312; t116 * t61 - t342 * t74 + 0.2e1 * (t2 / 0.2e1 + t11 * t309 + t19 * t311) * m(7) + t41;];
tauc  = t6(:);
