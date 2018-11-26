% Calculate vector of centrifugal and coriolis load on the joints for
% S6RPRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
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
% Datum: 2018-11-23 16:17
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RPRRPR3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR3_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR3_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR3_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR3_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:16:53
% EndTime: 2018-11-23 16:17:01
% DurationCPUTime: 7.55s
% Computational Cost: add. (5672->586), mult. (13072->757), div. (0->0), fcn. (7875->8), ass. (0->274)
t367 = Ifges(5,1) + Ifges(6,1);
t359 = Ifges(6,4) + Ifges(5,5);
t191 = sin(qJ(4));
t258 = qJD(4) * t191;
t195 = cos(qJ(3));
t263 = qJD(1) * t195;
t318 = pkin(8) - pkin(9);
t192 = sin(qJ(3));
t264 = qJD(1) * t192;
t181 = sin(pkin(10)) * pkin(1) + pkin(7);
t168 = t181 * qJD(1);
t137 = t195 * qJD(2) - t192 * t168;
t236 = pkin(3) * t192 - pkin(8) * t195;
t161 = t236 * qJD(1);
t194 = cos(qJ(4));
t79 = t194 * t137 + t191 * t161;
t68 = qJ(5) * t264 + t79;
t370 = pkin(9) * t191 * t263 + t318 * t258 + t68;
t173 = t318 * t194;
t268 = t194 * t195;
t252 = pkin(9) * t268;
t319 = pkin(4) + pkin(5);
t78 = -t191 * t137 + t161 * t194;
t369 = qJD(4) * t173 - (-t192 * t319 - t252) * qJD(1) + t78;
t262 = qJD(2) * t192;
t138 = t168 * t195 + t262;
t127 = qJD(3) * pkin(8) + t138;
t246 = -cos(pkin(10)) * pkin(1) - pkin(2);
t149 = -pkin(3) * t195 - pkin(8) * t192 + t246;
t130 = t149 * qJD(1);
t62 = -t191 * t127 + t194 * t130;
t337 = qJD(5) - t62;
t368 = qJD(3) / 0.2e1;
t366 = Ifges(5,6) - Ifges(6,6);
t254 = t194 * qJD(3);
t155 = t191 * t264 - t254;
t261 = qJD(3) * t191;
t156 = t194 * t264 + t261;
t190 = sin(qJ(6));
t193 = cos(qJ(6));
t212 = t155 * t190 + t156 * t193;
t89 = t155 * t193 - t156 * t190;
t83 = Ifges(7,4) * t89;
t365 = Ifges(7,2) * t212 - t83;
t364 = -pkin(9) * t156 + t337;
t184 = Ifges(4,4) * t263;
t244 = Ifges(4,5) * t368;
t63 = t194 * t127 + t191 * t130;
t213 = t191 * t63 + t194 * t62;
t180 = qJD(4) - t263;
t52 = -pkin(4) * t180 + t337;
t175 = t180 * qJ(5);
t55 = t175 + t63;
t214 = t191 * t55 - t194 * t52;
t287 = Ifges(6,5) * t194;
t221 = Ifges(6,3) * t191 + t287;
t292 = Ifges(5,4) * t194;
t225 = -Ifges(5,2) * t191 + t292;
t230 = mrSges(6,1) * t191 - mrSges(6,3) * t194;
t232 = mrSges(5,1) * t191 + mrSges(5,2) * t194;
t237 = qJD(3) * pkin(3) + t137;
t285 = Ifges(6,6) * t191;
t286 = Ifges(5,6) * t191;
t290 = Ifges(5,5) * t194;
t291 = Ifges(6,4) * t194;
t303 = t194 / 0.2e1;
t305 = t191 / 0.2e1;
t306 = -t191 / 0.2e1;
t311 = t156 / 0.2e1;
t313 = t155 / 0.2e1;
t314 = -t155 / 0.2e1;
t152 = Ifges(5,4) * t155;
t289 = Ifges(6,5) * t155;
t343 = t156 * t367 + t180 * t359 - t152 + t289;
t347 = t180 / 0.2e1;
t288 = Ifges(6,5) * t191;
t293 = Ifges(5,4) * t191;
t353 = t194 * t367 + t288 - t293;
t204 = qJ(5) * t156 + t237;
t61 = pkin(4) * t155 - t204;
t151 = Ifges(6,5) * t156;
t72 = Ifges(6,6) * t180 + Ifges(6,3) * t155 + t151;
t294 = Ifges(5,4) * t156;
t75 = -Ifges(5,2) * t155 + Ifges(5,6) * t180 + t294;
t329 = t214 * mrSges(6,2) + t213 * mrSges(5,3) + t237 * t232 - t221 * t313 - t225 * t314 - t61 * t230 - t72 * t305 - t75 * t306 - t353 * t311 - (-t286 + t290 + t285 + t291) * t347 - t343 * t303;
t361 = t264 / 0.2e1;
t363 = Ifges(4,1) * t361 + t184 / 0.2e1 + t244 - t137 * mrSges(4,3) - t329;
t260 = qJD(3) * t192;
t241 = qJD(1) * t260;
t176 = qJD(6) - t180;
t310 = -t176 / 0.2e1;
t245 = t192 * t258;
t253 = qJD(3) * qJD(4);
t114 = t194 * t253 + (t195 * t254 - t245) * qJD(1);
t257 = qJD(4) * t194;
t259 = qJD(3) * t195;
t115 = t191 * t253 + (t191 * t259 + t192 * t257) * qJD(1);
t27 = qJD(6) * t89 + t114 * t193 + t115 * t190;
t28 = -qJD(6) * t212 - t114 * t190 + t115 * t193;
t344 = Ifges(7,5) * t27 + Ifges(7,6) * t28;
t41 = -t155 * t319 + t204;
t32 = -t180 * t319 + t364;
t43 = pkin(9) * t155 + t63;
t37 = t175 + t43;
t8 = -t190 * t37 + t193 * t32;
t9 = t190 * t32 + t193 * t37;
t362 = (t212 * t9 + t8 * t89) * mrSges(7,3) - Ifges(7,3) * t241 + (Ifges(7,5) * t89 - Ifges(7,6) * t212) * t310 - t41 * (mrSges(7,1) * t212 + mrSges(7,2) * t89) + t344;
t360 = -qJD(3) / 0.2e1;
t172 = t318 * t191;
t107 = t172 * t190 + t173 * t193;
t358 = -qJD(6) * t107 + t190 * t370 + t193 * t369;
t106 = t172 * t193 - t173 * t190;
t357 = qJD(6) * t106 + t190 * t369 - t193 * t370;
t356 = t359 * t241 + (-Ifges(5,4) + Ifges(6,5)) * t115 + t367 * t114;
t270 = qJ(5) * t194;
t207 = -t191 * t319 + t270;
t256 = qJD(5) * t191;
t355 = t262 - (qJD(1) * t207 - t168) * t195 + qJD(4) * t207 + t256;
t354 = t191 * t367 - t287 + t292;
t302 = Ifges(7,4) * t212;
t350 = Ifges(7,1) * t89 - t302;
t325 = t27 / 0.2e1;
t324 = t28 / 0.2e1;
t35 = Ifges(7,1) * t212 + Ifges(7,5) * t176 + t83;
t348 = t35 / 0.2e1;
t321 = -t212 / 0.2e1;
t346 = -t241 / 0.2e1;
t243 = Ifges(4,6) * t360;
t166 = t193 * qJ(5) - t190 * t319;
t339 = -qJD(6) * t166 - t190 * t364 - t193 * t43;
t165 = -t190 * qJ(5) - t193 * t319;
t338 = qJD(6) * t165 - t190 * t43 + t193 * t364;
t211 = t190 * t194 - t191 * t193;
t140 = t211 * t192;
t129 = qJD(3) * t138;
t336 = -qJ(5) * t114 - qJD(5) * t156 + t129;
t335 = t191 * t359 + t194 * t366;
t128 = t137 * qJD(3);
t164 = t236 * qJD(3);
t148 = qJD(1) * t164;
t23 = -t127 * t258 + t194 * t128 + t130 * t257 + t191 * t148;
t24 = -t127 * t257 - t191 * t128 - t130 * t258 + t148 * t194;
t215 = -t191 * t24 + t194 * t23;
t17 = qJ(5) * t241 + t180 * qJD(5) + t23;
t18 = -pkin(4) * t241 - t24;
t217 = t17 * t194 + t18 * t191;
t334 = qJD(4) - qJD(6);
t271 = qJ(5) * t191;
t331 = -t194 * t319 - t271;
t248 = mrSges(4,3) * t264;
t275 = -qJD(3) * mrSges(4,1) + mrSges(5,1) * t155 + mrSges(5,2) * t156 + t248;
t330 = -m(4) * t137 - m(5) * t237 + t275;
t170 = t246 * qJD(1);
t249 = Ifges(6,6) / 0.2e1 - Ifges(5,6) / 0.2e1;
t250 = Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1;
t251 = Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1;
t328 = t249 * t155 + t251 * t156 + t250 * t180 + t170 * mrSges(4,1) + t55 * mrSges(6,3) + t62 * mrSges(5,1) + t9 * mrSges(7,2) + t243 - (Ifges(4,4) * t192 + Ifges(4,2) * t195) * qJD(1) / 0.2e1 - t176 * Ifges(7,3) - t212 * Ifges(7,5) - t89 * Ifges(7,6) + Ifges(5,6) * t314 + Ifges(6,6) * t313 - t52 * mrSges(6,1) - t63 * mrSges(5,2) - t8 * mrSges(7,1) + (Ifges(5,3) + Ifges(6,2)) * t347 + t359 * t311;
t327 = Ifges(7,4) * t325 + Ifges(7,2) * t324 + Ifges(7,6) * t346;
t326 = Ifges(7,1) * t325 + Ifges(7,4) * t324 + Ifges(7,5) * t346;
t323 = -t89 / 0.2e1;
t322 = t89 / 0.2e1;
t320 = t212 / 0.2e1;
t317 = t114 / 0.2e1;
t316 = -t115 / 0.2e1;
t315 = t115 / 0.2e1;
t312 = -t156 / 0.2e1;
t309 = t176 / 0.2e1;
t308 = -t180 / 0.2e1;
t304 = -t194 / 0.2e1;
t301 = pkin(9) * t192;
t38 = -mrSges(7,1) * t89 + mrSges(7,2) * t212;
t98 = mrSges(6,1) * t155 - mrSges(6,3) * t156;
t297 = -t38 + t98;
t296 = mrSges(5,3) * t155;
t295 = mrSges(5,3) * t156;
t282 = t170 * mrSges(4,2);
t206 = t211 * t195;
t131 = qJD(1) * t206;
t96 = t334 * t211;
t277 = -t131 + t96;
t210 = t190 * t191 + t193 * t194;
t205 = t210 * t195;
t132 = qJD(1) * t205;
t95 = t334 * t210;
t276 = t132 - t95;
t272 = qJ(5) * t155;
t269 = t181 * t191;
t119 = -mrSges(5,2) * t180 - t296;
t122 = -mrSges(6,2) * t155 + mrSges(6,3) * t180;
t267 = t119 + t122;
t120 = mrSges(5,1) * t180 - t295;
t121 = -mrSges(6,1) * t180 + mrSges(6,2) * t156;
t266 = -t120 + t121;
t265 = t149 * t257 + t191 * t164;
t160 = t181 * t268;
t94 = t191 * t149 + t160;
t255 = qJD(5) * t194;
t247 = mrSges(4,3) * t263;
t242 = m(4) * t181 + mrSges(4,3);
t240 = -pkin(4) - t269;
t159 = t195 * t269;
t93 = t149 * t194 - t159;
t10 = pkin(9) * t115 + t17;
t13 = -pkin(9) * t114 - t241 * t319 - t24;
t1 = qJD(6) * t8 + t10 * t193 + t13 * t190;
t2 = -qJD(6) * t9 - t10 * t190 + t13 * t193;
t238 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t81 = -qJ(5) * t195 + t94;
t234 = qJD(4) * t160 + t149 * t258 - t164 * t194;
t233 = mrSges(5,1) * t194 - mrSges(5,2) * t191;
t231 = mrSges(6,1) * t194 + mrSges(6,3) * t191;
t224 = Ifges(5,2) * t194 + t293;
t220 = -Ifges(6,3) * t194 + t288;
t219 = pkin(4) * t194 + t271;
t218 = pkin(4) * t191 - t270;
t187 = t195 * pkin(4);
t64 = pkin(5) * t195 + t159 + t187 + (-t149 - t301) * t194;
t71 = t191 * t301 + t81;
t21 = -t190 * t71 + t193 * t64;
t22 = t190 * t64 + t193 * t71;
t65 = -mrSges(7,2) * t176 + mrSges(7,3) * t89;
t66 = mrSges(7,1) * t176 - mrSges(7,3) * t212;
t216 = -t190 * t66 + t193 * t65;
t86 = -mrSges(6,1) * t241 + t114 * mrSges(6,2);
t208 = t181 + t218;
t203 = -t181 + t207;
t84 = -mrSges(6,2) * t115 + mrSges(6,3) * t241;
t85 = mrSges(5,1) * t241 - mrSges(5,3) * t114;
t87 = -mrSges(5,2) * t241 - mrSges(5,3) * t115;
t202 = (t84 + t87) * t194 + (-t85 + t86) * t191;
t201 = -t191 * t267 + t194 * t266;
t44 = (-t192 * t254 - t195 * t258) * t181 + t265;
t200 = -t24 * mrSges(5,1) + t18 * mrSges(6,1) + t23 * mrSges(5,2) - t17 * mrSges(6,3) + t238;
t183 = qJ(5) * t260;
t179 = Ifges(6,2) * t241;
t178 = Ifges(5,3) * t241;
t171 = -qJD(3) * mrSges(4,2) + t247;
t167 = -pkin(3) - t219;
t150 = pkin(3) - t331;
t143 = qJD(4) * t218 - t256;
t141 = t210 * t192;
t116 = t208 * t192;
t105 = Ifges(6,4) * t114;
t104 = Ifges(5,5) * t114;
t103 = Ifges(5,6) * t115;
t102 = Ifges(6,6) * t115;
t97 = pkin(4) * t156 + t272;
t92 = t203 * t192;
t82 = t187 - t93;
t80 = t262 + (qJD(1) * t218 + t168) * t195;
t70 = -pkin(4) * t264 - t78;
t67 = -t156 * t319 - t272;
t59 = (qJD(4) * t219 - t255) * t192 + t208 * t259;
t58 = mrSges(5,1) * t115 + mrSges(5,2) * t114;
t57 = mrSges(6,1) * t115 - mrSges(6,3) * t114;
t49 = t114 * Ifges(5,4) - t115 * Ifges(5,2) + Ifges(5,6) * t241;
t48 = t114 * Ifges(6,5) + Ifges(6,6) * t241 + t115 * Ifges(6,3);
t47 = qJD(3) * t205 + t140 * t334;
t46 = -qJD(3) * t206 + t192 * t95;
t45 = t260 * t269 - t234;
t40 = t240 * t260 + t234;
t39 = (qJD(4) * t331 + t255) * t192 + t203 * t259;
t36 = -qJD(5) * t195 + t183 + t44;
t34 = Ifges(7,2) * t89 + Ifges(7,6) * t176 + t302;
t31 = pkin(4) * t115 + t336;
t30 = t183 + (pkin(9) * qJD(4) - qJD(3) * t181) * t194 * t192 + (-qJD(5) + (pkin(9) * qJD(3) - qJD(4) * t181) * t191) * t195 + t265;
t29 = pkin(9) * t245 + (-t252 + (-pkin(5) + t240) * t192) * qJD(3) + t234;
t20 = mrSges(7,2) * t241 + mrSges(7,3) * t28;
t19 = -mrSges(7,1) * t241 - mrSges(7,3) * t27;
t16 = -t115 * t319 - t336;
t7 = -mrSges(7,1) * t28 + mrSges(7,2) * t27;
t4 = -qJD(6) * t22 - t190 * t30 + t193 * t29;
t3 = qJD(6) * t21 + t190 * t29 + t193 * t30;
t5 = [((t244 + t330 * t181 + 0.2e1 * t282 + 0.3e1 / 0.2e1 * t184 + t363) * qJD(3) - t179 / 0.2e1 - t178 / 0.2e1 - t249 * t115 - t105 / 0.2e1 - t102 / 0.2e1 + t103 / 0.2e1 - t104 / 0.2e1 - t251 * t114 + t200 + t242 * t128 + t344) * t195 + (Ifges(7,1) * t141 - Ifges(7,4) * t140) * t325 + (Ifges(7,4) * t141 - Ifges(7,2) * t140) * t324 + (-t1 * t140 - t141 * t2 + t46 * t9 - t47 * t8) * mrSges(7,3) + t16 * (mrSges(7,1) * t140 + mrSges(7,2) * t141) + t141 * t326 - t140 * t327 + (Ifges(7,1) * t47 + Ifges(7,4) * t46) * t320 + (Ifges(7,4) * t47 + Ifges(7,2) * t46) * t322 + (Ifges(7,5) * t47 + Ifges(7,6) * t46) * t309 + m(5) * (t23 * t94 + t24 * t93 + t44 * t63 + t45 * t62) + m(7) * (t1 * t22 + t16 * t92 + t2 * t21 + t3 * t9 + t39 * t41 + t4 * t8) + m(6) * (t116 * t31 + t17 * t81 + t18 * t82 + t36 * t55 + t40 * t52 + t59 * t61) + t47 * t348 + (t31 * t230 + t221 * t315 + t225 * t316 + t48 * t305 + t49 * t306 + (-t191 * t23 - t194 * t24) * mrSges(5,3) + (-t17 * t191 + t18 * t194) * mrSges(6,2) + (mrSges(4,3) + t232) * t129 + (t61 * t231 - t237 * t233 + t220 * t314 + t224 * t313 + t75 * t304 + (t191 * t62 - t194 * t63) * mrSges(5,3) + (-t191 * t52 - t194 * t55) * mrSges(6,2) + t354 * t312 + t335 * t308 + t343 * t306) * qJD(4) + ((t246 * mrSges(4,1) + (-0.3e1 / 0.2e1 * Ifges(4,4) + t291 / 0.2e1 + t285 / 0.2e1 + t290 / 0.2e1 - t286 / 0.2e1) * t192 - Ifges(7,5) * t141 / 0.2e1 + Ifges(7,6) * t140 / 0.2e1 + (0.3e1 / 0.2e1 * Ifges(4,1) - 0.3e1 / 0.2e1 * Ifges(4,2) - Ifges(7,3) - t250) * t195) * qJD(1) + t243 - t242 * t138 + t328) * qJD(3) + t353 * t317 + (qJD(4) * t72 + t356) * t303 + (t58 + (m(4) + m(5)) * t129 - qJD(3) * t171) * t181) * t192 + t21 * t19 + t22 * t20 + t39 * t38 + t46 * t34 / 0.2e1 + t41 * (-mrSges(7,1) * t46 + mrSges(7,2) * t47) + t3 * t65 + t4 * t66 + t81 * t84 + t82 * t86 + t92 * t7 + t93 * t85 + t94 * t87 + t59 * t98 + t116 * t57 + t44 * t119 + t45 * t120 + t40 * t121 + t36 * t122; t47 * t65 + t46 * t66 - t140 * t19 + t141 * t20 + m(7) * (t1 * t141 - t140 * t2 + t46 * t8 + t47 * t9) + (-t57 - t58 + t7 + (t191 * t266 + t194 * t267 + t171 - t247) * qJD(3) + m(7) * t16 + m(6) * (t254 * t55 + t261 * t52 - t31) + m(5) * (t254 * t63 - t261 * t62 - t129)) * t195 + (t201 * qJD(4) + m(4) * t128 + m(6) * (t257 * t52 - t258 * t55 + t217) + m(5) * (-t257 * t62 - t258 * t63 + t215) + t202 + (m(6) * t61 - m(7) * t41 - t248 + t297 + t330) * qJD(3)) * t192; (((-Ifges(7,5) * t211 - Ifges(7,6) * t210) * t360 + t138 * mrSges(4,3) + t243 + Ifges(4,4) * t361 + t335 * t368 - t328) * t192 + ((-Ifges(4,1) / 0.2e1 + Ifges(4,2) / 0.2e1) * t264 - t184 / 0.2e1 - t282 + t244 - t363) * t195) * qJD(1) + (-t329 + (-m(5) * t213 - m(6) * t214 + t201) * pkin(8)) * qJD(4) + (Ifges(7,1) * t132 - Ifges(7,4) * t131) * t321 + (Ifges(7,4) * t132 - Ifges(7,2) * t131) * t323 + (Ifges(7,5) * t132 - Ifges(7,6) * t131) * t310 - m(5) * (-t138 * t237 + t62 * t78 + t63 * t79) - t210 * t327 - t211 * t326 + (-Ifges(7,4) * t211 - Ifges(7,2) * t210) * t324 + (-Ifges(7,1) * t211 - Ifges(7,4) * t210) * t325 + (-t1 * t210 + t2 * t211 + t276 * t8 - t277 * t9) * mrSges(7,3) + t16 * (mrSges(7,1) * t210 - mrSges(7,2) * t211) + (Ifges(7,1) * t95 - Ifges(7,4) * t96) * t320 + (Ifges(7,4) * t95 - Ifges(7,2) * t96) * t322 + (Ifges(7,5) * t95 - Ifges(7,6) * t96) * t309 + (-t96 / 0.2e1 + t131 / 0.2e1) * t34 + t202 * pkin(8) + t220 * t315 + t224 * t316 + t49 * t303 + t48 * t304 - t31 * t231 + (-mrSges(4,1) - t233) * t129 + t215 * mrSges(5,3) + t217 * mrSges(6,2) - t275 * t138 + (mrSges(7,1) * t277 - mrSges(7,2) * t276) * t41 - m(6) * (t52 * t70 + t55 * t68 + t61 * t80) + (t95 / 0.2e1 - t132 / 0.2e1) * t35 + (-t80 + t143) * t98 + t354 * t317 + t355 * t38 + t356 * t305 + t357 * t65 + t358 * t66 + (t1 * t107 + t106 * t2 + t150 * t16 + t355 * t41 + t357 * t9 + t358 * t8) * m(7) - pkin(3) * t58 + t106 * t19 + t107 * t20 + m(6) * (pkin(8) * t217 + t143 * t61 + t167 * t31) + m(5) * (-pkin(3) * t129 + pkin(8) * t215) - t79 * t119 - t78 * t120 - t70 * t121 - t68 * t122 - t128 * mrSges(4,2) + t150 * t7 + t167 * t57 - t137 * t171; (-t359 * t155 - t156 * t366) * t308 + (-t155 * t367 + t151 - t294 + t72) * t312 + (t155 * t52 + t156 * t55) * mrSges(6,2) + t89 * t348 + t237 * (mrSges(5,1) * t156 - mrSges(5,2) * t155) + (Ifges(6,3) * t156 - t289) * t314 + t179 + t178 + t75 * t311 + t105 + t102 - t103 + t104 - t200 + (-t266 + t295) * t63 + (-t267 - t296) * t62 + t365 * t323 + (t34 - t350) * t321 - t67 * t38 + qJ(5) * t84 - pkin(4) * t86 - t97 * t98 + (-pkin(4) * t18 + qJ(5) * t17 + t337 * t55 - t52 * t63 - t61 * t97) * m(6) + qJD(5) * t122 + t338 * t65 + t339 * t66 + (t1 * t166 + t165 * t2 + t338 * t9 + t339 * t8 - t41 * t67) * m(7) - t61 * (mrSges(6,1) * t156 + mrSges(6,3) * t155) + (-Ifges(5,2) * t156 - t152 + t343) * t313 + t165 * t19 + t166 * t20 - t362; t193 * t19 + t190 * t20 + t297 * t156 + t216 * qJD(6) + (-t122 - t216) * t180 + t86 + (t1 * t190 - t156 * t41 + t193 * t2 + t176 * (-t190 * t8 + t193 * t9)) * m(7) + (t156 * t61 - t180 * t55 + t18) * m(6); t350 * t321 + t34 * t320 - t8 * t65 + t9 * t66 + t238 + (t35 - t365) * t323 + t362;];
tauc  = t5(:);
