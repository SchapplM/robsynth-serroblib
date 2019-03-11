% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
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
% Datum: 2019-03-09 03:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRPRR4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR4_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR4_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR4_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR4_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:44:28
% EndTime: 2019-03-09 03:44:40
% DurationCPUTime: 5.59s
% Computational Cost: add. (5505->504), mult. (12202->690), div. (0->0), fcn. (7124->8), ass. (0->247)
t186 = cos(qJ(3));
t182 = sin(qJ(5));
t183 = sin(qJ(3));
t245 = t182 * t183;
t198 = pkin(5) * t186 - pkin(9) * t245;
t234 = qJD(5) * t182;
t187 = -pkin(3) - pkin(8);
t268 = pkin(9) - t187;
t172 = t183 * qJD(1);
t168 = pkin(3) * t172;
t206 = pkin(8) * t183 - qJ(4) * t186;
t117 = qJD(1) * t206 + t168;
t185 = cos(qJ(5));
t162 = sin(pkin(10)) * pkin(1) + pkin(7);
t146 = t162 * qJD(1);
t171 = t183 * qJD(2);
t116 = t186 * t146 + t171;
t239 = qJD(1) * t186;
t96 = pkin(4) * t239 + t116;
t51 = -t117 * t182 + t185 * t96;
t328 = -qJD(1) * t198 + t268 * t234 - t51;
t144 = t268 * t185;
t222 = t185 * t172;
t52 = t185 * t117 + t182 * t96;
t327 = pkin(9) * t222 + qJD(5) * t144 + t52;
t217 = pkin(4) * qJD(1) + t146;
t201 = t217 * t183;
t173 = t186 * qJD(2);
t321 = qJD(4) - t173;
t77 = qJD(3) * t187 + t201 + t321;
t224 = -cos(pkin(10)) * pkin(1) - pkin(2);
t197 = -qJ(4) * t183 + t224;
t113 = t186 * t187 + t197;
t88 = t113 * qJD(1);
t38 = -t182 * t88 + t185 * t77;
t39 = t182 * t77 + t185 * t88;
t204 = t182 * t38 - t185 * t39;
t266 = Ifges(6,4) * t182;
t208 = Ifges(6,2) * t185 + t266;
t265 = Ifges(6,4) * t185;
t210 = Ifges(6,1) * t182 + t265;
t213 = mrSges(6,1) * t185 - mrSges(6,2) * t182;
t261 = Ifges(6,6) * t185;
t264 = Ifges(6,5) * t182;
t273 = -t185 / 0.2e1;
t274 = -t182 / 0.2e1;
t161 = t172 + qJD(5);
t275 = -t161 / 0.2e1;
t237 = qJD(3) * t185;
t138 = -t182 * t239 + t237;
t280 = -t138 / 0.2e1;
t137 = -qJD(3) * t182 - t185 * t239;
t281 = -t137 / 0.2e1;
t267 = Ifges(6,4) * t138;
t60 = Ifges(6,2) * t137 + Ifges(6,6) * t161 + t267;
t131 = Ifges(6,4) * t137;
t61 = Ifges(6,1) * t138 + Ifges(6,5) * t161 + t131;
t176 = qJD(3) * qJ(4);
t85 = t176 + t96;
t326 = t204 * mrSges(6,3) + (t261 + t264) * t275 + t208 * t281 + t210 * t280 + t85 * t213 + t273 * t60 + t274 * t61;
t181 = sin(qJ(6));
t184 = cos(qJ(6));
t36 = pkin(9) * t137 + t39;
t253 = t184 * t36;
t35 = -pkin(9) * t138 + t38;
t31 = pkin(5) * t161 + t35;
t10 = t181 * t31 + t253;
t236 = qJD(3) * t186;
t219 = qJD(1) * t236;
t158 = Ifges(7,3) * t219;
t215 = t184 * t137 - t138 * t181;
t71 = t137 * t181 + t138 * t184;
t272 = Ifges(7,4) * t71;
t229 = qJD(5) + qJD(6);
t153 = t172 + t229;
t277 = -t153 / 0.2e1;
t287 = -t71 / 0.2e1;
t289 = -t215 / 0.2e1;
t64 = Ifges(7,4) * t215;
t34 = Ifges(7,1) * t71 + Ifges(7,5) * t153 + t64;
t56 = -pkin(5) * t137 + t85;
t254 = t181 * t36;
t9 = t184 * t31 - t254;
t325 = t158 + (Ifges(7,5) * t215 - Ifges(7,6) * t71) * t277 + (t10 * t71 + t215 * t9) * mrSges(7,3) + (-Ifges(7,2) * t71 + t34 + t64) * t289 - t56 * (mrSges(7,1) * t71 + mrSges(7,2) * t215) + (Ifges(7,1) * t215 - t272) * t287;
t312 = mrSges(5,1) + mrSges(4,3);
t324 = mrSges(5,2) - mrSges(4,1);
t269 = pkin(4) + t162;
t102 = -t176 - t116;
t130 = -pkin(3) * t186 + t197;
t103 = t130 * qJD(1);
t148 = t224 * qJD(1);
t313 = qJD(3) / 0.2e1;
t314 = -qJD(3) / 0.2e1;
t315 = -qJD(1) / 0.2e1;
t323 = t103 * mrSges(5,2) + t116 * mrSges(4,3) + Ifges(4,6) * t313 + (Ifges(4,4) * t183 + t186 * Ifges(4,2)) * qJD(1) / 0.2e1 + Ifges(5,5) * t314 + (-Ifges(5,6) * t183 - t186 * Ifges(5,3)) * t315 - t102 * mrSges(5,1) - t148 * mrSges(4,1) + t326;
t160 = qJD(3) * t168;
t235 = qJD(4) * t183;
t190 = qJD(3) * t206 - t235;
t83 = qJD(1) * t190 + t160;
t84 = (t186 * t217 + t171) * qJD(3);
t16 = -qJD(5) * t39 - t182 * t83 + t185 * t84;
t232 = qJD(5) * t186;
t220 = t185 * t232;
t230 = qJD(3) * qJD(5);
t238 = qJD(3) * t183;
t93 = -t182 * t230 + (t182 * t238 - t220) * qJD(1);
t13 = pkin(5) * t219 - pkin(9) * t93 + t16;
t233 = qJD(5) * t185;
t15 = t182 * t84 + t185 * t83 + t77 * t233 - t234 * t88;
t221 = t182 * t232;
t191 = t183 * t237 + t221;
t94 = qJD(1) * t191 - t185 * t230;
t14 = pkin(9) * t94 + t15;
t2 = qJD(6) * t9 + t13 * t181 + t14 * t184;
t27 = qJD(6) * t215 + t181 * t94 + t184 * t93;
t28 = -qJD(6) * t71 - t181 * t93 + t184 * t94;
t304 = qJD(6) * t10;
t3 = t13 * t184 - t14 * t181 - t304;
t322 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,5) * t27 + Ifges(7,6) * t28;
t244 = t184 * t185;
t246 = t181 * t182;
t199 = -t244 + t246;
t200 = t181 * t185 + t184 * t182;
t193 = t200 * t183;
t110 = qJD(1) * t193;
t79 = t229 * t200;
t250 = -t79 - t110;
t109 = -t172 * t246 + t184 * t222;
t231 = qJD(6) * t181;
t78 = -t181 * t234 - t182 * t231 + t229 * t244;
t251 = t78 + t109;
t320 = -t10 * t251 + t199 * t3 - t2 * t200 - t250 * t9;
t33 = Ifges(7,2) * t215 + Ifges(7,6) * t153 + t272;
t318 = t33 / 0.2e1;
t143 = t268 * t182;
t82 = -t143 * t184 - t144 * t181;
t311 = -qJD(6) * t82 + t327 * t181 + t184 * t328;
t81 = t143 * t181 - t144 * t184;
t310 = qJD(6) * t81 + t181 * t328 - t327 * t184;
t223 = -pkin(5) * t185 - pkin(4);
t305 = -(qJD(1) * t223 - t146) * t183 + pkin(5) * t233 + t321;
t132 = t269 * t183;
t120 = t182 * t132;
t58 = t185 * t113 + t120;
t258 = t15 * t182;
t303 = t16 * t185 + t258;
t115 = t146 * t183 - t173;
t302 = -qJD(4) - t115;
t97 = -mrSges(6,2) * t161 + mrSges(6,3) * t137;
t98 = mrSges(6,1) * t161 - mrSges(6,3) * t138;
t202 = t182 * t98 - t185 * t97;
t65 = mrSges(6,1) * t219 - mrSges(6,3) * t93;
t66 = -mrSges(6,2) * t219 + mrSges(6,3) * t94;
t203 = t182 * t66 + t185 * t65;
t300 = t202 * qJD(5) - t203;
t299 = t186 * t229;
t298 = -m(4) * t116 + m(5) * t102;
t297 = t16 * mrSges(6,1) - t15 * mrSges(6,2) + Ifges(6,5) * t93 + Ifges(6,6) * t94 + t322;
t241 = qJD(3) * t324 + t172 * t312;
t99 = -qJD(3) * pkin(3) - t302;
t296 = m(4) * t115 + m(5) * t99 + t241;
t294 = 0.2e1 * t162;
t293 = m(4) / 0.2e1;
t292 = t27 / 0.2e1;
t291 = t28 / 0.2e1;
t290 = t60 / 0.2e1;
t288 = t215 / 0.2e1;
t286 = t71 / 0.2e1;
t285 = m(7) * t56;
t118 = t199 * t186;
t283 = t118 / 0.2e1;
t119 = t200 * t186;
t282 = -t119 / 0.2e1;
t279 = -t200 / 0.2e1;
t278 = -t199 / 0.2e1;
t276 = t153 / 0.2e1;
t263 = Ifges(6,5) * t185;
t262 = Ifges(6,6) * t182;
t252 = t185 * t38;
t150 = -mrSges(5,1) * t239 - qJD(3) * mrSges(5,3);
t80 = -mrSges(6,1) * t137 + mrSges(6,2) * t138;
t249 = t80 - t150;
t107 = t116 * qJD(3);
t248 = t107 * t186;
t243 = t185 * t186;
t149 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t239;
t242 = t149 - t150;
t133 = t269 * t186;
t37 = -mrSges(7,1) * t215 + mrSges(7,2) * t71;
t228 = -t37 - t249;
t227 = -Ifges(5,4) / 0.2e1 + Ifges(4,5) / 0.2e1;
t226 = Ifges(5,5) / 0.2e1 - Ifges(4,6) / 0.2e1;
t225 = -0.3e1 / 0.2e1 * Ifges(5,6) - 0.3e1 / 0.2e1 * Ifges(4,4);
t218 = pkin(9) * t186 - t113;
t170 = pkin(3) * t238;
t100 = t170 + t190;
t123 = t269 * t236;
t216 = -t100 * t182 + t185 * t123;
t164 = qJD(3) * t173;
t106 = -t146 * t238 + t164;
t212 = mrSges(6,1) * t182 + mrSges(6,2) * t185;
t211 = Ifges(6,1) * t185 - t266;
t209 = -Ifges(6,2) * t182 + t265;
t121 = t185 * t132;
t50 = pkin(5) * t183 + t182 * t218 + t121;
t53 = -pkin(9) * t243 + t58;
t20 = -t181 * t53 + t184 * t50;
t21 = t181 * t50 + t184 * t53;
t194 = -qJ(4) * t236 - t235;
t29 = t185 * t100 - t113 * t234 + t182 * t123 + t132 * t233;
t175 = qJD(3) * qJD(4);
t73 = -qJD(3) * t201 + t164 + t175;
t167 = Ifges(4,4) * t239;
t189 = t115 * mrSges(4,3) + t148 * mrSges(4,2) + t38 * mrSges(6,1) + t9 * mrSges(7,1) + t99 * mrSges(5,1) + Ifges(4,1) * t172 / 0.2e1 + Ifges(4,5) * t313 + t167 / 0.2e1 + Ifges(5,4) * t314 + (-Ifges(5,2) * t183 - Ifges(5,6) * t186) * t315 + t153 * Ifges(7,3) + t71 * Ifges(7,5) + t215 * Ifges(7,6) + t161 * Ifges(6,3) + t138 * Ifges(6,5) + t137 * Ifges(6,6) - t10 * mrSges(7,2) - t103 * mrSges(5,3) - t39 * mrSges(6,2);
t163 = pkin(5) * t182 + qJ(4);
t159 = Ifges(6,3) * t219;
t142 = -qJ(4) * t239 + t168;
t141 = (mrSges(5,2) * t186 - mrSges(5,3) * t183) * qJD(1);
t124 = t170 + t194;
t122 = t269 * t238;
t105 = pkin(5) * t243 + t133;
t104 = qJD(1) * t194 + t160;
t95 = t173 - t201;
t92 = -t106 - t175;
t67 = -pkin(5) * t221 + (-t162 + t223) * t238;
t57 = -t113 * t182 + t121;
t55 = mrSges(7,1) * t153 - mrSges(7,3) * t71;
t54 = -mrSges(7,2) * t153 + mrSges(7,3) * t215;
t49 = -mrSges(6,1) * t94 + mrSges(6,2) * t93;
t47 = -pkin(5) * t94 + t73;
t46 = t93 * Ifges(6,1) + t94 * Ifges(6,4) + Ifges(6,5) * t219;
t45 = t93 * Ifges(6,4) + t94 * Ifges(6,2) + Ifges(6,6) * t219;
t44 = -t199 * t238 + t200 * t299;
t43 = qJD(3) * t193 + t199 * t299;
t30 = -qJD(5) * t58 + t216;
t24 = -mrSges(7,2) * t219 + mrSges(7,3) * t28;
t23 = mrSges(7,1) * t219 - mrSges(7,3) * t27;
t22 = pkin(9) * t191 + t29;
t19 = t198 * qJD(3) + (t185 * t218 - t120) * qJD(5) + t216;
t12 = t184 * t35 - t254;
t11 = -t181 * t35 - t253;
t8 = -mrSges(7,1) * t28 + mrSges(7,2) * t27;
t7 = t27 * Ifges(7,1) + t28 * Ifges(7,4) + Ifges(7,5) * t219;
t6 = t27 * Ifges(7,4) + t28 * Ifges(7,2) + Ifges(7,6) * t219;
t5 = -qJD(6) * t21 - t181 * t22 + t184 * t19;
t4 = qJD(6) * t20 + t181 * t19 + t184 * t22;
t1 = [(t159 / 0.2e1 + t158 / 0.2e1 - t104 * mrSges(5,3) + ((t293 + m(5) / 0.2e1) * t294 + t312) * t107 + ((mrSges(4,1) * t224 - t130 * mrSges(5,2) + t183 * t225) * qJD(1) + t226 * qJD(3) + (-t242 + t298) * t162 - t323) * qJD(3) + t297) * t183 + m(7) * (t10 * t4 + t105 * t47 + t2 * t21 + t20 * t3 + t5 * t9 + t56 * t67) + (t45 * t273 + t104 * mrSges(5,2) + t73 * t213 - t93 * t210 / 0.2e1 - t94 * t208 / 0.2e1 - t92 * mrSges(5,1) + t106 * mrSges(4,3) + t46 * t274 + (-t15 * t185 + t16 * t182) * mrSges(6,3) + (t106 * t293 - m(5) * t92 / 0.2e1) * t294 + (-t85 * t212 + t209 * t281 + t211 * t280 + t161 * (t262 - t263) / 0.2e1 + t61 * t273 + t182 * t290 + (t182 * t39 + t252) * mrSges(6,3)) * qJD(5) + (((-t264 / 0.2e1 - t261 / 0.2e1 - t225) * t186 + Ifges(7,5) * t282 + Ifges(7,6) * t283 - t130 * mrSges(5,3) + t224 * mrSges(4,2) + (Ifges(7,3) / 0.2e1 + Ifges(6,3) / 0.2e1 + 0.3e1 / 0.2e1 * Ifges(5,2) + 0.3e1 / 0.2e1 * Ifges(4,1) - 0.3e1 / 0.2e1 * Ifges(5,3) - 0.3e1 / 0.2e1 * Ifges(4,2)) * t183) * qJD(1) + t189 + t227 * qJD(3) + t296 * t162) * qJD(3)) * t186 + (t10 * t44 + t118 * t2 + t119 * t3 - t43 * t9) * mrSges(7,3) + (-Ifges(7,4) * t119 + Ifges(7,2) * t118) * t291 + (-Ifges(7,1) * t119 + Ifges(7,4) * t118) * t292 + t47 * (-mrSges(7,1) * t118 - mrSges(7,2) * t119) + m(6) * (-t122 * t85 + t133 * t73 + t15 * t58 + t16 * t57 + t29 * t39 + t30 * t38) + m(5) * (t103 * t124 + t104 * t130) + t7 * t282 + t6 * t283 + (Ifges(7,1) * t43 + Ifges(7,4) * t44) * t286 + (Ifges(7,4) * t43 + Ifges(7,2) * t44) * t288 + (Ifges(7,5) * t43 + Ifges(7,6) * t44) * t276 + t44 * t318 + t20 * t23 + t21 * t24 + t43 * t34 / 0.2e1 + t4 * t54 + t5 * t55 + t56 * (-mrSges(7,1) * t44 + mrSges(7,2) * t43) + t57 * t65 + t58 * t66 + t67 * t37 + t29 * t97 + t30 * t98 + t105 * t8 - t122 * t80 + t133 * t49 + t124 * t141; t118 * t23 - t119 * t24 + t43 * t54 + t44 * t55 + (t49 + t8) * t183 + t300 * t186 + m(5) * (-t183 * t92 - t248) + m(4) * (t106 * t183 - t248) + m(7) * (t10 * t43 + t118 * t3 - t119 * t2 + t183 * t47 + t44 * t9) + m(6) * (-t16 * t243 + t183 * t73 - t186 * t258 - t39 * t220 + t38 * t221) + (m(6) * t39 * t245 + t312 * qJD(1) * (-t183 ^ 2 - t186 ^ 2) + (m(6) * t252 + t182 * t97 + t185 * t98 + t296) * t183 + (m(6) * t85 + t149 - t228 + t285 - t298) * t186) * qJD(3); t326 * qJD(5) + t310 * t54 + (t10 * t310 + t163 * t47 + t2 * t82 + t3 * t81 + t305 * t56 + t311 * t9) * m(7) + t311 * t55 + (-pkin(3) * t107 - qJ(4) * t92 + t102 * t302 - t103 * t142 - t116 * t99) * m(5) + (t203 + m(6) * t303 + (-m(6) * t204 - t202) * qJD(5)) * t187 + (((Ifges(4,4) / 0.2e1 + Ifges(5,6) / 0.2e1) * t172 + (-qJ(4) * mrSges(5,1) + t226) * qJD(3) + t323) * t183 + (-Ifges(5,6) * t239 / 0.2e1 + (-pkin(3) * mrSges(5,1) + t263 / 0.2e1 - t262 / 0.2e1 + Ifges(7,5) * t278 + Ifges(7,6) * t279 + t227) * qJD(3) - t189 - t167 / 0.2e1 + (-Ifges(4,1) / 0.2e1 - Ifges(5,2) / 0.2e1 + Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1) * t172) * t186) * qJD(1) + (-t78 / 0.2e1 - t109 / 0.2e1) * t33 + t320 * mrSges(7,3) + (-Ifges(7,1) * t199 - Ifges(7,4) * t200) * t292 + t47 * (mrSges(7,1) * t200 - mrSges(7,2) * t199) + (-Ifges(7,1) * t79 - Ifges(7,4) * t78) * t286 + (-Ifges(7,4) * t79 - Ifges(7,2) * t78) * t288 + (-Ifges(7,5) * t79 - Ifges(7,6) * t78) * t276 + (-t79 / 0.2e1 - t110 / 0.2e1) * t34 + t324 * t107 + t94 * t209 / 0.2e1 + t93 * t211 / 0.2e1 + t73 * t212 - t303 * mrSges(6,3) + (t73 * qJ(4) - t38 * t51 - t39 * t52 + (qJD(4) - t95) * t85) * m(6) + (-Ifges(7,4) * t199 - Ifges(7,2) * t200) * t291 + (Ifges(7,1) * t110 + Ifges(7,4) * t109) * t287 + (Ifges(7,4) * t110 + Ifges(7,2) * t109) * t289 + t45 * t274 + (Ifges(7,5) * t110 + Ifges(7,6) * t109) * t277 + t7 * t278 + t6 * t279 + t305 * t37 + qJ(4) * t49 + t81 * t23 + t82 * t24 - t92 * mrSges(5,3) - t95 * t80 - t52 * t97 - t51 * t98 - t106 * mrSges(4,2) - t241 * t116 + t242 * t115 - t142 * t141 + t249 * qJD(4) + (mrSges(7,1) * t251 + mrSges(7,2) * t250) * t56 + t163 * t8 + t185 * t46 / 0.2e1; t200 * t24 - t199 * t23 + t250 * t55 + t251 * t54 + t228 * qJD(3) + (mrSges(5,1) * t236 + (t141 - t202) * t183) * qJD(1) + (-qJD(3) * t56 - t320) * m(7) + (-qJD(3) * t85 - t161 * t204 + t303) * m(6) + (qJD(3) * t102 + t103 * t172 + t107) * m(5) - t300; (-Ifges(6,2) * t138 + t131 + t61) * t281 + (t137 * t38 + t138 * t39) * mrSges(6,3) + t159 + t297 + t71 * t318 + t138 * t290 + (Ifges(6,5) * t137 - Ifges(6,6) * t138) * t275 + (Ifges(6,1) * t137 - t267) * t280 - m(7) * (t10 * t12 + t11 * t9) + (-m(7) * t9 * t231 + (m(7) * t2 - qJD(6) * t55 + t24) * t181 + (t23 + t54 * qJD(6) + m(7) * (t3 + t304)) * t184 + (-t37 - t285) * t138) * pkin(5) - t12 * t54 - t11 * t55 - t38 * t97 + t39 * t98 - t85 * (mrSges(6,1) * t138 + mrSges(6,2) * t137) + t325; t10 * t55 + t33 * t286 - t9 * t54 + t322 + t325;];
tauc  = t1(:);
