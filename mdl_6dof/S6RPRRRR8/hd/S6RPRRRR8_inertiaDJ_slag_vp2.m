% Calculate time derivative of joint inertia matrix for
% S6RPRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:22
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRR8_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR8_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR8_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR8_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR8_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR8_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR8_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:19:47
% EndTime: 2019-03-09 07:19:56
% DurationCPUTime: 4.04s
% Computational Cost: add. (7129->450), mult. (14602->654), div. (0->0), fcn. (13516->8), ass. (0->218)
t182 = sin(qJ(5));
t184 = cos(qJ(5));
t158 = -mrSges(6,1) * t184 + mrSges(6,2) * t182;
t327 = -mrSges(5,1) + t158;
t301 = sin(qJ(4));
t302 = sin(qJ(3));
t303 = cos(qJ(4));
t304 = cos(qJ(3));
t146 = t301 * t302 - t303 * t304;
t325 = -qJD(3) - qJD(4);
t118 = t325 * t146;
t179 = t182 ^ 2;
t180 = t184 ^ 2;
t320 = t179 + t180;
t221 = t118 * t320;
t259 = qJD(5) * t184;
t181 = sin(qJ(6));
t183 = cos(qJ(6));
t200 = t181 * t182 - t183 * t184;
t316 = qJD(5) + qJD(6);
t115 = t316 * t200;
t148 = t181 * t184 + t182 * t183;
t116 = t316 * t148;
t261 = -Ifges(7,5) * t115 - Ifges(7,6) * t116;
t326 = Ifges(6,5) * t259 + t261;
t260 = qJD(5) * t182;
t234 = t146 * t260;
t147 = t301 * t304 + t303 * t302;
t117 = t325 * t147;
t270 = t117 * t184;
t194 = t234 + t270;
t324 = -Ifges(4,1) + Ifges(4,2);
t185 = -pkin(1) - pkin(7);
t323 = -pkin(8) + t185;
t226 = qJD(4) * t303;
t219 = pkin(3) * t226;
t202 = t184 * t219;
t172 = pkin(3) * t301 + pkin(9);
t296 = -pkin(10) - t172;
t222 = qJD(5) * t296;
t132 = t182 * t222 + t202;
t203 = t182 * t219;
t133 = t184 * t222 - t203;
t141 = t296 * t182;
t178 = t184 * pkin(10);
t263 = t172 * t184;
t142 = t178 + t263;
t97 = t141 * t183 - t142 * t181;
t44 = qJD(6) * t97 + t132 * t183 + t133 * t181;
t98 = t141 * t181 + t142 * t183;
t45 = -qJD(6) * t98 - t132 * t181 + t133 * t183;
t322 = t45 * mrSges(7,1) - t44 * mrSges(7,2);
t308 = -pkin(10) - pkin(9);
t161 = t308 * t182;
t162 = pkin(9) * t184 + t178;
t129 = t161 * t183 - t162 * t181;
t241 = qJD(5) * t308;
t152 = t182 * t241;
t153 = t184 * t241;
t75 = qJD(6) * t129 + t152 * t183 + t153 * t181;
t130 = t161 * t181 + t162 * t183;
t76 = -qJD(6) * t130 - t152 * t181 + t153 * t183;
t321 = t76 * mrSges(7,1) - t75 * mrSges(7,2);
t92 = t148 * t146;
t93 = t200 * t147;
t227 = qJD(3) * t302;
t228 = qJD(3) * t304;
t319 = -mrSges(4,1) * t227 - mrSges(4,2) * t228;
t266 = t146 * t182;
t105 = -mrSges(6,2) * t147 + mrSges(6,3) * t266;
t265 = t146 * t184;
t106 = mrSges(6,1) * t147 + mrSges(6,3) * t265;
t318 = -t105 * t260 - t106 * t259;
t271 = t117 * t182;
t195 = t146 * t259 - t271;
t317 = -t303 * t117 - t301 * t118;
t163 = pkin(3) * t228 + qJD(2);
t59 = pkin(4) * t118 - pkin(9) * t117 + t163;
t156 = t323 * t304;
t144 = t156 * qJD(3);
t155 = t323 * t302;
t192 = qJD(3) * t155;
t225 = qJD(4) * t301;
t65 = t144 * t303 - t155 * t225 + t156 * t226 - t192 * t301;
t223 = -t182 * t65 + t184 * t59;
t170 = t302 * pkin(3) + qJ(2);
t101 = pkin(4) * t147 + pkin(9) * t146 + t170;
t125 = t155 * t303 + t156 * t301;
t107 = t184 * t125;
t61 = t182 * t101 + t107;
t15 = -qJD(5) * t61 + t223;
t277 = t15 * t182;
t60 = t184 * t101 - t125 * t182;
t315 = -t60 * t259 - t61 * t260 - t277;
t121 = mrSges(7,1) * t200 + mrSges(7,2) * t148;
t208 = mrSges(6,1) * t182 + mrSges(6,2) * t184;
t149 = t208 * qJD(5);
t26 = -t116 * t147 - t118 * t200;
t28 = -t148 * t118 + t316 * t93;
t283 = t200 * mrSges(7,3);
t286 = t116 * mrSges(7,3);
t67 = mrSges(7,1) * t116 - mrSges(7,2) * t115;
t91 = t148 * t147;
t314 = -t118 * mrSges(5,2) - t26 * t283 + t93 * t286 + mrSges(6,3) * t221 + (t149 + t67) * t146 + (-t121 - t327) * t117 + (-t115 * t91 - t148 * t28) * mrSges(7,3);
t313 = (t129 * t115 - t130 * t116 - t76 * t148 - t200 * t75) * mrSges(7,3);
t312 = 2 * m(7);
t311 = -2 * mrSges(5,3);
t310 = 0.2e1 * qJD(2);
t309 = m(5) * pkin(3);
t300 = pkin(4) * t149;
t299 = t184 * pkin(5);
t294 = Ifges(6,4) * t182;
t293 = Ifges(6,4) * t184;
t292 = Ifges(6,6) * t182;
t291 = pkin(3) * qJD(4);
t290 = pkin(5) * qJD(6);
t289 = t115 * mrSges(7,3);
t42 = pkin(5) * t147 + pkin(10) * t265 + t60;
t48 = pkin(10) * t266 + t61;
t16 = -t181 * t48 + t183 * t42;
t288 = t115 * t16;
t124 = t155 * t301 - t303 * t156;
t66 = qJD(4) * t125 + t144 * t301 + t303 * t192;
t285 = t124 * t66;
t14 = t101 * t259 - t125 * t260 + t182 * t59 + t184 * t65;
t284 = t14 * t184;
t282 = t146 * mrSges(5,3);
t281 = t146 * t66;
t280 = t147 * t65;
t279 = t148 * mrSges(7,3);
t173 = -pkin(3) * t303 - pkin(4);
t157 = t173 - t299;
t276 = t157 * t67;
t174 = -pkin(4) - t299;
t275 = t174 * t67;
t49 = mrSges(6,1) * t118 - mrSges(6,3) * t194;
t274 = t182 * t49;
t269 = t124 * t117;
t86 = t146 * t117;
t218 = pkin(3) * t225;
t247 = pkin(5) * t260;
t154 = t218 + t247;
t264 = t154 * t121;
t262 = t173 * t149;
t258 = qJD(6) * t181;
t257 = qJD(6) * t183;
t256 = t118 * t311;
t255 = m(6) * t291;
t254 = mrSges(7,3) * t290;
t253 = t44 * t283;
t252 = t45 * t279;
t251 = t97 * t289;
t250 = t98 * t286;
t27 = -t200 * t117 + t316 * t92;
t29 = -t115 * t146 - t117 * t148;
t248 = Ifges(7,5) * t27 + Ifges(7,6) * t29 + Ifges(7,3) * t118;
t246 = t183 * t289;
t240 = t182 * t303;
t239 = t184 * t303;
t236 = t301 * t124;
t235 = t301 * t146;
t230 = t266 / 0.2e1;
t229 = t28 * mrSges(7,1) - t26 * mrSges(7,2);
t224 = t259 / 0.2e1;
t220 = t121 * t247;
t211 = mrSges(6,3) * t219;
t210 = mrSges(5,2) * t219;
t209 = mrSges(5,1) * t218;
t207 = Ifges(6,1) * t184 - t294;
t206 = -Ifges(6,2) * t182 + t293;
t17 = t181 * t42 + t183 * t48;
t50 = -mrSges(6,2) * t118 + mrSges(6,3) * t195;
t205 = t184 * t50 - t274;
t204 = t158 * t218;
t201 = -t269 + t281;
t10 = -pkin(10) * t270 + pkin(5) * t118 + (-t107 + (-pkin(10) * t146 - t101) * t182) * qJD(5) + t223;
t11 = pkin(10) * t195 + t14;
t3 = qJD(6) * t16 + t10 * t181 + t11 * t183;
t4 = -qJD(6) * t17 + t10 * t183 - t11 * t181;
t199 = t4 * mrSges(7,1) - t3 * mrSges(7,2) + t248;
t198 = t179 * t211;
t197 = t180 * t211;
t196 = -t183 * t200 * t254 + (-pkin(5) * t286 + t148 * t254) * t181 + t326;
t193 = t320 * t303;
t191 = -t277 + (-t182 * t61 - t184 * t60) * qJD(5);
t122 = Ifges(7,4) * t148 - Ifges(7,2) * t200;
t123 = Ifges(7,1) * t148 - Ifges(7,4) * t200;
t150 = t206 * qJD(5);
t151 = t207 * qJD(5);
t159 = Ifges(6,2) * t184 + t294;
t160 = Ifges(6,1) * t182 + t293;
t68 = -Ifges(7,4) * t115 - Ifges(7,2) * t116;
t69 = -Ifges(7,1) * t115 - Ifges(7,4) * t116;
t190 = -t115 * t123 - t116 * t122 + t148 * t69 + t184 * t150 + t182 * t151 - t159 * t260 + t160 * t259 - t200 * t68;
t188 = t191 + t284;
t187 = Ifges(6,5) * t194 + Ifges(6,6) * t195 + Ifges(6,3) * t118;
t33 = Ifges(6,4) * t194 + Ifges(6,2) * t195 + t118 * Ifges(6,6);
t34 = Ifges(6,1) * t194 + Ifges(6,4) * t195 + t118 * Ifges(6,5);
t37 = -pkin(5) * t195 + t66;
t94 = t200 * t146;
t46 = Ifges(7,4) * t94 + Ifges(7,2) * t92 + Ifges(7,6) * t147;
t47 = Ifges(7,1) * t94 + Ifges(7,4) * t92 + Ifges(7,5) * t147;
t7 = Ifges(7,4) * t27 + Ifges(7,2) * t29 + Ifges(7,6) * t118;
t8 = Ifges(7,1) * t27 + Ifges(7,4) * t29 + Ifges(7,5) * t118;
t81 = t147 * Ifges(6,6) - t146 * t206;
t82 = t147 * Ifges(6,5) - t146 * t207;
t83 = -pkin(5) * t266 + t124;
t186 = (qJD(5) * t230 + t270 / 0.2e1) * t160 + (Ifges(6,5) * t182 + Ifges(7,5) * t148 + Ifges(6,6) * t184 - Ifges(7,6) * t200) * t118 / 0.2e1 - t200 * t7 / 0.2e1 + t327 * t66 + (-Ifges(6,6) * t260 + t326) * t147 / 0.2e1 + t150 * t230 + t82 * t224 + (t146 * t224 - t271 / 0.2e1) * t159 - t3 * t283 + mrSges(6,3) * t284 - t17 * t286 - t65 * mrSges(5,2) + t83 * t67 + t92 * t68 / 0.2e1 + t94 * t69 / 0.2e1 - t115 * t47 / 0.2e1 - t116 * t46 / 0.2e1 + Ifges(5,5) * t117 - Ifges(5,6) * t118 + t37 * t121 + t29 * t122 / 0.2e1 + t27 * t123 / 0.2e1 + t148 * t8 / 0.2e1 + t124 * t149 + t182 * t34 / 0.2e1 + t184 * t33 / 0.2e1 - t81 * t260 / 0.2e1 - t151 * t265 / 0.2e1;
t140 = (-mrSges(7,1) * t181 - mrSges(7,2) * t183) * t290;
t100 = t208 * t146;
t80 = mrSges(7,1) * t147 - mrSges(7,3) * t94;
t79 = -mrSges(7,2) * t147 + mrSges(7,3) * t92;
t51 = -mrSges(7,1) * t92 + mrSges(7,2) * t94;
t39 = -mrSges(6,1) * t195 + mrSges(6,2) * t194;
t19 = -mrSges(7,2) * t118 + mrSges(7,3) * t29;
t18 = mrSges(7,1) * t118 - mrSges(7,3) * t27;
t9 = -mrSges(7,1) * t29 + mrSges(7,2) * t27;
t1 = [0.2e1 * (mrSges(5,2) * t170 - Ifges(5,1) * t146) * t117 + 0.2e1 * (-t117 * t147 + t118 * t146) * Ifges(5,4) + (mrSges(4,1) * t302 + mrSges(4,2) * t304 + mrSges(3,3)) * t310 + (0.2e1 * Ifges(4,4) * t302 + t304 * t324) * t227 + (-0.2e1 * Ifges(4,4) * t304 + t302 * t324) * t228 + t194 * t82 + t195 * t81 + (t281 + t280) * t311 + t125 * t256 + t147 * t187 + ((m(4) + m(3)) * t310 + 0.2e1 * (mrSges(4,1) * t304 - mrSges(4,2) * t302) * qJD(3)) * qJ(2) + (Ifges(7,5) * t94 + Ifges(7,6) * t92 + 0.2e1 * t170 * mrSges(5,1) + (-Ifges(6,5) * t184 + t292) * t146 + ((2 * Ifges(5,2)) + Ifges(7,3) + Ifges(6,3)) * t147) * t118 + t33 * t266 + (t16 * t4 + t17 * t3 + t37 * t83) * t312 + 0.2e1 * t16 * t18 + 0.2e1 * t17 * t19 + t29 * t46 + t27 * t47 + 0.2e1 * t37 * t51 + 0.2e1 * t60 * t49 + 0.2e1 * t61 * t50 + 0.2e1 * t3 * t79 + 0.2e1 * t4 * t80 + 0.2e1 * t83 * t9 + t92 * t7 + t94 * t8 - 0.2e1 * t66 * t100 + 0.2e1 * t14 * t105 + 0.2e1 * t15 * t106 + 0.2e1 * t124 * t39 + t147 * t248 + 0.2e1 * t163 * (mrSges(5,1) * t147 - mrSges(5,2) * t146) - t34 * t265 + 0.2e1 * mrSges(5,3) * t269 + 0.2e1 * m(5) * (t125 * t65 + t163 * t170 + t285) + 0.2e1 * m(6) * (t14 * t61 + t15 * t60 + t285); -t91 * t18 - t93 * t19 + t26 * t79 + t28 * t80 + (t39 + t9) * t146 + (t105 * t184 - t106 * t182) * t118 + (t100 - t51 + 0.2e1 * t282) * t117 + (t256 + (-t182 * t105 - t184 * t106) * qJD(5) + t205) * t147 + m(7) * (-t117 * t83 + t146 * t37 + t16 * t28 + t17 * t26 - t3 * t93 - t4 * t91) + m(6) * ((-t182 * t60 + t184 * t61) * t118 + t188 * t147 + t201) + m(5) * (t118 * t125 + t201 + t280); 0.2e1 * m(7) * (-t26 * t93 - t28 * t91 - t86) + 0.2e1 * m(6) * (t147 * t221 - t86) + 0.2e1 * m(5) * (t118 * t147 - t86); -t106 * t203 + m(7) * (t154 * t83 + t157 * t37 + t16 * t45 + t17 * t44 + t3 * t98 + t4 * t97) + t105 * t202 + m(6) * (t173 * t66 + (t239 * t61 - t240 * t60 + t236) * t291) + t186 + t50 * t263 + (t301 * t65 - t303 * t66 + (t125 * t303 + t236) * qJD(4)) * t309 - Ifges(4,5) * t227 - Ifges(4,6) * t228 + t44 * t79 + t45 * t80 + t97 * t18 + t98 * t19 + t154 * t51 + t157 * t9 + t173 * t39 - t4 * t279 + mrSges(7,3) * t288 + (-t100 - t282) * t218 + t319 * t185 + (m(6) * t188 - t274 + t318) * t172 + t315 * mrSges(6,3) + (pkin(3) * t317 - t147 * t219) * mrSges(5,3); m(6) * (-t173 * t117 + t172 * t221 + (t147 * t193 + t235) * t291) + ((t147 * t303 + t235) * qJD(4) - t317) * t309 + m(7) * (-t117 * t157 + t146 * t154 + t26 * t98 + t28 * t97 - t44 * t93 - t45 * t91) + t314 + t319; t190 + 0.2e1 * (t172 * t193 + t173 * t301) * t255 + (t154 * t157 + t44 * t98 + t45 * t97) * t312 + 0.2e1 * t197 + 0.2e1 * t198 - 0.2e1 * t210 - 0.2e1 * t209 - 0.2e1 * t252 - 0.2e1 * t253 - 0.2e1 * t250 + 0.2e1 * t251 + 0.2e1 * t204 + 0.2e1 * t264 + 0.2e1 * t276 + 0.2e1 * t262; (m(6) * (t284 + t315) + t205 + t318) * pkin(9) + t191 * mrSges(6,3) + m(7) * (t129 * t4 + t130 * t3 + t16 * t76 + t17 * t75 + t174 * t37 + t247 * t83) + t186 + (-t148 * t4 + t288) * mrSges(7,3) + t51 * t247 + t75 * t79 + t76 * t80 + t129 * t18 + t130 * t19 + t174 * t9 + (-m(6) * t66 - t39) * pkin(4); m(7) * (pkin(5) * t234 - t117 * t174 + t129 * t28 + t130 * t26 - t75 * t93 - t76 * t91) + m(6) * (pkin(4) * t117 + pkin(9) * t221) + t314; t190 + m(7) * (t129 * t45 + t130 * t44 + t154 * t174 + t157 * t247 + t75 * t98 + t76 * t97) + (-pkin(4) * t301 + pkin(9) * t193) * t255 + t197 + t198 - t210 - t209 - t252 - t253 - t250 + t251 + t204 + t220 - t300 + t264 + t276 + t262 + t275 + t313; 0.2e1 * t220 + 0.2e1 * t275 - 0.2e1 * t300 + (t129 * t76 + t130 * t75 + t174 * t247) * t312 + 0.2e1 * t313 + t190; t15 * mrSges(6,1) - t14 * mrSges(6,2) + (m(7) * (-t16 * t258 + t17 * t257 + t181 * t3 + t183 * t4) + t79 * t257 + t181 * t19 - t80 * t258 + t183 * t18) * pkin(5) + t187 + t199; (-t118 * t184 + t147 * t260) * mrSges(6,2) + (-t118 * t182 - t147 * t259) * mrSges(6,1) + m(7) * (t181 * t26 + t183 * t28 + (t181 * t91 - t183 * t93) * qJD(6)) * pkin(5) + t229; (-mrSges(6,1) * t240 - mrSges(6,2) * t239) * t291 + (t158 * t172 - t292) * qJD(5) + (m(7) * (t181 * t44 + t183 * t45 + t257 * t98 - t258 * t97) + t246) * pkin(5) + t196 + t322; (pkin(9) * t158 - t292) * qJD(5) + (m(7) * (-t129 * t258 + t130 * t257 + t181 * t75 + t183 * t76) + t246) * pkin(5) + t196 + t321; 0.2e1 * t140; t199; t229; t261 + t322; t261 + t321; t140; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
