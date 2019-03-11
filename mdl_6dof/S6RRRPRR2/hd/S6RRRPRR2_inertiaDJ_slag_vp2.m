% Calculate time derivative of joint inertia matrix for
% S6RRRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 18:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR2_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR2_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR2_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR2_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR2_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR2_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:06:30
% EndTime: 2019-03-09 18:06:39
% DurationCPUTime: 4.27s
% Computational Cost: add. (12881->451), mult. (27783->670), div. (0->0), fcn. (27976->10), ass. (0->210)
t193 = sin(qJ(5));
t197 = cos(qJ(5));
t174 = -mrSges(6,1) * t197 + mrSges(6,2) * t193;
t314 = t174 - mrSges(5,1);
t241 = qJD(5) * t197;
t192 = sin(qJ(6));
t196 = cos(qJ(6));
t210 = t192 * t193 - t196 * t197;
t310 = qJD(5) + qJD(6);
t129 = t310 * t210;
t165 = t192 * t197 + t193 * t196;
t130 = t310 * t165;
t246 = -Ifges(7,5) * t129 - Ifges(7,6) * t130;
t316 = Ifges(6,5) * t241 + t246;
t194 = sin(qJ(3));
t195 = sin(qJ(2));
t198 = cos(qJ(3));
t199 = cos(qJ(2));
t247 = t198 * t199;
t164 = -t194 * t195 + t247;
t166 = t194 * t199 + t198 * t195;
t191 = sin(pkin(11));
t263 = cos(pkin(11));
t123 = t191 * t164 + t263 * t166;
t311 = qJD(2) + qJD(3);
t131 = t311 * t164;
t132 = t311 * t166;
t97 = t263 * t131 - t191 * t132;
t206 = t123 * t241 + t193 * t97;
t264 = t197 * t97;
t96 = t131 * t191 + t263 * t132;
t315 = Ifges(6,5) * t264 + Ifges(6,3) * t96;
t184 = pkin(2) * t198 + pkin(3);
t218 = t263 * t194;
t150 = pkin(2) * t218 + t191 * t184;
t146 = pkin(9) + t150;
t285 = -pkin(10) - t146;
t136 = t285 * t193;
t188 = t197 * pkin(10);
t137 = t146 * t197 + t188;
t102 = t136 * t196 - t137 * t192;
t220 = qJD(5) * t285;
t250 = t191 * t194;
t280 = pkin(2) * qJD(3);
t148 = (t263 * t198 - t250) * t280;
t253 = t148 * t197;
t108 = t193 * t220 + t253;
t254 = t148 * t193;
t109 = t197 * t220 - t254;
t53 = t102 * qJD(6) + t108 * t196 + t109 * t192;
t103 = t136 * t192 + t137 * t196;
t54 = -t103 * qJD(6) - t108 * t192 + t109 * t196;
t313 = t54 * mrSges(7,1) - t53 * mrSges(7,2);
t292 = pkin(3) * t191;
t182 = pkin(9) + t292;
t284 = -pkin(10) - t182;
t158 = t284 * t193;
t252 = t182 * t197;
t159 = t188 + t252;
t118 = t158 * t196 - t159 * t192;
t219 = qJD(5) * t284;
t151 = t193 * t219;
t152 = t197 * t219;
t80 = t118 * qJD(6) + t151 * t196 + t152 * t192;
t119 = t158 * t192 + t159 * t196;
t81 = -t119 * qJD(6) - t151 * t192 + t152 * t196;
t312 = t81 * mrSges(7,1) - t80 * mrSges(7,2);
t84 = t210 * t123;
t217 = (t193 ^ 2 + t197 ^ 2) * t148;
t295 = -pkin(8) - pkin(7);
t178 = t295 * t195;
t167 = t194 * t178;
t179 = t295 * t199;
t139 = -t198 * t179 + t167;
t215 = mrSges(6,1) * t193 + mrSges(6,2) * t197;
t170 = t215 * qJD(5);
t242 = qJD(5) * t193;
t201 = (t295 * t247 - t167) * qJD(2);
t209 = -t131 * qJ(4) - t166 * qJD(4);
t243 = qJD(3) * t198;
t244 = qJD(3) * t194;
t106 = qJD(2) * t166 * t295 + t178 * t243 + t179 * t244;
t66 = -qJ(4) * t132 + qJD(4) * t164 + t106;
t41 = t263 * t66 + (-t178 * t244 + t179 * t243 + t201 + t209) * t191;
t117 = qJD(2) * t195 * pkin(2) + pkin(3) * t132;
t51 = pkin(4) * t96 - pkin(9) * t97 + t117;
t221 = -t193 * t41 + t197 * t51;
t116 = qJ(4) * t164 + t139;
t138 = t198 * t178 + t179 * t194;
t204 = -qJ(4) * t166 + t138;
t76 = t263 * t116 + t191 * t204;
t73 = t197 * t76;
t122 = -t263 * t164 + t166 * t191;
t185 = -pkin(2) * t199 - pkin(1);
t142 = -t164 * pkin(3) + t185;
t74 = t122 * pkin(4) - t123 * pkin(9) + t142;
t46 = t193 * t74 + t73;
t14 = -t46 * qJD(5) + t221;
t276 = t14 * t193;
t45 = -t193 * t76 + t197 * t74;
t309 = -t45 * t241 - t46 * t242 - t276;
t226 = t123 * t242;
t205 = t226 - t264;
t47 = mrSges(6,1) * t96 + t205 * mrSges(6,3);
t258 = t123 * t193;
t88 = -mrSges(6,2) * t122 - mrSges(6,3) * t258;
t257 = t123 * t197;
t89 = mrSges(6,1) * t122 - mrSges(6,3) * t257;
t308 = -t193 * t47 - t89 * t241 - t88 * t242;
t307 = t118 * t129 - t119 * t130 - t165 * t81 - t210 * t80;
t306 = t102 * t129 - t103 * t130 - t165 * t54 - t210 * t53;
t26 = -t130 * t123 - t210 * t97;
t27 = -t165 * t97 + t310 * t84;
t10 = Ifges(7,1) * t26 + Ifges(7,4) * t27 + Ifges(7,5) * t96;
t100 = -Ifges(7,1) * t129 - Ifges(7,4) * t130;
t107 = -qJD(3) * t139 + t201;
t133 = mrSges(7,1) * t210 + mrSges(7,2) * t165;
t134 = Ifges(7,4) * t165 - Ifges(7,2) * t210;
t135 = Ifges(7,1) * t165 - Ifges(7,4) * t210;
t32 = pkin(5) * t122 - pkin(10) * t257 + t45;
t37 = -pkin(10) * t258 + t46;
t16 = -t192 * t37 + t196 * t32;
t17 = t192 * t32 + t196 * t37;
t282 = Ifges(6,4) * t197;
t213 = -Ifges(6,2) * t193 + t282;
t171 = t213 * qJD(5);
t283 = Ifges(6,4) * t193;
t214 = Ifges(6,1) * t197 - t283;
t172 = t214 * qJD(5);
t175 = Ifges(6,2) * t197 + t283;
t176 = Ifges(6,1) * t193 + t282;
t212 = Ifges(6,5) * t193 + Ifges(6,6) * t197;
t223 = -t242 / 0.2e1;
t40 = t191 * t66 - t263 * (t107 + t209);
t23 = t206 * pkin(5) + t40;
t277 = t130 * mrSges(7,3);
t13 = t193 * t51 + t197 * t41 + t74 * t241 - t76 * t242;
t278 = t13 * t197;
t296 = t96 / 0.2e1;
t5 = -pkin(10) * t264 + pkin(5) * t96 + (-t73 + (pkin(10) * t123 - t74) * t193) * qJD(5) + t221;
t8 = -t206 * pkin(10) + t13;
t3 = t16 * qJD(6) + t192 * t5 + t196 * t8;
t30 = -t205 * Ifges(6,4) - t206 * Ifges(6,2) + Ifges(6,6) * t96;
t31 = -t205 * Ifges(6,1) - t206 * Ifges(6,4) + Ifges(6,5) * t96;
t4 = -t17 * qJD(6) - t192 * t8 + t196 * t5;
t83 = t165 * t123;
t42 = -Ifges(7,4) * t84 - Ifges(7,2) * t83 + Ifges(7,6) * t122;
t43 = -Ifges(7,1) * t84 - Ifges(7,4) * t83 + Ifges(7,5) * t122;
t75 = t116 * t191 - t263 * t204;
t58 = pkin(5) * t258 + t75;
t63 = Ifges(6,6) * t122 + t213 * t123;
t64 = Ifges(6,5) * t122 + t214 * t123;
t9 = Ifges(7,4) * t26 + Ifges(7,2) * t27 + Ifges(7,6) * t96;
t98 = t130 * mrSges(7,1) - t129 * mrSges(7,2);
t99 = -Ifges(7,4) * t129 - Ifges(7,2) * t130;
t305 = t314 * t40 - t206 * t175 / 0.2e1 + (-Ifges(6,6) * t242 + t316) * t122 / 0.2e1 + (-Ifges(7,6) * t296 - t9 / 0.2e1 - t3 * mrSges(7,3)) * t210 + (t16 * t129 - t4 * t165) * mrSges(7,3) + t212 * t296 - t17 * t277 + mrSges(6,3) * t278 + t197 * t30 / 0.2e1 + t193 * t31 / 0.2e1 + t75 * t170 - Ifges(4,6) * t132 + t23 * t133 + t27 * t134 / 0.2e1 + t26 * t135 / 0.2e1 - t129 * t43 / 0.2e1 - t130 * t42 / 0.2e1 + Ifges(4,5) * t131 - t84 * t100 / 0.2e1 - t106 * mrSges(4,2) + t107 * mrSges(4,1) - Ifges(5,6) * t96 + Ifges(5,5) * t97 + t58 * t98 - t83 * t99 / 0.2e1 + (t123 * t223 + t264 / 0.2e1) * t176 - t41 * mrSges(5,2) + (Ifges(7,5) * t296 + t10 / 0.2e1) * t165 + t64 * t241 / 0.2e1 + t172 * t257 / 0.2e1 - t171 * t258 / 0.2e1 + t63 * t223;
t304 = 2 * m(5);
t303 = 2 * m(6);
t302 = 2 * m(7);
t301 = 0.2e1 * t40;
t300 = 0.2e1 * t75;
t299 = 0.2e1 * t117;
t298 = 0.2e1 * t185;
t297 = m(5) * pkin(3);
t291 = t197 * pkin(5);
t289 = t40 * t75;
t286 = t96 * mrSges(5,3);
t281 = Ifges(6,6) * t193;
t279 = pkin(5) * qJD(6);
t149 = -pkin(2) * t250 + t263 * t184;
t145 = -pkin(4) - t149;
t141 = t145 - t291;
t275 = t141 * t98;
t147 = (t191 * t198 + t218) * t280;
t274 = t147 * t75;
t227 = t263 * pkin(3);
t183 = -t227 - pkin(4);
t173 = t183 - t291;
t267 = t173 * t98;
t234 = pkin(5) * t242;
t140 = t147 + t234;
t256 = t140 * t133;
t255 = t145 * t170;
t251 = t183 * t170;
t240 = qJD(6) * t192;
t239 = qJD(6) * t196;
t238 = 0.2e1 * mrSges(7,3);
t237 = 0.2e1 * t199;
t236 = Ifges(7,5) * t26 + Ifges(7,6) * t27 + Ifges(7,3) * t96;
t235 = mrSges(7,3) * t279;
t233 = t196 * t129 * mrSges(7,3);
t224 = t96 * mrSges(5,1) + t97 * mrSges(5,2);
t222 = -(2 * Ifges(5,4)) - t281;
t216 = t133 * t234;
t211 = -t193 * t89 + t197 * t88;
t208 = t4 * mrSges(7,1) - t3 * mrSges(7,2) + t236;
t207 = -t196 * t210 * t235 + (-pkin(5) * t277 + t165 * t235) * t192 + t316;
t203 = -t276 + (-t193 * t46 - t197 * t45) * qJD(5);
t202 = t165 * t100 - t129 * t135 - t130 * t134 + t197 * t171 + t193 * t172 - t175 * t242 + t176 * t241 - t210 * t99;
t161 = (-mrSges(7,1) * t192 - mrSges(7,2) * t196) * t279;
t87 = t215 * t123;
t62 = mrSges(7,1) * t122 + mrSges(7,3) * t84;
t61 = -mrSges(7,2) * t122 - mrSges(7,3) * t83;
t55 = mrSges(7,1) * t83 - mrSges(7,2) * t84;
t48 = -mrSges(6,2) * t96 - t206 * mrSges(6,3);
t36 = t206 * mrSges(6,1) - t205 * mrSges(6,2);
t19 = -mrSges(7,2) * t96 + mrSges(7,3) * t27;
t18 = mrSges(7,1) * t96 - mrSges(7,3) * t26;
t11 = -mrSges(7,1) * t27 + mrSges(7,2) * t26;
t1 = [0.2e1 * m(4) * (t106 * t139 + t107 * t138) + 0.2e1 * t142 * t224 + 0.2e1 * (t131 * t164 - t132 * t166) * Ifges(4,4) + 0.2e1 * (t106 * t164 - t107 * t166 - t131 * t138 - t132 * t139) * mrSges(4,3) + (mrSges(4,1) * t132 + mrSges(4,2) * t131) * t298 + (t16 * t4 + t17 * t3 + t23 * t58) * t302 + (t13 * t46 + t14 * t45 + t289) * t303 + (t117 * t142 + t41 * t76 + t289) * t304 + t36 * t300 + t87 * t301 + (mrSges(5,1) * t299 - 0.2e1 * t41 * mrSges(5,3) + t222 * t97 + ((2 * Ifges(5,2)) + Ifges(6,3) + Ifges(7,3)) * t96 + t236 + t315) * t122 + t96 * (-Ifges(7,5) * t84 - Ifges(7,6) * t83) - 0.2e1 * t132 * Ifges(4,2) * t164 + 0.2e1 * t131 * t166 * Ifges(4,1) + 0.2e1 * t14 * t89 - t83 * t9 - t84 * t10 + 0.2e1 * t13 * t88 + 0.2e1 * t3 * t61 + 0.2e1 * t4 * t62 + 0.2e1 * t23 * t55 + 0.2e1 * t58 * t11 + t27 * t42 + t26 * t43 + 0.2e1 * t45 * t47 + 0.2e1 * t46 * t48 + 0.2e1 * t17 * t19 + 0.2e1 * t16 * t18 - 0.2e1 * t76 * t286 + ((-mrSges(3,2) * pkin(1) + Ifges(3,4) * t199) * t237 + (0.2e1 * pkin(2) * (-mrSges(4,1) * t164 + mrSges(4,2) * t166) - 0.2e1 * pkin(1) * mrSges(3,1) + m(4) * pkin(2) * t298 - 0.2e1 * Ifges(3,4) * t195 + (Ifges(3,1) - Ifges(3,2)) * t237) * t195) * qJD(2) + (mrSges(5,3) * t300 - t193 * t63 + t197 * t64) * t97 + (mrSges(5,2) * t299 + mrSges(5,3) * t301 + 0.2e1 * Ifges(5,1) * t97 - t193 * t30 + t197 * t31 + (Ifges(6,5) * t197 + t222) * t96 + (-t122 * t212 - t193 * t64 - t197 * t63) * qJD(5)) * t123; t203 * mrSges(6,3) + (t197 * t48 + m(6) * (t278 + t309) + t308) * t146 + (m(4) * (t106 * t194 + t107 * t198 + (-t138 * t194 + t139 * t198) * qJD(3)) + (-t198 * t131 - t194 * t132 + (t164 * t198 + t166 * t194) * qJD(3)) * mrSges(4,3)) * pkin(2) + (-t122 * t148 + t123 * t147 - t149 * t97 - t150 * t96) * mrSges(5,3) + m(7) * (t102 * t4 + t103 * t3 + t140 * t58 + t141 * t23 + t16 * t54 + t17 * t53) + m(5) * (t148 * t76 - t149 * t40 + t150 * t41 + t274) + (Ifges(3,5) * t199 - Ifges(3,6) * t195 + (-mrSges(3,1) * t199 + mrSges(3,2) * t195) * pkin(7)) * qJD(2) + t211 * t148 + t145 * t36 + t147 * t87 + t140 * t55 + t141 * t11 + t102 * t18 + t103 * t19 + t53 * t61 + t54 * t62 + m(6) * (t145 * t40 + t46 * t253 - t45 * t254 + t274) + t305; t202 + t306 * t238 - 0.2e1 * t148 * mrSges(5,2) + (t102 * t54 + t103 * t53 + t140 * t141) * t302 + (-t147 * t149 + t148 * t150) * t304 + (t145 * t147 + t146 * t217) * t303 + 0.2e1 * t255 + 0.2e1 * t256 + 0.2e1 * t275 + 0.2e1 * t314 * t147 + 0.2e1 * mrSges(6,3) * t217 + 0.2e1 * (-mrSges(4,1) * t194 - mrSges(4,2) * t198) * t280; m(7) * (t118 * t4 + t119 * t3 + t16 * t81 + t17 * t80 + t173 * t23 + t58 * t234) + (t191 * t41 - t263 * t40) * t297 + t55 * t234 - t286 * t292 - t97 * mrSges(5,3) * t227 + t173 * t11 + t118 * t18 + t119 * t19 + t81 * t62 + t80 * t61 + t48 * t252 + (m(6) * t40 + t36) * t183 + (m(6) * (t203 + t278) + t308) * t182 + t309 * mrSges(6,3) + t305; t216 + t202 + m(7) * (t102 * t81 + t103 * t80 + t118 * t54 + t119 * t53 + t140 * t173 + t141 * t234) + t267 + t275 + t255 + t256 + t251 + (t191 * t297 - mrSges(5,2)) * t148 + (m(6) * t183 - t263 * t297 + t314) * t147 + (-mrSges(4,1) * t244 - mrSges(4,2) * t243) * pkin(2) + (t306 + t307) * mrSges(7,3) + (m(6) * t182 + mrSges(6,3)) * t217; 0.2e1 * t216 + 0.2e1 * t267 + 0.2e1 * t251 + (t118 * t81 + t119 * t80 + t173 * t234) * t302 + t307 * t238 + t202; -t129 * t61 - t130 * t62 - t210 * t18 + t165 * t19 + t193 * t48 + t197 * t47 + t211 * qJD(5) + m(7) * (-t129 * t17 - t130 * t16 + t165 * t3 - t210 * t4) + m(6) * (t13 * t193 + t14 * t197 + (-t193 * t45 + t197 * t46) * qJD(5)) + m(5) * t117 + t224; m(7) * (-t102 * t130 - t103 * t129 + t165 * t53 - t210 * t54); m(7) * (-t118 * t130 - t119 * t129 + t165 * t80 - t210 * t81); (-t129 * t165 + t130 * t210) * t302; -Ifges(6,5) * t226 + t14 * mrSges(6,1) - t13 * mrSges(6,2) - t206 * Ifges(6,6) + (-t62 * t240 + t196 * t18 + m(7) * (-t16 * t240 + t17 * t239 + t192 * t3 + t196 * t4) + t61 * t239 + t192 * t19) * pkin(5) + t208 + t315; -t215 * t148 + (t174 * t146 - t281) * qJD(5) + (m(7) * (-t102 * t240 + t103 * t239 + t192 * t53 + t196 * t54) + t233) * pkin(5) + t207 + t313; (t174 * t182 - t281) * qJD(5) + (m(7) * (-t118 * t240 + t119 * t239 + t192 * t80 + t196 * t81) + t233) * pkin(5) + t207 + t312; -t170 + m(7) * (-t129 * t192 - t130 * t196 + (t165 * t196 + t192 * t210) * qJD(6)) * pkin(5) - t98; 0.2e1 * t161; t208; t246 + t313; t246 + t312; -t98; t161; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
