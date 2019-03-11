% Calculate time derivative of joint inertia matrix for
% S6RRRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 01:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRP2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP2_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP2_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP2_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP2_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP2_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP2_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:00:43
% EndTime: 2019-03-10 01:00:55
% DurationCPUTime: 5.32s
% Computational Cost: add. (10585->417), mult. (23174->598), div. (0->0), fcn. (22510->8), ass. (0->211)
t330 = Ifges(6,1) + Ifges(7,1);
t195 = sin(qJ(3));
t196 = sin(qJ(2));
t199 = cos(qJ(3));
t200 = cos(qJ(2));
t255 = t199 * t200;
t151 = -t195 * t196 + t255;
t152 = t195 * t200 + t199 * t196;
t194 = sin(qJ(4));
t198 = cos(qJ(4));
t215 = t198 * t151 - t152 * t194;
t287 = Ifges(7,4) + Ifges(6,5);
t329 = t215 * t287;
t197 = cos(qJ(5));
t193 = sin(qJ(5));
t278 = Ifges(7,5) * t193;
t280 = Ifges(6,4) * t193;
t328 = t330 * t197 + t278 - t280;
t277 = Ifges(7,5) * t197;
t279 = Ifges(6,4) * t197;
t327 = t193 * t330 - t277 + t279;
t326 = Ifges(6,6) * t197 + t287 * t193;
t325 = Ifges(7,2) + Ifges(6,3);
t115 = t151 * t194 + t152 * t198;
t246 = qJD(5) * t193;
t314 = qJD(2) + qJD(3);
t119 = t314 * t151;
t120 = t314 * t152;
t60 = qJD(4) * t215 + t119 * t198 - t120 * t194;
t269 = t197 * t60;
t212 = t115 * t246 - t269;
t245 = qJD(5) * t197;
t234 = t115 * t245;
t271 = t193 * t60;
t213 = t234 + t271;
t61 = qJD(4) * t115 + t119 * t194 + t198 * t120;
t324 = t287 * t61 + (-Ifges(6,4) + Ifges(7,5)) * t213 - t330 * t212;
t285 = t328 * t115 - t329;
t323 = t328 * qJD(5);
t322 = Ifges(7,6) * t246 + t287 * t245;
t192 = t197 ^ 2;
t247 = qJD(4) * t198;
t240 = pkin(3) * t247;
t228 = t192 * t240;
t191 = t193 ^ 2;
t229 = t191 * t240;
t320 = t228 + t229;
t182 = pkin(2) * t199 + pkin(3);
t248 = qJD(4) * t194;
t260 = t194 * t195;
t110 = t182 * t247 + (-t195 * t248 + (t198 * t199 - t260) * qJD(3)) * pkin(2);
t266 = t110 * t192;
t267 = t110 * t191;
t319 = t266 + t267;
t301 = -pkin(8) - pkin(7);
t173 = t301 * t196;
t174 = t301 * t200;
t153 = t195 * t173;
t208 = (t255 * t301 - t153) * qJD(2);
t249 = qJD(3) * t199;
t250 = qJD(3) * t195;
t294 = t119 * pkin(9);
t126 = t199 * t173 + t174 * t195;
t104 = -pkin(9) * t152 + t126;
t127 = -t199 * t174 + t153;
t105 = pkin(9) * t151 + t127;
t316 = t198 * t104 - t105 * t194;
t79 = qJD(2) * t152 * t301 + t173 * t249 + t174 * t250;
t66 = -pkin(9) * t120 + t79;
t23 = t198 * t66 + t316 * qJD(4) + (-t173 * t250 + t174 * t249 + t208 - t294) * t194;
t108 = qJD(2) * t196 * pkin(2) + pkin(3) * t120;
t30 = pkin(4) * t61 - pkin(10) * t60 + t108;
t183 = -pkin(2) * t200 - pkin(1);
t130 = -t151 * pkin(3) + t183;
t68 = -pkin(4) * t215 - t115 * pkin(10) + t130;
t70 = t104 * t194 + t105 * t198;
t6 = t193 * t30 + t197 * t23 + t68 * t245 - t246 * t70;
t2 = qJ(6) * t61 - qJD(6) * t215 + t6;
t34 = -t193 * t70 + t197 * t68;
t33 = pkin(5) * t215 - t34;
t284 = t193 * t68 + t197 * t70;
t315 = qJD(5) * t284;
t7 = -t193 * t23 + t197 * t30 - t315;
t4 = -pkin(5) * t61 - t7;
t318 = t193 * t4 + t197 * t2 + t33 * t245;
t241 = pkin(3) * t248;
t19 = mrSges(6,1) * t213 - mrSges(6,2) * t212;
t80 = -qJD(3) * t127 + t208;
t24 = qJD(4) * t70 + t194 * t66 - t198 * (t80 - t294);
t317 = m(6) * t24 + t19;
t134 = pkin(5) * t246 - qJ(6) * t245 - qJD(6) * t193;
t163 = -t197 * mrSges(7,1) - t193 * mrSges(7,3);
t122 = t134 * t163;
t224 = t193 * mrSges(7,1) - t197 * mrSges(7,3);
t156 = t224 * qJD(5);
t219 = pkin(5) * t197 + qJ(6) * t193;
t162 = -pkin(4) - t219;
t128 = t162 * t156;
t225 = t193 * mrSges(6,1) + t197 * mrSges(6,2);
t157 = t225 * qJD(5);
t295 = pkin(4) * t157;
t313 = t122 + t128 - t295;
t312 = t80 * mrSges(4,1) - t79 * mrSges(4,2) + Ifges(4,5) * t119 - Ifges(4,6) * t120;
t129 = t134 + t241;
t113 = t129 * t163;
t296 = pkin(3) * t198;
t144 = t162 - t296;
t121 = t144 * t156;
t181 = -pkin(4) - t296;
t131 = t181 * t157;
t164 = -t197 * mrSges(6,1) + t193 * mrSges(6,2);
t143 = t164 * t241;
t169 = mrSges(6,3) * t229;
t170 = mrSges(7,2) * t229;
t171 = mrSges(6,3) * t228;
t172 = mrSges(7,2) * t228;
t311 = t113 + t121 + t131 + t143 + t169 + t170 + t171 + t172;
t310 = 2 * m(5);
t309 = 0.2e1 * m(6);
t308 = 2 * m(7);
t307 = -2 * Ifges(5,4);
t306 = 0.2e1 * t24;
t305 = -0.2e1 * t316;
t304 = 0.2e1 * t108;
t303 = 0.2e1 * t183;
t302 = t61 / 0.2e1;
t166 = Ifges(6,2) * t197 + t280;
t298 = -t166 / 0.2e1;
t291 = t197 * t6;
t290 = t24 * t316;
t289 = t7 * t193;
t220 = Ifges(7,3) * t193 + t277;
t48 = -Ifges(7,6) * t215 + t115 * t220;
t221 = -Ifges(6,2) * t193 + t279;
t276 = Ifges(6,6) * t215;
t49 = t115 * t221 - t276;
t286 = t48 - t49;
t263 = t115 * t193;
t73 = mrSges(6,2) * t215 - mrSges(6,3) * t263;
t76 = -mrSges(7,2) * t263 - mrSges(7,3) * t215;
t283 = t73 + t76;
t262 = t115 * t197;
t74 = -mrSges(6,1) * t215 - mrSges(6,3) * t262;
t75 = mrSges(7,1) * t215 + mrSges(7,2) * t262;
t282 = -t74 + t75;
t259 = t195 * t198;
t136 = pkin(2) * t259 + t194 * t182;
t133 = pkin(10) + t136;
t281 = t319 * t133;
t275 = Ifges(6,6) * t193;
t273 = t110 * mrSges(5,2);
t111 = t182 * t248 + (t195 * t247 + (t194 * t199 + t259) * qJD(3)) * pkin(2);
t272 = t111 * t316;
t270 = t194 * t316;
t265 = t110 * t193;
t264 = t110 * t197;
t261 = t193 * t198;
t257 = t197 * t198;
t254 = t319 * pkin(10);
t180 = pkin(3) * t194 + pkin(10);
t253 = t320 * t180;
t252 = t320 * pkin(10);
t244 = qJD(6) * t197;
t243 = 0.2e1 * t200;
t242 = m(7) * t244;
t32 = -qJ(6) * t215 + t284;
t239 = t32 * t246;
t232 = -t246 / 0.2e1;
t135 = -pkin(2) * t260 + t182 * t198;
t230 = t133 * t320 + t180 * t319;
t227 = mrSges(7,2) * t244 + t322;
t218 = pkin(5) * t193 - qJ(6) * t197;
t214 = t213 * Ifges(7,6) + t269 * t287 + t325 * t61;
t211 = (-mrSges(4,1) * t195 - mrSges(4,2) * t199) * qJD(3) * pkin(2);
t210 = (-mrSges(5,1) * t194 - mrSges(5,2) * t198) * qJD(4) * pkin(3);
t26 = -t61 * mrSges(7,1) - mrSges(7,2) * t212;
t209 = -mrSges(7,2) * t219 - t275;
t158 = t220 * qJD(5);
t159 = t221 * qJD(5);
t165 = -Ifges(7,3) * t197 + t278;
t207 = (t165 - t166) * t246 + t327 * t245 + (-t158 + t159) * t197 + t323 * t193;
t206 = -m(7) * t219 + t163 + t164;
t205 = -m(7) * t218 - t224 - t225;
t100 = mrSges(6,3) * t266;
t101 = mrSges(7,2) * t266;
t123 = -t135 + t162;
t106 = t123 * t156;
t109 = t111 * mrSges(5,1);
t132 = -pkin(4) - t135;
t118 = t132 * t157;
t83 = t111 + t134;
t77 = t83 * t163;
t94 = t111 * t164;
t98 = mrSges(6,3) * t267;
t99 = mrSges(7,2) * t267;
t204 = t100 + t101 + t106 - t109 + t118 + t207 + t77 + t94 + t98 + t99;
t14 = -Ifges(7,5) * t212 + Ifges(7,6) * t61 + Ifges(7,3) * t213;
t15 = -Ifges(6,4) * t212 - Ifges(6,2) * t213 + Ifges(6,6) * t61;
t37 = t115 * t218 - t316;
t9 = t218 * t60 + (qJD(5) * t219 - t244) * t115 + t24;
t203 = t37 * t156 - t316 * t157 + t9 * t163 + Ifges(5,5) * t60 - Ifges(5,6) * t61 - t23 * mrSges(5,2) + mrSges(6,3) * t291 + t234 * t298 + t48 * t246 / 0.2e1 + t49 * t232 + t326 * t302 - (-Ifges(6,6) * t246 + t322) * t215 / 0.2e1 + t324 * t193 / 0.2e1 + (t298 + t165 / 0.2e1) * t271 + (-t159 / 0.2e1 + t158 / 0.2e1) * t263 + (t164 - mrSges(5,1)) * t24 + t323 * t262 / 0.2e1 + (t15 / 0.2e1 - t14 / 0.2e1 - Ifges(7,6) * t302) * t197 + (t115 * t165 + t285) * t245 / 0.2e1 + t318 * mrSges(7,2) + t327 * (t115 * t232 + t269 / 0.2e1);
t25 = mrSges(6,1) * t61 + mrSges(6,3) * t212;
t27 = -mrSges(6,2) * t61 - mrSges(6,3) * t213;
t28 = -mrSges(7,2) * t213 + mrSges(7,3) * t61;
t202 = (t27 + t28) * t197 + (-t25 + t26) * t193 + (-t193 * t283 + t197 * t282) * qJD(5) + m(7) * (-t239 + t318) + m(6) * (-t245 * t34 - t246 * t284 - t289 + t291);
t201 = t203 - mrSges(7,2) * t239 + (-t289 + (-t193 * t284 - t197 * t34) * qJD(5)) * mrSges(6,3);
t185 = mrSges(7,2) * t245;
t72 = t225 * t115;
t71 = t224 * t115;
t18 = mrSges(7,1) * t213 + mrSges(7,3) * t212;
t1 = [(t108 * t130 + t23 * t70 - t290) * t310 + (t284 * t6 + t34 * t7 - t290) * t309 + (mrSges(4,1) * t120 + mrSges(4,2) * t119) * t303 + 0.2e1 * (t119 * t151 - t120 * t152) * Ifges(4,4) + 0.2e1 * (-t119 * t126 - t120 * t127 + t151 * t79 - t152 * t80) * mrSges(4,3) + 0.2e1 * t2 * t76 + 0.2e1 * t9 * t71 + 0.2e1 * t6 * t73 + 0.2e1 * t7 * t74 + 0.2e1 * t4 * t75 + 0.2e1 * t32 * t28 + 0.2e1 * t33 * t26 + 0.2e1 * t34 * t25 + 0.2e1 * t37 * t18 + (mrSges(5,2) * t304 + mrSges(5,3) * t306 + 0.2e1 * Ifges(5,1) * t60 + t324 * t197 + (t14 - t15) * t193 + (t307 + t287 * t197 + (-Ifges(6,6) + Ifges(7,6)) * t193) * t61 + ((t276 + t286) * t197 + (-t285 + t329) * t193) * qJD(5)) * t115 + (0.2e1 * t130 * mrSges(5,2) + mrSges(5,3) * t305 + t286 * t193 + t285 * t197) * t60 + ((-mrSges(3,2) * pkin(1) + Ifges(3,4) * t200) * t243 + (0.2e1 * pkin(2) * (-mrSges(4,1) * t151 + mrSges(4,2) * t152) - 0.2e1 * pkin(1) * mrSges(3,1) + m(4) * pkin(2) * t303 - 0.2e1 * Ifges(3,4) * t196 + (-Ifges(3,2) + Ifges(3,1)) * t243) * t196) * qJD(2) + t19 * t305 + t72 * t306 + (t2 * t32 + t33 * t4 + t37 * t9) * t308 - (mrSges(5,1) * t304 - 0.2e1 * mrSges(5,3) * t23 + (t307 - t275) * t60 + ((2 * Ifges(5,2)) + t325) * t61 + t214) * t215 + 0.2e1 * (mrSges(5,1) * t130 - mrSges(5,3) * t70) * t61 + 0.2e1 * t284 * t27 + 0.2e1 * m(4) * (t126 * t80 + t127 * t79) + 0.2e1 * t119 * t152 * Ifges(4,1) - 0.2e1 * t120 * Ifges(4,2) * t151; t203 + t132 * t19 + t123 * t18 + t111 * t72 + t83 * t71 + (m(4) * (t195 * t79 + t199 * t80 + (-t126 * t195 + t127 * t199) * qJD(3)) + (-t199 * t119 - t195 * t120 + (t151 * t199 + t152 * t195) * qJD(3)) * mrSges(4,3)) * pkin(2) + (t110 * t215 + t111 * t115 - t135 * t60 - t136 * t61) * mrSges(5,3) + m(6) * (t132 * t24 + t264 * t284 - t265 * t34 - t272) + m(7) * (t123 * t9 + t264 * t32 + t265 * t33 + t37 * t83) + t202 * t133 + (Ifges(3,5) * t200 - Ifges(3,6) * t196 + (-mrSges(3,1) * t200 + mrSges(3,2) * t196) * pkin(7)) * qJD(2) + (-t32 * qJD(5) * mrSges(7,2) + t282 * t110 + (-t7 - t315) * mrSges(6,3)) * t193 + (-t34 * qJD(5) * mrSges(6,3) + t110 * t283) * t197 + m(5) * (t110 * t70 - t135 * t24 + t136 * t23 - t272) + t312; 0.2e1 * t118 - 0.2e1 * t109 + 0.2e1 * t106 + 0.2e1 * t101 + 0.2e1 * t99 + 0.2e1 * t100 + 0.2e1 * t98 + 0.2e1 * t94 + 0.2e1 * t77 + t207 - 0.2e1 * t273 + (t111 * t132 + t281) * t309 + (t110 * t136 - t111 * t135) * t310 + (t123 * t83 + t281) * t308 + 0.2e1 * t211; t201 + m(7) * (t129 * t37 + t144 * t9) + t144 * t18 + t129 * t71 + t202 * t180 + (m(5) * (t194 * t23 - t198 * t24) + (-t194 * t61 - t198 * t60) * mrSges(5,3) + ((t115 * mrSges(5,3) + t72) * t194 + (mrSges(5,3) * t215 + t193 * t282 + t197 * t283) * t198 + m(6) * (t257 * t284 - t261 * t34 - t270) + m(7) * (t257 * t32 + t261 * t33) + m(5) * (t198 * t70 - t270)) * qJD(4)) * pkin(3) + t317 * t181 + t312; t204 + t211 - mrSges(5,1) * t241 + m(5) * (t110 * t194 - t111 * t198 - t135 * t248 + t136 * t247) * pkin(3) + m(7) * (t123 * t129 + t144 * t83 + t230) + (-t110 - t240) * mrSges(5,2) + (t111 * t181 + t132 * t241 + t230) * m(6) + t311; 0.2e1 * t143 + 0.2e1 * t131 + 0.2e1 * t121 + 0.2e1 * t113 + 0.2e1 * t169 + 0.2e1 * t170 + 0.2e1 * t171 + 0.2e1 * t172 + t207 + (t129 * t144 + t253) * t308 + (t181 * t241 + t253) * t309 + 0.2e1 * t210; t201 + t162 * t18 + t134 * t71 + t202 * pkin(10) + m(7) * (t134 * t37 + t162 * t9) - t317 * pkin(4); t204 - t273 + m(7) * (t123 * t134 + t162 * t83 + t254) + m(6) * (-pkin(4) * t111 + t254) + t313; t210 + t207 + m(7) * (t129 * t162 + t134 * t144 + t252) + m(6) * (-pkin(4) * t241 + t252) + t311 + t313; t134 * t162 * t308 + 0.2e1 * t122 + 0.2e1 * t128 + t207 - 0.2e1 * t295; qJD(6) * t76 + qJ(6) * t28 + t2 * mrSges(7,3) - t6 * mrSges(6,2) - t4 * mrSges(7,1) - pkin(5) * t26 + m(7) * (-pkin(5) * t4 + qJ(6) * t2 + qJD(6) * t32) - Ifges(6,6) * t271 + t7 * mrSges(6,1) - t326 * t115 * qJD(5) + t214; t133 * t242 + t205 * t110 + (t133 * t206 + t209) * qJD(5) + t227; t180 * t242 + t205 * t240 + (t180 * t206 + t209) * qJD(5) + t227; t209 * qJD(5) + (qJD(5) * t206 + t242) * pkin(10) + t227; 0.2e1 * (m(7) * qJ(6) + mrSges(7,3)) * qJD(6); m(7) * t4 + t26; t185 + (t133 * t245 + t265) * m(7); t185 + (t180 * t245 + t193 * t240) * m(7); m(7) * pkin(10) * t245 + t185; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
