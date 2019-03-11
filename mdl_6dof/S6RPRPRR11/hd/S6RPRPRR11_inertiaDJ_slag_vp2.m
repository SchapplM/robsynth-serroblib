% Calculate time derivative of joint inertia matrix for
% S6RPRPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-03-09 04:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRR11_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR11_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR11_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR11_inertiaDJ_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR11_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR11_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR11_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:11:28
% EndTime: 2019-03-09 04:11:43
% DurationCPUTime: 6.56s
% Computational Cost: add. (16238->554), mult. (46958->836), div. (0->0), fcn. (51063->14), ass. (0->243)
t188 = sin(pkin(13));
t192 = cos(pkin(13));
t241 = t188 ^ 2 + t192 ^ 2;
t310 = qJ(4) * t241;
t197 = sin(qJ(5));
t200 = cos(qJ(5));
t164 = t188 * t197 - t200 * t192;
t161 = t164 * qJD(5);
t165 = t188 * t200 + t192 * t197;
t199 = cos(qJ(6));
t196 = sin(qJ(6));
t235 = qJD(6) * t196;
t203 = t199 * t161 + t165 * t235;
t198 = sin(qJ(3));
t239 = qJD(3) * t198;
t189 = sin(pkin(12));
t191 = sin(pkin(6));
t193 = cos(pkin(12));
t195 = cos(pkin(6));
t194 = cos(pkin(7));
t201 = cos(qJ(3));
t246 = t194 * t201;
t190 = sin(pkin(7));
t249 = t190 * t201;
t309 = t191 * (-t189 * t198 + t193 * t246) + t195 * t249;
t308 = -2 * Ifges(4,4);
t307 = t196 / 0.2e1;
t283 = t199 / 0.2e1;
t135 = t309 * qJD(3);
t247 = t194 * t198;
t250 = t190 * t198;
t141 = t195 * t250 + (t189 * t201 + t193 * t247) * t191;
t248 = t191 * t193;
t158 = -t190 * t248 + t194 * t195;
t106 = -t141 * t188 + t158 * t192;
t107 = t141 * t192 + t158 * t188;
t209 = t200 * t106 - t107 * t197;
t56 = qJD(5) * t209 - t135 * t164;
t81 = t106 * t197 + t107 * t200;
t57 = qJD(5) * t81 + t135 * t165;
t26 = t57 * mrSges(6,1) + t56 * mrSges(6,2);
t257 = t135 * t192;
t258 = t135 * t188;
t93 = mrSges(5,1) * t258 + mrSges(5,2) * t257;
t306 = -t26 - t93;
t159 = -t188 * t250 + t192 * t194;
t160 = t188 * t194 + t192 * t250;
t304 = -t159 * t188 + t160 * t192;
t282 = pkin(1) * t195;
t181 = t193 * t282;
t252 = t189 * t191;
t142 = pkin(2) * t195 + t181 + (-pkin(9) * t194 - qJ(2)) * t252;
t151 = (-pkin(9) * t189 * t190 - pkin(2) * t193 - pkin(1)) * t191;
t102 = -t142 * t190 + t194 * t151;
t74 = -pkin(3) * t309 - qJ(4) * t141 + t102;
t242 = qJ(2) * t248 + t189 * t282;
t137 = (t190 * t195 + t194 * t248) * pkin(9) + t242;
t125 = t201 * t137;
t87 = t142 * t247 + t151 * t250 + t125;
t76 = qJ(4) * t158 + t87;
t47 = -t188 * t76 + t192 * t74;
t48 = t188 * t74 + t192 * t76;
t303 = -t188 * t47 + t192 * t48;
t256 = t142 * t194;
t208 = t151 * t190 + t256;
t86 = -t198 * t137 + t201 * t208;
t35 = -pkin(4) * t309 - pkin(10) * t107 + t47;
t40 = pkin(10) * t106 + t48;
t278 = t197 * t35 + t200 * t40;
t136 = t141 * qJD(3);
t240 = qJD(2) * t191;
t228 = t189 * t240;
t219 = t194 * t228;
t238 = qJD(3) * t201;
t225 = t190 * t238;
t227 = t193 * t240;
t70 = t238 * t256 + t151 * t225 + t201 * t227 + (-qJD(3) * t137 - t219) * t198;
t67 = qJD(4) * t158 + t70;
t220 = t190 * t228;
t82 = pkin(3) * t136 - qJ(4) * t135 - qJD(4) * t141 + t220;
t43 = -t188 * t67 + t192 * t82;
t32 = pkin(4) * t136 - pkin(10) * t257 + t43;
t44 = t188 * t82 + t192 * t67;
t39 = -pkin(10) * t258 + t44;
t6 = -qJD(5) * t278 - t197 * t39 + t200 * t32;
t174 = -mrSges(7,1) * t199 + mrSges(7,2) * t196;
t302 = -m(7) * pkin(5) - mrSges(6,1) + t174;
t301 = 2 * m(5);
t300 = 2 * m(6);
t299 = 0.2e1 * m(7);
t298 = -2 * mrSges(6,3);
t281 = pkin(10) + qJ(4);
t173 = t281 * t192;
t222 = qJD(5) * t281;
t232 = t200 * qJD(4);
t233 = t197 * qJD(4);
t236 = qJD(5) * t200;
t116 = t192 * t233 + t173 * t236 + (-t197 * t222 + t232) * t188;
t297 = 0.2e1 * t116;
t223 = t281 * t188;
t148 = t173 * t197 + t200 * t223;
t296 = 0.2e1 * t148;
t63 = -t196 * t309 + t199 * t81;
t30 = -qJD(6) * t63 + t136 * t199 - t196 * t56;
t295 = t30 / 0.2e1;
t162 = t165 * qJD(5);
t234 = qJD(6) * t199;
t204 = -t196 * t161 + t165 * t234;
t89 = -Ifges(7,4) * t203 - Ifges(7,2) * t204 + Ifges(7,6) * t162;
t294 = t89 / 0.2e1;
t90 = -Ifges(7,1) * t203 - Ifges(7,4) * t204 + Ifges(7,5) * t162;
t293 = t90 / 0.2e1;
t272 = Ifges(7,4) * t196;
t216 = Ifges(7,1) * t199 - t272;
t114 = Ifges(7,5) * t164 + t165 * t216;
t292 = t114 / 0.2e1;
t291 = -t165 / 0.2e1;
t184 = Ifges(7,5) * t234;
t290 = -Ifges(7,6) * t235 / 0.2e1 + t184 / 0.2e1;
t169 = t216 * qJD(6);
t289 = t169 / 0.2e1;
t288 = Ifges(7,5) * t307 + Ifges(7,6) * t283;
t176 = Ifges(7,2) * t199 + t272;
t287 = -t176 / 0.2e1;
t271 = Ifges(7,4) * t199;
t177 = Ifges(7,1) * t196 + t271;
t286 = t177 / 0.2e1;
t285 = t192 / 0.2e1;
t284 = -t196 / 0.2e1;
t62 = -t196 * t81 - t199 * t309;
t29 = qJD(6) * t62 + t136 * t196 + t199 * t56;
t12 = -mrSges(7,1) * t30 + mrSges(7,2) * t29;
t49 = mrSges(6,1) * t136 - mrSges(6,3) * t56;
t280 = t12 - t49;
t37 = -mrSges(7,1) * t62 + mrSges(7,2) * t63;
t66 = -mrSges(6,1) * t309 - mrSges(6,3) * t81;
t279 = t37 - t66;
t277 = mrSges(4,3) * t135;
t276 = mrSges(4,3) * t136;
t275 = mrSges(7,3) * t165;
t274 = Ifges(5,4) * t188;
t273 = Ifges(5,4) * t192;
t270 = Ifges(7,6) * t196;
t269 = t136 * Ifges(6,5);
t268 = t136 * Ifges(6,6);
t71 = (t189 * t246 + t193 * t198) * t240 + (t198 * t208 + t125) * qJD(3);
t265 = t201 * t71;
t264 = -mrSges(5,1) * t192 + mrSges(5,2) * t188 - mrSges(4,1);
t183 = -pkin(4) * t192 - pkin(3);
t134 = pkin(5) * t164 - pkin(11) * t165 + t183;
t149 = t200 * t173 - t197 * t223;
t98 = t134 * t199 - t149 * t196;
t262 = qJD(6) * t98;
t99 = t134 * t196 + t149 * t199;
t261 = qJD(6) * t99;
t260 = t116 * t148;
t118 = t159 * t197 + t160 * t200;
t101 = qJD(5) * t118 + t165 * t225;
t207 = t200 * t159 - t160 * t197;
t259 = t207 * t101;
t253 = t190 ^ 2 * t198;
t244 = Ifges(4,5) * t135 - Ifges(4,6) * t136;
t243 = -Ifges(6,5) * t161 - Ifges(6,6) * t162;
t237 = qJD(5) * t197;
t9 = Ifges(7,5) * t29 + Ifges(7,6) * t30 + Ifges(7,3) * t57;
t231 = Ifges(6,5) * t56 - Ifges(6,6) * t57 + Ifges(6,3) * t136;
t226 = t190 * t239;
t221 = t241 * mrSges(5,3);
t129 = t162 * mrSges(6,1) - t161 * mrSges(6,2);
t59 = pkin(4) * t258 + t71;
t19 = t57 * pkin(5) - t56 * pkin(11) + t59;
t5 = t197 * t32 + t200 * t39 + t35 * t236 - t237 * t40;
t3 = pkin(11) * t136 + t5;
t14 = -pkin(11) * t309 + t278;
t79 = -t158 * pkin(3) - t86;
t58 = -t106 * pkin(4) + t79;
t25 = -pkin(5) * t209 - t81 * pkin(11) + t58;
t7 = -t14 * t196 + t199 * t25;
t1 = qJD(6) * t7 + t19 * t196 + t199 * t3;
t8 = t14 * t199 + t196 * t25;
t2 = -qJD(6) * t8 + t19 * t199 - t196 * t3;
t218 = t1 * t199 - t2 * t196;
t217 = mrSges(7,1) * t196 + mrSges(7,2) * t199;
t215 = -Ifges(7,2) * t196 + t271;
t214 = -t188 * t43 + t192 * t44;
t91 = mrSges(5,2) * t309 + mrSges(5,3) * t106;
t92 = -mrSges(5,1) * t309 - mrSges(5,3) * t107;
t213 = -t188 * t92 + t192 * t91;
t15 = -t197 * t40 + t200 * t35;
t210 = t101 * t148 - t116 * t207;
t23 = Ifges(7,4) * t63 + Ifges(7,2) * t62 - Ifges(7,6) * t209;
t24 = Ifges(7,1) * t63 + Ifges(7,4) * t62 - Ifges(7,5) * t209;
t206 = t23 * t284 + t24 * t283;
t108 = -t196 * t118 - t199 * t249;
t205 = -t199 * t118 + t196 * t249;
t88 = -Ifges(7,5) * t203 - Ifges(7,6) * t204 + Ifges(7,3) * t162;
t168 = t215 * qJD(6);
t166 = t217 * qJD(6);
t146 = Ifges(6,1) * t165 - Ifges(6,4) * t164;
t145 = Ifges(6,4) * t165 - Ifges(6,2) * t164;
t144 = mrSges(6,1) * t164 + mrSges(6,2) * t165;
t139 = mrSges(7,1) * t164 - t199 * t275;
t138 = -mrSges(7,2) * t164 - t196 * t275;
t133 = pkin(5) * t162 + pkin(11) * t161;
t132 = t217 * t165;
t131 = -Ifges(6,1) * t161 - Ifges(6,4) * t162;
t130 = -Ifges(6,4) * t161 - Ifges(6,2) * t162;
t115 = t192 * t232 - t173 * t237 + (-t200 * t222 - t233) * t188;
t113 = Ifges(7,6) * t164 + t165 * t215;
t112 = Ifges(7,3) * t164 + (Ifges(7,5) * t199 - t270) * t165;
t111 = mrSges(4,1) * t158 - mrSges(4,3) * t141;
t110 = -mrSges(4,2) * t158 + mrSges(4,3) * t309;
t105 = -mrSges(7,2) * t162 - mrSges(7,3) * t204;
t104 = mrSges(7,1) * t162 + mrSges(7,3) * t203;
t100 = qJD(5) * t207 - t164 * t225;
t97 = mrSges(4,1) * t136 + mrSges(4,2) * t135;
t96 = mrSges(5,1) * t136 - mrSges(5,3) * t257;
t95 = -mrSges(5,2) * t136 - mrSges(5,3) * t258;
t94 = mrSges(7,1) * t204 - mrSges(7,2) * t203;
t85 = -mrSges(5,1) * t106 + mrSges(5,2) * t107;
t84 = t136 * Ifges(5,5) + (t192 * Ifges(5,1) - t274) * t135;
t83 = t136 * Ifges(5,6) + (-t188 * Ifges(5,2) + t273) * t135;
t69 = qJD(6) * t205 - t196 * t100 + t199 * t226;
t68 = qJD(6) * t108 + t199 * t100 + t196 * t226;
t65 = mrSges(6,2) * t309 + mrSges(6,3) * t209;
t61 = -t115 * t196 + t133 * t199 - t261;
t60 = t115 * t199 + t133 * t196 + t262;
t51 = -mrSges(6,1) * t209 + mrSges(6,2) * t81;
t50 = -mrSges(6,2) * t136 - mrSges(6,3) * t57;
t46 = Ifges(6,1) * t81 + Ifges(6,4) * t209 - Ifges(6,5) * t309;
t45 = Ifges(6,4) * t81 + Ifges(6,2) * t209 - Ifges(6,6) * t309;
t42 = -mrSges(7,1) * t209 - mrSges(7,3) * t63;
t41 = mrSges(7,2) * t209 + mrSges(7,3) * t62;
t22 = Ifges(7,5) * t63 + Ifges(7,6) * t62 - Ifges(7,3) * t209;
t21 = Ifges(6,1) * t56 - Ifges(6,4) * t57 + t269;
t20 = Ifges(6,4) * t56 - Ifges(6,2) * t57 + t268;
t18 = -mrSges(7,2) * t57 + mrSges(7,3) * t30;
t17 = mrSges(7,1) * t57 - mrSges(7,3) * t29;
t13 = pkin(5) * t309 - t15;
t11 = Ifges(7,1) * t29 + Ifges(7,4) * t30 + Ifges(7,5) * t57;
t10 = Ifges(7,4) * t29 + Ifges(7,2) * t30 + Ifges(7,6) * t57;
t4 = -pkin(5) * t136 - t6;
t16 = [-(t9 - t20) * t209 + (t22 - t45) * t57 + (t1 * t8 + t13 * t4 + t2 * t7) * t299 + (t43 * t47 + t44 * t48 + t71 * t79) * t301 + (t15 * t6 + t278 * t5 + t58 * t59) * t300 + 0.2e1 * t278 * t50 + 0.2e1 * m(3) * (t242 * t193 + (qJ(2) * t252 - t181) * t189) * t240 + 0.2e1 * m(4) * (t102 * t220 + t70 * t87 - t71 * t86) + 0.2e1 * t13 * t12 + 0.2e1 * t7 * t17 + 0.2e1 * t8 * t18 + t29 * t24 + t30 * t23 + 0.2e1 * t4 * t37 + t158 * t244 + 0.2e1 * t1 * t41 + 0.2e1 * t2 * t42 + 0.2e1 * (-mrSges(3,2) * t195 + mrSges(3,3) * t248) * t227 + 0.2e1 * t15 * t49 - 0.2e1 * (mrSges(3,1) * t195 - mrSges(3,3) * t252) * t228 + t56 * t46 + 0.2e1 * t58 * t26 + 0.2e1 * t59 * t51 + t62 * t10 + t63 * t11 + 0.2e1 * t5 * t65 + 0.2e1 * t6 * t66 + t81 * t21 + 0.2e1 * t71 * t85 + 0.2e1 * t44 * t91 - 0.2e1 * t87 * t276 + 0.2e1 * t43 * t92 - 0.2e1 * t86 * t277 + 0.2e1 * t79 * t93 + 0.2e1 * t48 * t95 + 0.2e1 * t47 * t96 + 0.2e1 * t102 * t97 + t106 * t83 + t107 * t84 + 0.2e1 * t70 * t110 - 0.2e1 * t71 * t111 + (0.2e1 * Ifges(4,1) * t141 + Ifges(4,5) * t158 - (Ifges(5,5) * t192 - Ifges(5,6) * t188 + t308) * t309) * t135 + 0.2e1 * (-mrSges(4,1) * t309 + mrSges(4,2) * t141) * t220 - (Ifges(5,4) * t107 + Ifges(5,2) * t106 - Ifges(5,6) * t309) * t258 + (t141 * t308 + Ifges(5,5) * t107 + Ifges(6,5) * t81 - Ifges(4,6) * t158 + Ifges(5,6) * t106 + Ifges(6,6) * t209 - ((2 * Ifges(4,2)) + (2 * Ifges(5,3)) + Ifges(6,3)) * t309) * t136 + (Ifges(5,1) * t107 + Ifges(5,4) * t106 - Ifges(5,5) * t309) * t257 - t309 * t231; t100 * t65 + t108 * t17 - t205 * t18 + t118 * t50 + t159 * t96 + t160 * t95 + t194 * t97 + t68 * t41 + t69 * t42 - t280 * t207 + t279 * t101 + m(7) * (-t1 * t205 + t101 * t13 + t108 * t2 - t207 * t4 + t68 * t8 + t69 * t7) + m(6) * (t100 * t278 - t101 * t15 + t118 * t5 + t207 * t6) + m(5) * (t159 * t43 + t160 * t44) + (-t198 * t276 + (-t277 + t306) * t201 + ((t110 + t213) * t201 + (-t111 + t51 + t85) * t198) * qJD(3) + m(6) * (-t201 * t59 + t239 * t58) + m(5) * (t238 * t303 + t239 * t79 - t265) + m(4) * (t198 * t70 + t238 * t87 - t239 * t86 + t219 - t265)) * t190; 0.2e1 * m(7) * (t108 * t69 - t205 * t68 - t259) + 0.2e1 * m(5) * (t190 * t304 - t253) * t238 + 0.2e1 * (t118 * t100 - t238 * t253 - t259) * m(6); -(t88 / 0.2e1 - t130 / 0.2e1) * t209 + t213 * qJD(4) + t214 * mrSges(5,3) + m(7) * (t1 * t99 + t116 * t13 + t148 * t4 + t2 * t98 + t60 * t8 + t61 * t7) - (-t15 * mrSges(6,3) + t46 / 0.2e1 + t206) * t161 + (t112 / 0.2e1 - t145 / 0.2e1) * t57 + t63 * t293 + t62 * t294 + t113 * t295 + t83 * t285 + t29 * t292 + m(6) * (t115 * t278 - t116 * t15 - t148 * t6 + t149 * t5 + t183 * t59) + (-t278 * mrSges(6,3) + t22 / 0.2e1 - t45 / 0.2e1) * t162 + m(5) * (-pkin(3) * t71 + t214 * qJ(4) + qJD(4) * t303) + (-t188 * t96 + t192 * t95) * qJ(4) + t244 + t60 * t41 + t61 * t42 + t264 * t71 - t70 * mrSges(4,2) + (-t5 * mrSges(6,3) + t9 / 0.2e1 - t20 / 0.2e1 - t268 / 0.2e1) * t164 - pkin(3) * t93 + t13 * t94 + t279 * t116 + t98 * t17 + t280 * t148 + t99 * t18 + t7 * t104 + t8 * t105 + t115 * t65 + (-t6 * mrSges(6,3) + t10 * t284 + t11 * t283 + t21 / 0.2e1 + t269 / 0.2e1 + (t24 * t284 - t199 * t23 / 0.2e1) * qJD(6)) * t165 + (-t188 * (Ifges(5,2) * t192 + t274) / 0.2e1 + (Ifges(5,1) * t188 + t273) * t285) * t135 + t58 * t129 + t81 * t131 / 0.2e1 + t4 * t132 + t1 * t138 + t2 * t139 + t59 * t144 + t56 * t146 / 0.2e1 + t149 * t50 + t183 * t26 + t188 * t84 / 0.2e1 + t136 * (Ifges(5,5) * t188 + Ifges(5,6) * t192) / 0.2e1 - t309 * t243 / 0.2e1; t101 * t132 + t108 * t104 - t205 * t105 - t207 * t94 + t68 * t138 + t69 * t139 + m(7) * (t108 * t61 - t205 * t60 + t68 * t99 + t69 * t98 + t210) + m(6) * (t100 * t149 + t115 * t118 + t210) + m(5) * t304 * qJD(4) + ((-m(5) * pkin(3) + m(6) * t183 + t144 + t264) * t239 + (-t129 + (m(5) * t310 - mrSges(4,2) + t221) * qJD(3)) * t201) * t190 + (-t100 * t164 + t101 * t165 - t118 * t162 + t161 * t207) * mrSges(6,3); 0.2e1 * t98 * t104 + 0.2e1 * t99 * t105 + t132 * t297 + 0.2e1 * t183 * t129 + 0.2e1 * t60 * t138 + 0.2e1 * t61 * t139 + t94 * t296 + (t60 * t99 + t61 * t98 + t260) * t299 + (t115 * t149 + t260) * t300 + (t115 * t298 - t130 + t88) * t164 + (t149 * t298 + t112 - t145) * t162 - (mrSges(6,3) * t296 - t113 * t196 + t114 * t199 + t146) * t161 + (t301 * t310 + 0.2e1 * t221) * qJD(4) + (mrSges(6,3) * t297 - t196 * t89 + t199 * t90 + t131 + (-t113 * t199 - t114 * t196) * qJD(6)) * t165; t199 * t17 + t196 * t18 + (-t196 * t42 + t199 * t41) * qJD(6) + m(7) * (t1 * t196 + t199 * t2 + (-t196 * t7 + t199 * t8) * qJD(6)) + m(6) * t59 + m(5) * t71 - t306; m(7) * (t196 * t68 + t199 * t69 + (-t108 * t196 - t199 * t205) * qJD(6)) + (m(5) + m(6)) * t226; m(7) * (t196 * t60 + t199 * t61 + (-t196 * t98 + t199 * t99) * qJD(6)) + t138 * t234 + t196 * t105 - t139 * t235 + t199 * t104 + t129; 0; -t5 * mrSges(6,2) + t6 * mrSges(6,1) + t13 * t166 - t209 * t290 + t62 * t168 / 0.2e1 + t63 * t289 + t4 * t174 + t57 * t288 + t176 * t295 + t29 * t286 + t11 * t307 + t10 * t283 + t206 * qJD(6) + (-m(7) * t4 - t12) * pkin(5) + ((-t196 * t8 - t199 * t7) * qJD(6) + t218) * mrSges(7,3) + (-t41 * t235 - t42 * t234 - t196 * t17 + t199 * t18 + m(7) * (-t234 * t7 - t235 * t8 + t218)) * pkin(11) + t231; -t100 * mrSges(6,2) - t207 * t166 + (m(7) * pkin(11) + mrSges(7,3)) * (-t69 * t196 + t68 * t199 + (-t108 * t199 + t196 * t205) * qJD(6)) + t302 * t101; -pkin(5) * t94 - t115 * mrSges(6,2) + t148 * t166 + t164 * t290 + t162 * t288 + t302 * t116 + (t60 * mrSges(7,3) + t165 * t289 - t161 * t286 + t294 + (-t98 * mrSges(7,3) + t165 * t287 + t292) * qJD(6) + (m(7) * (t60 - t262) + t105 - qJD(6) * t139) * pkin(11)) * t199 + (-t61 * mrSges(7,3) + t168 * t291 - t161 * t287 + t293 + (-t99 * mrSges(7,3) + t177 * t291 - t113 / 0.2e1) * qJD(6) + (-qJD(6) * t138 - t104 + m(7) * (-t61 - t261)) * pkin(11)) * t196 + t243; 0; -0.2e1 * pkin(5) * t166 + t168 * t199 + t169 * t196 + (-t176 * t196 + t177 * t199) * qJD(6); mrSges(7,1) * t2 - mrSges(7,2) * t1 + t9; mrSges(7,1) * t69 - mrSges(7,2) * t68; mrSges(7,1) * t61 - mrSges(7,2) * t60 + t88; -t166; t184 + (pkin(11) * t174 - t270) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t16(1) t16(2) t16(4) t16(7) t16(11) t16(16); t16(2) t16(3) t16(5) t16(8) t16(12) t16(17); t16(4) t16(5) t16(6) t16(9) t16(13) t16(18); t16(7) t16(8) t16(9) t16(10) t16(14) t16(19); t16(11) t16(12) t16(13) t16(14) t16(15) t16(20); t16(16) t16(17) t16(18) t16(19) t16(20) t16(21);];
Mq  = res;
