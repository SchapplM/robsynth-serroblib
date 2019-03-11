% Calculate time derivative of joint inertia matrix for
% S6RRPRRP13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
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
% Datum: 2019-03-09 13:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRP13_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP13_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP13_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP13_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP13_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP13_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP13_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:55:37
% EndTime: 2019-03-09 12:55:52
% DurationCPUTime: 6.94s
% Computational Cost: add. (5452->660), mult. (14059->898), div. (0->0), fcn. (12352->8), ass. (0->265)
t204 = cos(qJ(5));
t321 = Ifges(6,6) + Ifges(7,6);
t324 = t321 * t204;
t201 = sin(qJ(5));
t322 = Ifges(6,5) + Ifges(7,5);
t323 = t322 * t201;
t202 = sin(qJ(4));
t260 = qJD(4) * t202;
t232 = t201 * t260;
t205 = cos(qJ(4));
t253 = qJD(5) * t205;
t233 = t204 * t253;
t209 = t232 - t233;
t200 = cos(pkin(6));
t199 = sin(pkin(6));
t206 = cos(qJ(2));
t277 = t199 * t206;
t242 = t202 * t277;
t132 = t200 * t205 - t242;
t203 = sin(qJ(2));
t278 = t199 * t203;
t100 = t132 * t204 + t201 * t278;
t263 = qJD(2) * t206;
t238 = t199 * t263;
t131 = t200 * t202 + t205 * t277;
t264 = qJD(2) * t199;
t239 = t203 * t264;
t97 = -qJD(4) * t131 + t202 * t239;
t45 = -qJD(5) * t100 - t201 * t97 + t204 * t238;
t99 = -t132 * t201 + t204 * t278;
t46 = qJD(5) * t99 + t201 * t238 + t204 * t97;
t258 = qJD(4) * t205;
t98 = -qJD(4) * t242 + t200 * t258 - t205 * t239;
t7 = Ifges(7,5) * t46 + Ifges(7,6) * t45 + Ifges(7,3) * t98;
t8 = Ifges(6,5) * t46 + Ifges(6,6) * t45 + Ifges(6,3) * t98;
t320 = t7 + t8;
t316 = -m(6) * pkin(10) - mrSges(6,3);
t207 = -pkin(2) - pkin(9);
t257 = qJD(4) * t207;
t108 = -pkin(5) * t209 + t202 * t257;
t234 = t201 * t253;
t259 = qJD(4) * t204;
t237 = t202 * t259;
t210 = t234 + t237;
t83 = -mrSges(7,1) * t209 - mrSges(7,2) * t210;
t315 = m(7) * t108 + t83;
t298 = pkin(10) * t205;
t299 = pkin(4) * t202;
t158 = qJ(3) - t298 + t299;
t275 = t202 * t207;
t115 = t201 * t158 + t204 * t275;
t265 = t201 ^ 2 + t204 ^ 2;
t163 = -mrSges(7,1) * t204 + mrSges(7,2) * t201;
t188 = -pkin(5) * t204 - pkin(4);
t314 = m(7) * t188 + t163;
t313 = 2 * m(6);
t312 = 0.2e1 * m(7);
t311 = -2 * pkin(1);
t310 = 2 * mrSges(4,1);
t309 = -2 * mrSges(3,3);
t308 = -2 * mrSges(7,3);
t280 = qJ(3) * t203;
t211 = -pkin(2) * t206 - t280;
t117 = (-pkin(1) + t211) * t199;
t307 = -0.2e1 * t117;
t306 = m(6) * pkin(4);
t304 = m(7) * pkin(5);
t303 = pkin(3) + pkin(8);
t302 = m(6) * t202;
t300 = pkin(1) * t206;
t297 = mrSges(6,2) + mrSges(7,2);
t295 = Ifges(4,6) + Ifges(3,4);
t293 = -qJ(6) - pkin(10);
t103 = (t206 * t207 - pkin(1) - t280) * t199;
t184 = pkin(8) * t278;
t240 = -pkin(2) - t300;
t88 = pkin(3) * t278 + t184 + (-pkin(9) + t240) * t200;
t53 = t205 * t103 + t202 * t88;
t48 = pkin(10) * t278 + t53;
t187 = t200 * t203 * pkin(1);
t134 = pkin(8) * t277 + t187;
t116 = -t200 * qJ(3) - t134;
t102 = pkin(3) * t277 - t116;
t57 = pkin(4) * t131 - pkin(10) * t132 + t102;
t20 = t201 * t57 + t204 * t48;
t62 = -mrSges(7,2) * t131 + mrSges(7,3) * t99;
t63 = -mrSges(6,2) * t131 + mrSges(6,3) * t99;
t292 = t62 + t63;
t64 = mrSges(7,1) * t131 - mrSges(7,3) * t100;
t65 = mrSges(6,1) * t131 - mrSges(6,3) * t100;
t291 = -t64 - t65;
t15 = -mrSges(6,1) * t45 + mrSges(6,2) * t46;
t68 = mrSges(5,1) * t238 - mrSges(5,3) * t97;
t290 = t68 - t15;
t289 = Ifges(5,4) * t202;
t288 = Ifges(5,4) * t205;
t287 = Ifges(6,4) * t201;
t286 = Ifges(6,4) * t204;
t285 = Ifges(7,4) * t201;
t284 = Ifges(7,4) * t204;
t251 = t200 * t300;
t179 = qJD(2) * t251;
t125 = -pkin(8) * t239 + t179;
t283 = t125 * mrSges(3,2);
t107 = mrSges(5,1) * t278 - mrSges(5,3) * t132;
t59 = -mrSges(6,1) * t99 + mrSges(6,2) * t100;
t282 = -t107 + t59;
t164 = -mrSges(6,1) * t204 + mrSges(6,2) * t201;
t281 = t164 - mrSges(5,1);
t126 = t134 * qJD(2);
t279 = t126 * t203;
t276 = t201 * t205;
t274 = t204 * t205;
t273 = t205 * t207;
t212 = -Ifges(7,2) * t201 + t284;
t121 = Ifges(7,6) * t202 + t205 * t212;
t213 = -Ifges(6,2) * t201 + t286;
t122 = Ifges(6,6) * t202 + t205 * t213;
t272 = t121 + t122;
t214 = Ifges(7,1) * t204 - t285;
t123 = Ifges(7,5) * t202 + t205 * t214;
t215 = Ifges(6,1) * t204 - t287;
t124 = Ifges(6,5) * t202 + t205 * t215;
t271 = -t123 - t124;
t154 = -mrSges(7,2) * t202 - mrSges(7,3) * t276;
t155 = -mrSges(6,2) * t202 - mrSges(6,3) * t276;
t270 = t154 + t155;
t156 = mrSges(7,1) * t202 - mrSges(7,3) * t274;
t157 = mrSges(6,1) * t202 - mrSges(6,3) * t274;
t269 = -t156 - t157;
t268 = Ifges(3,5) * t238 + Ifges(4,5) * t239;
t267 = Ifges(7,6) * t232 + Ifges(7,3) * t258;
t266 = Ifges(6,6) * t232 + Ifges(6,3) * t258;
t254 = qJD(5) * t204;
t256 = qJD(5) * t201;
t143 = mrSges(7,1) * t256 + mrSges(7,2) * t254;
t262 = qJD(3) * t203;
t261 = qJD(4) * t201;
t255 = qJD(5) * t202;
t252 = qJD(6) * t204;
t10 = Ifges(6,4) * t46 + Ifges(6,2) * t45 + Ifges(6,6) * t98;
t9 = Ifges(7,4) * t46 + Ifges(7,2) * t45 + Ifges(7,6) * t98;
t250 = -t9 / 0.2e1 - t10 / 0.2e1;
t249 = Ifges(5,5) * t97 - Ifges(5,6) * t98 + Ifges(5,3) * t238;
t11 = Ifges(7,1) * t46 + Ifges(7,4) * t45 + Ifges(7,5) * t98;
t12 = Ifges(6,1) * t46 + Ifges(6,4) * t45 + Ifges(6,5) * t98;
t248 = t11 / 0.2e1 + t12 / 0.2e1;
t34 = Ifges(7,4) * t100 + Ifges(7,2) * t99 + Ifges(7,6) * t131;
t35 = Ifges(6,4) * t100 + Ifges(6,2) * t99 + Ifges(6,6) * t131;
t247 = t34 / 0.2e1 + t35 / 0.2e1;
t36 = Ifges(7,1) * t100 + Ifges(7,4) * t99 + Ifges(7,5) * t131;
t37 = Ifges(6,1) * t100 + Ifges(6,4) * t99 + Ifges(6,5) * t131;
t246 = -t36 / 0.2e1 - t37 / 0.2e1;
t169 = Ifges(7,2) * t204 + t285;
t72 = -t169 * t253 + (Ifges(7,6) * t205 - t202 * t212) * qJD(4);
t170 = Ifges(6,2) * t204 + t287;
t73 = -t170 * t253 + (Ifges(6,6) * t205 - t202 * t213) * qJD(4);
t245 = t72 / 0.2e1 + t73 / 0.2e1;
t172 = Ifges(7,1) * t201 + t284;
t74 = -t172 * t253 + (Ifges(7,5) * t205 - t202 * t214) * qJD(4);
t173 = Ifges(6,1) * t201 + t286;
t75 = -t173 * t253 + (Ifges(6,5) * t205 - t202 * t215) * qJD(4);
t244 = t74 / 0.2e1 + t75 / 0.2e1;
t243 = mrSges(7,1) + t304;
t140 = qJD(3) + (pkin(4) * t205 + pkin(10) * t202) * qJD(4);
t236 = t205 * t257;
t241 = t201 * t140 + t158 * t254 + t204 * t236;
t235 = t201 * t255;
t231 = -t278 / 0.2e1;
t230 = t121 / 0.2e1 + t122 / 0.2e1;
t229 = t123 / 0.2e1 + t124 / 0.2e1;
t193 = Ifges(7,5) * t254;
t194 = Ifges(6,5) * t254;
t228 = t193 / 0.2e1 + t194 / 0.2e1 - t321 * t256 / 0.2e1;
t148 = t212 * qJD(5);
t149 = t213 * qJD(5);
t227 = t148 / 0.2e1 + t149 / 0.2e1;
t151 = t214 * qJD(5);
t152 = t215 * qJD(5);
t226 = t151 / 0.2e1 + t152 / 0.2e1;
t225 = t323 / 0.2e1 + t324 / 0.2e1;
t224 = t170 / 0.2e1 + t169 / 0.2e1;
t223 = t172 / 0.2e1 + t173 / 0.2e1;
t14 = -t45 * mrSges(7,1) + t46 * mrSges(7,2);
t222 = (mrSges(4,2) - mrSges(3,1)) * t126;
t221 = -t201 * t207 + pkin(5);
t19 = -t201 * t48 + t204 * t57;
t52 = -t202 * t103 + t205 * t88;
t220 = qJD(5) * t293;
t219 = -mrSges(6,1) - t243;
t218 = Ifges(5,5) * t238;
t217 = Ifges(5,6) * t238;
t216 = mrSges(6,1) * t201 + mrSges(6,2) * t204;
t104 = (t277 * t303 + t187) * qJD(2);
t178 = pkin(2) * t239;
t82 = t178 + (-t262 + (pkin(9) * t203 - qJ(3) * t206) * qJD(2)) * t199;
t22 = -t103 * t258 + t104 * t205 - t202 * t82 - t88 * t260;
t21 = -t103 * t260 + t202 * t104 + t205 * t82 + t88 * t258;
t17 = pkin(10) * t238 + t21;
t195 = t200 * qJD(3);
t87 = -t239 * t303 + t179 + t195;
t30 = pkin(4) * t98 - pkin(10) * t97 + t87;
t3 = t204 * t17 + t201 * t30 + t57 * t254 - t256 * t48;
t47 = -pkin(4) * t278 - t52;
t208 = -t115 * qJD(5) + t204 * t140;
t18 = -pkin(4) * t238 - t22;
t4 = -qJD(5) * t20 - t17 * t201 + t204 * t30;
t174 = Ifges(5,1) * t205 - t289;
t171 = -Ifges(5,2) * t202 + t288;
t166 = t202 * mrSges(5,1) + t205 * mrSges(5,2);
t165 = t293 * t204;
t162 = t293 * t201;
t153 = (-Ifges(5,1) * t202 - t288) * qJD(4);
t150 = (-Ifges(5,2) * t205 - t289) * qJD(4);
t145 = (mrSges(5,1) * t205 - mrSges(5,2) * t202) * qJD(4);
t144 = t216 * qJD(5);
t142 = -mrSges(4,1) * t277 - mrSges(4,3) * t200;
t141 = (pkin(5) * t201 - t207) * t205;
t139 = t204 * t158;
t137 = t216 * t205;
t136 = (mrSges(7,1) * t201 + mrSges(7,2) * t204) * t205;
t133 = -t184 + t251;
t128 = -qJD(6) * t201 + t204 * t220;
t127 = t201 * t220 + t252;
t120 = Ifges(6,3) * t202 + (Ifges(6,5) * t204 - Ifges(6,6) * t201) * t205;
t119 = Ifges(7,3) * t202 + (Ifges(7,5) * t204 - Ifges(7,6) * t201) * t205;
t118 = t200 * t240 + t184;
t114 = -t201 * t275 + t139;
t113 = -t125 - t195;
t112 = -mrSges(6,2) * t258 + mrSges(6,3) * t209;
t111 = -mrSges(7,2) * t258 + mrSges(7,3) * t209;
t110 = mrSges(6,1) * t258 + mrSges(6,3) * t210;
t109 = mrSges(7,1) * t258 + mrSges(7,3) * t210;
t106 = -mrSges(5,2) * t278 - mrSges(5,3) * t131;
t105 = t178 + (-qJ(3) * t263 - t262) * t199;
t95 = -qJ(6) * t276 + t115;
t85 = -qJ(6) * t274 + t202 * t221 + t139;
t84 = -mrSges(6,1) * t209 - mrSges(6,2) * t210;
t78 = mrSges(5,1) * t131 + mrSges(5,2) * t132;
t71 = -Ifges(6,5) * t210 - Ifges(6,6) * t233 + t266;
t70 = -Ifges(7,5) * t210 - Ifges(7,6) * t233 + t267;
t69 = -mrSges(5,2) * t238 - mrSges(5,3) * t98;
t67 = Ifges(5,1) * t132 - Ifges(5,4) * t131 + Ifges(5,5) * t278;
t66 = Ifges(5,4) * t132 - Ifges(5,2) * t131 + Ifges(5,6) * t278;
t61 = -t201 * t236 + t208;
t60 = -t207 * t235 + t241;
t58 = -mrSges(7,1) * t99 + mrSges(7,2) * t100;
t56 = mrSges(5,1) * t98 + mrSges(5,2) * t97;
t50 = Ifges(5,1) * t97 - Ifges(5,4) * t98 + t218;
t49 = Ifges(5,4) * t97 - Ifges(5,2) * t98 + t217;
t38 = -qJ(6) * t233 + (-qJD(6) * t205 + (qJ(6) * qJD(4) - qJD(5) * t207) * t202) * t201 + t241;
t33 = Ifges(6,5) * t100 + Ifges(6,6) * t99 + Ifges(6,3) * t131;
t32 = Ifges(7,5) * t100 + Ifges(7,6) * t99 + Ifges(7,3) * t131;
t31 = qJ(6) * t237 + (qJ(6) * t256 + qJD(4) * t221 - t252) * t205 + t208;
t27 = -pkin(5) * t99 + t47;
t26 = mrSges(6,1) * t98 - mrSges(6,3) * t46;
t25 = mrSges(7,1) * t98 - mrSges(7,3) * t46;
t24 = -mrSges(6,2) * t98 + mrSges(6,3) * t45;
t23 = -mrSges(7,2) * t98 + mrSges(7,3) * t45;
t13 = qJ(6) * t99 + t20;
t6 = pkin(5) * t131 - qJ(6) * t100 + t19;
t5 = -pkin(5) * t45 + t18;
t2 = qJ(6) * t45 + qJD(6) * t99 + t3;
t1 = pkin(5) * t98 - qJ(6) * t46 - qJD(6) * t100 + t4;
t16 = [(0.2e1 * t222 + t268 - 0.2e1 * t283) * t200 + (0.2e1 * t105 * (mrSges(4,2) * t206 - mrSges(4,3) * t203) + t279 * t310 + t203 * t249 + 0.2e1 * (t125 * t206 + t279) * mrSges(3,3) + ((t116 * t310 + mrSges(4,2) * t307 + t134 * t309 + (Ifges(4,5) - (2 * Ifges(3,6))) * t200 + (mrSges(3,1) * t311 - 0.2e1 * t203 * t295) * t199) * t203 + (t118 * t310 + t133 * t309 + mrSges(4,3) * t307 + Ifges(5,5) * t132 - Ifges(5,6) * t131 + (-(2 * Ifges(4,4)) + Ifges(3,5)) * t200 + (mrSges(3,2) * t311 + 0.2e1 * t206 * t295) * t199 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) + (2 * Ifges(4,2)) - (2 * Ifges(4,3)) + Ifges(5,3)) * t278) * t206) * qJD(2)) * t199 + (t34 + t35) * t45 + 0.2e1 * t113 * t142 + t132 * t50 + 0.2e1 * t21 * t106 + 0.2e1 * t22 * t107 + 0.2e1 * t102 * t56 + t97 * t67 + 0.2e1 * t87 * t78 + 0.2e1 * t52 * t68 + 0.2e1 * t53 * t69 + 0.2e1 * t2 * t62 + 0.2e1 * t3 * t63 + 0.2e1 * t1 * t64 + 0.2e1 * t4 * t65 + 0.2e1 * t5 * t58 + 0.2e1 * t18 * t59 + 0.2e1 * t47 * t15 + 0.2e1 * t13 * t23 + 0.2e1 * t20 * t24 + 0.2e1 * t6 * t25 + 0.2e1 * t19 * t26 + 0.2e1 * t27 * t14 + (t1 * t6 + t13 * t2 + t27 * t5) * t312 + (t18 * t47 + t19 * t4 + t20 * t3) * t313 + (-t49 + t320) * t131 + (t36 + t37) * t46 + (t33 + t32 - t66) * t98 + (t9 + t10) * t99 + 0.2e1 * m(4) * (t105 * t117 + t113 * t116 + t118 * t126) + 0.2e1 * m(3) * (t125 * t134 - t126 * t133) + 0.2e1 * m(5) * (t102 * t87 + t21 * t53 + t22 * t52) + (t11 + t12) * t100; t222 + (-t49 / 0.2e1 + t7 / 0.2e1 + t8 / 0.2e1 + t207 * t69 - t21 * mrSges(5,3) - t217 / 0.2e1) * t202 - t283 + t97 * t174 / 0.2e1 + t3 * t155 + t1 * t156 + t4 * t157 + t87 * t166 + t18 * t137 + t141 * t14 + t102 * t145 + t132 * t153 / 0.2e1 + t2 * t154 + t5 * t136 + t114 * t26 + t115 * t24 + t108 * t58 + t6 * t109 + t19 * t110 + t13 * t111 + t20 * t112 - t113 * mrSges(4,3) + t95 * t23 + t27 * t83 + t47 * t84 + t85 * t25 + t38 * t62 + t60 * t63 + t31 * t64 + t61 * t65 + qJ(3) * t56 + (mrSges(4,1) * t211 - Ifges(4,4) * t206 - Ifges(3,6) * t203) * t264 + (-t150 / 0.2e1 + t70 / 0.2e1 + t71 / 0.2e1) * t131 + m(7) * (t1 * t85 + t108 * t27 + t13 * t38 + t141 * t5 + t2 * t95 + t31 * t6) + m(4) * (-pkin(2) * t126 - qJ(3) * t113 - qJD(3) * t116) + (-t171 / 0.2e1 + t119 / 0.2e1 + t120 / 0.2e1) * t98 + (-t142 + t78) * qJD(3) + t268 + t229 * t46 + t230 * t45 + t244 * t100 + t245 * t99 + m(6) * (t114 * t4 + t115 * t3 - t18 * t273 + t61 * t19 + t60 * t20) + m(5) * (qJ(3) * t87 + qJD(3) * t102 + t21 * t275 + t22 * t273) + (t50 / 0.2e1 - t22 * mrSges(5,3) + t218 / 0.2e1 + t290 * t207 + t248 * t204 + t250 * t201 + (t201 * t246 - t204 * t247) * qJD(5)) * t205 + ((Ifges(5,6) * t231 - t66 / 0.2e1 + t32 / 0.2e1 + t33 / 0.2e1 - t53 * mrSges(5,3)) * t205 + (t52 * mrSges(5,3) + Ifges(5,5) * t231 - t67 / 0.2e1 + t246 * t204 + t247 * t201) * t202 + (t205 * t106 + t282 * t202 + t47 * t302 + m(5) * (-t202 * t52 + t205 * t53)) * t207) * qJD(4); 0.2e1 * qJ(3) * t145 + 0.2e1 * t108 * t136 + 0.2e1 * t85 * t109 + 0.2e1 * t114 * t110 + 0.2e1 * t95 * t111 + 0.2e1 * t115 * t112 + 0.2e1 * t141 * t83 + 0.2e1 * t38 * t154 + 0.2e1 * t60 * t155 + 0.2e1 * t31 * t156 + 0.2e1 * t61 * t157 + (t114 * t61 + t115 * t60) * t313 + (t108 * t141 + t31 * t85 + t38 * t95) * t312 + 0.2e1 * (mrSges(4,3) + t166 + (m(4) + m(5)) * qJ(3)) * qJD(3) + (-t150 + t70 + t71 + (0.2e1 * t137 * t207 + t201 * t272 + t204 * t271 - t174) * qJD(4)) * t202 + (-0.2e1 * t207 * t84 + t153 + (t74 + t75) * t204 + (-t72 - t73) * t201 + (t201 * t271 - t204 * t272) * qJD(5) + (-0.2e1 * t207 ^ 2 * t302 + t119 + t120 - t171) * qJD(4)) * t205; m(4) * t126 + mrSges(4,1) * t238 + (-t14 + (t201 * t291 + t204 * t292 + t106) * qJD(4) + m(6) * (-t19 * t261 + t20 * t259 - t18) + m(7) * (t13 * t259 - t261 * t6 - t5) + m(5) * (qJD(4) * t53 + t22) + t290) * t205 + (t69 + (t23 + t24) * t204 + (-t25 - t26) * t201 + (t58 + t282) * qJD(4) + (-t201 * t292 + t204 * t291) * qJD(5) + m(6) * (qJD(4) * t47 - t19 * t254 - t20 * t256 - t201 * t4 + t204 * t3) + m(7) * (qJD(4) * t27 - t1 * t201 - t13 * t256 + t2 * t204 - t254 * t6) + m(5) * (-qJD(4) * t52 + t21)) * t202; (-t84 + (t270 * t204 + t269 * t201 + m(6) * (-t114 * t201 + t115 * t204) + m(7) * (-t201 * t85 + t204 * t95)) * qJD(4) - t315) * t205 + ((t111 + t112) * t204 + (-t109 - t110) * t201 + (t136 + t137) * qJD(4) + (-t201 * t270 + t204 * t269) * qJD(5) + m(6) * (-t114 * t254 - t115 * t256 - t201 * t61 + t204 * t60 - 0.2e1 * t236) + m(7) * (qJD(4) * t141 - t201 * t31 + t204 * t38 - t254 * t85 - t256 * t95)) * t202; 0.4e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * (-0.1e1 + t265) * t202 * t258; t188 * t14 + t162 * t25 + t5 * t163 - t165 * t23 + t27 * t143 + t47 * t144 + t127 * t62 + t128 * t64 - t21 * mrSges(5,2) + t22 * mrSges(5,1) - pkin(4) * t15 + t249 + t223 * t46 + t224 * t45 + t225 * t98 + t226 * t100 + t227 * t99 + t228 * t131 + (t2 * mrSges(7,3) + t3 * mrSges(6,3) + (-t19 * mrSges(6,3) - t6 * mrSges(7,3) - t246) * qJD(5) + (m(6) * (-qJD(5) * t19 + t3) + t24 - qJD(5) * t65) * pkin(10) - t250) * t204 + (-t4 * mrSges(6,3) - t1 * mrSges(7,3) + (-m(6) * t4 - t26) * pkin(10) + (-t13 * mrSges(7,3) + pkin(5) * t58 - pkin(10) * t63 + t316 * t20 + t27 * t304 - t247) * qJD(5) + t248) * t201 + m(7) * (t1 * t162 + t127 * t13 + t128 * t6 - t165 * t2 + t188 * t5) + (t164 - t306) * t18; t188 * t83 + t128 * t156 + t162 * t109 + t108 * t163 - t165 * t111 + t141 * t143 + t127 * t154 - pkin(4) * t84 - t144 * t273 + m(7) * (t108 * t188 + t127 * t95 + t128 * t85 + t162 * t31 - t165 * t38) + t228 * t202 + ((-t207 * mrSges(5,2) - Ifges(5,6) + t225) * t205 + (-Ifges(5,5) + (t281 - t306) * t207) * t202) * qJD(4) + (-t61 * mrSges(6,3) - t31 * mrSges(7,3) - t227 * t205 + t224 * t260 + (-m(6) * t61 - t110) * pkin(10) + (-t95 * mrSges(7,3) + pkin(5) * t136 - pkin(10) * t155 + t316 * t115 + t141 * t304 - t205 * t223 - t230) * qJD(5) + t244) * t201 + (t60 * mrSges(6,3) + t38 * mrSges(7,3) + t226 * t205 - t223 * t260 + (m(6) * t60 + t112) * pkin(10) + (-t114 * mrSges(6,3) - t85 * mrSges(7,3) - t224 * t205 + (-m(6) * t114 - t157) * pkin(10) + t229) * qJD(5) + t245) * t204; m(7) * (-pkin(5) * t234 + t165 * t235) + (m(6) * (t265 * t298 - t299) + m(7) * (-t162 * t276 - t165 * t274)) * qJD(4) + (-t144 - t143 + (-mrSges(5,2) + (mrSges(6,3) + mrSges(7,3)) * t265) * qJD(4)) * t205 + (m(7) * (t127 * t204 - t128 * t201 - t162 * t254) + (t281 + t314) * qJD(4)) * t202; -0.2e1 * pkin(4) * t144 + 0.2e1 * t188 * t143 + (-t127 * t165 + t128 * t162) * t312 + (t128 * t308 + t151 + t152 + (0.2e1 * t314 * pkin(5) - t165 * t308 - t169 - t170) * qJD(5)) * t201 + (0.2e1 * t127 * mrSges(7,3) + t148 + t149 + (t162 * t308 + t172 + t173) * qJD(5)) * t204; mrSges(6,1) * t4 + mrSges(7,1) * t1 - mrSges(6,2) * t3 - mrSges(7,2) * t2 + (m(7) * t1 + t25) * pkin(5) + t320; mrSges(6,1) * t61 + mrSges(7,1) * t31 - mrSges(6,2) * t60 - mrSges(7,2) * t38 - t322 * t237 + (m(7) * t31 + t109) * pkin(5) + (-t323 - t324) * t253 + t266 + t267; (t201 * t297 + t204 * t219) * t255 + (t201 * t219 - t204 * t297) * t258; -mrSges(7,2) * t127 + t193 + t194 + t243 * t128 + ((-mrSges(6,1) * pkin(10) - mrSges(7,3) * pkin(5)) * t204 + (mrSges(6,2) * pkin(10) - t321) * t201) * qJD(5); 0; m(7) * t5 + t14; t315; m(7) * t260; t256 * t304 + t143; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t16(1) t16(2) t16(4) t16(7) t16(11) t16(16); t16(2) t16(3) t16(5) t16(8) t16(12) t16(17); t16(4) t16(5) t16(6) t16(9) t16(13) t16(18); t16(7) t16(8) t16(9) t16(10) t16(14) t16(19); t16(11) t16(12) t16(13) t16(14) t16(15) t16(20); t16(16) t16(17) t16(18) t16(19) t16(20) t16(21);];
Mq  = res;
