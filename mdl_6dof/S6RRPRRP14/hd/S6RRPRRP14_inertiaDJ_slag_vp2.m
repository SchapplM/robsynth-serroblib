% Calculate time derivative of joint inertia matrix for
% S6RRPRRP14
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
% Datum: 2019-03-09 13:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRP14_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP14_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP14_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP14_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP14_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP14_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP14_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:03:49
% EndTime: 2019-03-09 13:04:04
% DurationCPUTime: 6.70s
% Computational Cost: add. (5469->639), mult. (14037->867), div. (0->0), fcn. (12281->8), ass. (0->259)
t196 = sin(qJ(5));
t319 = Ifges(7,4) + Ifges(6,5);
t320 = t319 * t196;
t199 = cos(qJ(5));
t265 = t196 ^ 2 + t199 ^ 2;
t197 = sin(qJ(4));
t200 = cos(qJ(4));
t153 = pkin(4) * t197 - pkin(10) * t200 + qJ(3);
t202 = -pkin(2) - pkin(9);
t280 = t197 * t202;
t315 = t196 * t153 + t199 * t280;
t318 = qJD(5) * t315;
t194 = sin(pkin(6));
t201 = cos(qJ(2));
t263 = qJD(2) * t201;
t240 = t194 * t263;
t195 = cos(pkin(6));
t282 = t194 * t201;
t128 = t195 * t197 + t200 * t282;
t198 = sin(qJ(2));
t264 = qJD(2) * t194;
t241 = t198 * t264;
t94 = -qJD(4) * t128 + t197 * t241;
t244 = t197 * t282;
t129 = t195 * t200 - t244;
t283 = t194 * t198;
t97 = t129 * t199 + t196 * t283;
t42 = qJD(5) * t97 + t196 * t94 - t199 * t240;
t210 = -t129 * t196 + t199 * t283;
t43 = qJD(5) * t210 + t196 * t240 + t199 * t94;
t258 = qJD(4) * t200;
t95 = -qJD(4) * t244 + t195 * t258 - t200 * t241;
t7 = Ifges(6,5) * t43 - Ifges(6,6) * t42 + Ifges(6,3) * t95;
t8 = Ifges(7,4) * t43 + Ifges(7,2) * t95 + Ifges(7,6) * t42;
t317 = t7 + t8;
t136 = qJD(3) + (pkin(4) * t200 + pkin(10) * t197) * qJD(4);
t314 = t199 * t136 - t318;
t285 = qJ(3) * t198;
t101 = (t201 * t202 - pkin(1) - t285) * t194;
t184 = t195 * t198 * pkin(1);
t305 = pkin(3) + pkin(8);
t103 = (t282 * t305 + t184) * qJD(2);
t260 = qJD(4) * t197;
t174 = pkin(2) * t241;
t262 = qJD(3) * t198;
t81 = t174 + (-t262 + (pkin(9) * t198 - qJ(3) * t201) * qJD(2)) * t194;
t181 = pkin(8) * t283;
t303 = pkin(1) * t201;
t242 = -pkin(2) - t303;
t86 = pkin(3) * t283 + t181 + (-pkin(9) + t242) * t195;
t21 = -t101 * t260 + t197 * t103 + t200 * t81 + t86 * t258;
t17 = pkin(10) * t240 + t21;
t252 = t195 * t303;
t175 = qJD(2) * t252;
t190 = t195 * qJD(3);
t85 = -t241 * t305 + t175 + t190;
t29 = pkin(4) * t95 - pkin(10) * t94 + t85;
t50 = t200 * t101 + t197 * t86;
t45 = pkin(10) * t283 + t50;
t131 = pkin(8) * t282 + t184;
t115 = -t195 * qJ(3) - t131;
t100 = pkin(3) * t282 - t115;
t54 = pkin(4) * t128 - pkin(10) * t129 + t100;
t300 = t196 * t54 + t199 * t45;
t4 = -qJD(5) * t300 - t17 * t196 + t199 * t29;
t215 = pkin(5) * t199 + qJ(6) * t196;
t254 = qJD(6) * t199;
t313 = qJD(5) * t215 - t254;
t312 = 2 * m(6);
t311 = 2 * m(7);
t310 = -0.2e1 * pkin(1);
t309 = 2 * mrSges(4,1);
t308 = -2 * mrSges(3,3);
t216 = -pkin(2) * t201 - t285;
t116 = (-pkin(1) + t216) * t194;
t307 = -0.2e1 * t116;
t306 = m(6) * pkin(4);
t304 = m(6) * t197;
t301 = Ifges(4,6) + Ifges(3,4);
t60 = mrSges(7,2) * t210 + mrSges(7,3) * t128;
t61 = -mrSges(6,2) * t128 + mrSges(6,3) * t210;
t299 = t60 + t61;
t62 = mrSges(6,1) * t128 - mrSges(6,3) * t97;
t63 = -mrSges(7,1) * t128 + mrSges(7,2) * t97;
t298 = -t62 + t63;
t15 = mrSges(6,1) * t42 + mrSges(6,2) * t43;
t67 = mrSges(5,1) * t240 - mrSges(5,3) * t94;
t297 = t67 - t15;
t296 = Ifges(5,4) * t197;
t295 = Ifges(5,4) * t200;
t294 = Ifges(6,4) * t196;
t293 = Ifges(6,4) * t199;
t292 = Ifges(7,5) * t196;
t291 = Ifges(7,5) * t199;
t290 = Ifges(6,6) * t199;
t289 = Ifges(7,6) * t196;
t124 = -pkin(8) * t241 + t175;
t288 = t124 * mrSges(3,2);
t106 = mrSges(5,1) * t283 - mrSges(5,3) * t129;
t56 = -mrSges(6,1) * t210 + mrSges(6,2) * t97;
t287 = -t106 + t56;
t160 = -t199 * mrSges(6,1) + t196 * mrSges(6,2);
t286 = t160 - mrSges(5,1);
t125 = t131 * qJD(2);
t284 = t125 * t198;
t281 = t196 * t200;
t278 = t199 * t153;
t277 = t199 * t200;
t276 = t200 * t202;
t255 = qJD(5) * t200;
t259 = qJD(4) * t199;
t208 = t196 * t255 + t197 * t259;
t107 = mrSges(6,1) * t258 + mrSges(6,3) * t208;
t108 = -mrSges(7,1) * t258 - mrSges(7,2) * t208;
t275 = -t107 + t108;
t236 = t196 * t260;
t237 = t199 * t255;
t207 = t236 - t237;
t109 = -mrSges(6,2) * t258 + mrSges(6,3) * t207;
t110 = mrSges(7,2) * t207 + mrSges(7,3) * t258;
t274 = t109 + t110;
t217 = Ifges(7,3) * t196 + t291;
t118 = Ifges(7,6) * t197 + t200 * t217;
t218 = -Ifges(6,2) * t196 + t293;
t121 = Ifges(6,6) * t197 + t200 * t218;
t273 = t118 - t121;
t219 = Ifges(7,1) * t199 + t292;
t122 = Ifges(7,4) * t197 + t200 * t219;
t220 = Ifges(6,1) * t199 - t294;
t123 = Ifges(6,5) * t197 + t200 * t220;
t272 = -t122 - t123;
t149 = -mrSges(6,2) * t197 - mrSges(6,3) * t281;
t152 = -mrSges(7,2) * t281 + mrSges(7,3) * t197;
t271 = t149 + t152;
t150 = mrSges(6,1) * t197 - mrSges(6,3) * t277;
t151 = -mrSges(7,1) * t197 + mrSges(7,2) * t277;
t270 = -t150 + t151;
t269 = Ifges(3,5) * t240 + Ifges(4,5) * t241;
t268 = Ifges(7,2) * t258 + Ifges(7,6) * t237;
t267 = Ifges(6,6) * t236 + Ifges(6,3) * t258;
t266 = t265 * pkin(10) * t258;
t256 = qJD(5) * t199;
t257 = qJD(5) * t196;
t143 = Ifges(7,4) * t256 + Ifges(7,6) * t257;
t261 = qJD(4) * t196;
t6 = Ifges(7,5) * t43 + Ifges(7,6) * t95 + Ifges(7,3) * t42;
t9 = Ifges(6,4) * t43 - Ifges(6,2) * t42 + Ifges(6,6) * t95;
t253 = t6 / 0.2e1 - t9 / 0.2e1;
t251 = Ifges(5,5) * t94 - Ifges(5,6) * t95 + Ifges(5,3) * t240;
t10 = Ifges(7,1) * t43 + Ifges(7,4) * t95 + Ifges(7,5) * t42;
t11 = Ifges(6,1) * t43 - Ifges(6,4) * t42 + Ifges(6,5) * t95;
t249 = t10 / 0.2e1 + t11 / 0.2e1;
t30 = Ifges(7,5) * t97 + Ifges(7,6) * t128 - Ifges(7,3) * t210;
t33 = Ifges(6,4) * t97 + Ifges(6,2) * t210 + Ifges(6,6) * t128;
t248 = t30 / 0.2e1 - t33 / 0.2e1;
t34 = Ifges(7,1) * t97 + Ifges(7,4) * t128 - Ifges(7,5) * t210;
t35 = Ifges(6,1) * t97 + Ifges(6,4) * t210 + Ifges(6,5) * t128;
t247 = -t34 / 0.2e1 - t35 / 0.2e1;
t162 = -Ifges(7,3) * t199 + t292;
t69 = -t162 * t255 + (Ifges(7,6) * t200 - t197 * t217) * qJD(4);
t165 = Ifges(6,2) * t199 + t294;
t72 = -t165 * t255 + (Ifges(6,6) * t200 - t197 * t218) * qJD(4);
t246 = -t72 / 0.2e1 + t69 / 0.2e1;
t167 = Ifges(7,1) * t196 - t291;
t73 = -t167 * t255 + (Ifges(7,4) * t200 - t197 * t219) * qJD(4);
t168 = Ifges(6,1) * t196 + t293;
t74 = -t168 * t255 + (Ifges(6,5) * t200 - t197 * t220) * qJD(4);
t245 = t73 / 0.2e1 + t74 / 0.2e1;
t239 = t202 * t258;
t243 = t196 * t136 + t153 * t256 + t199 * t239;
t238 = t202 * t257;
t235 = -t283 / 0.2e1;
t234 = t118 / 0.2e1 - t121 / 0.2e1;
t233 = t122 / 0.2e1 + t123 / 0.2e1;
t141 = t217 * qJD(5);
t144 = t218 * qJD(5);
t232 = t141 / 0.2e1 - t144 / 0.2e1;
t142 = Ifges(6,5) * t256 - Ifges(6,6) * t257;
t231 = t142 / 0.2e1 + t143 / 0.2e1;
t146 = t219 * qJD(5);
t147 = t220 * qJD(5);
t230 = t146 / 0.2e1 + t147 / 0.2e1;
t229 = t162 / 0.2e1 - t165 / 0.2e1;
t228 = t290 / 0.2e1 - Ifges(7,6) * t199 / 0.2e1 + t320 / 0.2e1;
t227 = t167 / 0.2e1 + t168 / 0.2e1;
t26 = -t95 * mrSges(7,1) + t43 * mrSges(7,2);
t226 = (mrSges(4,2) - mrSges(3,1)) * t125;
t225 = t196 * t202 - pkin(5);
t49 = -t197 * t101 + t200 * t86;
t224 = Ifges(5,5) * t240;
t223 = Ifges(5,6) * t240;
t222 = t196 * mrSges(6,1) + t199 * mrSges(6,2);
t159 = -t199 * mrSges(7,1) - t196 * mrSges(7,3);
t221 = t196 * mrSges(7,1) - t199 * mrSges(7,3);
t214 = pkin(5) * t196 - qJ(6) * t199;
t19 = -t196 * t45 + t199 * t54;
t22 = -t101 * t258 + t103 * t200 - t197 * t81 - t86 * t260;
t211 = -t202 + t214;
t3 = t199 * t17 + t196 * t29 + t54 * t256 - t257 * t45;
t44 = -pkin(4) * t283 - t49;
t18 = -pkin(4) * t240 - t22;
t1 = qJ(6) * t95 + qJD(6) * t128 + t3;
t12 = qJ(6) * t128 + t300;
t13 = -pkin(5) * t128 - t19;
t2 = -pkin(5) * t95 - t4;
t206 = t1 * t199 - t12 * t257 + t13 * t256 + t196 * t2;
t205 = -t19 * t256 - t196 * t4 + t199 * t3 - t257 * t300;
t24 = -mrSges(6,2) * t95 - mrSges(6,3) * t42;
t25 = mrSges(6,1) * t95 - mrSges(6,3) * t43;
t27 = -mrSges(7,2) * t42 + mrSges(7,3) * t95;
t204 = (t24 + t27) * t199 + (-t25 + t26) * t196 + (-t196 * t299 + t199 * t298) * qJD(5);
t203 = m(7) * t254 + (-m(7) * t215 + t159 + t160) * qJD(5);
t169 = Ifges(5,1) * t200 - t296;
t166 = -Ifges(5,2) * t197 + t295;
t161 = mrSges(5,1) * t197 + mrSges(5,2) * t200;
t154 = -pkin(4) - t215;
t148 = (-Ifges(5,1) * t197 - t295) * qJD(4);
t145 = (-Ifges(5,2) * t200 - t296) * qJD(4);
t140 = (mrSges(5,1) * t200 - mrSges(5,2) * t197) * qJD(4);
t139 = t222 * qJD(5);
t138 = t221 * qJD(5);
t137 = -mrSges(4,1) * t282 - mrSges(4,3) * t195;
t134 = t222 * t200;
t133 = t221 * t200;
t130 = -t181 + t252;
t126 = qJD(5) * t214 - qJD(6) * t196;
t120 = Ifges(7,2) * t197 + (Ifges(7,4) * t199 + t289) * t200;
t119 = Ifges(6,3) * t197 + (Ifges(6,5) * t199 - Ifges(6,6) * t196) * t200;
t117 = t195 * t242 + t181;
t114 = t211 * t200;
t112 = -t196 * t280 + t278;
t111 = -t124 - t190;
t105 = -mrSges(5,2) * t283 - mrSges(5,3) * t128;
t104 = t174 + (-qJ(3) * t263 - t262) * t194;
t102 = t197 * t225 - t278;
t99 = qJ(6) * t197 + t315;
t83 = -mrSges(6,1) * t207 - mrSges(6,2) * t208;
t82 = -mrSges(7,1) * t207 + mrSges(7,3) * t208;
t77 = mrSges(5,1) * t128 + mrSges(5,2) * t129;
t71 = -Ifges(7,4) * t208 - Ifges(7,6) * t236 + t268;
t70 = -Ifges(6,5) * t208 - Ifges(6,6) * t237 + t267;
t68 = -mrSges(5,2) * t240 - mrSges(5,3) * t95;
t66 = Ifges(5,1) * t129 - Ifges(5,4) * t128 + Ifges(5,5) * t283;
t65 = Ifges(5,4) * t129 - Ifges(5,2) * t128 + Ifges(5,6) * t283;
t64 = t313 * t200 - t211 * t260;
t59 = -t196 * t239 + t314;
t58 = -t197 * t238 + t243;
t57 = t225 * t258 - t314;
t55 = -mrSges(7,1) * t210 - mrSges(7,3) * t97;
t53 = mrSges(5,1) * t95 + mrSges(5,2) * t94;
t52 = qJ(6) * t258 + (qJD(6) - t238) * t197 + t243;
t47 = Ifges(5,1) * t94 - Ifges(5,4) * t95 + t224;
t46 = Ifges(5,4) * t94 - Ifges(5,2) * t95 + t223;
t32 = Ifges(7,4) * t97 + Ifges(7,2) * t128 - Ifges(7,6) * t210;
t31 = Ifges(6,5) * t97 + Ifges(6,6) * t210 + Ifges(6,3) * t128;
t23 = -pkin(5) * t210 - qJ(6) * t97 + t44;
t14 = mrSges(7,1) * t42 - mrSges(7,3) * t43;
t5 = pkin(5) * t42 - qJ(6) * t43 - qJD(6) * t97 + t18;
t16 = [(t30 - t33) * t42 - (-t9 + t6) * t210 + 0.2e1 * t111 * t137 + t129 * t47 + 0.2e1 * t21 * t105 + 0.2e1 * t22 * t106 + 0.2e1 * t100 * t53 + t94 * t66 + 0.2e1 * t85 * t77 + 0.2e1 * t49 * t67 + 0.2e1 * t50 * t68 + 0.2e1 * t1 * t60 + 0.2e1 * t3 * t61 + 0.2e1 * t4 * t62 + 0.2e1 * t2 * t63 + 0.2e1 * t5 * t55 + 0.2e1 * t18 * t56 + 0.2e1 * t44 * t15 + 0.2e1 * t19 * t25 + 0.2e1 * t13 * t26 + 0.2e1 * t12 * t27 + 0.2e1 * t23 * t14 + (t34 + t35) * t43 + (t1 * t12 + t13 * t2 + t23 * t5) * t311 + (t198 * t251 + 0.2e1 * t104 * (mrSges(4,2) * t201 - mrSges(4,3) * t198) + t284 * t309 + 0.2e1 * (t124 * t201 + t284) * mrSges(3,3) + ((t115 * t309 + mrSges(4,2) * t307 + t131 * t308 + (Ifges(4,5) - (2 * Ifges(3,6))) * t195 + (mrSges(3,1) * t310 - 0.2e1 * t198 * t301) * t194) * t198 + (t117 * t309 + t130 * t308 + mrSges(4,3) * t307 + Ifges(5,5) * t129 - Ifges(5,6) * t128 + (-(2 * Ifges(4,4)) + Ifges(3,5)) * t195 + (mrSges(3,2) * t310 + 0.2e1 * t201 * t301) * t194 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) + (2 * Ifges(4,2)) - (2 * Ifges(4,3)) + Ifges(5,3)) * t283) * t201) * qJD(2)) * t194 + (0.2e1 * t226 + t269 - 0.2e1 * t288) * t195 + 0.2e1 * t300 * t24 + (t18 * t44 + t19 * t4 + t3 * t300) * t312 + (t31 + t32 - t65) * t95 + (t11 + t10) * t97 + 0.2e1 * m(3) * (t124 * t131 - t125 * t130) + 0.2e1 * m(4) * (t104 * t116 + t111 * t115 + t117 * t125) + 0.2e1 * m(5) * (t100 * t85 + t21 * t50 + t22 * t49) + (-t46 + t317) * t128; m(4) * (-pkin(2) * t125 - qJ(3) * t111 - qJD(3) * t115) + m(7) * (t1 * t99 + t102 * t2 + t114 * t5 + t12 * t52 + t13 * t57 + t23 * t64) - t246 * t210 + t315 * t24 + m(6) * (t112 * t4 - t18 * t276 + t59 * t19 + t3 * t315 + t300 * t58) - t288 + t1 * t152 + t85 * t161 + t94 * t169 / 0.2e1 + t5 * t133 + t18 * t134 + m(5) * (qJ(3) * t85 + qJD(3) * t100 + t21 * t280 + t22 * t276) + t100 * t140 + t129 * t148 / 0.2e1 + t3 * t149 + t4 * t150 + t2 * t151 + t19 * t107 + t13 * t108 + t12 * t110 - t111 * mrSges(4,3) + t112 * t25 + t114 * t14 + t102 * t26 + t99 * t27 + t23 * t82 + t44 * t83 + t52 * t60 + t58 * t61 + t59 * t62 + t57 * t63 + t64 * t55 + qJ(3) * t53 + t245 * t97 + t233 * t43 + t234 * t42 + (t202 * t68 - t21 * mrSges(5,3) - t223 / 0.2e1 - t46 / 0.2e1 + t7 / 0.2e1 + t8 / 0.2e1) * t197 + (-t145 / 0.2e1 + t70 / 0.2e1 + t71 / 0.2e1) * t128 + (-t166 / 0.2e1 + t119 / 0.2e1 + t120 / 0.2e1) * t95 + (mrSges(4,1) * t216 - Ifges(4,4) * t201 - Ifges(3,6) * t198) * t264 + t300 * t109 + (-t137 + t77) * qJD(3) + (-t22 * mrSges(5,3) + t224 / 0.2e1 + t47 / 0.2e1 + t297 * t202 + t249 * t199 + t253 * t196 + (t196 * t247 + t199 * t248) * qJD(5)) * t200 + t226 + ((Ifges(5,6) * t235 - t65 / 0.2e1 + t31 / 0.2e1 + t32 / 0.2e1 - t50 * mrSges(5,3)) * t200 + (Ifges(5,5) * t235 - t66 / 0.2e1 + t49 * mrSges(5,3) + t247 * t199 - t248 * t196) * t197 + (t200 * t105 + t287 * t197 + m(5) * (-t197 * t49 + t200 * t50) + t44 * t304) * t202) * qJD(4) + t269; 0.2e1 * qJ(3) * t140 + 0.2e1 * t102 * t108 + 0.2e1 * t112 * t107 + 0.2e1 * t315 * t109 + 0.2e1 * t99 * t110 + 0.2e1 * t114 * t82 + 0.2e1 * t64 * t133 + 0.2e1 * t58 * t149 + 0.2e1 * t59 * t150 + 0.2e1 * t57 * t151 + 0.2e1 * t52 * t152 + (t112 * t59 + t315 * t58) * t312 + (t102 * t57 + t114 * t64 + t52 * t99) * t311 + 0.2e1 * (mrSges(4,3) + t161 + (m(4) + m(5)) * qJ(3)) * qJD(3) + (-t145 + t70 + t71 + (0.2e1 * t134 * t202 - t196 * t273 + t199 * t272 - t169) * qJD(4)) * t197 + (-0.2e1 * t202 * t83 + t148 + (t73 + t74) * t199 + (t69 - t72) * t196 + (t196 * t272 + t199 * t273) * qJD(5) + (-0.2e1 * t202 ^ 2 * t304 + t119 + t120 - t166) * qJD(4)) * t200; m(4) * t125 + mrSges(4,1) * t240 + (-t14 + (t196 * t298 + t199 * t299 + t105) * qJD(4) + m(6) * (-t19 * t261 + t259 * t300 - t18) + m(7) * (t12 * t259 + t13 * t261 - t5) + m(5) * (qJD(4) * t50 + t22) + t297) * t200 + (t68 + (t55 + t287) * qJD(4) + m(6) * (qJD(4) * t44 + t205) + m(7) * (qJD(4) * t23 + t206) + m(5) * (-qJD(4) * t49 + t21) + t204) * t197; (-m(7) * t64 - t82 - t83 + (t271 * t199 + t270 * t196 + m(7) * (t102 * t196 + t199 * t99) + m(6) * (-t112 * t196 + t199 * t315)) * qJD(4)) * t200 + (t274 * t199 + t275 * t196 + (t133 + t134) * qJD(4) + (-t196 * t271 + t199 * t270) * qJD(5) + m(7) * (qJD(4) * t114 + t102 * t256 + t196 * t57 + t199 * t52 - t257 * t99) + m(6) * (-t112 * t256 - t196 * t59 + t199 * t58 - t257 * t315 - 0.2e1 * t239)) * t197; 0.4e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * (-0.1e1 + t265) * t197 * t258; t154 * t14 + t5 * t159 + t23 * t138 + t44 * t139 + t126 * t55 + ((t13 * mrSges(7,2) - t19 * mrSges(6,3) - t247) * t199 + (-t12 * mrSges(7,2) - mrSges(6,3) * t300 + t248) * t196) * qJD(5) + (t2 * mrSges(7,2) - t4 * mrSges(6,3) + t249) * t196 - t21 * mrSges(5,2) + t22 * mrSges(5,1) - pkin(4) * t15 + t228 * t95 + t229 * t42 + t230 * t97 + t231 * t128 - t232 * t210 + t227 * t43 + (m(6) * t205 + m(7) * t206 + t204) * pkin(10) + t251 + m(7) * (t126 * t23 + t154 * t5) + (t1 * mrSges(7,2) + t3 * mrSges(6,3) - t253) * t199 + (t160 - t306) * t18; t154 * t82 + t64 * t159 + t126 * t133 + t114 * t138 - pkin(4) * t83 + m(7) * (t114 * t126 + t154 * t64) - t139 * t276 + t231 * t197 + ((-t202 * mrSges(5,2) - Ifges(5,6) + t228) * t200 + (-Ifges(5,5) + (t286 - t306) * t202) * t197) * qJD(4) + (t58 * mrSges(6,3) + t52 * mrSges(7,2) + t230 * t200 - t227 * t260 + (t102 * mrSges(7,2) - t112 * mrSges(6,3) + t229 * t200 + t233) * qJD(5) + (t270 * qJD(5) + m(7) * (qJD(5) * t102 + t52) + m(6) * (-qJD(5) * t112 + t58) + t274) * pkin(10) - t246) * t199 + (-t59 * mrSges(6,3) + t57 * mrSges(7,2) + t232 * t200 - t229 * t260 + (-t99 * mrSges(7,2) - mrSges(6,3) * t315 - t227 * t200 + t234) * qJD(5) + (-t271 * qJD(5) + m(7) * (-qJD(5) * t99 + t57) + m(6) * (-t59 - t318) + t275) * pkin(10) + t245) * t196; (t159 + t286) * t260 + m(6) * (-pkin(4) * t260 + t266) + m(7) * (t154 * t260 + t266) + (-t139 - t138 - m(7) * t126 + (-mrSges(5,2) + (mrSges(7,2) + mrSges(6,3)) * t265) * qJD(4)) * t200; -0.2e1 * pkin(4) * t139 + 0.2e1 * t138 * t154 + (-t141 + t144) * t199 + (t146 + t147) * t196 + 0.2e1 * (m(7) * t154 + t159) * t126 + ((t167 + t168) * t199 + (t162 - t165) * t196) * qJD(5); -pkin(5) * t26 + m(7) * (-pkin(5) * t2 + qJ(6) * t1 + qJD(6) * t12) + qJD(6) * t60 + qJ(6) * t27 + t1 * mrSges(7,3) - t2 * mrSges(7,1) - t3 * mrSges(6,2) + t4 * mrSges(6,1) + t317; -pkin(5) * t108 + m(7) * (-pkin(5) * t57 + qJ(6) * t52 + qJD(6) * t99) + qJD(6) * t152 + qJ(6) * t110 + t52 * mrSges(7,3) - t58 * mrSges(6,2) + t59 * mrSges(6,1) - t57 * mrSges(7,1) + (-t290 - t320) * t255 + (-t199 * t319 - t289) * t260 + t267 + t268; (-m(7) * t214 - t221 - t222) * t258 + t203 * t197; -t313 * mrSges(7,2) + t203 * pkin(10) + t142 + t143; 0.2e1 * (m(7) * qJ(6) + mrSges(7,3)) * qJD(6); m(7) * t2 + t26; m(7) * t57 + t108; (t196 * t258 + t197 * t256) * m(7); (m(7) * pkin(10) + mrSges(7,2)) * t256; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t16(1) t16(2) t16(4) t16(7) t16(11) t16(16); t16(2) t16(3) t16(5) t16(8) t16(12) t16(17); t16(4) t16(5) t16(6) t16(9) t16(13) t16(18); t16(7) t16(8) t16(9) t16(10) t16(14) t16(19); t16(11) t16(12) t16(13) t16(14) t16(15) t16(20); t16(16) t16(17) t16(18) t16(19) t16(20) t16(21);];
Mq  = res;
