% Calculate time derivative of joint inertia matrix for
% S6RPPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3,theta4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPPRR3_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR3_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR3_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR3_inertiaDJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR3_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPPRR3_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPPRR3_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:33:28
% EndTime: 2019-03-09 01:33:38
% DurationCPUTime: 6.71s
% Computational Cost: add. (13935->472), mult. (24269->700), div. (0->0), fcn. (27213->10), ass. (0->220)
t264 = sin(pkin(9));
t265 = cos(pkin(9));
t276 = sin(qJ(1));
t277 = cos(qJ(1));
t154 = t277 * t264 - t276 * t265;
t153 = -t264 * t276 - t265 * t277;
t170 = pkin(10) + qJ(5);
t162 = sin(t170);
t163 = cos(t170);
t202 = -Icges(6,5) * t163 + Icges(6,6) * t162;
t90 = -Icges(6,3) * t153 + t154 * t202;
t298 = t154 * t90;
t91 = Icges(6,3) * t154 + t153 * t202;
t299 = t153 * t91;
t300 = -t298 + t299;
t285 = t153 ^ 2;
t284 = t154 ^ 2;
t143 = t153 * qJD(1);
t144 = t154 * qJD(1);
t243 = qJD(5) * t162;
t192 = t144 * t163 + t153 * t243;
t191 = t192 * rSges(6,1) + t143 * rSges(6,3);
t242 = qJD(5) * t163;
t256 = t144 * t162;
t193 = t153 * t242 - t256;
t194 = -t143 * t163 + t154 * t243;
t258 = t143 * t162;
t195 = t154 * t242 + t258;
t220 = -rSges(6,1) * t163 + rSges(6,2) * t162;
t96 = -t153 * rSges(6,3) + t154 * t220;
t253 = t153 * t163;
t254 = t153 * t162;
t97 = -rSges(6,1) * t253 + rSges(6,2) * t254 + t154 * rSges(6,3);
t27 = t143 * t96 + t154 * (rSges(6,1) * t194 + t144 * rSges(6,3)) - t144 * t97 + t153 * t191 + (t153 * t193 + t154 * t195) * rSges(6,2);
t287 = 2 * m(6);
t296 = t27 * t287;
t245 = t277 * pkin(1) + t276 * qJ(2);
t235 = t277 * pkin(2) + t245;
t295 = -rSges(3,1) * t277 - rSges(3,3) * t276 - t245;
t174 = sin(qJ(6));
t175 = cos(qJ(6));
t182 = qJD(6) * t154 + t192;
t240 = qJD(6) * t163;
t222 = t153 * t240 + t143;
t61 = -t174 * t182 + t175 * t222;
t62 = t174 * t222 + t175 * t182;
t39 = t62 * rSges(7,1) + t61 * rSges(7,2) - t193 * rSges(7,3);
t294 = t192 * pkin(5) - t193 * pkin(8) + t39;
t272 = pkin(5) * t163;
t224 = pkin(8) * t162 + t272;
t293 = t154 * t224;
t172 = cos(pkin(10));
t160 = pkin(4) * t172 + pkin(3);
t278 = rSges(7,3) + pkin(8);
t288 = t162 * t278 + t160 + t272;
t286 = 2 * m(7);
t283 = t144 / 0.2e1;
t282 = -t153 / 0.2e1;
t281 = -t154 / 0.2e1;
t280 = t154 / 0.2e1;
t279 = t163 / 0.2e1;
t140 = t220 * qJD(5);
t275 = m(6) * t140;
t151 = -rSges(6,1) * t162 - rSges(6,2) * t163;
t274 = m(6) * t151;
t273 = pkin(5) * t162;
t250 = t163 * t174;
t108 = -t153 * t175 + t154 * t250;
t249 = t163 * t175;
t109 = -t153 * t174 - t154 * t249;
t252 = t154 * t162;
t73 = Icges(7,4) * t109 + Icges(7,2) * t108 - Icges(7,6) * t252;
t75 = Icges(7,1) * t109 + Icges(7,4) * t108 - Icges(7,5) * t252;
t211 = t174 * t73 - t175 * t75;
t183 = -qJD(6) * t153 + t194;
t221 = t154 * t240 + t144;
t59 = -t174 * t183 + t175 * t221;
t60 = t174 * t221 + t175 * t183;
t32 = Icges(7,5) * t60 + Icges(7,6) * t59 - Icges(7,3) * t195;
t34 = Icges(7,4) * t60 + Icges(7,2) * t59 - Icges(7,6) * t195;
t36 = Icges(7,1) * t60 + Icges(7,4) * t59 - Icges(7,5) * t195;
t71 = Icges(7,5) * t109 + Icges(7,6) * t108 - Icges(7,3) * t252;
t9 = (qJD(5) * t211 + t32) * t163 + (-qJD(5) * t71 + t174 * t34 - t175 * t36 + (t174 * t75 + t175 * t73) * qJD(6)) * t162;
t271 = t9 * t153;
t110 = t153 * t250 + t154 * t175;
t111 = -t153 * t249 + t154 * t174;
t74 = Icges(7,4) * t111 + Icges(7,2) * t110 - Icges(7,6) * t254;
t76 = Icges(7,1) * t111 + Icges(7,4) * t110 - Icges(7,5) * t254;
t210 = t174 * t74 - t175 * t76;
t33 = Icges(7,5) * t62 + Icges(7,6) * t61 - Icges(7,3) * t193;
t35 = Icges(7,4) * t62 + Icges(7,2) * t61 - Icges(7,6) * t193;
t37 = Icges(7,1) * t62 + Icges(7,4) * t61 - Icges(7,5) * t193;
t72 = Icges(7,5) * t111 + Icges(7,6) * t110 - Icges(7,3) * t254;
t10 = (qJD(5) * t210 + t33) * t163 + (-qJD(5) * t72 + t174 * t35 - t175 * t37 + (t174 * t76 + t175 * t74) * qJD(6)) * t162;
t270 = t10 * t154;
t269 = rSges(5,3) + qJ(4);
t173 = -pkin(7) - qJ(4);
t268 = rSges(6,3) - t173;
t219 = -t109 * rSges(7,1) - t108 * rSges(7,2);
t78 = -rSges(7,3) * t252 - t219;
t267 = t78 - t293;
t79 = t111 * rSges(7,1) + t110 * rSges(7,2) - rSges(7,3) * t254;
t266 = -pkin(5) * t253 - pkin(8) * t254 + t79;
t263 = Icges(6,4) * t162;
t262 = Icges(6,4) * t163;
t261 = Icges(7,4) * t174;
t260 = Icges(7,4) * t175;
t203 = Icges(7,2) * t174 - t260;
t118 = Icges(7,6) * t163 + t162 * t203;
t259 = t118 * t174;
t201 = -Icges(7,5) * t175 + Icges(7,6) * t174;
t117 = Icges(7,3) * t163 + t162 * t201;
t251 = t162 * t117;
t218 = -rSges(7,1) * t175 + rSges(7,2) * t174;
t241 = qJD(6) * t162;
t101 = (rSges(7,1) * t174 + rSges(7,2) * t175) * t241 + (-rSges(7,3) * t162 + t163 * t218) * qJD(5);
t248 = t224 * qJD(5) - t101;
t120 = t163 * rSges(7,3) + t162 * t218;
t223 = -pkin(8) * t163 + t273;
t247 = t120 - t223;
t167 = t277 * qJ(2);
t246 = qJD(1) * t167 + qJD(2) * t276;
t244 = qJD(5) * t154;
t239 = t276 * pkin(1);
t30 = t162 * t211 + t163 * t71;
t205 = -Icges(7,1) * t175 + t261;
t119 = Icges(7,5) * t163 + t162 * t205;
t49 = t108 * t118 + t109 * t119 - t154 * t251;
t238 = t30 / 0.2e1 + t49 / 0.2e1;
t31 = t162 * t210 + t163 * t72;
t50 = t110 * t118 + t111 * t119 - t153 * t251;
t237 = t31 / 0.2e1 + t50 / 0.2e1;
t236 = t119 * t249;
t228 = m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1;
t225 = t60 * rSges(7,1) + t59 * rSges(7,2);
t23 = t108 * t73 + t109 * t75 - t252 * t71;
t24 = t108 * t74 + t109 * t76 - t252 * t72;
t217 = -t153 * t24 - t154 * t23;
t25 = t110 * t73 + t111 * t75 - t254 * t71;
t26 = t110 * t74 + t111 * t76 - t254 * t72;
t216 = -t153 * t26 - t154 * t25;
t215 = -t153 * t31 - t154 * t30;
t214 = -t153 * t78 + t154 * t79;
t207 = -t153 * t160 - t154 * t173 + t235;
t206 = -Icges(6,1) * t163 + t263;
t204 = Icges(6,2) * t162 - t262;
t200 = t143 * t154 - t144 * t153;
t199 = rSges(5,1) * t172 - rSges(5,2) * sin(pkin(10)) + pkin(3);
t198 = t160 - t220;
t197 = -t162 * t32 - t242 * t71;
t196 = -t162 * t33 - t242 * t72;
t190 = -pkin(2) * t276 - t239;
t187 = -t153 * t276 + t154 * t277;
t185 = t167 + t190;
t100 = (Icges(7,1) * t174 + t260) * t241 + (-Icges(7,5) * t162 + t163 * t205) * qJD(5);
t98 = (Icges(7,5) * t174 + Icges(7,6) * t175) * t241 + (-Icges(7,3) * t162 + t163 * t201) * qJD(5);
t99 = (Icges(7,2) * t175 + t261) * t241 + (-Icges(7,6) * t162 + t163 * t203) * qJD(5);
t184 = t163 * t98 + t242 * t259 + (t118 * t175 + t119 * t174) * t241 + (-t100 * t175 + t174 * t99) * t162;
t181 = -rSges(3,1) * t276 + rSges(3,3) * t277 - t239;
t180 = qJD(1) * t190 + t246;
t165 = qJD(2) * t277;
t179 = -t235 * qJD(1) + t165;
t178 = t154 * qJD(4) + t180;
t177 = t153 * qJD(4) + t179;
t176 = -t143 * t173 + t144 * t160 + t178;
t134 = t167 + t181;
t116 = t295 * qJD(1) + t165;
t115 = qJD(1) * t181 + t246;
t103 = -rSges(4,1) * t153 - rSges(4,2) * t154 + t235;
t102 = t154 * rSges(4,1) - t153 * rSges(4,2) + t185;
t95 = Icges(6,5) * t154 + t153 * t206;
t94 = -Icges(6,5) * t153 + t154 * t206;
t93 = Icges(6,6) * t154 + t153 * t204;
t92 = -Icges(6,6) * t153 + t154 * t204;
t88 = t143 * rSges(4,1) + t144 * rSges(4,2) + t179;
t87 = t144 * rSges(4,1) - t143 * rSges(4,2) + t180;
t84 = t247 * t153;
t83 = t247 * t154;
t81 = -t153 * t199 + t154 * t269 + t235;
t80 = t153 * t269 + t154 * t199 + t185;
t77 = t163 * t117 + (-t119 * t175 + t259) * t162;
t70 = t97 + t207;
t69 = t153 * t268 + t154 * t198 + t185;
t64 = Icges(6,5) * t192 + Icges(6,6) * t193 + Icges(6,3) * t143;
t63 = Icges(6,5) * t194 + Icges(6,6) * t195 + Icges(6,3) * t144;
t56 = t143 * t199 - t144 * t269 + t177;
t55 = t143 * t269 + t144 * t199 + t178;
t54 = t120 * t254 + t163 * t79;
t53 = -t120 * t252 - t163 * t78;
t52 = t144 * t247 + t153 * t248;
t51 = -t143 * t247 + t154 * t248;
t48 = t143 * t198 - t144 * t268 + t151 * t244 + t177;
t47 = rSges(6,2) * t193 + t176 + t191;
t46 = t207 + t266;
t45 = -t153 * t173 + t154 * t288 + t185 + t219;
t40 = t214 * t162;
t38 = -rSges(7,3) * t195 + t225;
t29 = t153 * t266 + t154 * t267;
t28 = ((-t236 - t251) * qJD(5) + t184) * t163;
t22 = t144 * t173 + (t163 * t278 - t273) * t244 + t288 * t143 + t177 - t225;
t21 = t176 + t294;
t20 = (qJD(5) * t120 * t153 + t39) * t163 + (-qJD(5) * t79 + t101 * t153 - t120 * t144) * t162;
t19 = (-t120 * t244 - t38) * t163 + (qJD(5) * t78 - t101 * t154 - t120 * t143) * t162;
t18 = t100 * t111 + t110 * t99 - t117 * t193 + t118 * t61 + t119 * t62 - t254 * t98;
t17 = t100 * t109 + t108 * t99 - t117 * t195 + t118 * t59 + t119 * t60 - t252 * t98;
t16 = -t153 * t25 + t154 * t26;
t15 = -t153 * t23 + t154 * t24;
t14 = t162 * t216 + t163 * t50;
t13 = t162 * t217 + t163 * t49;
t12 = t214 * t242 + (t143 * t79 + t144 * t78 - t153 * t38 + t154 * t39) * t162;
t11 = t294 * t153 - t266 * t144 + (t223 * t244 + t38) * t154 - (-t267 + t293) * t143;
t8 = t110 * t35 + t111 * t37 + t153 * t196 + t256 * t72 + t61 * t74 + t62 * t76;
t7 = t110 * t34 + t111 * t36 + t153 * t197 + t256 * t71 + t61 * t73 + t62 * t75;
t6 = t108 * t35 + t109 * t37 + t154 * t196 - t258 * t72 + t59 * t74 + t60 * t76;
t5 = t108 * t34 + t109 * t36 + t154 * t197 - t258 * t71 + t59 * t73 + t60 * t75;
t4 = t143 * t26 + t144 * t25 - t153 * t7 + t154 * t8;
t3 = t143 * t24 + t144 * t23 - t153 * t5 + t154 * t6;
t2 = (qJD(5) * t216 + t18) * t163 + (-qJD(5) * t50 - t143 * t25 + t144 * t26 - t153 * t8 - t154 * t7) * t162;
t1 = (qJD(5) * t217 + t17) * t163 + (-qJD(5) * t49 - t143 * t23 + t144 * t24 - t153 * t6 - t154 * t5) * t162;
t41 = [(t21 * t46 + t22 * t45) * t286 - qJD(5) * t236 + (t47 * t70 + t48 * t69) * t287 + 0.2e1 * m(4) * (t102 * t88 + t103 * t87) + 0.2e1 * m(5) * (t55 * t81 + t56 * t80) + 0.2e1 * m(3) * (-t115 * t295 + t116 * t134) + t184 + (Icges(6,1) * t162 - t204 + t262) * t242 + (-Icges(6,2) * t163 - t117 - t206 - t263) * t243; m(7) * (t276 * t22 - t277 * t21 + (t276 * t46 + t277 * t45) * qJD(1)) + m(6) * (t276 * t48 - t277 * t47 + (t276 * t70 + t277 * t69) * qJD(1)) + m(4) * (t276 * t88 - t277 * t87 + (t102 * t277 + t103 * t276) * qJD(1)) + m(5) * (t276 * t56 - t277 * t55 + (t276 * t81 + t277 * t80) * qJD(1)) + m(3) * (t276 * t116 - t277 * t115 + (t134 * t277 - t276 * t295) * qJD(1)); 0; 0; 0; 0; m(7) * (t143 * t45 + t144 * t46 - t153 * t21 + t154 * t22) + m(6) * (t143 * t69 + t144 * t70 - t153 * t47 + t154 * t48) + m(5) * (t143 * t80 + t144 * t81 - t153 * t55 + t154 * t56); 0.2e1 * t228 * (qJD(1) * t187 + t143 * t276 - t144 * t277); 0; 0.4e1 * t228 * t200; -t271 / 0.2e1 + t270 / 0.2e1 + (t285 / 0.2e1 + t284 / 0.2e1) * t202 * qJD(5) + (-t162 * t94 / 0.2e1 - t163 * t92 / 0.2e1 + t238) * t144 - (t162 * t95 / 0.2e1 + t93 * t279 - t237) * t143 + m(7) * (-t21 * t83 - t22 * t84 + t52 * t45 + t51 * t46) + (-t153 * t69 - t154 * t70) * t275 + (-t143 * t70 + t144 * t69 - t153 * t48 - t154 * t47) * t274 + t200 * (-Icges(6,5) * t162 - Icges(6,6) * t163) + (qJD(5) * (t162 * t92 - t163 * t94) - t162 * (Icges(6,1) * t194 + Icges(6,4) * t195 + Icges(6,5) * t144) - t163 * (Icges(6,4) * t194 + Icges(6,2) * t195 + Icges(6,6) * t144) + t17) * t282 + (qJD(5) * (t162 * t93 - t163 * t95) - t162 * (Icges(6,1) * t192 + Icges(6,4) * t193 + Icges(6,5) * t143) - t163 * (Icges(6,4) * t192 + Icges(6,2) * t193 + Icges(6,6) * t143) + t18) * t280; m(7) * (t52 * t276 - t51 * t277 + (-t276 * t83 - t277 * t84) * qJD(1)) + t187 * t275 + (t276 * t144 + t277 * t143 + (-t153 * t277 - t154 * t276) * qJD(1)) * t274; -m(6) * t27 - m(7) * t11; m(7) * (-t143 * t84 - t144 * t83 - t153 * t51 + t154 * t52); ((t284 + t285) * t140 + t200 * t151) * t151 * t287 + (t29 * t11 - t51 * t83 - t52 * t84) * t286 + (t284 * t64 + t96 * t296 + t4) * t154 + (-t284 * t63 + t97 * t296 - t3 + (-t153 * t63 + t154 * t64) * t153) * t153 + (0.3e1 * t90 * t285 + t15 + (-t299 - t300) * t154) * t144 + (0.3e1 * t284 * t91 + t16 + (-t298 + t300) * t153) * t143; m(7) * (t19 * t45 + t20 * t46 + t21 * t54 + t22 * t53) + t28 + (-t153 * t237 - t154 * t238) * t242 + (-t77 * qJD(5) + (-t9 / 0.2e1 - t17 / 0.2e1) * t154 + (-t10 / 0.2e1 - t18 / 0.2e1) * t153 + t237 * t144 - t238 * t143) * t162; m(7) * (t19 * t276 - t20 * t277 + (t276 * t54 + t277 * t53) * qJD(1)); -m(7) * t12; m(7) * (t143 * t53 + t144 * t54 - t153 * t20 + t154 * t19); t143 * t14 / 0.2e1 + t2 * t280 + t13 * t283 + t1 * t282 + (t31 * t143 + t30 * t144 + t270 - t271) * t279 + m(7) * (t40 * t11 + t12 * t29 - t19 * t84 - t20 * t83 + t54 * t51 + t53 * t52) + (t15 * t281 + t16 * t282) * t242 + (t16 * t283 + t4 * t282 - t143 * t15 / 0.2e1 + t3 * t281 - qJD(5) * (-t153 * t30 + t154 * t31) / 0.2e1) * t162; (t12 * t40 + t19 * t53 + t20 * t54) * t286 + (t28 + (-t154 * t13 - t153 * t14 + t163 * t215) * qJD(5)) * t163 + (t144 * t14 - t153 * t2 - t143 * t13 - t154 * t1 + t163 * (-t10 * t153 - t30 * t143 + t31 * t144 - t9 * t154) + (-t162 * t215 - 0.2e1 * t77 * t163) * qJD(5)) * t162;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t41(1) t41(2) t41(4) t41(7) t41(11) t41(16); t41(2) t41(3) t41(5) t41(8) t41(12) t41(17); t41(4) t41(5) t41(6) t41(9) t41(13) t41(18); t41(7) t41(8) t41(9) t41(10) t41(14) t41(19); t41(11) t41(12) t41(13) t41(14) t41(15) t41(20); t41(16) t41(17) t41(18) t41(19) t41(20) t41(21);];
Mq  = res;
