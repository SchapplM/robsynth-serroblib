% Calculate time derivative of joint inertia matrix for
% S6RPPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3]';
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
% Datum: 2019-03-09 01:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPPRR4_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR4_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR4_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR4_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR4_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPPRR4_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPPRR4_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:35:18
% EndTime: 2019-03-09 01:35:27
% DurationCPUTime: 5.53s
% Computational Cost: add. (9959->487), mult. (24364->714), div. (0->0), fcn. (27301->8), ass. (0->232)
t274 = sin(pkin(9));
t275 = cos(pkin(9));
t286 = sin(qJ(1));
t287 = cos(qJ(1));
t147 = -t274 * t286 - t275 * t287;
t146 = t147 ^ 2;
t139 = t147 * qJD(1);
t148 = t287 * t274 - t286 * t275;
t140 = t148 * qJD(1);
t174 = sin(qJ(5));
t176 = cos(qJ(5));
t224 = rSges(6,1) * t174 + rSges(6,2) * t176;
t200 = t147 * t224;
t225 = rSges(6,1) * t176 - t174 * rSges(6,2);
t251 = qJD(5) * t176;
t205 = t139 * t174 + t148 * t251;
t267 = t139 * t176;
t234 = t205 * rSges(6,1) + rSges(6,2) * t267 + t140 * rSges(6,3);
t252 = qJD(5) * t174;
t237 = t148 * t252;
t253 = qJD(5) * t147;
t262 = t148 * t176;
t263 = t148 * t174;
t97 = rSges(6,1) * t263 + rSges(6,2) * t262 - t147 * rSges(6,3);
t98 = -t148 * rSges(6,3) - t200;
t27 = t139 * t97 + t148 * (-rSges(6,2) * t237 + t234) + (t139 * rSges(6,3) + t225 * t253) * t147 + (t98 - t200) * t140;
t298 = 2 * m(6);
t308 = t27 * t298;
t307 = t225 ^ 2 * t298;
t288 = rSges(7,3) + pkin(8);
t300 = t288 * t176;
t306 = qJ(4) - t300;
t254 = t287 * pkin(1) + t286 * qJ(2);
t241 = t287 * pkin(2) + t254;
t305 = -rSges(3,1) * t287 - rSges(3,3) * t286 - t254;
t283 = pkin(5) * t174;
t304 = t300 - t283;
t273 = Icges(6,4) * t174;
t212 = Icges(6,2) * t176 + t273;
t93 = -Icges(6,6) * t147 + t148 * t212;
t272 = Icges(6,4) * t176;
t214 = Icges(6,1) * t174 + t272;
t95 = -Icges(6,5) * t147 + t148 * t214;
t216 = t174 * t95 + t176 * t93;
t94 = -Icges(6,6) * t148 - t147 * t212;
t96 = -Icges(6,5) * t148 - t147 * t214;
t215 = t174 * t96 + t176 * t94;
t198 = t215 * t148;
t199 = t216 * t147;
t210 = Icges(6,5) * t174 + Icges(6,6) * t176;
t91 = -Icges(6,3) * t147 + t148 * t210;
t92 = -Icges(6,3) * t148 - t147 * t210;
t303 = t147 * t92 + t148 * t91 - t198 + t199;
t228 = -pkin(8) * t176 + t283;
t302 = t147 * t228;
t154 = t224 * qJD(5);
t301 = t154 * t225 * t298;
t265 = t140 * t176;
t202 = t147 * t252 + t265;
t299 = qJD(1) * (t147 * t286 - t148 * t287) - t139 * t286 + t140 * t287;
t297 = 2 * m(7);
t296 = -0.2e1 * t225;
t295 = -t139 / 0.2e1;
t294 = -t147 / 0.2e1;
t293 = -t148 / 0.2e1;
t292 = t174 / 0.2e1;
t291 = -t176 / 0.2e1;
t290 = rSges(5,2) - pkin(3);
t289 = rSges(6,3) + pkin(3);
t285 = m(6) * t154;
t284 = m(6) * t225;
t282 = pkin(5) * t176;
t281 = t140 * pkin(7);
t173 = sin(qJ(6));
t175 = cos(qJ(6));
t260 = t173 * t174;
t108 = -t147 * t175 - t148 * t260;
t259 = t174 * t175;
t109 = -t147 * t173 + t148 * t259;
t76 = Icges(7,4) * t109 + Icges(7,2) * t108 - Icges(7,6) * t262;
t78 = Icges(7,1) * t109 + Icges(7,4) * t108 - Icges(7,5) * t262;
t218 = t173 * t76 - t175 * t78;
t204 = t237 - t267;
t190 = -qJD(6) * t147 + t205;
t250 = qJD(6) * t174;
t226 = -t148 * t250 + t140;
t64 = -t173 * t190 + t175 * t226;
t65 = t173 * t226 + t175 * t190;
t33 = Icges(7,5) * t65 + Icges(7,6) * t64 + Icges(7,3) * t204;
t35 = Icges(7,4) * t65 + Icges(7,2) * t64 + Icges(7,6) * t204;
t37 = Icges(7,1) * t65 + Icges(7,4) * t64 + Icges(7,5) * t204;
t74 = Icges(7,5) * t109 + Icges(7,6) * t108 - Icges(7,3) * t262;
t9 = (-qJD(5) * t218 - t33) * t174 + (-qJD(5) * t74 + t173 * t35 - t175 * t37 + (t173 * t78 + t175 * t76) * qJD(6)) * t176;
t280 = t9 * t147;
t110 = t147 * t260 - t148 * t175;
t111 = -t147 * t259 - t148 * t173;
t264 = t147 * t176;
t77 = Icges(7,4) * t111 + Icges(7,2) * t110 + Icges(7,6) * t264;
t79 = Icges(7,1) * t111 + Icges(7,4) * t110 + Icges(7,5) * t264;
t217 = t173 * t77 - t175 * t79;
t203 = -t140 * t174 + t147 * t251;
t189 = qJD(6) * t148 + t203;
t227 = t147 * t250 - t139;
t66 = t173 * t189 + t175 * t227;
t67 = t173 * t227 - t175 * t189;
t34 = Icges(7,5) * t67 + Icges(7,6) * t66 - Icges(7,3) * t202;
t36 = Icges(7,4) * t67 + Icges(7,2) * t66 - Icges(7,6) * t202;
t38 = Icges(7,1) * t67 + Icges(7,4) * t66 - Icges(7,5) * t202;
t75 = Icges(7,5) * t111 + Icges(7,6) * t110 + Icges(7,3) * t264;
t10 = (-qJD(5) * t217 - t34) * t174 + (-qJD(5) * t75 + t173 * t36 - t175 * t38 + (t173 * t79 + t175 * t77) * qJD(6)) * t176;
t279 = t10 * t148;
t278 = rSges(5,3) + qJ(4);
t138 = pkin(5) * t263;
t248 = pkin(8) * t262;
t257 = t109 * rSges(7,1) + t108 * rSges(7,2);
t80 = -rSges(7,3) * t262 + t257;
t277 = t138 - t248 + t80;
t223 = -t111 * rSges(7,1) - t110 * rSges(7,2);
t81 = rSges(7,3) * t264 - t223;
t276 = -t81 + t302;
t271 = Icges(7,4) * t173;
t270 = Icges(7,4) * t175;
t213 = Icges(7,1) * t175 - t271;
t124 = -Icges(7,5) * t174 - t176 * t213;
t269 = t124 * t175;
t266 = t140 * t147;
t222 = rSges(7,1) * t175 - rSges(7,2) * t173;
t249 = qJD(6) * t176;
t105 = (rSges(7,1) * t173 + rSges(7,2) * t175) * t249 + (-rSges(7,3) * t176 + t174 * t222) * qJD(5);
t258 = t228 * qJD(5) + t105;
t125 = -t174 * rSges(7,3) - t176 * t222;
t229 = t174 * pkin(8) + t282;
t256 = t125 - t229;
t170 = t287 * qJ(2);
t255 = qJD(1) * t170 + qJD(2) * t286;
t247 = t65 * rSges(7,1) + t64 * rSges(7,2) + rSges(7,3) * t237;
t246 = t286 * pkin(1);
t29 = -t174 * t74 + t176 * t218;
t209 = Icges(7,5) * t175 - Icges(7,6) * t173;
t122 = -Icges(7,3) * t174 - t176 * t209;
t211 = -Icges(7,2) * t173 + t270;
t123 = -Icges(7,6) * t174 - t176 * t211;
t50 = t108 * t123 + t109 * t124 - t122 * t262;
t245 = t29 / 0.2e1 + t50 / 0.2e1;
t30 = -t174 * t75 + t176 * t217;
t51 = t110 * t123 + t111 * t124 + t122 * t264;
t244 = -t30 / 0.2e1 - t51 / 0.2e1;
t243 = t123 * t260;
t242 = t205 * pkin(5) + pkin(8) * t237;
t235 = m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1;
t233 = -t147 * pkin(3) + t241;
t230 = t67 * rSges(7,1) + t66 * rSges(7,2);
t23 = t108 * t76 + t109 * t78 - t262 * t74;
t24 = t108 * t77 + t109 * t79 - t262 * t75;
t221 = t147 * t24 - t148 * t23;
t25 = t110 * t76 + t111 * t78 + t264 * t74;
t26 = t110 * t77 + t111 * t79 + t264 * t75;
t220 = t147 * t26 - t148 * t25;
t219 = t30 * t147 - t29 * t148;
t208 = t139 * t148 - t266;
t207 = -t176 * t33 + t252 * t74;
t206 = -t176 * t34 + t252 * t75;
t201 = -t147 * pkin(7) + t233;
t197 = -pkin(2) * t286 - t246;
t192 = t170 + t197;
t191 = t147 * qJ(4) + t192;
t182 = qJD(1) * t197 + t255;
t181 = t140 * pkin(3) + t182;
t46 = t281 + t139 * qJ(4) + (-rSges(6,2) * t252 + qJD(4)) * t148 + t181 + t234;
t168 = qJD(2) * t287;
t180 = -t241 * qJD(1) + t168;
t178 = -t140 * qJ(4) + t147 * qJD(4) + t180;
t177 = t139 * pkin(7) + t178;
t47 = rSges(6,1) * t203 - rSges(6,2) * t202 + t139 * t289 + t177;
t187 = t148 * pkin(7) + t191;
t60 = t148 * t289 + t187 + t200;
t61 = t148 * qJ(4) + t201 + t97;
t188 = -t139 * t60 - t140 * t61 + t147 * t46 - t148 * t47;
t102 = (Icges(7,5) * t173 + Icges(7,6) * t175) * t249 + (-Icges(7,3) * t176 + t174 * t209) * qJD(5);
t103 = (Icges(7,2) * t175 + t271) * t249 + (-Icges(7,6) * t176 + t174 * t211) * qJD(5);
t104 = (Icges(7,1) * t173 + t270) * t249 + (-Icges(7,5) * t176 + t174 * t213) * qJD(5);
t186 = -t174 * t102 + t252 * t269 + (-t104 * t176 + t123 * t249) * t175 + (t103 * t176 + t124 * t249) * t173;
t185 = -rSges(3,1) * t286 + rSges(3,3) * t287 - t246;
t179 = t148 * qJD(4) + t181;
t130 = t170 + t185;
t115 = t305 * qJD(1) + t168;
t114 = qJD(1) * t185 + t255;
t101 = -rSges(4,1) * t147 - rSges(4,2) * t148 + t241;
t100 = t148 * rSges(4,1) - t147 * rSges(4,2) + t192;
t90 = t256 * t147;
t89 = t256 * t148;
t88 = t139 * rSges(4,1) + t140 * rSges(4,2) + t180;
t87 = t140 * rSges(4,1) - t139 * rSges(4,2) + t182;
t85 = rSges(5,2) * t147 + t148 * t278 + t233;
t84 = t147 * rSges(5,3) - t148 * t290 + t191;
t83 = -t174 * t122 + (t123 * t173 - t269) * t176;
t69 = -Icges(6,5) * t203 + Icges(6,6) * t202 - Icges(6,3) * t139;
t68 = Icges(6,5) * t205 - Icges(6,6) * t204 + Icges(6,3) * t140;
t59 = -t140 * rSges(5,3) - t139 * t290 + t178;
t58 = -t140 * rSges(5,2) + t139 * t278 + t179;
t55 = t125 * t264 + t174 * t81;
t54 = t125 * t262 - t174 * t80;
t53 = -t140 * t256 + t147 * t258;
t52 = -t139 * t256 - t148 * t258;
t49 = t306 * t148 + t138 + t201 + t257;
t48 = t148 * pkin(3) - t304 * t147 + t187 + t223;
t41 = (-t147 * t80 - t148 * t81) * t176;
t40 = -rSges(7,3) * t202 + t230;
t39 = -rSges(7,3) * t267 + t247;
t31 = (-t122 * t176 - t243) * qJD(5) + t186;
t28 = t147 * t276 + t148 * t277;
t22 = t139 * pkin(3) + t304 * t140 + (t174 * t288 + t282) * t253 + t177 - t230;
t21 = t306 * t139 + t179 + t242 + t247 + t281;
t20 = (-t125 * t253 + t40) * t174 + (qJD(5) * t81 + t105 * t147 - t125 * t140) * t176;
t19 = (-qJD(5) * t125 * t148 - t39) * t174 + (-qJD(5) * t80 + t105 * t148 + t125 * t139) * t176;
t18 = t102 * t264 + t110 * t103 + t111 * t104 - t122 * t202 + t66 * t123 + t67 * t124;
t17 = -t102 * t262 + t108 * t103 + t109 * t104 + t122 * t204 + t64 * t123 + t65 * t124;
t16 = -t147 * t25 - t148 * t26;
t15 = -t147 * t23 - t148 * t24;
t14 = -t51 * t174 + t176 * t220;
t13 = -t50 * t174 + t176 * t221;
t12 = t39 * t264 - t81 * t237 + (t139 * t81 + t148 * t40) * t176 - t202 * t80;
t11 = (t39 + t242) * t148 - (t248 - t277) * t139 + (t229 * t253 - t40) * t147 + (-t276 - t302) * t140;
t8 = t110 * t36 + t111 * t38 - t147 * t206 - t265 * t75 + t66 * t77 + t67 * t79;
t7 = t110 * t35 + t111 * t37 - t147 * t207 - t265 * t74 + t66 * t76 + t67 * t78;
t6 = t108 * t36 + t109 * t38 + t148 * t206 - t267 * t75 + t64 * t77 + t65 * t79;
t5 = t108 * t35 + t109 * t37 + t148 * t207 - t267 * t74 + t64 * t76 + t65 * t78;
t4 = -t139 * t26 + t140 * t25 - t147 * t7 - t148 * t8;
t3 = -t139 * t24 + t140 * t23 - t147 * t5 - t148 * t6;
t2 = (-qJD(5) * t220 - t18) * t174 + (-qJD(5) * t51 - t139 * t25 - t140 * t26 + t147 * t8 - t148 * t7) * t176;
t1 = (-qJD(5) * t221 - t17) * t174 + (-qJD(5) * t50 - t139 * t23 - t140 * t24 + t147 * t6 - t148 * t5) * t176;
t32 = [-qJD(5) * t243 + (t21 * t49 + t22 * t48) * t297 + (t46 * t61 + t47 * t60) * t298 + 0.2e1 * m(5) * (t58 * t85 + t59 * t84) + 0.2e1 * m(3) * (-t114 * t305 + t115 * t130) + 0.2e1 * m(4) * (t100 * t88 + t101 * t87) + t186 + (-Icges(6,1) * t176 + t212 + t273) * t252 + (Icges(6,2) * t174 - t122 - t214 - t272) * t251; m(7) * (t286 * t22 - t287 * t21 + (t286 * t49 + t287 * t48) * qJD(1)) + m(6) * (t286 * t47 - t287 * t46 + (t286 * t61 + t287 * t60) * qJD(1)) + m(5) * (t286 * t59 - t287 * t58 + (t286 * t85 + t287 * t84) * qJD(1)) + m(3) * (t286 * t115 - t287 * t114 + (t130 * t287 - t286 * t305) * qJD(1)) + m(4) * (t286 * t88 - t287 * t87 + (t100 * t287 + t101 * t286) * qJD(1)); 0; 0; 0; 0; m(7) * (t139 * t48 + t140 * t49 - t147 * t21 + t148 * t22) - m(6) * t188 + m(5) * (t139 * t84 + t140 * t85 - t147 * t58 + t148 * t59); -0.2e1 * t235 * t299; 0; 0.4e1 * t235 * t208; -t279 / 0.2e1 - t280 / 0.2e1 + (t148 ^ 2 / 0.2e1 + t146 / 0.2e1) * t210 * qJD(5) + (t291 * t95 + t292 * t93 + t245) * t140 - (t291 * t96 + t292 * t94 - t244) * t139 + m(7) * (t21 * t90 - t22 * t89 + t52 * t48 + t53 * t49) + (t147 * t61 - t148 * t60) * t285 - t188 * t284 + t208 * (-Icges(6,5) * t176 + Icges(6,6) * t174) + (qJD(5) * t216 + t174 * (Icges(6,4) * t205 - Icges(6,2) * t204 + Icges(6,6) * t140) - t176 * (Icges(6,1) * t205 - Icges(6,4) * t204 + Icges(6,5) * t140) + t17) * t294 + (qJD(5) * t215 + t174 * (-Icges(6,4) * t203 + Icges(6,2) * t202 - Icges(6,6) * t139) - t176 * (-Icges(6,1) * t203 + Icges(6,4) * t202 - Icges(6,5) * t139) + t18) * t293; m(7) * (t52 * t286 - t53 * t287 + (t286 * t90 - t287 * t89) * qJD(1)) + (-t147 * t287 - t148 * t286) * t285 - t299 * t284; -m(6) * t27 - m(7) * t11; m(6) * (t266 * t296 - t146 * t154 + (-t139 * t296 - t148 * t154) * t148) + m(7) * (-t139 * t89 + t140 * t90 - t147 * t53 + t148 * t52); -t146 * t301 + t140 * t15 - t139 * t16 + (t28 * t11 - t52 * t89 + t53 * t90) * t297 + (t97 * t308 - t4 + (-t148 * t69 - t301) * t148 + (t198 + t303) * t140 + (-0.2e1 * t147 * t215 - 0.3e1 * t148 * t92 + t307) * t139) * t148 + (-t98 * t308 - t3 - t146 * t68 + (-t147 * t69 - t148 * t68) * t148 + (t199 - t303 + (t215 - t91) * t148) * t139 + (0.3e1 * t147 * t91 - t307 + (-t216 + t92) * t148) * t140) * t147; m(7) * (t19 * t49 + t20 * t48 + t21 * t54 + t22 * t55) + (-t31 + (t147 * t244 + t148 * t245) * qJD(5)) * t174 + (-t83 * qJD(5) + (-t9 / 0.2e1 - t17 / 0.2e1) * t148 + (t10 / 0.2e1 + t18 / 0.2e1) * t147 + t244 * t140 - t245 * t139) * t176; m(7) * (t20 * t286 - t19 * t287 + (t286 * t54 + t287 * t55) * qJD(1)); m(7) * t12; m(7) * (t139 * t55 + t140 * t54 - t147 * t19 + t148 * t20); t140 * t13 / 0.2e1 + t1 * t294 + t14 * t295 + t2 * t293 - t174 * (-t30 * t139 + t29 * t140 - t279 - t280) / 0.2e1 + m(7) * (t41 * t11 - t12 * t28 + t19 * t90 - t20 * t89 + t55 * t52 + t54 * t53) + (t148 * t15 / 0.2e1 + t16 * t294) * t252 + (t15 * t295 + t3 * t293 - t140 * t16 / 0.2e1 + t147 * t4 / 0.2e1 - qJD(5) * (-t147 * t29 - t148 * t30) / 0.2e1) * t176; (-t12 * t41 + t19 * t54 + t20 * t55) * t297 + (t31 * t174 + (t148 * t13 - t147 * t14 + t174 * t219) * qJD(5)) * t174 + (-t139 * t13 - t148 * t1 - t140 * t14 + t147 * t2 - t174 * (t10 * t147 - t29 * t139 - t30 * t140 - t9 * t148) + (0.2e1 * t83 * t174 - t176 * t219) * qJD(5)) * t176;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t32(1) t32(2) t32(4) t32(7) t32(11) t32(16); t32(2) t32(3) t32(5) t32(8) t32(12) t32(17); t32(4) t32(5) t32(6) t32(9) t32(13) t32(18); t32(7) t32(8) t32(9) t32(10) t32(14) t32(19); t32(11) t32(12) t32(13) t32(14) t32(15) t32(20); t32(16) t32(17) t32(18) t32(19) t32(20) t32(21);];
Mq  = res;
