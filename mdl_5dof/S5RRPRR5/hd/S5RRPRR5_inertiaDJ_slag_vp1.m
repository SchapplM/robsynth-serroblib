% Calculate time derivative of joint inertia matrix for
% S5RRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% m [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:03
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR5_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR5_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR5_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR5_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR5_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:02:12
% EndTime: 2022-01-20 11:02:28
% DurationCPUTime: 5.47s
% Computational Cost: add. (9454->375), mult. (6258->532), div. (0->0), fcn. (4706->10), ass. (0->228)
t185 = cos(pkin(9));
t314 = rSges(4,1) * t185 + pkin(2) - rSges(4,2) * sin(pkin(9));
t182 = qJD(1) + qJD(2);
t183 = qJ(1) + qJ(2);
t175 = sin(t183);
t176 = cos(t183);
t310 = rSges(4,3) + qJ(3);
t75 = t310 * t175 + t314 * t176;
t313 = t182 * t75;
t309 = t314 * t175;
t180 = pkin(9) + qJ(4);
t172 = sin(t180);
t173 = cos(t180);
t257 = qJD(4) * t175;
t269 = t176 * t182;
t308 = -t172 * t269 - t173 * t257;
t285 = rSges(5,2) * t172;
t289 = rSges(5,1) * t173;
t307 = -t285 + t289;
t277 = Icges(5,4) * t173;
t212 = -Icges(5,2) * t172 + t277;
t202 = t212 * t176;
t93 = Icges(5,6) * t175 + t202;
t278 = Icges(5,4) * t172;
t214 = Icges(5,1) * t173 - t278;
t204 = t214 * t176;
t95 = Icges(5,5) * t175 + t204;
t217 = t172 * t93 - t173 * t95;
t306 = t175 * t217;
t174 = qJ(5) + t180;
t166 = sin(t174);
t167 = cos(t174);
t275 = Icges(6,4) * t167;
t211 = -Icges(6,2) * t166 + t275;
t201 = t211 * t176;
t83 = Icges(6,6) * t175 + t201;
t276 = Icges(6,4) * t166;
t213 = Icges(6,1) * t167 - t276;
t203 = t213 * t176;
t85 = Icges(6,5) * t175 + t203;
t219 = t166 * t83 - t167 * t85;
t305 = t175 * t219;
t92 = -Icges(5,6) * t176 + t175 * t212;
t94 = -Icges(5,5) * t176 + t175 * t214;
t218 = t172 * t92 - t173 * t94;
t304 = t176 * t218;
t82 = -Icges(6,6) * t176 + t175 * t211;
t84 = -Icges(6,5) * t176 + t175 * t213;
t220 = t166 * t82 - t167 * t84;
t303 = t176 * t220;
t159 = t175 * rSges(6,3);
t288 = rSges(6,1) * t167;
t302 = t176 * t288 + t159;
t181 = qJD(4) + qJD(5);
t118 = Icges(6,2) * t167 + t276;
t119 = Icges(6,1) * t166 + t275;
t208 = t118 * t166 - t119 * t167;
t209 = Icges(6,5) * t167 - Icges(6,6) * t166;
t301 = t209 * t181 + t182 * t208;
t128 = Icges(5,2) * t173 + t278;
t129 = Icges(5,1) * t172 + t277;
t207 = t128 * t172 - t129 * t173;
t210 = Icges(5,5) * t173 - Icges(5,6) * t172;
t300 = t210 * qJD(4) + t182 * t207;
t299 = 2 * m(3);
t298 = 2 * m(4);
t297 = 2 * m(5);
t296 = 2 * m(6);
t295 = t175 / 0.2e1;
t294 = -t176 / 0.2e1;
t116 = t307 * qJD(4);
t293 = m(5) * t116;
t130 = rSges(5,1) * t172 + rSges(5,2) * t173;
t292 = m(5) * t130;
t187 = sin(qJ(1));
t291 = pkin(1) * t187;
t186 = -pkin(7) - qJ(3);
t168 = t185 * pkin(3) + pkin(2);
t284 = rSges(6,2) * t166;
t139 = t175 * t284;
t263 = t176 * rSges(6,3) + t139;
t287 = rSges(6,1) * t175;
t86 = t167 * t287 - t263;
t254 = t176 * t284;
t87 = -t254 + t302;
t50 = t175 * t86 + t176 * t87;
t283 = pkin(1) * qJD(1);
t160 = t175 * rSges(5,3);
t274 = t118 * t181;
t273 = t119 * t181;
t272 = t166 * t181;
t271 = t167 * t181;
t270 = t175 * t182;
t179 = pkin(8) - t186;
t268 = t179 * t182;
t267 = t182 * t186;
t135 = pkin(4) * t173 + t168;
t266 = t176 * t135 + t179 * t175;
t265 = rSges(6,3) * t269 + t182 * t139;
t249 = t172 * t270;
t264 = rSges(5,2) * t249 + rSges(5,3) * t269;
t245 = t172 * t257;
t262 = -pkin(4) * t245 + t175 * t268;
t261 = t176 * rSges(5,3) + t175 * t285;
t260 = t175 ^ 2 + t176 ^ 2;
t259 = qJD(4) * t172;
t258 = qJD(4) * t173;
t253 = rSges(6,2) * t271;
t247 = -t175 * t253 - t182 * t254 - t272 * t287;
t256 = t175 * (t302 * t182 + t247) + t176 * (-t176 * t253 + (-t167 * t270 - t176 * t272) * rSges(6,1) + t265) + t86 * t269;
t252 = t187 * t283;
t188 = cos(qJ(1));
t251 = t188 * t283;
t250 = pkin(4) * t259;
t246 = -rSges(5,1) * t245 + t308 * rSges(5,2);
t243 = t270 / 0.2e1;
t242 = t269 / 0.2e1;
t120 = rSges(6,1) * t166 + rSges(6,2) * t167;
t240 = -pkin(4) * t172 - t120;
t141 = t176 * t268;
t155 = qJD(3) * t175;
t224 = -t135 - t288;
t215 = t224 * t175;
t18 = t141 + t155 + t182 * t215 + (-t120 * t181 - t250) * t176 + t265;
t16 = t18 - t252;
t63 = t176 * t179 + t215 + t263;
t59 = t63 - t291;
t239 = t182 * t59 - t16;
t156 = qJD(3) * t176;
t19 = t156 + (t176 * t224 - t159) * t182 - t247 - t262;
t17 = t19 - t251;
t178 = t188 * pkin(1);
t64 = t87 + t266;
t60 = t178 + t64;
t238 = t182 * t60 + t17;
t237 = t182 * t63 - t18;
t197 = Icges(6,6) * t182 - t274;
t46 = t176 * t197 - t211 * t270;
t236 = t181 * t85 + t46;
t47 = t175 * t197 + t182 * t201;
t235 = t181 * t84 + t47;
t198 = Icges(6,5) * t182 - t273;
t48 = t176 * t198 - t213 * t270;
t234 = -t181 * t83 + t48;
t49 = t175 * t198 + t182 * t203;
t233 = -t181 * t82 + t49;
t232 = t182 * t64 + t19;
t205 = t130 * qJD(4);
t225 = -t168 - t289;
t216 = t225 * t175;
t30 = t155 + t182 * t216 + (-t205 - t267) * t176 + t264;
t28 = t30 - t252;
t67 = -t176 * t186 + t216 + t261;
t65 = t67 - t291;
t231 = t182 * t65 - t28;
t143 = t175 * t267;
t31 = t143 + t156 + (t176 * t225 - t160) * t182 - t246;
t29 = t31 - t251;
t223 = -t176 * t168 + t175 * t186;
t97 = t307 * t176 + t160;
t68 = -t223 + t97;
t66 = t178 + t68;
t230 = t182 * t66 + t29;
t229 = t182 * t67 - t30;
t228 = t182 * t68 + t31;
t104 = (-t284 + t288) * t181;
t38 = t120 * t270 - t104 * t176 + (-t258 * t176 + t249) * pkin(4);
t78 = t240 * t175;
t227 = t182 * t78 + t38;
t39 = t308 * pkin(4) - t104 * t175 - t120 * t269;
t79 = t240 * t176;
t226 = t182 * t79 - t39;
t132 = t176 * rSges(3,1) - rSges(3,2) * t175;
t109 = -rSges(3,1) * t269 + rSges(3,2) * t270;
t80 = -Icges(6,3) * t176 + t175 * t209;
t20 = -t175 * t220 - t176 * t80;
t199 = t209 * t176;
t81 = Icges(6,3) * t175 + t199;
t21 = -t176 * t81 - t305;
t22 = t175 * t80 - t303;
t23 = t175 * t81 - t176 * t219;
t117 = Icges(6,5) * t166 + Icges(6,6) * t167;
t196 = Icges(6,3) * t182 - t117 * t181;
t44 = t176 * t196 - t209 * t270;
t45 = t175 * t196 + t182 * t199;
t222 = -t176 * ((t176 * t45 + (t21 + t303) * t182) * t176 + (t20 * t182 + (-t166 * t46 + t167 * t48 - t271 * t83 - t272 * t85) * t175 + (-t44 + (t182 * t85 - t233) * t167 + (-t182 * t83 + t235) * t166) * t176) * t175) + t175 * ((t175 * t44 + (t22 + t305) * t182) * t175 + (t23 * t182 + (t166 * t47 - t167 * t49 + t271 * t82 + t272 * t84) * t176 + (-t45 + (t182 * t84 + t234) * t167 + (-t182 * t82 - t236) * t166) * t175) * t176) + (t175 * t21 - t176 * t20) * t270 + (t175 * t23 - t176 * t22) * t269;
t131 = -rSges(3,1) * t175 - rSges(3,2) * t176;
t127 = Icges(5,5) * t172 + Icges(5,6) * t173;
t102 = t211 * t181;
t103 = t213 * t181;
t190 = t117 * t182 + (t103 - t274) * t167 + (-t102 - t273) * t166;
t206 = (t166 * t234 + t167 * t236 + t301 * t175 + t190 * t176) * t295 + (t166 * t233 + t167 * t235 + t190 * t175 - t301 * t176) * t294 + (-t117 * t176 + t166 * t84 + t167 * t82 - t175 * t208) * t243 + (t117 * t175 + t166 * t85 + t167 * t83 - t176 * t208) * t242;
t108 = t131 * t182;
t200 = t210 * t176;
t74 = t310 * t176 - t309;
t195 = Icges(5,5) * t182 - qJD(4) * t129;
t194 = Icges(5,6) * t182 - qJD(4) * t128;
t193 = Icges(5,3) * t182 - qJD(4) * t127;
t61 = -t309 * t182 + t310 * t269 + t155;
t114 = t212 * qJD(4);
t115 = t214 * qJD(4);
t192 = t167 * t102 + t166 * t103 + t173 * t114 + t172 * t115 - t118 * t272 + t119 * t271 - t128 * t259 + t129 * t258;
t189 = -t114 * t172 + t115 * t173 + t127 * t182 + (-t128 * t173 - t129 * t172) * qJD(4);
t191 = t206 + (t300 * t175 + t189 * t176 - qJD(4) * t217 + t172 * (t176 * t195 - t214 * t270) + t173 * (t176 * t194 - t212 * t270)) * t295 + (t189 * t175 - t300 * t176 - qJD(4) * t218 + t172 * (t175 * t195 + t182 * t204) + t173 * (t175 * t194 + t182 * t202)) * t294 + (-t127 * t176 + t172 * t94 + t173 * t92 - t175 * t207) * t243 + (t127 * t175 + t172 * t95 + t173 * t93 - t176 * t207) * t242;
t62 = t156 - t313;
t112 = t132 + t178;
t111 = t131 - t291;
t100 = t109 - t251;
t99 = t108 - t252;
t96 = t175 * t289 - t261;
t91 = Icges(5,3) * t175 + t200;
t90 = -Icges(5,3) * t176 + t175 * t210;
t72 = t178 + t75;
t71 = t74 - t291;
t70 = t223 + t266;
t69 = (-t179 - t186) * t176 + (t135 - t168) * t175;
t54 = t175 * t193 + t182 * t200;
t53 = t176 * t193 - t210 * t270;
t52 = t62 - t251;
t51 = t61 - t252;
t27 = t175 * t91 - t176 * t217;
t26 = t175 * t90 - t304;
t25 = -t176 * t91 - t306;
t24 = -t175 * t218 - t176 * t90;
t15 = t175 * t69 + t176 * t70 + t50;
t10 = -t270 * t87 + t256;
t3 = (t182 * t69 + t141 + (-t250 + t267) * t176) * t176 + (t143 + (-t70 - t87) * t182 + t262) * t175 + t256;
t1 = [(t16 * t60 + t17 * t59) * t296 + (t28 * t66 + t29 * t65) * t297 + (t51 * t72 + t52 * t71) * t298 + (t100 * t111 + t112 * t99) * t299 + t192; m(6) * (t16 * t64 + t17 * t63 + t18 * t60 + t19 * t59) + m(5) * (t28 * t68 + t29 * t67 + t30 * t66 + t31 * t65) + m(4) * (t51 * t75 + t52 * t74 + t61 * t72 + t62 * t71) + m(3) * (t100 * t131 + t108 * t112 + t109 * t111 + t132 * t99) + t192; (t18 * t64 + t19 * t63) * t296 + (t30 * t68 + t31 * t67) * t297 + (t61 * t75 + t62 * t74) * t298 + (t108 * t132 + t109 * t131) * t299 + t192; m(6) * (t175 * t238 + t176 * t239) + m(5) * (t175 * t230 + t176 * t231) + m(4) * ((t182 * t71 - t51) * t176 + (t182 * t72 + t52) * t175); m(6) * (t175 * t232 + t176 * t237) + m(5) * (t175 * t228 + t176 * t229) + m(4) * ((t182 * t74 - t61) * t176 + (t62 + t313) * t175); 0; t191 + m(6) * (t16 * t78 + t17 * t79 + t38 * t59 + t39 * t60) + (t175 * t231 - t176 * t230) * t292 + (-t175 * t66 - t176 * t65) * t293; t191 + m(6) * (t18 * t78 + t19 * t79 + t38 * t63 + t39 * t64) + (t175 * t229 - t176 * t228) * t292 + (-t175 * t68 - t176 * t67) * t293; m(6) * (t175 * t227 + t176 * t226); ((t175 * t96 + t176 * t97) * (((-t97 + t160) * t182 + t246) * t175 + (-t176 * t205 + t182 * t96 + t264) * t176) + t260 * t130 * t116) * t297 + (t175 * t27 - t176 * t26) * t269 + t175 * ((t175 * t53 + (t26 + t306) * t182) * t175 + (t27 * t182 + (t258 * t92 + t259 * t94) * t176 + (-t54 + (-qJD(4) * t93 + t182 * t94) * t173 + (-qJD(4) * t95 - t182 * t92) * t172) * t175) * t176) + (t175 * t25 - t176 * t24) * t270 - t176 * ((t176 * t54 + (t25 + t304) * t182) * t176 + (t24 * t182 + (-t258 * t93 - t259 * t95) * t175 + (-t53 + (qJD(4) * t92 + t182 * t95) * t173 + (qJD(4) * t94 - t182 * t93) * t172) * t176) * t175) + (t15 * t3 + t38 * t79 + t39 * t78) * t296 + t222; m(6) * ((-t175 * t60 - t176 * t59) * t104 + (t175 * t239 - t176 * t238) * t120) + t206; m(6) * ((-t175 * t64 - t176 * t63) * t104 + (t175 * t237 - t176 * t232) * t120) + t206; 0; m(6) * (t10 * t15 + t3 * t50 + (-t175 * t78 - t176 * t79) * t104 + (t175 * t226 - t176 * t227) * t120) + t222; (t104 * t120 * t260 + t10 * t50) * t296 + t222;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
