% Calculate time derivative of joint inertia matrix for
% S5RPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2022-01-23 09:23
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR4_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR4_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR4_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR4_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR4_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:22:25
% EndTime: 2022-01-23 09:22:36
% DurationCPUTime: 5.87s
% Computational Cost: add. (7630->394), mult. (5952->532), div. (0->0), fcn. (4556->10), ass. (0->234)
t344 = Icges(4,3) + Icges(5,3);
t170 = qJ(3) + pkin(9);
t161 = sin(t170);
t163 = cos(t170);
t173 = sin(qJ(3));
t175 = cos(qJ(3));
t343 = Icges(4,5) * t175 + Icges(5,5) * t163 - Icges(4,6) * t173 - Icges(5,6) * t161;
t171 = qJ(1) + pkin(8);
t162 = sin(t171);
t164 = cos(t171);
t341 = t344 * t162 + t343 * t164;
t342 = t343 * t162 - t344 * t164;
t288 = Icges(4,4) * t175;
t209 = -Icges(4,2) * t173 + t288;
t102 = Icges(4,6) * t162 + t164 * t209;
t289 = Icges(4,4) * t173;
t214 = Icges(4,1) * t175 - t289;
t104 = Icges(4,5) * t162 + t164 * t214;
t286 = Icges(5,4) * t163;
t207 = -Icges(5,2) * t161 + t286;
t92 = Icges(5,6) * t162 + t164 * t207;
t287 = Icges(5,4) * t161;
t212 = Icges(5,1) * t163 - t287;
t94 = Icges(5,5) * t162 + t164 * t212;
t319 = t102 * t173 - t104 * t175 + t161 * t92 - t163 * t94;
t340 = t319 * t162;
t339 = t319 * t164;
t165 = qJ(5) + t170;
t156 = sin(t165);
t297 = rSges(6,2) * t156;
t157 = cos(t165);
t300 = rSges(6,1) * t157;
t225 = -t297 + t300;
t88 = t162 * rSges(6,3) + t164 * t225;
t338 = (-Icges(4,5) * t173 - Icges(5,5) * t161 - Icges(4,6) * t175 - Icges(5,6) * t163) * qJD(3);
t101 = -Icges(4,6) * t164 + t162 * t209;
t103 = -Icges(4,5) * t164 + t162 * t214;
t91 = -Icges(5,6) * t164 + t162 * t207;
t93 = -Icges(5,5) * t164 + t162 * t212;
t337 = t101 * t173 - t103 * t175 + t161 * t91 - t163 * t93;
t335 = t341 * qJD(1);
t334 = t162 * t342 - t340 + (-t337 - t341) * t164;
t159 = t162 ^ 2;
t160 = t164 ^ 2;
t315 = 2 * m(5);
t126 = rSges(5,1) * t161 + rSges(5,2) * t163;
t191 = qJD(3) * t126;
t150 = qJD(4) * t162;
t248 = qJD(3) * t173;
t242 = pkin(3) * t248;
t231 = t164 * t242;
t151 = qJD(4) * t164;
t172 = -qJ(4) - pkin(6);
t254 = qJD(1) * t162;
t258 = t162 * t242 + t172 * t254;
t241 = -t151 - t258;
t253 = qJD(1) * t164;
t304 = -pkin(6) - t172;
t166 = t175 * pkin(3);
t158 = t166 + pkin(2);
t305 = pkin(2) - t158;
t306 = pkin(6) * t162;
t323 = t162 * t305;
t155 = t164 * pkin(6);
t263 = t164 * t172;
t79 = t155 + t263 - t323;
t245 = t162 * ((-t164 * t305 - t306) * qJD(1) + t241) + t164 * (-t231 + t150 + (t164 * t304 + t323) * qJD(1)) + t79 * t253;
t298 = rSges(5,2) * t161;
t259 = rSges(5,3) * t253 + t254 * t298;
t153 = t162 * rSges(5,3);
t321 = -t164 * t298 + t153;
t135 = t164 * t158;
t80 = -pkin(2) * t164 + t162 * t304 + t135;
t301 = rSges(5,1) * t163;
t226 = -t298 + t301;
t95 = -rSges(5,3) * t164 + t162 * t226;
t96 = t164 * t301 + t321;
t4 = (qJD(1) * t95 - t164 * t191 + t259) * t164 + (-t162 * t191 + (-t80 - t96 + t321) * qJD(1)) * t162 + t245;
t333 = t315 * t4;
t332 = t337 * t162 + t164 * t342;
t329 = t341 * t162 - t339;
t328 = -qJD(1) * t342 + t338 * t164;
t327 = -t338 * t162 - t335;
t299 = rSges(4,2) * t173;
t302 = rSges(4,1) * t175;
t227 = -t299 + t302;
t295 = rSges(4,3) * t164;
t105 = t162 * t227 - t295;
t244 = t164 * t299;
t154 = t162 * rSges(4,3);
t257 = t164 * t302 + t154;
t106 = -t244 + t257;
t144 = rSges(4,1) * t173 + rSges(4,2) * t175;
t192 = qJD(3) * t144;
t252 = qJD(1) * t173;
t240 = t162 * t252;
t178 = rSges(4,2) * t240 + rSges(4,3) * t253 - t164 * t192;
t13 = (qJD(1) * t105 + t178) * t164 + (-t162 * t192 + (-t106 - t244 + t154) * qJD(1)) * t162;
t316 = 2 * m(4);
t326 = t13 * t316;
t202 = Icges(6,5) * t157 - Icges(6,6) * t156;
t81 = -Icges(6,3) * t164 + t162 * t202;
t325 = qJD(1) * t81;
t296 = rSges(6,2) * t157;
t120 = rSges(6,1) * t156 + t296;
t169 = qJD(3) + qJD(5);
t193 = t120 * t169;
t284 = Icges(6,4) * t157;
t205 = -Icges(6,2) * t156 + t284;
t84 = Icges(6,6) * t162 + t164 * t205;
t285 = Icges(6,4) * t156;
t210 = Icges(6,1) * t157 - t285;
t86 = Icges(6,5) * t162 + t164 * t210;
t223 = t156 * t84 - t157 * t86;
t324 = t162 * t223;
t83 = -Icges(6,6) * t164 + t162 * t205;
t85 = -Icges(6,5) * t164 + t162 * t210;
t224 = t156 * t83 - t157 * t85;
t322 = t164 * t224;
t310 = sin(qJ(1)) * pkin(1);
t320 = t155 - t310;
t118 = Icges(6,2) * t157 + t285;
t119 = Icges(6,1) * t156 + t284;
t199 = t118 * t156 - t119 * t157;
t317 = qJD(1) * t199 + t202 * t169;
t314 = 2 * m(6);
t313 = t162 / 0.2e1;
t312 = -t164 / 0.2e1;
t311 = m(4) * t144;
t309 = pkin(3) * t173;
t308 = pkin(4) * t161;
t307 = pkin(4) * t163;
t167 = cos(qJ(1)) * pkin(1);
t303 = t162 * t79 + t164 * t80;
t87 = -rSges(6,3) * t164 + t162 * t225;
t36 = t162 * t87 + t164 * t88;
t294 = t161 * t93;
t293 = t161 * t94;
t292 = t163 * t91;
t291 = t163 * t92;
t168 = pkin(7) - t172;
t290 = rSges(6,3) + t168;
t82 = Icges(6,3) * t162 + t164 * t202;
t274 = qJD(1) * t82;
t272 = t101 * t175;
t271 = t102 * t175;
t270 = t103 * t173;
t269 = t104 * t173;
t268 = t118 * t169;
t267 = t119 * t169;
t266 = t156 * t169;
t265 = t157 * t169;
t264 = t164 * t169;
t146 = t168 * t162;
t229 = -t308 - t309;
t124 = t229 * qJD(3);
t262 = t164 * t124 + t168 * t253;
t132 = t158 + t307;
t261 = t164 * t132 + t146;
t260 = rSges(6,3) * t253 + t254 * t297;
t256 = t159 + t160;
t251 = qJD(3) * t161;
t250 = qJD(3) * t162;
t249 = qJD(3) * t163;
t247 = qJD(3) * t175;
t246 = t162 * (qJD(1) * t88 - t162 * t193) + t164 * (-t264 * t296 + (-t156 * t264 - t157 * t254) * rSges(6,1) + t260) + t87 * t253;
t237 = -t126 - t309;
t43 = -qJD(1) * t83 - t164 * t268;
t236 = t169 * t86 + t43;
t44 = qJD(1) * t84 - t162 * t268;
t235 = t169 * t85 + t44;
t45 = -qJD(1) * t85 - t164 * t267;
t234 = -t169 * t84 + t45;
t46 = qJD(1) * t86 - t162 * t267;
t233 = t169 * t83 - t46;
t232 = t162 * t172 - t135;
t16 = -t162 * t224 - t164 * t81;
t17 = -t164 * t82 - t324;
t18 = t162 * t81 - t322;
t19 = t162 * t82 - t164 * t223;
t117 = Icges(6,5) * t156 + Icges(6,6) * t157;
t187 = t117 * t169;
t41 = -t164 * t187 - t325;
t42 = -t162 * t187 + t274;
t230 = -t164 * ((t164 * t42 + (t17 + t322) * qJD(1)) * t164 + (t16 * qJD(1) + (-t156 * t43 + t157 * t45 - t265 * t84 - t266 * t86 + t274) * t162 + (-t41 + t233 * t157 + t235 * t156 + (-t223 - t81) * qJD(1)) * t164) * t162) + t162 * ((t162 * t41 + (t18 + t324) * qJD(1)) * t162 + (t19 * qJD(1) + (t156 * t44 - t157 * t46 + t265 * t83 + t266 * t85 - t325) * t164 + (-t42 + t234 * t157 - t236 * t156 + (-t224 + t82) * qJD(1)) * t162) * t164) + (-t16 * t164 + t162 * t17) * t254 + (t162 * t19 - t164 * t18) * t253;
t98 = t237 * t164;
t195 = -t132 - t225;
t49 = t162 * t195 + t164 * t290 - t310;
t50 = t167 + t88 + t261;
t218 = t162 * t50 + t164 * t49;
t197 = -t120 + t229;
t69 = t197 * t162;
t70 = t197 * t164;
t217 = t162 * t69 + t164 * t70;
t213 = Icges(4,1) * t173 + t288;
t211 = Icges(5,1) * t161 + t286;
t208 = Icges(4,2) * t175 + t289;
t206 = Icges(5,2) * t163 + t287;
t198 = -pkin(2) - t227;
t196 = -t158 - t226;
t108 = t205 * t169;
t109 = t210 * t169;
t177 = qJD(1) * t117 + (t109 - t268) * t157 + (-t108 - t267) * t156;
t190 = (t156 * t234 + t157 * t236 + t317 * t162 + t177 * t164) * t313 + (-t156 * t233 + t157 * t235 + t177 * t162 - t164 * t317) * t312 + (-t117 * t164 + t156 * t85 + t157 * t83 - t162 * t199) * t254 / 0.2e1 + (t117 * t162 + t156 * t86 + t157 * t84 - t164 * t199) * t253 / 0.2e1;
t186 = qJD(3) * t213;
t185 = qJD(3) * t211;
t184 = qJD(3) * t208;
t183 = qJD(3) * t206;
t110 = t225 * t169;
t179 = -t110 + (-t166 - t307) * qJD(3);
t139 = pkin(3) * t240;
t133 = t227 * qJD(3);
t116 = t226 * qJD(3);
t97 = t237 * t162;
t72 = t306 + t167 + (pkin(2) - t299) * t164 + t257;
t71 = t162 * t198 + t295 + t320;
t68 = t232 + t261;
t67 = (-t168 - t172) * t164 + t162 * (t132 - t158);
t60 = t167 - t232 + t96;
t59 = -t310 + (rSges(5,3) - t172) * t164 + t196 * t162;
t48 = -t126 * t253 - t116 * t162 + (-t162 * t247 - t164 * t252) * pkin(3);
t47 = t126 * t254 + t139 + (-pkin(3) * t247 - t116) * t164;
t40 = t144 * t250 + (-t167 + (-rSges(4,3) - pkin(6)) * t162 + t198 * t164) * qJD(1);
t39 = ((-pkin(2) - t302) * t162 + t320) * qJD(1) + t178;
t31 = qJD(1) * t70 + t162 * t179;
t30 = t139 + (t120 + t308) * t254 + t179 * t164;
t25 = t126 * t250 + (t164 * t196 - t153 - t167) * qJD(1) - t241;
t24 = t150 + qJD(3) * t98 + (-t310 - t263 + (-t158 - t301) * t162) * qJD(1) + t259;
t15 = t151 + (-t124 + t193) * t162 + (-t162 * t290 + t164 * t195 - t167) * qJD(1);
t14 = t150 - t164 * t193 + (-t310 + (-t132 - t300) * t162) * qJD(1) + t260 + t262;
t12 = -t254 * t88 + t246;
t7 = t162 * t67 + t164 * t68 + t303 + t36;
t3 = t162 * (t162 * t124 + t258) + t164 * (t231 + t262) + ((t67 + t263) * t164 + (-t68 - t80 - t88 + t146) * t162) * qJD(1) + t245 + t246;
t1 = [(t14 * t50 + t15 * t49) * t314 + t119 * t265 + t156 * t109 - t118 * t266 + t157 * t108 + (t24 * t60 + t59 * t25) * t315 + (t39 * t72 + t40 * t71) * t316 + (t212 - t206) * t251 + (t211 + t207) * t249 + (t214 - t208) * t248 + (t213 + t209) * t247; 0; 0; m(4) * ((-t162 * t39 - t164 * t40) * t144 + (-t162 * t72 - t164 * t71) * t133) + m(5) * (t24 * t97 + t25 * t98 + t47 * t59 + t48 * t60) + m(6) * (t14 * t69 + t15 * t70 + t30 * t49 + t31 * t50) + ((-t72 * t311 + t293 / 0.2e1 + t291 / 0.2e1 + t271 / 0.2e1 + t269 / 0.2e1) * t164 + (t294 / 0.2e1 + t292 / 0.2e1 + t71 * t311 + t272 / 0.2e1 + t270 / 0.2e1) * t162) * qJD(1) + t190 + (-qJD(3) * t319 + t161 * (-qJD(1) * t93 - t164 * t185) + t163 * (-qJD(1) * t91 - t164 * t183) + t173 * (-qJD(1) * t103 - t164 * t186) + t175 * (-qJD(1) * t101 - t164 * t184)) * t313 + (-qJD(3) * t337 + t161 * (qJD(1) * t94 - t162 * t185) + t163 * (qJD(1) * t92 - t162 * t183) + t173 * (qJD(1) * t104 - t162 * t186) + t175 * (qJD(1) * t102 - t162 * t184)) * t312 + t343 * qJD(3) * (t159 / 0.2e1 + t160 / 0.2e1); m(4) * t13 + m(5) * t4 + m(6) * t3; (t7 * t3 + t30 * t70 + t31 * t69) * t314 + (t303 * t4 + t98 * t47 + t97 * t48) * t315 + t256 * t144 * t133 * t316 + t230 + (t106 * t326 + t327 * t160 + t332 * t254 + t96 * t333 + (-t164 * t337 - t334) * t253) * t164 + (t95 * t333 + t105 * t326 + t328 * t159 + t329 * t253 + ((t102 * t247 + t104 * t248 + t249 * t92 + t251 * t94 + t327 - t335) * t162 + (t101 * t247 + t103 * t248 + t249 * t91 + t251 * t93 + t328) * t164 + ((-t269 - t271 - t291 - t293) * t162 + (-t292 - t294 - t270 - t272) * t164) * qJD(3) + ((-t337 + t341) * t162 + t339 + t329 + t332) * qJD(1)) * t164 + (t334 + t340) * t254) * t162; m(6) * (qJD(1) * t218 - t14 * t164 + t15 * t162) + m(5) * (t162 * t25 - t164 * t24 + (t162 * t60 + t164 * t59) * qJD(1)); 0; m(6) * (qJD(1) * t217 + t162 * t30 - t164 * t31) + m(5) * (t162 * t47 - t164 * t48 + (t162 * t97 + t164 * t98) * qJD(1)); 0; m(6) * (-t218 * t110 + (-t14 * t162 - t15 * t164 + (t162 * t49 - t164 * t50) * qJD(1)) * t120) + t190; m(6) * t12; m(6) * (t12 * t7 + t36 * t3 - t217 * t110 + (-t162 * t31 - t164 * t30 + (t162 * t70 - t164 * t69) * qJD(1)) * t120) + t230; 0; (t110 * t120 * t256 + t36 * t12) * t314 + t230;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
