% Calculate time derivative of joint inertia matrix for
% S5RPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% m_mdh [6x1]
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
% Datum: 2019-12-05 18:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR5_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR5_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR5_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR5_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR5_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:16:07
% EndTime: 2019-12-05 18:16:16
% DurationCPUTime: 3.55s
% Computational Cost: add. (9111->334), mult. (5758->481), div. (0->0), fcn. (4324->10), ass. (0->206)
t161 = qJ(1) + pkin(9);
t156 = qJ(3) + t161;
t150 = sin(t156);
t162 = qJ(4) + qJ(5);
t158 = cos(t162);
t268 = rSges(6,1) * t158;
t151 = cos(t156);
t157 = sin(t162);
t266 = rSges(6,2) * t157;
t293 = t151 * rSges(6,3) + t150 * t266;
t77 = -t150 * t268 + t293;
t298 = sin(qJ(1)) * pkin(1) + pkin(2) * sin(t161);
t165 = cos(qJ(4));
t235 = qJD(4) * t165;
t160 = qJD(1) + qJD(3);
t163 = sin(qJ(4));
t243 = t160 * t163;
t297 = t150 * t235 + t151 * t243;
t159 = qJD(4) + qJD(5);
t244 = t158 * t159;
t245 = t157 * t159;
t296 = rSges(6,1) * t245 + rSges(6,2) * t244;
t260 = Icges(5,4) * t165;
t195 = -Icges(5,2) * t163 + t260;
t83 = Icges(5,6) * t151 - t195 * t150;
t261 = Icges(5,4) * t163;
t197 = Icges(5,1) * t165 - t261;
t85 = Icges(5,5) * t151 - t197 * t150;
t200 = t163 * t83 - t165 * t85;
t295 = t200 * t151;
t258 = Icges(6,4) * t158;
t194 = -Icges(6,2) * t157 + t258;
t73 = Icges(6,6) * t151 - t194 * t150;
t259 = Icges(6,4) * t157;
t196 = Icges(6,1) * t158 - t259;
t75 = Icges(6,5) * t151 - t196 * t150;
t205 = t157 * t73 - t158 * t75;
t294 = t205 * t151;
t116 = Icges(6,1) * t157 + t258;
t249 = t116 * t159;
t292 = -Icges(6,5) * t160 + t249;
t115 = Icges(6,2) * t158 + t259;
t250 = t115 * t159;
t291 = -Icges(6,6) * t160 + t250;
t114 = Icges(6,5) * t157 + Icges(6,6) * t158;
t290 = -Icges(6,3) * t160 + t114 * t159;
t289 = -t296 * t151 + t77 * t160;
t138 = Icges(5,5) * t163 + Icges(5,6) * t165;
t288 = -Icges(5,3) * t160 + t138 * qJD(4);
t139 = Icges(5,2) * t165 + t261;
t287 = -Icges(5,6) * t160 + t139 * qJD(4);
t140 = Icges(5,1) * t163 + t260;
t286 = -Icges(5,5) * t160 + t140 * qJD(4);
t121 = t195 * qJD(4);
t122 = t197 * qJD(4);
t285 = t121 * t163 - t122 * t165 - t138 * t160 + (t139 * t165 + t140 * t163) * qJD(4);
t97 = t194 * t159;
t98 = t196 * t159;
t284 = (t97 + t249) * t157 + (-t98 + t250) * t158 - t114 * t160;
t283 = 2 * m(4);
t282 = 2 * m(5);
t281 = 2 * m(6);
t280 = t150 / 0.2e1;
t279 = t151 / 0.2e1;
t278 = -rSges(5,3) - pkin(7);
t267 = rSges(5,2) * t163;
t269 = rSges(5,1) * t165;
t125 = (-t267 + t269) * qJD(4);
t277 = m(5) * t125;
t141 = rSges(5,1) * t163 + rSges(5,2) * t165;
t276 = m(5) * t141;
t274 = t150 * pkin(7);
t147 = t151 * pkin(7);
t152 = pkin(4) * t165 + pkin(3);
t271 = pkin(3) - t152;
t217 = t150 * t271;
t167 = -pkin(8) - pkin(7);
t246 = t151 * t167;
t270 = t147 - t217 + t246 - t77;
t265 = t150 * rSges(5,3);
t264 = t150 * rSges(6,3);
t146 = t151 * rSges(5,3);
t248 = t150 * t160;
t247 = t151 * t160;
t242 = t160 * t167;
t236 = qJD(4) * t163;
t224 = t150 * t236;
t241 = -pkin(4) * t224 - t150 * t242;
t93 = rSges(4,1) * t248 + rSges(4,2) * t247;
t136 = t150 * t267;
t239 = t136 + t146;
t238 = t298 * qJD(1);
t237 = t150 ^ 2 + t151 ^ 2;
t234 = t150 * t269;
t124 = t151 * t266;
t222 = t151 * t236;
t229 = -pkin(4) * t222 - t151 * t242 - t152 * t248;
t228 = t160 * t124 + t296 * t150;
t221 = t151 * t235;
t226 = -rSges(5,1) * t222 - rSges(5,2) * t221 - t160 * t234;
t225 = rSges(5,1) * t224 + t297 * rSges(5,2);
t220 = -t248 / 0.2e1;
t219 = t247 / 0.2e1;
t218 = -pkin(3) - t269;
t117 = rSges(6,1) * t157 + rSges(6,2) * t158;
t216 = pkin(4) * t163 + t117;
t178 = t194 * t160;
t45 = -t150 * t178 - t291 * t151;
t76 = Icges(6,5) * t150 + t196 * t151;
t215 = t159 * t76 + t45;
t46 = t291 * t150 - t151 * t178;
t214 = t159 * t75 + t46;
t180 = t196 * t160;
t47 = -t150 * t180 - t292 * t151;
t74 = Icges(6,6) * t150 + t194 * t151;
t213 = -t159 * t74 + t47;
t48 = t292 * t150 - t151 * t180;
t212 = -t159 * t73 + t48;
t101 = -rSges(4,1) * t151 + t150 * rSges(4,2);
t211 = -t152 - t268;
t94 = -rSges(4,1) * t247 + rSges(4,2) * t248;
t207 = -pkin(2) * cos(t161) - cos(qJ(1)) * pkin(1);
t100 = -rSges(4,1) * t150 - rSges(4,2) * t151;
t204 = t157 * t74 - t158 * t76;
t201 = t163 * t85 + t165 * t83;
t84 = Icges(5,6) * t150 + t195 * t151;
t86 = Icges(5,5) * t150 + t197 * t151;
t199 = t163 * t86 + t165 * t84;
t198 = t163 * t84 - t165 * t86;
t193 = Icges(5,5) * t165 - Icges(5,6) * t163;
t192 = Icges(6,5) * t158 - Icges(6,6) * t157;
t191 = t115 * t157 - t116 * t158;
t189 = t139 * t163 - t140 * t165;
t71 = Icges(6,3) * t151 - t192 * t150;
t17 = t205 * t150 + t151 * t71;
t184 = t204 * t150;
t72 = Icges(6,3) * t150 + t192 * t151;
t18 = t151 * t72 + t184;
t19 = t150 * t71 - t294;
t20 = t150 * t72 - t204 * t151;
t176 = t192 * t160;
t43 = -t150 * t176 - t290 * t151;
t44 = t290 * t150 - t151 * t176;
t188 = -(t150 * t18 + t151 * t17) * t248 + t150 * ((t150 * t43 + (-t19 + t184) * t160) * t150 + (t20 * t160 + (-t157 * t46 + t158 * t48 - t73 * t244 - t75 * t245) * t151 + (t44 + (-t160 * t75 + t213) * t158 + (t160 * t73 - t215) * t157) * t150) * t151) + t151 * ((t151 * t44 + (t18 + t294) * t160) * t151 + (-t17 * t160 + (t157 * t45 - t158 * t47 + t74 * t244 + t76 * t245) * t150 + (t43 + (-t160 * t76 - t212) * t158 + (t160 * t74 + t214) * t157) * t151) * t150) + (t150 * t20 + t151 * t19) * t247;
t187 = -t151 * t268 - t264;
t185 = t207 * qJD(1);
t183 = t198 * t150;
t172 = t192 * t159 + t191 * t160;
t182 = (t172 * t150 - t284 * t151 + t213 * t157 + t215 * t158) * t280 + (t150 * t284 + t172 * t151 + t212 * t157 + t214 * t158) * t279 + (t151 * t114 + t191 * t150 + t157 * t75 + t158 * t73) * t220 + (t150 * t114 - t191 * t151 + t157 * t76 + t158 * t74) * t219;
t181 = t197 * t160;
t179 = t195 * t160;
t177 = t193 * t160;
t173 = t211 * t151 - t264;
t171 = t193 * qJD(4) + t189 * t160;
t64 = t218 * t150 + t147 + t239;
t170 = t278 * t150 + t218 * t151;
t137 = t151 * t267;
t65 = t137 + t170;
t142 = t150 * t167;
t63 = t124 + t142 + t173;
t62 = t211 * t150 - t246 + t293;
t169 = -t115 * t245 + t116 * t244 + t165 * t121 + t163 * t122 - t139 * t236 + t140 * t235 + t157 * t98 + t158 * t97;
t135 = pkin(3) * t248;
t33 = t135 + (t278 * t151 - t136) * t160 - t226;
t168 = t182 + (-t198 * qJD(4) + t171 * t150 - t285 * t151 + t163 * (-t150 * t181 - t286 * t151) + t165 * (-t150 * t179 - t287 * t151)) * t280 + (-t200 * qJD(4) + t285 * t150 + t171 * t151 + t163 * (t286 * t150 - t151 * t181) + t165 * (t287 * t150 - t151 * t179)) * t279 + (t151 * t138 + t189 * t150 + t201) * t220 + (t150 * t138 - t189 * t151 + t199) * t219;
t27 = -t229 - t289;
t34 = t170 * t160 + t225;
t28 = t173 * t160 + t228 - t241;
t99 = (-t266 + t268) * t159;
t90 = t101 + t207;
t89 = t100 - t298;
t88 = t151 * t269 - t137 + t265;
t87 = -t234 + t239;
t82 = Icges(5,3) * t150 + t193 * t151;
t81 = Icges(5,3) * t151 - t193 * t150;
t80 = t216 * t151;
t79 = t216 * t150;
t78 = -t124 - t187;
t70 = t185 + t94;
t69 = t238 + t93;
t68 = -t271 * t151 - t142 - t274;
t66 = t151 * t78;
t61 = t207 + t65;
t60 = -t298 + t64;
t55 = t288 * t150 - t151 * t177;
t54 = -t150 * t177 - t288 * t151;
t53 = t207 + t63;
t52 = -t298 + t62;
t49 = t187 * t160 + t228;
t42 = t297 * pkin(4) + t117 * t247 + t150 * t99;
t41 = t117 * t248 - t151 * t99 + (t150 * t243 - t221) * pkin(4);
t38 = -t150 * t77 + t66;
t37 = t151 * t289;
t30 = t185 + t34;
t29 = t33 + t238;
t26 = t150 * t82 - t198 * t151;
t25 = t150 * t81 - t295;
t24 = t151 * t82 + t183;
t23 = t200 * t150 + t151 * t81;
t22 = t185 + t28;
t21 = t27 + t238;
t16 = t270 * t150 + t151 * t68 + t66;
t11 = t151 * t226 - t150 * t225 + ((-t87 + t146) * t151 + (t265 - t88 + (t267 + t269) * t151) * t150) * t160;
t8 = -t77 * t247 + t37 + (-t160 * t78 - t49) * t150;
t3 = t151 * (t135 + t229) + t37 + (-t49 + t241) * t150 + ((-t68 - t78 - t274) * t150 + (-t217 + t270 - t147) * t151) * t160;
t1 = [(t69 * t90 + t70 * t89) * t283 + (t29 * t61 + t30 * t60) * t282 + (t21 * t53 + t22 * t52) * t281 + t169; 0; 0; m(4) * (t100 * t70 + t101 * t69 + t89 * t94 + t90 * t93) + m(5) * (t29 * t65 + t30 * t64 + t33 * t61 + t34 * t60) + m(6) * (t21 * t63 + t22 * t62 + t27 * t53 + t28 * t52) + t169; 0; (t27 * t63 + t28 * t62) * t281 + (t33 * t65 + t34 * t64) * t282 + (t100 * t94 + t101 * t93) * t283 + t169; ((t160 * t61 - t30) * t151 + (t160 * t60 + t29) * t150) * t276 + (t150 * t61 - t151 * t60) * t277 + m(6) * (t21 * t79 - t22 * t80 + t41 * t52 + t42 * t53) + t168; m(5) * t11 + m(6) * t3; ((t160 * t65 - t34) * t151 + (t160 * t64 + t33) * t150) * t276 + (t150 * t65 - t151 * t64) * t277 + m(6) * (t27 * t79 - t28 * t80 + t41 * t62 + t42 * t63) + t168; ((-t150 * t87 + t151 * t88) * t11 + t237 * t141 * t125) * t282 - (t150 * t24 + t151 * t23) * t248 + t151 * ((t151 * t55 + (t24 + t295) * t160) * t151 + (-t23 * t160 + (t84 * t235 + t86 * t236) * t150 + (t201 * qJD(4) + t198 * t160 + t54) * t151) * t150) + (t150 * t26 + t151 * t25) * t247 + t150 * ((t150 * t54 + (-t25 + t183) * t160) * t150 + (t26 * t160 + (-t83 * t235 - t85 * t236) * t151 + (-t199 * qJD(4) + t200 * t160 + t55) * t150) * t151) + (t16 * t3 - t41 * t80 + t42 * t79) * t281 + t188; m(6) * ((t150 * t53 - t151 * t52) * t99 + ((t160 * t53 - t22) * t151 + (t160 * t52 + t21) * t150) * t117) + t182; m(6) * t8; m(6) * ((t150 * t63 - t151 * t62) * t99 + ((t160 * t63 - t28) * t151 + (t160 * t62 + t27) * t150) * t117) + t182; m(6) * (t16 * t8 + t3 * t38 + (t150 * t79 + t151 * t80) * t99 + ((t160 * t79 - t41) * t151 + (-t160 * t80 + t42) * t150) * t117) + t188; (t237 * t99 * t117 + t38 * t8) * t281 + t188;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
