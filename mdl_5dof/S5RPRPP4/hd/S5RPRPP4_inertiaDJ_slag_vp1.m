% Calculate time derivative of joint inertia matrix for
% S5RPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta4]';
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
% Datum: 2019-12-31 18:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPP4_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP4_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP4_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP4_inertiaDJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP4_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPP4_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPP4_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:14:18
% EndTime: 2019-12-31 18:14:32
% DurationCPUTime: 8.71s
% Computational Cost: add. (2975->400), mult. (4800->571), div. (0->0), fcn. (3614->6), ass. (0->226)
t163 = cos(qJ(1));
t156 = qJ(3) + pkin(7);
t146 = cos(t156);
t284 = rSges(6,3) + qJ(5);
t227 = t284 * t146;
t145 = sin(t156);
t307 = rSges(6,1) + pkin(4);
t236 = t307 * t145;
t332 = t227 - t236;
t341 = t163 * t332;
t276 = Icges(6,5) * t145;
t281 = Icges(5,4) * t145;
t340 = t276 - t281 + (Icges(5,1) + Icges(6,1)) * t146;
t339 = t145 * t284;
t246 = qJD(5) * t146;
t261 = t146 * t307 + t339;
t338 = (qJD(3) * t261 - t246) * t163;
t337 = t340 * qJD(3);
t161 = sin(qJ(1));
t336 = t161 / 0.2e1;
t247 = qJD(3) * t163;
t335 = -t247 / 0.2e1;
t334 = -qJD(1) / 0.2e1;
t333 = qJD(1) / 0.2e1;
t160 = sin(qJ(3));
t162 = cos(qJ(3));
t248 = qJD(3) * t162;
t253 = qJD(1) * t163;
t330 = t160 * t253 + t161 * t248;
t249 = qJD(3) * t161;
t329 = -t145 * t253 - t146 * t249;
t294 = rSges(5,2) * t146;
t219 = rSges(5,1) * t145 + t294;
t177 = t163 * t219;
t220 = rSges(4,1) * t160 + rSges(4,2) * t162;
t178 = t163 * t220;
t283 = Icges(4,4) * t160;
t193 = Icges(4,2) * t162 + t283;
t77 = Icges(4,6) * t163 + t161 * t193;
t282 = Icges(4,4) * t162;
t199 = Icges(4,1) * t160 + t282;
t79 = Icges(4,5) * t163 + t161 * t199;
t206 = t160 * t79 + t162 * t77;
t172 = t206 * t163;
t191 = Icges(5,2) * t146 + t281;
t63 = Icges(5,6) * t163 + t161 * t191;
t280 = Icges(5,4) * t146;
t197 = Icges(5,1) * t145 + t280;
t67 = Icges(5,5) * t163 + t161 * t197;
t210 = t145 * t67 + t146 * t63;
t174 = t210 * t163;
t183 = -Icges(6,3) * t146 + t276;
t57 = Icges(6,6) * t163 + t161 * t183;
t275 = Icges(6,5) * t146;
t195 = Icges(6,1) * t145 - t275;
t65 = Icges(6,4) * t163 + t161 * t195;
t212 = t145 * t65 - t146 * t57;
t176 = t212 * t163;
t154 = t163 * rSges(6,2);
t265 = t145 * t161;
t326 = t265 * t307 + t154;
t232 = t162 * t253;
t238 = rSges(4,1) * t330 + rSges(4,2) * t232;
t245 = -rSges(4,3) - pkin(1) - pkin(6);
t250 = qJD(3) * t160;
t257 = qJ(2) * t253 + qJD(2) * t161;
t23 = (-rSges(4,2) * t250 + qJD(1) * t245) * t161 + t238 + t257;
t296 = rSges(4,2) * t160;
t122 = rSges(4,1) * t162 - t296;
t149 = qJD(2) * t163;
t24 = t149 + t122 * t247 + (t245 * t163 + (-qJ(2) - t220) * t161) * qJD(1);
t325 = t161 * t24 - t163 * t23;
t324 = -Icges(6,6) * t161 + t163 * t183;
t185 = Icges(5,5) * t145 + Icges(5,6) * t146;
t323 = -Icges(5,3) * t161 + t163 * t185;
t187 = Icges(4,5) * t160 + Icges(4,6) * t162;
t322 = -Icges(4,3) * t161 + t163 * t187;
t189 = Icges(6,4) * t145 - Icges(6,6) * t146;
t321 = -Icges(6,2) * t161 + t163 * t189;
t320 = -Icges(5,6) * t161 + t163 * t191;
t319 = -Icges(4,6) * t161 + t163 * t193;
t318 = -Icges(6,4) * t161 + t163 * t195;
t317 = -Icges(5,5) * t161 + t163 * t197;
t316 = -Icges(4,5) * t161 + t163 * t199;
t262 = t161 * t162;
t141 = pkin(3) * t262;
t41 = t161 * t261 + t141;
t301 = pkin(3) * t162;
t42 = (-t261 - t301) * t163;
t225 = qJD(1) * t261;
t241 = pkin(3) * t247;
t259 = qJD(1) * t141 + t160 * t241;
t299 = -qJD(3) * t332 - qJD(5) * t145;
t5 = t161 * t225 + t163 * t299 + t259;
t135 = pkin(3) * t232;
t242 = pkin(3) * t250;
t6 = t135 + t163 * t225 + (-t242 - t299) * t161;
t315 = qJD(1) * (t161 * t42 + t163 * t41) + t6 * t161 - t163 * t5;
t164 = t161 * t246 - t249 * t339 + t307 * t329;
t159 = -qJ(4) - pkin(6);
t254 = qJD(1) * t161;
t224 = pkin(3) * t330 + qJD(4) * t163 + t159 * t254;
t180 = t224 + t257;
t305 = -rSges(6,2) - pkin(1);
t235 = t305 * t161;
t2 = (-t163 * t227 + t235) * qJD(1) - t164 + t180;
t151 = t163 * qJ(2);
t302 = pkin(3) * t160;
t258 = t161 * t159 + t163 * t302;
t237 = t151 + t258;
t21 = t235 + t237 - t341;
t263 = t160 * t161;
t140 = pkin(3) * t263;
t256 = t163 * pkin(1) + t161 * qJ(2);
t182 = -t159 * t163 + t140 + t256;
t22 = -t161 * t227 + t182 + t326;
t221 = qJD(4) * t161 - t159 * t253 - t162 * t241;
t181 = t149 - t221;
t228 = -qJ(2) - t302;
t3 = t338 + (t305 * t163 + (t228 + t332) * t161) * qJD(1) + t181;
t314 = qJD(1) * (t161 * t22 + t163 * t21) + t161 * t3 - t163 * t2;
t312 = 2 * m(4);
t311 = 2 * m(5);
t310 = 2 * m(6);
t157 = t161 ^ 2;
t158 = t163 ^ 2;
t309 = -t160 / 0.2e1;
t308 = t162 / 0.2e1;
t306 = rSges(3,2) - pkin(1);
t304 = -rSges(5,3) - pkin(1);
t303 = m(4) * t122;
t300 = pkin(6) * t163;
t152 = t163 * rSges(5,3);
t264 = t146 * t161;
t71 = rSges(5,1) * t265 + rSges(5,2) * t264 + t152;
t88 = t140 + (-pkin(6) - t159) * t163;
t298 = -t71 - t88;
t291 = t161 * rSges(6,2);
t297 = t291 + t341;
t295 = rSges(5,2) * t145;
t293 = t160 * t77;
t292 = t160 * t319;
t290 = t161 * rSges(4,3);
t289 = t161 * rSges(5,3);
t287 = t162 * t79;
t286 = t162 * t316;
t153 = t163 * rSges(4,3);
t59 = Icges(5,3) * t163 + t161 * t185;
t268 = qJD(1) * t59;
t61 = Icges(6,2) * t163 + t161 * t189;
t267 = qJD(1) * t61;
t75 = Icges(4,3) * t163 + t161 * t187;
t266 = qJD(1) * t75;
t255 = t157 + t158;
t252 = qJD(3) * t145;
t251 = qJD(3) * t146;
t244 = t264 * t284 - t326 - t88;
t243 = m(6) * t252;
t240 = -(t185 + t189) * qJD(3) / 0.2e1;
t239 = rSges(5,1) * t329 - t253 * t294;
t83 = rSges(4,1) * t263 + rSges(4,2) * t262 + t153;
t109 = t220 * qJD(3);
t226 = t255 * t109;
t105 = rSges(5,1) * t146 - t295;
t211 = -t145 * t318 + t146 * t324;
t209 = -t145 * t317 - t146 * t320;
t205 = -t160 * t316 - t162 * t319;
t203 = t161 * t21 - t163 * t22;
t201 = t161 * t41 - t163 * t42;
t200 = Icges(4,1) * t162 - t283;
t194 = -Icges(4,2) * t160 + t282;
t192 = -Icges(5,2) * t145 + t280;
t190 = Icges(6,4) * t146 + Icges(6,6) * t145;
t188 = Icges(4,5) * t162 - Icges(4,6) * t160;
t186 = Icges(5,5) * t146 - Icges(5,6) * t145;
t184 = Icges(6,3) * t145 + t275;
t179 = rSges(3,3) * t163 + t161 * t306;
t175 = t211 * t161;
t173 = t209 * t161;
t171 = t205 * t161;
t170 = qJD(3) * t200;
t167 = qJD(3) * t194;
t166 = qJD(3) * t192;
t165 = qJD(3) * t184;
t96 = t219 * qJD(3);
t87 = -t161 * pkin(6) - t258;
t86 = -rSges(3,2) * t163 + t161 * rSges(3,3) + t256;
t85 = t151 + t179;
t84 = t290 - t178;
t74 = t163 * t87;
t73 = t289 - t177;
t56 = (-t105 - t301) * t163;
t55 = t105 * t161 + t141;
t54 = t149 + (t306 * t163 + (-rSges(3,3) - qJ(2)) * t161) * qJD(1);
t53 = qJD(1) * t179 + t257;
t52 = t256 + t83 + t300;
t51 = t161 * t245 + t151 + t178;
t50 = pkin(6) * t254 + t224;
t49 = t163 * ((t140 - t300) * qJD(1) + t221);
t44 = qJD(1) * t322 + t188 * t249;
t43 = -t188 * t247 + t266;
t40 = t182 + t71;
t39 = t161 * t304 + t177 + t237;
t32 = qJD(1) * t321 + t190 * t249;
t31 = -t190 * t247 + t267;
t30 = qJD(1) * t323 + t186 * t249;
t29 = -t186 * t247 + t268;
t26 = t105 * t253 + t135 + (-t96 - t242) * t161;
t25 = t105 * t254 + t163 * t96 + t259;
t20 = -t161 * t322 - t163 * t205;
t19 = t161 * t75 - t172;
t18 = -t163 * t322 + t171;
t17 = t206 * t161 + t163 * t75;
t16 = -t161 * t323 - t163 * t209;
t15 = t161 * t59 - t174;
t14 = -t161 * t321 - t163 * t211;
t13 = t161 * t61 - t176;
t12 = -t163 * t323 + t173;
t11 = t210 * t161 + t163 * t59;
t10 = -t163 * t321 + t175;
t9 = t212 * t161 + t163 * t61;
t8 = t105 * t247 + (t304 * t163 + (-t219 + t228) * t161) * qJD(1) + t181;
t7 = (-rSges(5,2) * t252 + qJD(1) * t304) * t161 + t180 - t239;
t4 = t161 * t244 + t163 * t297 + t74;
t1 = t49 + (-t50 + (-t87 + t291 - t297) * qJD(1) + t164) * t161 + (-t338 + (t161 * t236 + t154 + t244) * qJD(1)) * t163;
t27 = [(t2 * t22 + t21 * t3) * t310 - t160 * t170 - t199 * t248 - t162 * t167 + t193 * t250 + (t39 * t8 + t40 * t7) * t311 + (t23 * t52 + t24 * t51) * t312 + 0.2e1 * m(3) * (t53 * t86 + t54 * t85) + (t191 - t183) * t252 + (-t197 - t195) * t251 + (-t166 + t165) * t146 - t337 * t145; m(6) * t314 + m(5) * (t161 * t8 - t163 * t7 + (t161 * t40 + t163 * t39) * qJD(1)) + m(4) * ((t161 * t52 + t163 * t51) * qJD(1) + t325) + m(3) * (t161 * t54 - t163 * t53 + (t161 * t86 + t163 * t85) * qJD(1)); 0; (t240 * t163 + (qJD(1) * t316 + t161 * t170) * t308 + (qJD(1) * t319 + t161 * t167) * t309) * t163 + (t240 * t161 + (qJD(1) * t79 - t200 * t247) * t308 + (qJD(1) * t77 - t194 * t247) * t309) * t161 + m(4) * (t325 * t122 - (t161 * t51 - t163 * t52) * t109) + m(5) * (t25 * t40 + t26 * t39 + t55 * t8 + t56 * t7) + m(6) * (t2 * t42 + t21 * t6 + t22 * t5 + t3 * t41) + (t337 * t163 * t336 + t340 * t161 * t335 + ((t317 + t318) * t163 + (t65 + t67) * t161) * t333) * t146 + ((t324 * t333 + t165 * t336 + t320 * t334 - t161 * t166 / 0.2e1) * t163 + (t57 * t333 + t184 * t335 + t63 * t334 + t192 * t247 / 0.2e1) * t161) * t145 + ((t52 * t303 + t293 / 0.2e1 - t287 / 0.2e1 + (-t65 / 0.2e1 - t67 / 0.2e1) * t146 + (-t57 / 0.2e1 + t63 / 0.2e1) * t145) * t161 + (t292 / 0.2e1 - t286 / 0.2e1 + t51 * t303 + (-t318 / 0.2e1 - t317 / 0.2e1) * t146 + (-t324 / 0.2e1 + t320 / 0.2e1) * t145) * t163) * qJD(1) + (-(t157 / 0.2e1 + t158 / 0.2e1) * t187 - t172 / 0.2e1 - t176 / 0.2e1 - t174 / 0.2e1 + (-t205 / 0.2e1 - t211 / 0.2e1 - t209 / 0.2e1) * t161) * qJD(3); m(5) * (t26 * t161 - t163 * t25 + (t161 * t56 + t163 * t55) * qJD(1)) + m(6) * t315 - m(4) * t226; (t1 * t4 + t41 * t6 + t42 * t5) * t310 + (t55 * t26 + t56 * t25 + (t161 * t298 + t163 * t73 + t74) * (t49 + (-t50 + t239) * t161 + (-t105 * t158 + t157 * t295) * qJD(3) + ((t152 + t298) * t163 + (t177 - t73 - t87 + t289) * t161) * qJD(1))) * t311 + t163 * ((t163 * t44 + (t18 + t172) * qJD(1)) * t163 + (-t17 * qJD(1) + (-t248 * t316 + t250 * t319) * t161 + (t43 + (t287 - t293) * qJD(3) + (t205 - t75) * qJD(1)) * t163) * t161) + t163 * ((t163 * t32 + (t10 + t176) * qJD(1)) * t163 + (-t9 * qJD(1) + (-t251 * t318 - t252 * t324) * t161 + (t31 + (t145 * t57 + t146 * t65) * qJD(3) + (t211 - t61) * qJD(1)) * t163) * t161) + t163 * ((t163 * t30 + (t12 + t174) * qJD(1)) * t163 + (-t11 * qJD(1) + (-t251 * t317 + t252 * t320) * t161 + (t29 + (-t145 * t63 + t146 * t67) * qJD(3) + (t209 - t59) * qJD(1)) * t163) * t161) + t161 * ((t161 * t29 + (-t15 + t173) * qJD(1)) * t161 + (t16 * qJD(1) + (-t251 * t67 + t252 * t63 + t268) * t163 + (t30 + (-t145 * t320 + t146 * t317) * qJD(3) + t210 * qJD(1)) * t161) * t163) + ((-t161 * t83 + t163 * t84) * (-t161 * t238 + (-t122 * t158 + t157 * t296) * qJD(3) + ((-t83 + t153) * t163 + (t178 - t84 + t290) * t161) * qJD(1)) - t122 * t226) * t312 + t161 * ((t161 * t43 + (-t19 + t171) * qJD(1)) * t161 + (t20 * qJD(1) + (-t248 * t79 + t250 * t77 + t266) * t163 + (t44 + (t286 - t292) * qJD(3) + t206 * qJD(1)) * t161) * t163) + t161 * ((t161 * t31 + (-t13 + t175) * qJD(1)) * t161 + (t14 * qJD(1) + (-t251 * t65 - t252 * t57 + t267) * t163 + (t32 + (t145 * t324 + t146 * t318) * qJD(3) + t212 * qJD(1)) * t161) * t163) + ((-t11 - t17 - t9) * t163 + (-t10 - t12 - t18) * t161) * t254 + ((t13 + t15 + t19) * t163 + (t14 + t16 + t20) * t161) * t253; m(6) * (-qJD(1) * t203 + t161 * t2 + t163 * t3) + m(5) * (t161 * t7 + t163 * t8 + (-t161 * t39 + t163 * t40) * qJD(1)); 0; m(6) * (-qJD(1) * t201 + t161 * t5 + t163 * t6) + m(5) * (t161 * t25 + t163 * t26 + (-t161 * t55 + t163 * t56) * qJD(1)); 0; m(6) * (-t146 * t314 + t203 * t252); t255 * t243; m(6) * ((qJD(3) * t201 + t1) * t145 + (qJD(3) * t4 - t315) * t146); 0; 0.2e1 * (0.1e1 - t255) * t146 * t243;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t27(1), t27(2), t27(4), t27(7), t27(11); t27(2), t27(3), t27(5), t27(8), t27(12); t27(4), t27(5), t27(6), t27(9), t27(13); t27(7), t27(8), t27(9), t27(10), t27(14); t27(11), t27(12), t27(13), t27(14), t27(15);];
Mq = res;
