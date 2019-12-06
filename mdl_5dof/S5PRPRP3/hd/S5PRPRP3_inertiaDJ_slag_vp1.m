% Calculate time derivative of joint inertia matrix for
% S5PRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRP3_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP3_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP3_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP3_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP3_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRP3_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRP3_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:32:40
% EndTime: 2019-12-05 15:32:58
% DurationCPUTime: 7.41s
% Computational Cost: add. (8712->408), mult. (13067->629), div. (0->0), fcn. (12237->8), ass. (0->220)
t169 = qJ(2) + pkin(8);
t166 = cos(t169);
t171 = cos(pkin(7));
t174 = sin(qJ(4));
t255 = t171 * t174;
t170 = sin(pkin(7));
t176 = cos(qJ(4));
t256 = t170 * t176;
t153 = -t166 * t255 + t256;
t254 = t171 * t176;
t257 = t170 * t174;
t154 = t166 * t254 + t257;
t165 = sin(t169);
t258 = t165 * t171;
t94 = Icges(6,5) * t154 + Icges(6,6) * t153 + Icges(6,3) * t258;
t96 = Icges(5,5) * t154 + Icges(5,6) * t153 + Icges(5,3) * t258;
t332 = -t96 - t94;
t100 = Icges(5,4) * t154 + Icges(5,2) * t153 + Icges(5,6) * t258;
t98 = Icges(6,4) * t154 + Icges(6,2) * t153 + Icges(6,6) * t258;
t341 = t100 + t98;
t102 = Icges(6,1) * t154 + Icges(6,4) * t153 + Icges(6,5) * t258;
t104 = Icges(5,1) * t154 + Icges(5,4) * t153 + Icges(5,5) * t258;
t340 = t102 + t104;
t151 = -t166 * t257 - t254;
t152 = t166 * t256 - t255;
t259 = t165 * t170;
t97 = Icges(6,4) * t152 + Icges(6,2) * t151 + Icges(6,6) * t259;
t99 = Icges(5,4) * t152 + Icges(5,2) * t151 + Icges(5,6) * t259;
t337 = t97 + t99;
t101 = Icges(6,1) * t152 + Icges(6,4) * t151 + Icges(6,5) * t259;
t103 = Icges(5,1) * t152 + Icges(5,4) * t151 + Icges(5,5) * t259;
t336 = t101 + t103;
t339 = t341 * t151 + t340 * t152 - t332 * t259;
t93 = Icges(6,5) * t152 + Icges(6,6) * t151 + Icges(6,3) * t259;
t95 = Icges(5,5) * t152 + Icges(5,6) * t151 + Icges(5,3) * t259;
t338 = t95 + t93;
t208 = Icges(6,5) * t176 - Icges(6,6) * t174;
t119 = -Icges(6,3) * t166 + t165 * t208;
t209 = Icges(5,5) * t176 - Icges(5,6) * t174;
t120 = -Icges(5,3) * t166 + t165 * t209;
t335 = t119 + t120;
t264 = Icges(6,4) * t176;
t210 = -Icges(6,2) * t174 + t264;
t121 = -Icges(6,6) * t166 + t165 * t210;
t266 = Icges(5,4) * t176;
t211 = -Icges(5,2) * t174 + t266;
t122 = -Icges(5,6) * t166 + t165 * t211;
t331 = t121 + t122;
t265 = Icges(6,4) * t174;
t214 = Icges(6,1) * t176 - t265;
t123 = -Icges(6,5) * t166 + t165 * t214;
t267 = Icges(5,4) * t174;
t215 = Icges(5,1) * t176 - t267;
t124 = -Icges(5,5) * t166 + t165 * t215;
t334 = t123 + t124;
t330 = -t341 * t174 + t340 * t176;
t329 = -t337 * t174 + t336 * t176;
t309 = t337 * t151 + t336 * t152 + t338 * t259;
t328 = t331 * t151 + t334 * t152 + t335 * t259;
t248 = qJD(4) * t165;
t87 = (-Icges(6,2) * t176 - t265) * t248 + (Icges(6,6) * t165 + t166 * t210) * qJD(2);
t88 = (-Icges(5,2) * t176 - t267) * t248 + (Icges(5,6) * t165 + t166 * t211) * qJD(2);
t327 = -t87 - t88;
t89 = (-Icges(6,1) * t174 - t264) * t248 + (Icges(6,5) * t165 + t166 * t214) * qJD(2);
t90 = (-Icges(5,1) * t174 - t266) * t248 + (Icges(5,5) * t165 + t166 * t215) * qJD(2);
t326 = t89 + t90;
t175 = sin(qJ(2));
t177 = cos(qJ(2));
t323 = (-Icges(3,5) * t175 - Icges(4,5) * t165 - Icges(3,6) * t177 - Icges(4,6) * t166) * qJD(2);
t262 = t120 * t166;
t263 = t119 * t166;
t322 = -t262 - t263;
t321 = t339 * t171;
t278 = t166 * ((-Icges(5,5) * t174 - Icges(5,6) * t176) * t248 + (Icges(5,3) * t165 + t166 * t209) * qJD(2));
t279 = t166 * ((-Icges(6,5) * t174 - Icges(6,6) * t176) * t248 + (Icges(6,3) * t165 + t166 * t208) * qJD(2));
t320 = -t278 - t279;
t319 = t329 * t165 - t166 * t338;
t318 = t330 * t165 + t332 * t166;
t167 = t170 ^ 2;
t168 = t171 ^ 2;
t317 = t167 + t168;
t316 = t331 * t174 - t176 * t334;
t299 = t317 * qJD(2);
t250 = qJD(2) * t165;
t240 = t174 * t250;
t110 = -qJD(4) * t152 + t170 * t240;
t239 = t176 * t250;
t111 = qJD(4) * t151 - t170 * t239;
t249 = qJD(2) * t166;
t238 = t170 * t249;
t64 = Icges(6,5) * t111 + Icges(6,6) * t110 + Icges(6,3) * t238;
t195 = t165 * t64 + t249 * t93;
t68 = Icges(6,4) * t111 + Icges(6,2) * t110 + Icges(6,6) * t238;
t72 = Icges(6,1) * t111 + Icges(6,4) * t110 + Icges(6,5) * t238;
t14 = t101 * t111 + t110 * t97 + t151 * t68 + t152 * t72 + t170 * t195;
t112 = -qJD(4) * t154 + t171 * t240;
t113 = qJD(4) * t153 - t171 * t239;
t237 = t171 * t249;
t65 = Icges(6,5) * t113 + Icges(6,6) * t112 + Icges(6,3) * t237;
t194 = t165 * t65 + t249 * t94;
t69 = Icges(6,4) * t113 + Icges(6,2) * t112 + Icges(6,6) * t237;
t73 = Icges(6,1) * t113 + Icges(6,4) * t112 + Icges(6,5) * t237;
t15 = t102 * t111 + t110 * t98 + t151 * t69 + t152 * t73 + t170 * t194;
t66 = Icges(5,5) * t111 + Icges(5,6) * t110 + Icges(5,3) * t238;
t193 = t165 * t66 + t249 * t95;
t70 = Icges(5,4) * t111 + Icges(5,2) * t110 + Icges(5,6) * t238;
t74 = Icges(5,1) * t111 + Icges(5,4) * t110 + Icges(5,5) * t238;
t16 = t103 * t111 + t110 * t99 + t151 * t70 + t152 * t74 + t170 * t193;
t67 = Icges(5,5) * t113 + Icges(5,6) * t112 + Icges(5,3) * t237;
t192 = t165 * t67 + t249 * t96;
t71 = Icges(5,4) * t113 + Icges(5,2) * t112 + Icges(5,6) * t237;
t75 = Icges(5,1) * t113 + Icges(5,4) * t112 + Icges(5,5) * t237;
t17 = t100 * t110 + t104 * t111 + t151 * t71 + t152 * t75 + t170 * t192;
t315 = (-t110 * t331 - t111 * t334 + t327 * t151 - t326 * t152) * t166 + ((t17 + t15) * t171 + (t16 + t14 + t320) * t170) * t165 + (((t322 + t309) * t170 + t321) * t166 + t328 * t165) * qJD(2);
t314 = (-t329 * qJD(2) + t64 + t66) * t166 + ((-t72 - t74) * t176 + (t68 + t70) * t174 + (t336 * t174 + t337 * t176) * qJD(4) - t338 * qJD(2)) * t165;
t313 = (-t330 * qJD(2) + t65 + t67) * t166 + ((-t73 - t75) * t176 + (t69 + t71) * t174 + (t340 * t174 + t341 * t176) * qJD(4) + t332 * qJD(2)) * t165;
t312 = t323 * t170;
t311 = t323 * t171;
t37 = t102 * t154 + t153 * t98 + t258 * t94;
t39 = t100 * t153 + t104 * t154 + t258 * t96;
t308 = t37 + t39;
t307 = t316 * t165 - t322;
t301 = t319 * t170 + t318 * t171;
t286 = pkin(4) * t176;
t298 = qJ(5) * t166 - t165 * t286;
t295 = 2 * m(5);
t294 = 2 * m(6);
t291 = t170 / 0.2e1;
t288 = pkin(2) * t175;
t247 = qJD(4) * t174;
t242 = pkin(4) * t247;
t178 = t298 * qJD(2) + qJD(5) * t165 - t166 * t242;
t246 = qJD(4) * t176;
t241 = pkin(4) * t246;
t283 = rSges(6,1) * t111 + rSges(6,2) * t110 + rSges(6,3) * t238 + t170 * t178 - t171 * t241;
t282 = rSges(6,1) * t113 + rSges(6,2) * t112 + rSges(6,3) * t237 + t170 * t241 + t171 * t178;
t180 = qJ(5) * t165 + t166 * t286;
t223 = rSges(6,1) * t176 - rSges(6,2) * t174;
t281 = -qJD(5) * t166 - t165 * t242 + (-rSges(6,1) * t174 - rSges(6,2) * t176) * t248 + (rSges(6,3) * t165 + t166 * t223 + t180) * qJD(2);
t280 = pkin(2) * qJD(2);
t36 = t101 * t154 + t153 * t97 + t258 * t93;
t277 = t170 * t36;
t38 = t103 * t154 + t153 * t99 + t258 * t95;
t276 = t170 * t38;
t273 = rSges(6,1) * t152 + rSges(6,2) * t151 + rSges(6,3) * t259 - pkin(4) * t255 + t170 * t180;
t272 = rSges(6,1) * t154 + rSges(6,2) * t153 + rSges(6,3) * t258 + pkin(4) * t257 + t171 * t180;
t224 = rSges(5,1) * t176 - rSges(5,2) * t174;
t128 = -rSges(5,3) * t166 + t165 * t224;
t261 = t128 * t166;
t253 = -rSges(6,3) * t166 + t165 * t223 - t298;
t252 = t317 * pkin(2) * t177;
t251 = t317 * t175 * t280;
t245 = m(6) * t250;
t243 = t177 * t280;
t159 = rSges(4,1) * t165 + rSges(4,2) * t166;
t236 = -t159 - t288;
t160 = pkin(3) * t165 - pkin(6) * t166;
t235 = -t160 - t288;
t234 = t253 * t171;
t233 = t253 * t170;
t227 = pkin(3) * t166 + pkin(6) * t165;
t232 = t227 * t317 + t252;
t231 = -t160 * t299 - t251;
t230 = -t128 + t235;
t225 = rSges(4,1) * t166 - rSges(4,2) * t165;
t229 = -t225 * qJD(2) - t243;
t228 = -t227 * qJD(2) - t243;
t162 = t175 * rSges(3,1) + rSges(3,2) * t177;
t106 = rSges(5,1) * t152 + rSges(5,2) * t151 + rSges(5,3) * t259;
t108 = rSges(5,1) * t154 + rSges(5,2) * t153 + rSges(5,3) * t258;
t206 = t106 * t171 - t108 * t170;
t92 = (-rSges(5,1) * t174 - rSges(5,2) * t176) * t248 + (rSges(5,3) * t165 + t166 * t224) * qJD(2);
t197 = t228 - t92;
t196 = t235 - t253;
t191 = t228 - t281;
t59 = -t159 * t299 - t251;
t189 = t59 * t225;
t179 = -t170 * t272 + t171 * t273;
t116 = t229 * t171;
t115 = t229 * t170;
t109 = t162 * t299;
t82 = t230 * t171;
t81 = t230 * t170;
t79 = rSges(5,1) * t113 + rSges(5,2) * t112 + rSges(5,3) * t237;
t77 = rSges(5,1) * t111 + rSges(5,2) * t110 + rSges(5,3) * t238;
t61 = t197 * t171;
t60 = t197 * t170;
t58 = -t108 * t166 - t128 * t258;
t57 = t106 * t166 + t128 * t259;
t56 = t196 * t171;
t55 = t196 * t170;
t52 = t206 * t165;
t51 = t120 * t258 + t122 * t153 + t124 * t154;
t50 = t119 * t258 + t121 * t153 + t123 * t154;
t47 = t191 * t171;
t46 = t191 * t170;
t41 = -t165 * t234 - t166 * t272;
t40 = t165 * t233 + t166 * t273;
t31 = t106 * t170 + t108 * t171 + t232;
t30 = t179 * t165;
t29 = -t92 * t258 - t166 * t79 + (t108 * t165 - t171 * t261) * qJD(2);
t28 = t92 * t259 + t166 * t77 + (-t106 * t165 + t170 * t261) * qJD(2);
t27 = t170 * t77 + t171 * t79 + t231;
t26 = (-t170 * t79 + t171 * t77) * t165 + t206 * t249;
t25 = t170 * t273 + t171 * t272 + t232;
t24 = t170 * t283 + t171 * t282 + t231;
t23 = -t282 * t166 - t281 * t258 + (t165 * t272 - t166 * t234) * qJD(2);
t22 = t283 * t166 + t281 * t259 + (-t165 * t273 + t166 * t233) * qJD(2);
t21 = t100 * t112 + t104 * t113 + t153 * t71 + t154 * t75 + t171 * t192;
t20 = t103 * t113 + t112 * t99 + t153 * t70 + t154 * t74 + t171 * t193;
t19 = t102 * t113 + t112 * t98 + t153 * t69 + t154 * t73 + t171 * t194;
t18 = t101 * t113 + t112 * t97 + t153 * t68 + t154 * t72 + t171 * t195;
t9 = (-t170 * t282 + t171 * t283) * t165 + t179 * t249;
t8 = t170 * t21 - t171 * t20;
t7 = t170 * t19 - t171 * t18;
t6 = -t16 * t171 + t17 * t170;
t5 = -t14 * t171 + t15 * t170;
t4 = -(t112 * t122 + t113 * t124 + t153 * t88 + t154 * t90) * t166 + (t20 * t170 + (t21 - t278) * t171) * t165 + (t51 * t165 + (t276 + (t39 - t262) * t171) * t166) * qJD(2);
t3 = -(t112 * t121 + t113 * t123 + t153 * t87 + t154 * t89) * t166 + (t18 * t170 + (t19 - t279) * t171) * t165 + (t50 * t165 + (t277 + (t37 - t263) * t171) * t166) * qJD(2);
t1 = [0; -m(3) * t109 + m(4) * t59 + m(5) * t27 + m(6) * t24; (t24 * t25 + t46 * t55 + t47 * t56) * t294 - t171 * t6 + (t27 * t31 + t60 * t81 + t61 * t82) * t295 - t171 * t5 + 0.2e1 * m(4) * (t252 * t59 + (t116 * t236 + t171 * t189) * t171) + 0.2e1 * m(3) * (qJD(2) * t162 - t109) * t317 * (rSges(3,1) * t177 - rSges(3,2) * t175) - t312 * t171 * t168 + (t7 + t8 + t311 * t167 + 0.2e1 * (t115 * t236 + t170 * t189) * m(4) + (-t312 * t170 + t311 * t171) * t171) * t170; 0; m(6) * (t170 * t47 - t171 * t46) + m(5) * (t170 * t61 - t171 * t60) + m(4) * (-t115 * t171 + t116 * t170); 0; m(5) * t26 + m(6) * t9; t4 * t291 + t3 * t291 + m(6) * (t22 * t56 + t23 * t55 + t24 * t30 + t25 * t9 + t40 * t47 + t41 * t46) + m(5) * (t26 * t31 + t27 * t52 + t28 * t82 + t29 * t81 + t57 * t61 + t58 * t60) + ((t8 / 0.2e1 + t7 / 0.2e1) * t171 + (t5 / 0.2e1 + t6 / 0.2e1) * t170) * t165 + (((t339 * t170 - t309 * t171) * t291 + ((-t36 - t38) * t171 + t308 * t170) * t171 / 0.2e1) * t166 + (t318 * t170 - t319 * t171) * t165 / 0.2e1) * qJD(2) - (-t313 * t170 + t314 * t171) * t166 / 0.2e1 - t315 * t171 / 0.2e1; m(5) * (t170 * t28 - t171 * t29) + m(6) * (t170 * t22 - t171 * t23); (t22 * t40 + t23 * t41 + t30 * t9) * t294 + (t26 * t52 + t28 * t57 + t29 * t58) * t295 + t315 * t259 + (t3 + t4) * t258 + (t301 * t250 + (t309 * t170 + t321) * t238 + (t308 * t171 + t276 + t277) * t237) * t165 + (t320 * t166 + t307 * t250 - t328 * t238 + (-t50 - t51) * t237 + ((t327 * t174 + t326 * t176 - t246 * t331 - t247 * t334) * t166 + t313 * t171 + t314 * t170) * t165 + ((-t316 * t166 - t301) * t166 + (t335 * t166 + t307) * t165) * qJD(2)) * t166; t245; m(6) * (-t166 * t24 + (t170 * t46 + t171 * t47) * t165 + (t165 * t25 + (t170 * t55 + t171 * t56) * t166) * qJD(2)); 0; m(6) * (-t166 * t9 + (t170 * t23 + t171 * t22) * t165 + (t165 * t30 + (t170 * t41 + t171 * t40) * t166) * qJD(2)); 0.2e1 * (-0.1e1 + t317) * t166 * t245;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
