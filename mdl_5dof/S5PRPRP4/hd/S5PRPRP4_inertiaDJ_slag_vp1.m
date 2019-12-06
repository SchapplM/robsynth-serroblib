% Calculate time derivative of joint inertia matrix for
% S5PRPRP4
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
% Datum: 2019-12-05 15:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRP4_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP4_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP4_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP4_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP4_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRP4_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRP4_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:35:04
% EndTime: 2019-12-05 15:35:23
% DurationCPUTime: 7.48s
% Computational Cost: add. (8429->409), mult. (13160->629), div. (0->0), fcn. (12546->8), ass. (0->220)
t170 = qJ(2) + pkin(8);
t169 = cos(t170);
t172 = cos(pkin(7));
t174 = sin(qJ(4));
t257 = t172 * t174;
t171 = sin(pkin(7));
t176 = cos(qJ(4));
t258 = t171 * t176;
t153 = t169 * t257 - t258;
t256 = t172 * t176;
t259 = t171 * t174;
t154 = t169 * t256 + t259;
t168 = sin(t170);
t260 = t168 * t172;
t92 = Icges(6,5) * t154 + Icges(6,6) * t260 + Icges(6,3) * t153;
t98 = Icges(5,4) * t154 - Icges(5,2) * t153 + Icges(5,6) * t260;
t341 = t92 - t98;
t94 = Icges(5,5) * t154 - Icges(5,6) * t153 + Icges(5,3) * t260;
t96 = Icges(6,4) * t154 + Icges(6,2) * t260 + Icges(6,6) * t153;
t332 = t96 + t94;
t100 = Icges(6,1) * t154 + Icges(6,4) * t260 + Icges(6,5) * t153;
t102 = Icges(5,1) * t154 - Icges(5,4) * t153 + Icges(5,5) * t260;
t340 = t100 + t102;
t151 = t169 * t259 + t256;
t152 = t169 * t258 - t257;
t261 = t168 * t171;
t91 = Icges(6,5) * t152 + Icges(6,6) * t261 + Icges(6,3) * t151;
t97 = Icges(5,4) * t152 - Icges(5,2) * t151 + Icges(5,6) * t261;
t338 = t91 - t97;
t101 = Icges(5,1) * t152 - Icges(5,4) * t151 + Icges(5,5) * t261;
t99 = Icges(6,1) * t152 + Icges(6,4) * t261 + Icges(6,5) * t151;
t336 = t101 + t99;
t339 = t341 * t151 + t340 * t152 + t332 * t261;
t93 = Icges(5,5) * t152 - Icges(5,6) * t151 + Icges(5,3) * t261;
t95 = Icges(6,4) * t152 + Icges(6,2) * t261 + Icges(6,6) * t151;
t337 = t95 + t93;
t266 = Icges(6,5) * t176;
t206 = Icges(6,3) * t174 + t266;
t118 = -Icges(6,6) * t169 + t168 * t206;
t268 = Icges(5,4) * t176;
t209 = -Icges(5,2) * t174 + t268;
t121 = -Icges(5,6) * t169 + t168 * t209;
t326 = t118 - t121;
t207 = Icges(5,5) * t176 - Icges(5,6) * t174;
t119 = -Icges(5,3) * t169 + t168 * t207;
t208 = Icges(6,4) * t176 + Icges(6,6) * t174;
t120 = -Icges(6,2) * t169 + t168 * t208;
t335 = t119 + t120;
t267 = Icges(6,5) * t174;
t212 = Icges(6,1) * t176 + t267;
t122 = -Icges(6,4) * t169 + t168 * t212;
t269 = Icges(5,4) * t174;
t213 = Icges(5,1) * t176 - t269;
t123 = -Icges(5,5) * t169 + t168 * t213;
t334 = t122 + t123;
t331 = t341 * t174 + t340 * t176;
t330 = t338 * t174 + t336 * t176;
t308 = t151 * t338 + t152 * t336 + t261 * t337;
t329 = t151 * t326 + t152 * t334 + t261 * t335;
t246 = qJD(4) * t168;
t83 = (Icges(6,3) * t176 - t267) * t246 + (Icges(6,6) * t168 + t169 * t206) * qJD(2);
t86 = (-Icges(5,2) * t176 - t269) * t246 + (Icges(5,6) * t168 + t169 * t209) * qJD(2);
t328 = t83 - t86;
t87 = (-Icges(6,1) * t174 + t266) * t246 + (Icges(6,4) * t168 + t169 * t212) * qJD(2);
t88 = (-Icges(5,1) * t174 - t268) * t246 + (Icges(5,5) * t168 + t169 * t213) * qJD(2);
t327 = t87 + t88;
t175 = sin(qJ(2));
t177 = cos(qJ(2));
t324 = (-Icges(3,5) * t175 - Icges(4,5) * t168 - Icges(3,6) * t177 - Icges(4,6) * t169) * qJD(2);
t264 = t120 * t169;
t265 = t119 * t169;
t323 = -t264 - t265;
t322 = t339 * t172;
t278 = t169 * ((-Icges(6,4) * t174 + Icges(6,6) * t176) * t246 + (Icges(6,2) * t168 + t169 * t208) * qJD(2));
t279 = t169 * ((-Icges(5,5) * t174 - Icges(5,6) * t176) * t246 + (Icges(5,3) * t168 + t169 * t207) * qJD(2));
t321 = -t278 - t279;
t320 = t330 * t168 - t169 * t337;
t319 = t331 * t168 - t332 * t169;
t318 = -t326 * t174 - t176 * t334;
t292 = t172 ^ 2;
t293 = t171 ^ 2;
t317 = t292 + t293;
t298 = qJD(2) * t317;
t244 = qJD(4) * t176;
t237 = t169 * t244;
t248 = qJD(2) * t171;
t110 = -t171 * t237 + (qJD(4) * t172 + t168 * t248) * t174;
t250 = qJD(2) * t168;
t241 = t176 * t250;
t111 = -qJD(4) * t151 - t171 * t241;
t249 = qJD(2) * t169;
t240 = t169 * t248;
t66 = Icges(6,4) * t111 + Icges(6,2) * t240 - Icges(6,6) * t110;
t192 = t168 * t66 + t249 * t95;
t62 = Icges(6,5) * t111 + Icges(6,6) * t240 - Icges(6,3) * t110;
t70 = Icges(6,1) * t111 + Icges(6,4) * t240 - Icges(6,5) * t110;
t14 = -t110 * t91 + t111 * t99 + t151 * t62 + t152 * t70 + t171 * t192;
t245 = qJD(4) * t174;
t247 = qJD(2) * t174;
t112 = -t171 * t245 - t172 * t237 + t247 * t260;
t113 = -qJD(4) * t153 - t172 * t241;
t239 = t172 * t249;
t67 = Icges(6,4) * t113 + Icges(6,2) * t239 - Icges(6,6) * t112;
t191 = t168 * t67 + t249 * t96;
t63 = Icges(6,5) * t113 + Icges(6,6) * t239 - Icges(6,3) * t112;
t71 = Icges(6,1) * t113 + Icges(6,4) * t239 - Icges(6,5) * t112;
t15 = t100 * t111 - t110 * t92 + t151 * t63 + t152 * t71 + t171 * t191;
t64 = Icges(5,5) * t111 + Icges(5,6) * t110 + Icges(5,3) * t240;
t194 = t168 * t64 + t249 * t93;
t68 = Icges(5,4) * t111 + Icges(5,2) * t110 + Icges(5,6) * t240;
t72 = Icges(5,1) * t111 + Icges(5,4) * t110 + Icges(5,5) * t240;
t16 = t101 * t111 + t110 * t97 - t151 * t68 + t152 * t72 + t171 * t194;
t65 = Icges(5,5) * t113 + Icges(5,6) * t112 + Icges(5,3) * t239;
t193 = t168 * t65 + t249 * t94;
t69 = Icges(5,4) * t113 + Icges(5,2) * t112 + Icges(5,6) * t239;
t73 = Icges(5,1) * t113 + Icges(5,4) * t112 + Icges(5,5) * t239;
t17 = t102 * t111 + t110 * t98 - t151 * t69 + t152 * t73 + t171 * t193;
t316 = (t326 * t110 - t111 * t334 - t328 * t151 - t327 * t152) * t169 + ((t15 + t17) * t172 + (t14 + t16 + t321) * t171) * t168 + (((t323 + t308) * t171 + t322) * t169 + t329 * t168) * qJD(2);
t315 = rSges(6,1) + pkin(4);
t314 = (-t330 * qJD(2) + t64 + t66) * t169 + ((-t70 - t72) * t176 + (-t62 + t68) * t174 + (t174 * t336 - t176 * t338) * qJD(4) - t337 * qJD(2)) * t168;
t313 = (t331 * qJD(2) - t65 - t67) * t169 + ((t71 + t73) * t176 + (t63 - t69) * t174 + (-t340 * t174 + t341 * t176) * qJD(4) + t332 * qJD(2)) * t168;
t312 = rSges(6,3) + qJ(5);
t311 = t324 * t171;
t310 = t324 * t172;
t37 = t100 * t154 + t153 * t92 + t260 * t96;
t39 = t102 * t154 - t153 * t98 + t260 * t94;
t307 = t37 + t39;
t306 = t318 * t168 - t323;
t301 = t320 * t171 + t319 * t172;
t295 = 2 * m(5);
t294 = 2 * m(6);
t289 = t171 / 0.2e1;
t286 = pkin(2) * t175;
t283 = rSges(6,2) * t240 + qJD(5) * t151 - t312 * t110 + t315 * t111;
t282 = rSges(6,2) * t239 + qJD(5) * t153 - t312 * t112 + t315 * t113;
t222 = pkin(4) * t176 + qJ(5) * t174;
t223 = rSges(6,1) * t176 + rSges(6,3) * t174;
t281 = t222 * t249 + (qJD(5) * t174 + (-pkin(4) * t174 + qJ(5) * t176) * qJD(4)) * t168 + (-rSges(6,1) * t174 + rSges(6,3) * t176) * t246 + (rSges(6,2) * t168 + t169 * t223) * qJD(2);
t280 = pkin(2) * qJD(2);
t36 = t153 * t91 + t154 * t99 + t260 * t95;
t277 = t171 * t36;
t38 = t101 * t154 - t153 * t97 + t260 * t93;
t276 = t171 * t38;
t224 = rSges(5,1) * t176 - rSges(5,2) * t174;
t127 = -rSges(5,3) * t169 + t168 * t224;
t263 = t127 * t169;
t255 = rSges(6,2) * t261 + t312 * t151 + t315 * t152;
t254 = rSges(6,2) * t260 + t312 * t153 + t315 * t154;
t253 = t317 * pkin(2) * t177;
t252 = -rSges(6,2) * t169 + (t222 + t223) * t168;
t251 = t317 * t175 * t280;
t242 = t177 * t280;
t238 = t169 * t247;
t159 = rSges(4,1) * t168 + rSges(4,2) * t169;
t236 = -t159 - t286;
t160 = pkin(3) * t168 - pkin(6) * t169;
t235 = -t160 - t286;
t234 = t252 * t172;
t233 = t252 * t171;
t227 = pkin(3) * t169 + pkin(6) * t168;
t232 = t227 * t317 + t253;
t231 = -t160 * t298 - t251;
t230 = -t127 + t235;
t225 = rSges(4,1) * t169 - rSges(4,2) * t168;
t229 = -t225 * qJD(2) - t242;
t228 = -t227 * qJD(2) - t242;
t164 = t175 * rSges(3,1) + rSges(3,2) * t177;
t104 = rSges(5,1) * t152 - rSges(5,2) * t151 + rSges(5,3) * t261;
t106 = rSges(5,1) * t154 - rSges(5,2) * t153 + rSges(5,3) * t260;
t205 = t104 * t172 - t106 * t171;
t90 = (-rSges(5,1) * t174 - rSges(5,2) * t176) * t246 + (rSges(5,3) * t168 + t169 * t224) * qJD(2);
t196 = t228 - t90;
t195 = t235 - t252;
t190 = t228 - t281;
t59 = -t159 * t298 - t251;
t188 = t59 * t225;
t179 = t168 * t244 + t238;
t178 = -t171 * t254 + t172 * t255;
t115 = t229 * t172;
t114 = t229 * t171;
t107 = t164 * t298;
t81 = t230 * t172;
t80 = t230 * t171;
t79 = t195 * t172;
t78 = t195 * t171;
t77 = rSges(5,1) * t113 + rSges(5,2) * t112 + rSges(5,3) * t239;
t75 = rSges(5,1) * t111 + rSges(5,2) * t110 + rSges(5,3) * t240;
t61 = t196 * t172;
t60 = t196 * t171;
t58 = -t106 * t169 - t127 * t260;
t57 = t104 * t169 + t127 * t261;
t52 = t205 * t168;
t51 = t119 * t260 - t121 * t153 + t123 * t154;
t50 = t118 * t153 + t120 * t260 + t122 * t154;
t47 = t190 * t172;
t46 = t190 * t171;
t45 = -t168 * t234 - t169 * t254;
t44 = t168 * t233 + t169 * t255;
t31 = t104 * t171 + t106 * t172 + t232;
t30 = t178 * t168;
t29 = -t90 * t260 - t169 * t77 + (t106 * t168 - t172 * t263) * qJD(2);
t28 = t90 * t261 + t169 * t75 + (-t104 * t168 + t171 * t263) * qJD(2);
t27 = t171 * t75 + t172 * t77 + t231;
t26 = t171 * t255 + t172 * t254 + t232;
t25 = (-t171 * t77 + t172 * t75) * t168 + t205 * t249;
t24 = t171 * t283 + t172 * t282 + t231;
t23 = -t282 * t169 - t281 * t260 + (t168 * t254 - t169 * t234) * qJD(2);
t22 = t283 * t169 + t281 * t261 + (-t168 * t255 + t169 * t233) * qJD(2);
t21 = t102 * t113 + t112 * t98 - t153 * t69 + t154 * t73 + t172 * t193;
t20 = t101 * t113 + t112 * t97 - t153 * t68 + t154 * t72 + t172 * t194;
t19 = t100 * t113 - t112 * t92 + t153 * t63 + t154 * t71 + t172 * t191;
t18 = -t112 * t91 + t113 * t99 + t153 * t62 + t154 * t70 + t172 * t192;
t9 = (-t171 * t282 + t172 * t283) * t168 + t178 * t249;
t8 = t171 * t21 - t172 * t20;
t7 = t171 * t19 - t172 * t18;
t6 = -t16 * t172 + t17 * t171;
t5 = -t14 * t172 + t15 * t171;
t4 = -(t112 * t121 + t113 * t123 - t153 * t86 + t154 * t88) * t169 + (t20 * t171 + (t21 - t279) * t172) * t168 + (t51 * t168 + (t276 + (t39 - t265) * t172) * t169) * qJD(2);
t3 = -(-t112 * t118 + t113 * t122 + t153 * t83 + t154 * t87) * t169 + (t18 * t171 + (t19 - t278) * t172) * t168 + (t50 * t168 + (t277 + (t37 - t264) * t172) * t169) * qJD(2);
t1 = [0; -m(3) * t107 + m(4) * t59 + m(5) * t27 + m(6) * t24; (t24 * t26 + t46 * t78 + t47 * t79) * t294 + (t27 * t31 + t60 * t80 + t61 * t81) * t295 - t172 * t6 - t172 * t5 + 0.2e1 * m(4) * (t253 * t59 + (t115 * t236 + t172 * t188) * t172) + 0.2e1 * m(3) * (qJD(2) * t164 - t107) * t317 * (rSges(3,1) * t177 - rSges(3,2) * t175) - t311 * t172 * t292 + (t7 + t8 + t310 * t293 + 0.2e1 * (t114 * t236 + t171 * t188) * m(4) + (-t311 * t171 + t310 * t172) * t172) * t171; 0; m(6) * (t171 * t47 - t172 * t46) + m(5) * (t171 * t61 - t172 * t60) + m(4) * (-t114 * t172 + t115 * t171); 0; m(5) * t25 + m(6) * t9; t3 * t289 + t4 * t289 + m(6) * (t22 * t79 + t23 * t78 + t24 * t30 + t26 * t9 + t44 * t47 + t45 * t46) + m(5) * (t25 * t31 + t27 * t52 + t28 * t81 + t29 * t80 + t57 * t61 + t58 * t60) + ((t8 / 0.2e1 + t7 / 0.2e1) * t172 + (t6 / 0.2e1 + t5 / 0.2e1) * t171) * t168 + (((t171 * t339 - t308 * t172) * t289 + ((-t36 - t38) * t172 + t307 * t171) * t172 / 0.2e1) * t169 + (t319 * t171 - t320 * t172) * t168 / 0.2e1) * qJD(2) - (t313 * t171 + t314 * t172) * t169 / 0.2e1 - t316 * t172 / 0.2e1; m(5) * (t171 * t28 - t172 * t29) + m(6) * (t171 * t22 - t172 * t23); (t22 * t44 + t23 * t45 + t30 * t9) * t294 + (t25 * t52 + t28 * t57 + t29 * t58) * t295 + t316 * t261 + (t3 + t4) * t260 + (t301 * t250 + (t308 * t171 + t322) * t240 + (t307 * t172 + t276 + t277) * t239) * t168 + (t321 * t169 + t306 * t250 - t329 * t240 + (-t50 - t51) * t239 + ((t328 * t174 + t327 * t176 + t326 * t244 - t245 * t334) * t169 - t313 * t172 + t314 * t171) * t168 + ((-t318 * t169 - t301) * t169 + (t169 * t335 + t306) * t168) * qJD(2)) * t169; t179 * m(6); m(6) * (t26 * t238 - t110 * t78 - t112 * t79 + t151 * t46 + t153 * t47 + (t174 * t24 + t244 * t26) * t168); m(6) * (t110 * t172 - t112 * t171); m(6) * (t30 * t238 - t110 * t45 - t112 * t44 + t151 * t23 + t153 * t22 + (t174 * t9 + t30 * t244) * t168); (t168 * t174 * t179 - t110 * t151 - t112 * t153) * t294;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
