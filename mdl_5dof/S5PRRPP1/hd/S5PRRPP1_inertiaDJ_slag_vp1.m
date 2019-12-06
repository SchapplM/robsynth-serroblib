% Calculate time derivative of joint inertia matrix for
% S5PRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
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
% Datum: 2019-12-05 16:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPP1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP1_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP1_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP1_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP1_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPP1_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPP1_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:06:12
% EndTime: 2019-12-05 16:06:26
% DurationCPUTime: 5.27s
% Computational Cost: add. (4702->299), mult. (4486->410), div. (0->0), fcn. (3401->6), ass. (0->187)
t320 = Icges(6,4) + Icges(5,5);
t319 = -Icges(5,6) + Icges(6,6);
t148 = sin(qJ(3));
t149 = cos(qJ(3));
t176 = Icges(4,5) * t149 - Icges(4,6) * t148;
t146 = qJ(3) + pkin(8);
t142 = sin(t146);
t144 = cos(t146);
t316 = t319 * t142 + t320 * t144;
t318 = t176 + t316;
t317 = Icges(6,2) + Icges(4,3) + Icges(5,3);
t247 = Icges(6,5) * t144;
t252 = Icges(5,4) * t144;
t315 = -t247 + t252 + (Icges(5,1) + Icges(6,1)) * t142;
t256 = rSges(6,3) + qJ(5);
t274 = rSges(6,1) + pkin(4);
t300 = t142 * t256 + t144 * t274;
t287 = t142 * t274;
t234 = -t144 * t256 + t287;
t145 = pkin(7) + qJ(2);
t141 = sin(t145);
t143 = cos(t145);
t307 = t317 * t141 + t318 * t143;
t248 = Icges(6,5) * t142;
t183 = Icges(6,1) * t144 + t248;
t68 = -Icges(6,4) * t143 + t141 * t183;
t253 = Icges(5,4) * t142;
t185 = Icges(5,1) * t144 - t253;
t70 = -Icges(5,5) * t143 + t141 * t185;
t314 = t70 + t68;
t69 = Icges(6,4) * t141 + t143 * t183;
t71 = Icges(5,5) * t141 + t143 * t185;
t313 = t71 + t69;
t224 = qJD(3) * t144;
t311 = t256 * t224;
t254 = Icges(4,4) * t149;
t181 = -Icges(4,2) * t148 + t254;
t81 = Icges(4,6) * t141 + t143 * t181;
t255 = Icges(4,4) * t148;
t187 = Icges(4,1) * t149 - t255;
t83 = Icges(4,5) * t141 + t143 * t187;
t188 = t148 * t81 - t149 * t83;
t179 = -Icges(5,2) * t142 + t252;
t67 = Icges(5,6) * t141 + t143 * t179;
t192 = t142 * t67 - t144 * t71;
t174 = Icges(6,3) * t142 + t247;
t61 = Icges(6,6) * t141 + t143 * t174;
t194 = t142 * t61 + t144 * t69;
t279 = t188 + t192 - t194;
t310 = t279 * t143;
t309 = t315 * qJD(3);
t308 = t318 * t141 - t317 * t143;
t306 = (-Icges(4,5) * t148 - Icges(4,6) * t149 - t320 * t142 + t319 * t144) * qJD(3);
t80 = -Icges(4,6) * t143 + t141 * t181;
t82 = -Icges(4,5) * t143 + t141 * t187;
t189 = t148 * t80 - t149 * t82;
t66 = -Icges(5,6) * t143 + t141 * t179;
t193 = t142 * t66 - t144 * t70;
t60 = -Icges(6,6) * t143 + t141 * t174;
t195 = t142 * t60 + t144 * t68;
t305 = t189 + t193 - t195;
t304 = t141 / 0.2e1;
t303 = -t143 / 0.2e1;
t302 = -qJD(2) / 0.2e1;
t301 = qJD(2) / 0.2e1;
t164 = t188 * t141;
t165 = t189 * t143;
t166 = t192 * t141;
t167 = t193 * t143;
t168 = t194 * t141;
t169 = t195 * t143;
t299 = t308 * t141 - t307 * t143 - t164 - t165 - t166 - t167 + t168 + t169;
t298 = t307 * qJD(2);
t297 = t313 * t141 + t314 * t143;
t296 = (t67 - t61) * t141 + (t66 - t60) * t143;
t139 = t141 ^ 2;
t140 = t143 ^ 2;
t105 = rSges(5,1) * t142 + rSges(5,2) * t144;
t162 = qJD(3) * t105;
t132 = qJD(4) * t141;
t223 = qJD(3) * t148;
t216 = pkin(3) * t223;
t147 = -qJ(4) - pkin(6);
t228 = qJD(2) * t147;
t213 = qJD(4) * t143 + (t216 + t228) * t141;
t229 = qJD(2) * t143;
t269 = -pkin(6) - t147;
t138 = pkin(3) * t149 + pkin(2);
t270 = pkin(2) - t138;
t271 = pkin(6) * t141;
t288 = t141 * t270;
t137 = t143 * pkin(6);
t235 = t143 * t147;
t58 = t137 + t235 - t288;
t220 = t141 * ((-t143 * t270 - t271) * qJD(2) - t213) + t143 * (-t143 * t216 + t132 + (t143 * t269 + t288) * qJD(2)) + t58 * t229;
t230 = qJD(2) * t141;
t262 = rSges(5,2) * t142;
t233 = rSges(5,3) * t229 + t230 * t262;
t134 = t141 * rSges(5,3);
t237 = t142 * t143;
t286 = -rSges(5,2) * t237 + t134;
t119 = t143 * t138;
t59 = -pkin(2) * t143 + t141 * t269 + t119;
t264 = rSges(5,1) * t144;
t204 = -t262 + t264;
t73 = -rSges(5,3) * t143 + t141 * t204;
t236 = t143 * t144;
t75 = rSges(5,1) * t236 + t286;
t2 = (qJD(2) * t73 - t143 * t162 + t233) * t143 + (-t141 * t162 + (-t59 - t75 + t286) * qJD(2)) * t141 + t220;
t276 = 2 * m(5);
t295 = t2 * t276;
t277 = 2 * m(4);
t135 = t141 * rSges(4,3);
t127 = t148 * rSges(4,1) + rSges(4,2) * t149;
t163 = qJD(3) * t127;
t227 = qJD(2) * t148;
t212 = t141 * t227;
t151 = rSges(4,2) * t212 + rSges(4,3) * t229 - t143 * t163;
t263 = rSges(4,2) * t148;
t218 = t143 * t263;
t265 = rSges(4,1) * t149;
t205 = -t263 + t265;
t261 = t143 * rSges(4,3);
t84 = t141 * t205 - t261;
t232 = t143 * t265 + t135;
t85 = -t218 + t232;
t4 = (qJD(2) * t84 + t151) * t143 + (-t141 * t163 + (-t218 - t85 + t135) * qJD(2)) * t141;
t294 = t277 * t4;
t285 = t305 * t141 + t308 * t143;
t282 = t307 * t141 - t310;
t281 = -t308 * qJD(2) + t306 * t143;
t280 = -t306 * t141 - t298;
t275 = 2 * m(6);
t273 = m(4) * t127;
t272 = pkin(3) * t148;
t268 = t141 * t58 + t143 * t59;
t267 = -rSges(6,2) * t143 + t141 * t300;
t136 = t141 * rSges(6,2);
t266 = t236 * t274 + t237 * t256 + t136;
t231 = t139 + t140;
t226 = qJD(3) * t141;
t225 = qJD(3) * t142;
t222 = qJD(3) * t149;
t221 = qJD(5) * t142;
t219 = m(6) * t225;
t215 = pkin(3) * t222;
t214 = t316 * qJD(3) / 0.2e1;
t210 = t274 * qJD(3);
t209 = -t105 - t272;
t208 = -t141 * t147 + t119;
t207 = rSges(6,2) * t229 + (t221 + t311) * t143;
t206 = -t234 - t272;
t152 = -t138 - t300;
t150 = t152 * t141;
t23 = (rSges(6,2) - t147) * t143 + t150;
t24 = t208 + t266;
t201 = t141 * t24 + t143 * t23;
t51 = t206 * t141;
t52 = t206 * t143;
t200 = t141 * t51 + t143 * t52;
t178 = Icges(5,2) * t144 + t253;
t172 = -qJD(3) * t300 + qJD(5) * t144 - t215;
t171 = -pkin(2) - t205;
t170 = -t138 - t204;
t157 = qJD(3) * t178;
t122 = pkin(3) * t212;
t112 = t205 * qJD(3);
t96 = t204 * qJD(3);
t77 = t209 * t143;
t76 = t209 * t141;
t57 = t271 + (pkin(2) - t263) * t143 + t232;
t56 = t141 * t171 + t137 + t261;
t50 = t208 + t75;
t49 = (rSges(5,3) - t147) * t143 + t170 * t141;
t28 = -t105 * t229 - t141 * t96 + (-t141 * t222 - t143 * t227) * pkin(3);
t27 = t105 * t230 + t122 + (-t96 - t215) * t143;
t26 = t127 * t226 + ((-rSges(4,3) - pkin(6)) * t141 + t171 * t143) * qJD(2);
t25 = (t137 + (-pkin(2) - t265) * t141) * qJD(2) + t151;
t22 = t105 * t226 + (t143 * t170 - t134) * qJD(2) + t213;
t21 = t132 + (-t138 - t264) * t230 + (qJD(3) * t209 - t228) * t143 + t233;
t16 = qJD(2) * t52 + t141 * t172;
t15 = t143 * t172 + t230 * t234 + t122;
t6 = (t234 * qJD(3) - t221) * t141 + (t143 * t152 - t136) * qJD(2) + t213;
t5 = t132 + (-t272 - t287) * t143 * qJD(3) + (t150 - t235) * qJD(2) + t207;
t3 = t141 * t267 + t143 * t266 + t268;
t1 = (qJD(2) * t267 - t210 * t237 + t207) * t143 + t220 + ((-t59 + t136 - t266) * qJD(2) + (t311 + (qJD(5) - t210) * t142) * t141) * t141;
t7 = [0; 0; (t25 * t57 + t26 * t56) * t277 + (t21 * t50 + t22 * t49) * t276 + (t23 * t6 + t24 * t5) * t275 + (-Icges(4,2) * t149 + t187 - t255) * t223 + (Icges(4,1) * t148 + t181 + t254) * t222 + (-Icges(6,3) * t144 - t178 + t183 + t185 + t248) * t225 + (t179 - t174 + t315) * t224; m(4) * t4 + m(5) * t2 + m(6) * t1; m(5) * (t21 * t76 + t22 * t77 + t27 * t49 + t28 * t50) + m(6) * (t15 * t23 + t16 * t24 + t5 * t51 + t52 * t6) + (m(4) * (-t112 * t56 - t127 * t26) + (t157 * t304 + t61 * t301 + t67 * t302) * t144 + (t313 * t302 + t309 * t304) * t142 + t214 * t143) * t143 + (m(4) * (-t112 * t57 - t127 * t25) + (t157 * t303 + t60 * t301 + t66 * t302) * t144 + (t314 * t302 + t309 * t303) * t142 + t214 * t141) * t141 + ((t139 / 0.2e1 + t140 / 0.2e1) * t176 + t165 / 0.2e1 - t164 / 0.2e1 - t169 / 0.2e1 + t167 / 0.2e1 + t168 / 0.2e1 - t166 / 0.2e1) * qJD(3) + ((-t57 * t273 + (t67 / 0.2e1 - t61 / 0.2e1) * t144 + (t71 / 0.2e1 + t69 / 0.2e1) * t142) * t143 + (t56 * t273 + (-t60 / 0.2e1 + t66 / 0.2e1) * t144 + (t68 / 0.2e1 + t70 / 0.2e1) * t142) * t141) * qJD(2); (t1 * t3 + t15 * t52 + t16 * t51) * t275 + (t268 * t2 + t77 * t27 + t76 * t28) * t276 + t231 * t127 * t112 * t277 + (t280 * t140 + t285 * t230 + t85 * t294 + t75 * t295 + (-t143 * t305 - t299) * t229) * t143 + (t73 * t295 + t84 * t294 + t281 * t139 + t282 * t229 + ((t222 * t81 + t223 * t83 + t280 - t298) * t141 + (t222 * t80 + t223 * t82 + t281) * t143 + t297 * t225 + t296 * t224 + ((-t148 * t83 - t149 * t81) * t141 + (-t148 * t82 - t149 * t80) * t143 - t296 * t144 - t297 * t142) * qJD(3) + ((-t305 + t307) * t141 + t310 + t282 + t285) * qJD(2)) * t143 + (t279 * t141 + t299) * t230) * t141; 0; m(5) * (t141 * t22 - t143 * t21 + (t141 * t50 + t143 * t49) * qJD(2)) + m(6) * (qJD(2) * t201 + t141 * t6 - t143 * t5); m(6) * (qJD(2) * t200 + t141 * t15 - t143 * t16) + m(5) * (t141 * t27 - t143 * t28 + (t141 * t76 + t143 * t77) * qJD(2)); 0; t219; m(6) * (t201 * t224 + (t141 * t5 + t143 * t6 + (-t141 * t23 + t143 * t24) * qJD(2)) * t142); m(6) * ((qJD(3) * t200 - t1) * t144 + (qJD(3) * t3 + t141 * t16 + t143 * t15 + (-t141 * t52 + t143 * t51) * qJD(2)) * t142); 0; 0.2e1 * (-0.1e1 + t231) * t144 * t219;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t7(1), t7(2), t7(4), t7(7), t7(11); t7(2), t7(3), t7(5), t7(8), t7(12); t7(4), t7(5), t7(6), t7(9), t7(13); t7(7), t7(8), t7(9), t7(10), t7(14); t7(11), t7(12), t7(13), t7(14), t7(15);];
Mq = res;
