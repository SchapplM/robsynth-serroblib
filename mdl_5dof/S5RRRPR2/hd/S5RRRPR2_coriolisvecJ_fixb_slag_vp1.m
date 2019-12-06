% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRPR2_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR2_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR2_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR2_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR2_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:40:44
% EndTime: 2019-12-05 18:40:50
% DurationCPUTime: 3.64s
% Computational Cost: add. (10345->327), mult. (5823->394), div. (0->0), fcn. (4286->10), ass. (0->216)
t157 = sin(qJ(5));
t159 = cos(qJ(5));
t133 = rSges(6,1) * t157 + rSges(6,2) * t159;
t318 = qJD(5) * t133;
t156 = qJ(1) + qJ(2);
t152 = qJ(3) + t156;
t143 = pkin(9) + t152;
t140 = sin(t143);
t155 = qJD(1) + qJD(2);
t148 = qJD(3) + t155;
t259 = t140 * t148;
t101 = rSges(5,2) * t259;
t141 = cos(t143);
t279 = rSges(5,1) * t141;
t145 = cos(t152);
t290 = pkin(3) * t145;
t217 = -t279 - t290;
t137 = t140 * rSges(5,2);
t94 = -t137 + t279;
t320 = -t101 + (-t217 - t94) * t148;
t151 = Icges(6,4) * t159;
t207 = -Icges(6,2) * t157 + t151;
t313 = Icges(6,1) * t157 + t151;
t246 = t313 + t207;
t267 = Icges(6,4) * t157;
t127 = Icges(6,2) * t159 + t267;
t130 = Icges(6,1) * t159 - t267;
t247 = t127 - t130;
t319 = (t157 * t246 + t159 * t247) * t148;
t317 = 0.2e1 * qJD(5);
t244 = qJD(5) * t140;
t95 = t133 * t244;
t287 = pkin(8) * t140;
t288 = pkin(4) * t141;
t96 = t287 + t288;
t316 = -t148 * t96 + t95;
t150 = cos(t156);
t250 = t150 * t155;
t241 = pkin(2) * t250;
t315 = t241 + t320;
t160 = cos(qJ(1));
t274 = pkin(1) * qJD(1);
t234 = t160 * t274;
t190 = t234 + t241;
t144 = sin(t152);
t100 = rSges(4,1) * t145 - t144 * rSges(4,2);
t91 = t148 * t100;
t64 = -t190 - t91;
t158 = sin(qJ(1));
t235 = t158 * t274;
t149 = sin(t156);
t214 = rSges(3,1) * t149 + rSges(3,2) * t150;
t92 = t214 * t155;
t88 = t92 + t235;
t254 = t141 * t159;
t255 = t141 * t157;
t68 = Icges(6,4) * t254 - Icges(6,2) * t255 + Icges(6,6) * t140;
t122 = Icges(6,4) * t255;
t70 = Icges(6,1) * t254 + Icges(6,5) * t140 - t122;
t209 = t157 * t68 - t159 * t70;
t314 = t141 * t209;
t258 = t140 * t157;
t123 = rSges(6,2) * t258;
t248 = t141 * rSges(6,3) + t123;
t253 = t144 * t148;
t240 = pkin(3) * t253;
t243 = qJD(5) * t141;
t286 = pkin(8) * t141;
t257 = t140 * t159;
t237 = rSges(6,1) * t257;
t204 = t237 - t248;
t61 = t148 * t204;
t177 = t133 * t243 + t240 + t61 - t148 * (-pkin(4) * t140 + t286);
t251 = t149 * t155;
t242 = pkin(2) * t251;
t166 = t177 + t242;
t28 = t166 + t235;
t124 = rSges(6,2) * t255;
t236 = rSges(6,1) * t254;
t275 = rSges(6,3) * t140;
t201 = -t236 - t275;
t71 = -t124 - t201;
t29 = t95 + (-t71 - t96 - t290) * t148 - t190;
t312 = t140 * t29 + t141 * t28;
t256 = t141 * t148;
t249 = -rSges(5,1) * t259 - rSges(5,2) * t256;
t52 = (-t94 - t290) * t148 - t190;
t212 = -rSges(5,1) * t140 - rSges(5,2) * t141;
t82 = t148 * t212;
t311 = (-t249 + t82) * t52;
t233 = t148 * t124 + t318 * t140;
t62 = t148 * t71;
t310 = -t233 - t62 + t316;
t125 = Icges(6,5) * t157 + Icges(6,6) * t159;
t307 = -Icges(6,3) * t148 + qJD(5) * t125;
t306 = -Icges(6,6) * t148 + qJD(5) * t127;
t111 = t207 * qJD(5);
t112 = t130 * qJD(5);
t305 = qJD(5) * (t127 * t159 + t157 * t313) + t111 * t157 - t112 * t159 - t125 * t148;
t304 = -Icges(6,5) * t148 + qJD(5) * t313;
t280 = -Icges(6,2) * t254 - t122 + t70;
t282 = t141 * t313 + t68;
t302 = t157 * t280 + t159 * t282;
t121 = Icges(6,4) * t258;
t69 = -Icges(6,1) * t257 + Icges(6,5) * t141 + t121;
t281 = Icges(6,2) * t257 + t121 + t69;
t67 = Icges(6,6) * t141 - t140 * t207;
t283 = -t140 * t313 + t67;
t301 = -t157 * t281 - t159 * t283;
t298 = -rSges(6,3) - pkin(8);
t297 = pkin(1) * t158;
t296 = pkin(1) * t160;
t295 = pkin(1) * qJD(1) ^ 2;
t294 = pkin(2) * t149;
t293 = pkin(2) * t150;
t292 = pkin(2) * t155 ^ 2;
t291 = pkin(3) * t144;
t289 = pkin(3) * t148 ^ 2;
t126 = Icges(6,5) * t159 - Icges(6,6) * t157;
t65 = Icges(6,3) * t141 - t126 * t140;
t285 = t141 * t65 + t67 * t258;
t284 = t140 * t65 + t69 * t254;
t277 = rSges(6,1) * t159;
t271 = t157 * t67;
t270 = t159 * t69;
t205 = t157 * t127 - t159 * t313;
t76 = t125 * t140;
t50 = -t141 * t205 + t76;
t269 = t50 * t148;
t262 = t125 * t141;
t261 = t126 * t148;
t84 = t133 * t140;
t260 = t133 * t141;
t252 = t145 * t148;
t74 = rSges(4,1) * t253 + rSges(4,2) * t252;
t146 = t158 * t295;
t245 = t149 * t292 + t146;
t239 = pkin(3) * t252;
t238 = t160 * t295;
t232 = t318 * t141 + t148 * t237;
t229 = t144 * t289 + t245;
t226 = -pkin(4) - t277;
t225 = -t244 / 0.2e1;
t223 = -t243 / 0.2e1;
t222 = t243 / 0.2e1;
t66 = Icges(6,5) * t254 - Icges(6,6) * t255 + Icges(6,3) * t140;
t221 = -t66 - t270;
t105 = rSges(3,1) * t150 - t149 * rSges(3,2);
t104 = pkin(4) * t259;
t220 = t104 + t232;
t93 = -rSges(3,1) * t250 + rSges(3,2) * t251;
t75 = -rSges(4,1) * t252 + rSges(4,2) * t253;
t213 = -rSges(4,1) * t144 - rSges(4,2) * t145;
t211 = -rSges(6,2) * t157 + t277;
t41 = t157 * t69 + t159 * t67;
t210 = -t270 + t271;
t42 = t157 * t70 + t159 * t68;
t25 = t141 * t66 - t257 * t70 + t68 * t258;
t203 = -t100 - t293;
t202 = t137 + t217;
t89 = -t105 * t155 - t234;
t24 = -t257 * t69 + t285;
t200 = (t140 * t25 + t141 * t24) * qJD(5);
t26 = -t255 * t67 + t284;
t27 = t140 * t66 - t314;
t199 = (t140 * t27 + t141 * t26) * qJD(5);
t198 = -t150 * t292 - t238;
t196 = -t239 - t241;
t195 = t248 + t286 - t291;
t194 = t130 * t148;
t193 = t207 * t148;
t191 = t235 + t242;
t189 = t213 - t294;
t188 = t212 - t291;
t187 = t220 + t242;
t183 = t75 - t241;
t181 = -t261 * t140 - t141 * t307 + t148 * t209;
t180 = t140 * t307 - t261 * t141 + t148 * t210;
t179 = t202 - t293;
t178 = t126 * qJD(5) + t148 * t205;
t176 = t195 - t294;
t175 = t124 - t236 - t288 - t290;
t174 = t188 - t294;
t173 = -t145 * t289 + t198;
t172 = t140 * t204 + t141 * t71;
t171 = t191 + t240;
t170 = t190 + t239;
t168 = t175 - t293;
t10 = t199 + t269;
t14 = -qJD(5) * t210 + t157 * (t140 * t304 - t141 * t194) + t159 * (t140 * t306 - t141 * t193);
t15 = -qJD(5) * t209 + t157 * (-t140 * t194 - t141 * t304) + t159 * (-t140 * t193 - t141 * t306);
t20 = t178 * t140 - t141 * t305;
t21 = t140 * t305 + t178 * t141;
t49 = t140 * t205 + t262;
t46 = t49 * t148;
t9 = t46 + t200;
t167 = (t46 + ((t27 + t285 + t314) * t141 + (-t26 + (t221 - t271) * t141 + t25 + t284) * t140) * qJD(5)) * t225 + (t10 - t269 + ((t25 + (-t66 + t271) * t141 - t284) * t141 + (t140 * t221 - t24 + t285) * t140) * qJD(5)) * t223 + (t14 + t21) * t222 + (t15 + t20 + t9) * t244 / 0.2e1 + (-qJD(5) * t205 + t111 * t159 + t112 * t157 + (t41 + t49) * t225 + (t42 + t50) * t222) * t148;
t43 = t148 * t248 - t232;
t44 = t148 * t201 + t233;
t163 = (-t44 - t62) * t140 + (t43 + t61) * t141;
t116 = t211 * qJD(5);
t18 = t116 * t244 + (t104 - t43 + (-pkin(8) * t148 + t318) * t141) * t148 + t229;
t19 = -t116 * t243 + (t44 + t316) * t148 + t173;
t162 = (t18 * t298 + t19 * t226) * t140 + (t29 * (-t123 + t291) - t28 * (-t275 - t287 - t290) + (-t226 * t28 + t29 * t298) * t141) * t148;
t90 = t148 * t213;
t73 = t155 * t93 - t238;
t72 = t155 * t92 + t146;
t63 = t191 - t90;
t54 = t148 * t75 + t198;
t53 = t148 * t74 + t245;
t51 = t171 - t82;
t48 = t148 * (-rSges(5,1) * t256 + t101) + t173;
t47 = -t148 * t249 + t229;
t32 = qJD(5) * t172 + qJD(4);
t11 = t163 * qJD(5);
t1 = [m(3) * (t72 * (-t105 - t296) + t73 * (-t214 - t297) + (t89 - t93 + t234) * t88) + t167 + (t18 * (t168 - t296) + t29 * (t187 + t235) + t19 * (t176 - t297) + t162 + (t190 - t170 - t29 + t310) * t28) * m(6) + (t47 * (t179 - t296) + t52 * (t171 - t249) + t48 * (t174 - t297) + (t234 - t170 - t52 + t315) * t51) * m(5) + (t53 * (t203 - t296) + t64 * (t191 + t74) + t54 * (t189 - t297) + (-t183 + t234) * t63) * m(4); t167 + (t18 * t168 + t19 * t176 + t162 + (t187 - t166) * t29 + (t241 + t196 + t310) * t28) * m(6) + (t48 * t174 + t47 * t179 + (t196 + t315) * t51 + t311) * m(5) + (t54 * t189 + t53 * t203 + (-t91 - t241 - t183) * t63 + (t90 + t74) * t64) * m(4) + (-t72 * t105 - t214 * t73 - t88 * t93 + t89 * t92 - (t88 * t105 + t214 * t89) * t155) * m(3); t167 + (t175 * t18 + t19 * t195 + t162 + (t220 - t177) * t29 + (-t239 + t310) * t28) * m(6) + (t48 * t188 + t47 * t202 + (-t239 + t320) * t51 + t311) * m(5) + (-t53 * t100 + t213 * t54 - t63 * t75 + t64 * t74 - (t63 * t100 - t213 * t64) * t148) * m(4); m(6) * t11; t148 * ((t148 * t42 + t14) * t141 + (-t148 * t41 + t15) * t140) / 0.2e1 - t148 * ((-t157 * t247 + t246 * t159) * t148 + ((t140 * t280 + t141 * t281) * t159 + (-t140 * t282 - t141 * t283) * t157) * qJD(5)) / 0.2e1 + ((t243 * t76 + t261) * t141 + (t319 + (t302 * t140 + (-t301 - t262) * t141) * qJD(5)) * t140) * t223 + ((-t244 * t262 + t261) * t140 + (-t319 + (t301 * t141 + (-t302 + t76) * t140) * qJD(5)) * t141) * t225 + (t148 * t20 + ((t180 * t140 + t148 * t27) * t141 + (t181 * t140 - t148 * t26) * t140) * t317) * t140 / 0.2e1 + (t148 * t21 + ((t180 * t141 + t148 * t25) * t141 + (t181 * t141 - t148 * t24) * t140) * t317) * t141 / 0.2e1 - (t9 + t200) * t259 / 0.2e1 + (t10 + t199) * t256 / 0.2e1 + (t11 * t172 + t163 * t32 + t18 * t84 - t19 * t260 - (t260 * t29 - t28 * t84) * t148 - (t32 * (-t140 * t84 - t141 * t260) + t312 * t211) * qJD(5) + (t29 * t256 - t28 * t259) * t133 + t312 * t116) * m(6);];
tauc = t1(:);
