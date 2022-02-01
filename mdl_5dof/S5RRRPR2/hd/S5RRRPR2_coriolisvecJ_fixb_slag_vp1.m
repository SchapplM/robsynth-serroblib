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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:31
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 11:30:33
% EndTime: 2022-01-20 11:30:39
% DurationCPUTime: 3.19s
% Computational Cost: add. (10345->319), mult. (5823->398), div. (0->0), fcn. (4286->10), ass. (0->210)
t155 = qJ(1) + qJ(2);
t151 = qJ(3) + t155;
t143 = pkin(9) + t151;
t139 = cos(t143);
t145 = cos(t151);
t140 = pkin(3) * t145;
t301 = -t139 * rSges(5,1) - t140;
t154 = qJD(1) + qJD(2);
t147 = qJD(3) + t154;
t253 = t139 * t147;
t101 = pkin(8) * t253;
t138 = sin(t143);
t156 = sin(qJ(5));
t255 = t138 * t156;
t115 = rSges(6,2) * t255;
t158 = cos(qJ(5));
t270 = rSges(6,2) * t158;
t231 = qJD(5) * t270;
t202 = rSges(6,3) * t253 + t147 * t115 - t139 * t231;
t242 = qJD(5) * t156;
t229 = t139 * t242;
t247 = t139 * rSges(6,3) + t115;
t254 = t138 * t158;
t69 = rSges(6,1) * t254 - t247;
t60 = t147 * t69;
t135 = t139 * pkin(8);
t93 = pkin(4) * t138 - t135;
t300 = -rSges(6,1) * t229 + t147 * t93 + t101 + t202 + t60;
t157 = sin(qJ(1));
t269 = pkin(1) * qJD(1);
t234 = t157 * t269;
t148 = sin(t155);
t149 = cos(t155);
t102 = rSges(3,1) * t148 + rSges(3,2) * t149;
t260 = t102 * t154;
t85 = -t234 - t260;
t299 = 0.2e1 * qJD(5);
t297 = -rSges(5,2) * t138 - t301;
t295 = t147 * t297;
t251 = t139 * t158;
t116 = rSges(6,1) * t251;
t294 = t138 * rSges(6,3) + t116;
t136 = t139 * pkin(4);
t94 = t138 * pkin(8) + t136;
t150 = Icges(6,4) * t158;
t205 = -Icges(6,2) * t156 + t150;
t121 = Icges(6,1) * t156 + t150;
t249 = t148 * t154;
t240 = pkin(2) * t249;
t193 = -t234 - t240;
t293 = t193 + t234;
t252 = t139 * t156;
t235 = rSges(6,2) * t252;
t70 = -t235 + t294;
t292 = t140 + t70 + t94;
t291 = -t240 + t300;
t118 = Icges(6,5) * t158 - Icges(6,6) * t156;
t117 = Icges(6,5) * t156 + Icges(6,6) * t158;
t177 = Icges(6,3) * t147 - qJD(5) * t117;
t195 = t205 * t139;
t66 = Icges(6,6) * t138 + t195;
t265 = t156 * t66;
t262 = Icges(6,4) * t156;
t122 = Icges(6,1) * t158 - t262;
t196 = t122 * t139;
t68 = Icges(6,5) * t138 + t196;
t207 = -t158 * t68 + t265;
t256 = t138 * t147;
t290 = -t118 * t256 + t139 * t177 + t147 * t207;
t194 = t118 * t139;
t65 = Icges(6,4) * t254 - Icges(6,2) * t255 - Icges(6,6) * t139;
t266 = t156 * t65;
t114 = Icges(6,4) * t255;
t67 = Icges(6,1) * t254 - Icges(6,5) * t139 - t114;
t208 = -t158 * t67 + t266;
t289 = t138 * t177 + (t194 + t208) * t147;
t119 = Icges(6,2) * t158 + t262;
t204 = t156 * t119 - t158 * t121;
t288 = t118 * qJD(5) + t147 * t204;
t63 = Icges(6,5) * t254 - Icges(6,6) * t255 - Icges(6,3) * t139;
t24 = -t138 * t208 - t139 * t63;
t273 = -Icges(6,2) * t254 - t114 + t67;
t275 = t121 * t138 + t65;
t287 = -t156 * t273 - t158 * t275;
t271 = rSges(6,1) * t158;
t226 = -pkin(4) - t271;
t144 = sin(t151);
t279 = pkin(3) * t144;
t159 = cos(qJ(1));
t233 = t159 * t269;
t248 = t149 * t154;
t239 = pkin(2) * t248;
t192 = t233 + t239;
t124 = rSges(6,1) * t156 + t270;
t244 = qJD(5) * t138;
t92 = t124 * t244;
t29 = t147 * t292 + t192 - t92;
t241 = t29 * t279;
t243 = qJD(5) * t139;
t230 = t124 * t243;
t191 = -t230 - t240;
t175 = t191 - t234;
t28 = (-t69 - t93 - t279) * t147 + t175;
t161 = (t28 * (-t116 - t136 - t140) - t241 + (t28 * (-rSges(6,3) - pkin(8)) + t29 * t226) * t138) * t147;
t232 = t147 * t235 + (rSges(6,1) * t242 + t231) * t138;
t286 = t161 + (-t92 + t232) * t28 - (-t28 * t292 - t241) * t147;
t146 = t147 ^ 2;
t283 = t147 / 0.2e1;
t282 = pkin(1) * t157;
t281 = pkin(2) * t148;
t280 = pkin(2) * t154 ^ 2;
t278 = pkin(3) * t146;
t152 = t159 * pkin(1);
t277 = -t138 * t63 - t67 * t251;
t64 = Icges(6,3) * t138 + t194;
t276 = t138 * t64 + t68 * t251;
t274 = -t121 * t139 - t66;
t272 = -t119 * t139 + t68;
t137 = t145 * rSges(4,1);
t97 = rSges(4,1) * t144 + rSges(4,2) * t145;
t87 = t147 * t97;
t98 = -rSges(4,2) * t144 + t137;
t267 = t147 * t98;
t258 = t117 * t139;
t49 = -t138 * t204 - t258;
t264 = t49 * t147;
t259 = t117 * t138;
t257 = t118 * t147;
t250 = t144 * t147;
t246 = -t119 + t122;
t245 = t121 + t205;
t238 = pkin(3) * t250;
t160 = qJD(1) ^ 2;
t237 = t160 * t282;
t236 = t160 * t152;
t90 = rSges(5,1) * t138 + rSges(5,2) * t139;
t189 = -t90 - t279;
t224 = -t244 / 0.2e1;
t221 = t243 / 0.2e1;
t55 = t68 * t254;
t219 = t139 * t64 - t55;
t218 = -t63 + t265;
t103 = t149 * rSges(3,1) - rSges(3,2) * t148;
t89 = rSges(3,1) * t248 - rSges(3,2) * t249;
t74 = -rSges(4,2) * t250 + t137 * t147;
t142 = pkin(2) * t149;
t214 = t142 + t98;
t211 = -rSges(6,2) * t156 + t271;
t210 = -t138 * t29 - t139 * t28;
t209 = t138 * t69 + t139 * t70;
t41 = t156 * t67 + t158 * t65;
t42 = t156 * t68 + t158 * t66;
t203 = t142 + t297;
t25 = -t255 * t66 - t219;
t200 = (t138 * t25 - t139 * t24) * qJD(5);
t26 = -t252 * t65 - t277;
t27 = -t252 * t66 + t276;
t199 = (t138 * t27 - t139 * t26) * qJD(5);
t198 = -t148 * t280 - t237;
t197 = -t149 * t280 - t236;
t190 = -t97 - t281;
t81 = t147 * t90;
t188 = -t238 - t81 - t240;
t187 = -t156 * t272 + t158 * t274;
t186 = -t74 - t239;
t184 = t189 - t281;
t183 = t142 + t292;
t182 = (-t156 * t245 + t158 * t246) * t147;
t181 = -t144 * t278 + t198;
t180 = -t145 * t278 + t197;
t179 = Icges(6,5) * t147 - qJD(5) * t121;
t178 = Icges(6,6) * t147 - qJD(5) * t119;
t174 = t138 * t226 + t135 + t247 - t279;
t43 = (-t147 * t254 - t229) * rSges(6,1) + t202;
t44 = t294 * t147 - t232;
t172 = (t43 + t60) * t139 + (-t147 * t70 + t44) * t138;
t50 = -t139 * t204 + t259;
t46 = t50 * t147;
t10 = t46 + t199;
t107 = t205 * qJD(5);
t108 = t122 * qJD(5);
t14 = -qJD(5) * t208 + t156 * (t138 * t179 + t147 * t196) + t158 * (t138 * t178 + t147 * t195);
t15 = -qJD(5) * t207 + t156 * (-t122 * t256 + t139 * t179) + t158 * (t139 * t178 - t205 * t256);
t164 = -t107 * t156 + t108 * t158 + t117 * t147 + (-t119 * t158 - t121 * t156) * qJD(5);
t20 = t288 * t138 + t164 * t139;
t21 = t164 * t138 - t288 * t139;
t9 = t200 + t264;
t171 = (t46 + ((t25 - t55 + (t64 + t266) * t139 + t277) * t139 + t276 * t138) * qJD(5)) * t221 + (-qJD(5) * t204 + t107 * t158 + t108 * t156) * t147 + (-t264 + ((t139 * t218 + t27 - t276) * t139 + (t138 * t218 + t219 + t26) * t138) * qJD(5) + t9) * t224 + (t15 + t20) * t244 / 0.2e1 - (t14 + t21 + t10) * t243 / 0.2e1 + ((t41 + t49) * t138 + (t42 + t50) * t139) * qJD(5) * t283;
t167 = t174 - t281;
t51 = t147 * t189 + t193;
t52 = t192 + t295;
t162 = (t52 * t189 + t51 * t301) * t147;
t110 = t211 * qJD(5);
t99 = rSges(5,2) * t256;
t86 = t103 * t154 + t233;
t83 = t124 * t139;
t82 = t124 * t138;
t72 = -t154 * t89 - t236;
t71 = -t154 * t260 - t237;
t62 = t192 + t267;
t61 = t193 - t87;
t54 = -t147 * t74 + t197;
t53 = -t147 * t87 + t198;
t48 = -t147 * (rSges(5,1) * t253 - t99) + t180;
t47 = -t146 * t90 + t181;
t32 = qJD(5) * t209 + qJD(4);
t19 = -t110 * t243 + (-t94 * t147 - t44 + t92) * t147 + t180;
t18 = -t110 * t244 + (-pkin(4) * t256 + t101 - t230 + t43) * t147 + t181;
t11 = t172 * qJD(5);
t1 = [m(3) * (t72 * (-t102 - t282) + t71 * (t103 + t152) + (-t89 - t233 + t86) * t85) + t171 + (t19 * (t167 - t282) + t28 * (-t192 + t232) + t18 * (t152 + t183) + t161 + (-t175 + t28 + t238 - t234 + t291) * t29) * m(6) + (t48 * (t184 - t282) + t51 * (-t192 + t99) + t47 * (t152 + t203) + t162 + (-t188 + t51 + t293) * t52) * m(5) + (t54 * (t190 - t282) + t61 * (t186 - t233) + t53 * (t152 + t214) + (t61 - t293 - t240) * t62) * m(4); t171 + (t19 * t167 + t18 * t183 + (-t191 + t291) * t29 + t286) * m(6) + (t48 * t184 + t47 * t203 + t162 + (-t188 - t240) * t52 + (t99 + t295) * t51) * m(5) + (t54 * t190 + t53 * t214 + (t186 + t239 + t267) * t61) * m(4) + (-(-t86 * t102 - t103 * t85) * t154 - t102 * t72 + t103 * t71 - t85 * t89 - t86 * t260) * m(3); t171 + (t19 * t174 + t18 * t292 + (t230 + t300) * t29 + t286) * m(6) + (t52 * t81 - (-t279 * t52 - t297 * t51) * t147 + t48 * t189 + t47 * t297 + t51 * t99 + t162) * m(5) + (t53 * t98 - t54 * t97 - t61 * t74 - t62 * t87 - (-t61 * t98 - t62 * t97) * t147) * m(4); m(6) * t11; ((t147 * t42 - t14) * t139 + (t147 * t41 + t15) * t138) * t283 + ((-t244 * t258 + t257) * t138 + (t182 + (-t287 * t139 + (t259 + t187) * t138) * qJD(5)) * t139) * t224 + ((-t243 * t259 - t257) * t139 + (t182 + (t187 * t138 + (-t287 + t258) * t139) * qJD(5)) * t138) * t221 - t147 * ((t246 * t156 + t245 * t158) * t147 + ((t138 * t272 - t139 * t273) * t158 + (t138 * t274 + t139 * t275) * t156) * qJD(5)) / 0.2e1 + (t20 * t147 + ((-t289 * t138 + t147 * t27) * t139 + (t290 * t138 + t147 * t26) * t138) * t299) * t138 / 0.2e1 - (t147 * t21 + ((t289 * t139 + t147 * t25) * t139 + (-t290 * t139 + t147 * t24) * t138) * t299) * t139 / 0.2e1 + (t9 + t200) * t256 / 0.2e1 + (t10 + t199) * t253 / 0.2e1 + (t11 * t209 + t32 * t172 + t210 * t110 + ((-t147 * t29 - t19) * t139 + (t147 * t28 - t18) * t138) * t124 - (t28 * t82 - t29 * t83) * t147 - (t32 * (-t138 * t82 - t139 * t83) + t210 * t211) * qJD(5)) * m(6);];
tauc = t1(:);
