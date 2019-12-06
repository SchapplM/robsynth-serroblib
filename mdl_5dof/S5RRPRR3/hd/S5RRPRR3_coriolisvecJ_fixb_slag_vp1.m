% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPRR3
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
% Datum: 2019-12-05 18:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRR3_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR3_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR3_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR3_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR3_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:30:18
% EndTime: 2019-12-05 18:30:24
% DurationCPUTime: 3.50s
% Computational Cost: add. (10079->309), mult. (5719->395), div. (0->0), fcn. (4238->10), ass. (0->206)
t156 = sin(qJ(5));
t158 = cos(qJ(5));
t133 = rSges(6,1) * t156 + rSges(6,2) * t158;
t315 = qJD(5) * t133;
t155 = qJ(1) + qJ(2);
t148 = pkin(9) + t155;
t145 = qJ(4) + t148;
t141 = cos(t145);
t252 = t141 * t156;
t120 = rSges(6,2) * t252;
t140 = sin(t145);
t154 = qJD(1) + qJD(2);
t147 = qJD(4) + t154;
t232 = t147 * t120 + t140 * t315;
t241 = qJD(5) * t140;
t95 = t133 * t241;
t286 = pkin(4) * t141;
t96 = pkin(8) * t140 + t286;
t313 = -t147 * t96 + t95;
t251 = t141 * t158;
t235 = rSges(6,1) * t251;
t194 = -rSges(6,3) * t140 - t235;
t71 = -t120 - t194;
t62 = t147 * t71;
t317 = -t232 - t62 + t313;
t151 = Icges(6,4) * t158;
t201 = -Icges(6,2) * t156 + t151;
t311 = Icges(6,1) * t156 + t151;
t243 = t311 + t201;
t264 = Icges(6,4) * t156;
t127 = Icges(6,2) * t158 + t264;
t130 = Icges(6,1) * t158 - t264;
t244 = t127 - t130;
t316 = (t243 * t156 + t244 * t158) * t147;
t314 = 0.2e1 * qJD(5);
t157 = sin(qJ(1));
t271 = pkin(1) * qJD(1);
t234 = t157 * t271;
t149 = sin(t155);
t150 = cos(t155);
t208 = rSges(3,1) * t149 + rSges(3,2) * t150;
t92 = t208 * t154;
t88 = t92 + t234;
t68 = Icges(6,4) * t251 - Icges(6,2) * t252 + Icges(6,6) * t140;
t118 = Icges(6,4) * t252;
t70 = Icges(6,1) * t251 + Icges(6,5) * t140 - t118;
t203 = t156 * t68 - t158 * t70;
t312 = t203 * t141;
t255 = t140 * t156;
t246 = rSges(6,2) * t255 + t141 * rSges(6,3);
t143 = sin(t148);
t139 = t143 * rSges(4,2);
t144 = cos(t148);
t276 = rSges(4,1) * t144;
t100 = -t139 + t276;
t250 = t143 * t154;
t121 = rSges(4,2) * t250;
t289 = pkin(2) * t150;
t211 = -t276 - t289;
t247 = t150 * t154;
t310 = -pkin(2) * t247 - t121 + (-t100 - t211) * t154;
t212 = pkin(3) * t144 + t289;
t159 = cos(qJ(1));
t233 = t159 * t271;
t171 = t212 * t154 + t233;
t309 = t233 - t171;
t287 = pkin(3) * t143;
t290 = pkin(2) * t149;
t213 = t287 + t290;
t173 = t213 * t154 + t234;
t240 = qJD(5) * t141;
t285 = pkin(8) * t141;
t254 = t140 * t158;
t236 = rSges(6,1) * t254;
t197 = t236 - t246;
t61 = t147 * t197;
t196 = t133 * t240 + t61 - t147 * (-pkin(4) * t140 + t285);
t28 = t173 + t196;
t29 = t95 + (-t71 - t96) * t147 - t171;
t308 = t140 * t29 + t141 * t28;
t206 = -rSges(5,1) * t140 - rSges(5,2) * t141;
t82 = t147 * t206;
t51 = -t82 + t173;
t94 = rSges(5,1) * t141 - t140 * rSges(5,2);
t83 = t147 * t94;
t52 = -t171 - t83;
t253 = t141 * t147;
t256 = t140 * t147;
t74 = rSges(5,1) * t256 + rSges(5,2) * t253;
t75 = -rSges(5,1) * t253 + rSges(5,2) * t256;
t307 = -t51 * t75 + t52 * t74;
t104 = pkin(4) * t256;
t231 = t141 * t315 + t147 * t236;
t219 = t104 + t231;
t306 = (-t196 + t219) * t29 + t317 * t28;
t125 = Icges(6,5) * t156 + Icges(6,6) * t158;
t303 = -Icges(6,3) * t147 + t125 * qJD(5);
t302 = -Icges(6,6) * t147 + t127 * qJD(5);
t111 = t201 * qJD(5);
t112 = t130 * qJD(5);
t301 = (t127 * t158 + t156 * t311) * qJD(5) + t111 * t156 - t112 * t158 - t125 * t147;
t300 = -Icges(6,5) * t147 + qJD(5) * t311;
t277 = -Icges(6,2) * t251 - t118 + t70;
t279 = t141 * t311 + t68;
t298 = t277 * t156 + t279 * t158;
t117 = Icges(6,4) * t255;
t69 = -Icges(6,1) * t254 + Icges(6,5) * t141 + t117;
t278 = Icges(6,2) * t254 + t117 + t69;
t67 = Icges(6,6) * t141 - t201 * t140;
t280 = -t140 * t311 + t67;
t297 = -t278 * t156 - t280 * t158;
t294 = -rSges(6,3) - pkin(8);
t293 = pkin(1) * t157;
t292 = pkin(1) * t159;
t291 = pkin(1) * qJD(1) ^ 2;
t153 = t154 ^ 2;
t288 = pkin(2) * t153;
t126 = Icges(6,5) * t158 - Icges(6,6) * t156;
t65 = Icges(6,3) * t141 - t126 * t140;
t282 = t141 * t65 + t67 * t255;
t281 = t140 * t65 + t69 * t251;
t274 = rSges(6,1) * t158;
t273 = rSges(6,2) * t156;
t268 = t156 * t67;
t267 = t158 * t69;
t199 = t127 * t156 - t158 * t311;
t76 = t125 * t140;
t50 = -t199 * t141 + t76;
t266 = t50 * t147;
t259 = t125 * t141;
t258 = t126 * t147;
t84 = t133 * t140;
t257 = t133 * t141;
t249 = t144 * t154;
t248 = t149 * t154;
t245 = -rSges(4,1) * t250 - rSges(4,2) * t249;
t146 = t157 * t291;
t242 = t149 * t288 + t146;
t237 = t159 * t291;
t228 = t153 * t287 + t242;
t225 = -pkin(4) - t274;
t224 = -t241 / 0.2e1;
t222 = -t240 / 0.2e1;
t221 = t240 / 0.2e1;
t66 = Icges(6,5) * t251 - Icges(6,6) * t252 + Icges(6,3) * t140;
t220 = -t66 - t267;
t105 = rSges(3,1) * t150 - t149 * rSges(3,2);
t214 = t246 + t285;
t93 = -rSges(3,1) * t247 + rSges(3,2) * t248;
t207 = -rSges(4,1) * t143 - rSges(4,2) * t144;
t205 = -t273 + t274;
t41 = t156 * t69 + t158 * t67;
t204 = -t267 + t268;
t42 = t156 * t70 + t158 * t68;
t25 = t141 * t66 - t70 * t254 + t68 * t255;
t195 = t139 + t211;
t89 = -t105 * t154 - t233;
t24 = -t69 * t254 + t282;
t193 = (t140 * t25 + t141 * t24) * qJD(5);
t26 = -t67 * t252 + t281;
t27 = t140 * t66 - t312;
t192 = (t140 * t27 + t141 * t26) * qJD(5);
t191 = t130 * t147;
t190 = t201 * t147;
t188 = pkin(2) * t248 + t234;
t187 = t120 - t235 - t286;
t186 = t207 - t290;
t181 = -t258 * t140 - t141 * t303 + t203 * t147;
t180 = t140 * t303 - t258 * t141 + t204 * t147;
t179 = -t212 - t94;
t178 = t126 * qJD(5) + t199 * t147;
t177 = -t213 + t214;
t176 = t206 - t213;
t175 = t140 * t197 + t141 * t71;
t174 = -t212 * t153 - t237;
t170 = t187 - t212;
t10 = t192 + t266;
t14 = -t204 * qJD(5) + t156 * (t140 * t300 - t141 * t191) + t158 * (t140 * t302 - t141 * t190);
t15 = -t203 * qJD(5) + t156 * (-t140 * t191 - t141 * t300) + t158 * (-t140 * t190 - t141 * t302);
t20 = t178 * t140 - t141 * t301;
t21 = t140 * t301 + t178 * t141;
t49 = t199 * t140 + t259;
t46 = t49 * t147;
t9 = t46 + t193;
t169 = (t46 + ((t27 + t282 + t312) * t141 + (-t26 + (t220 - t268) * t141 + t25 + t281) * t140) * qJD(5)) * t224 + (-t266 + ((t25 + (-t66 + t268) * t141 - t281) * t141 + (t220 * t140 - t24 + t282) * t140) * qJD(5) + t10) * t222 + (t14 + t21) * t221 + (t15 + t20 + t9) * t241 / 0.2e1 + (-t199 * qJD(5) + t111 * t158 + t112 * t156 + (t41 + t49) * t224 + (t42 + t50) * t221) * t147;
t165 = (t28 * t212 + t29 * t213) * t154;
t43 = t147 * t246 - t231;
t44 = t194 * t147 + t232;
t163 = (-t44 - t62) * t140 + (t43 + t61) * t141;
t113 = t205 * qJD(5);
t18 = t113 * t241 + (t104 - t43 + (-pkin(8) * t147 + t315) * t141) * t147 + t228;
t19 = -t113 * t240 + (t44 + t313) * t147 + t174;
t162 = (t18 * t294 + t19 * t225) * t140 + ((-t29 * t273 - t28 * t294) * t140 + (-t28 * t225 + t29 * t294) * t141) * t147;
t161 = t165 + t162;
t90 = t154 * t207;
t73 = t154 * t93 - t237;
t72 = t154 * t92 + t146;
t64 = -t233 + (-t100 - t289) * t154;
t63 = t188 - t90;
t57 = -t237 - t150 * t288 + t154 * (-rSges(4,1) * t249 + t121);
t56 = -t154 * t245 + t242;
t48 = t147 * t75 + t174;
t47 = t147 * t74 + t228;
t32 = t175 * qJD(5) + qJD(3);
t11 = t163 * qJD(5);
t1 = [t169 + m(3) * (t72 * (-t105 - t292) + t73 * (-t208 - t293) + (t89 - t93 + t233) * t88) + (t18 * (t170 - t292) + t29 * (t219 + t234) + t19 * (t177 - t293) + t161 + (-t29 + t309 + t317) * t28) * m(6) + (t47 * (t179 - t292) + t52 * (t234 + t74) + t48 * (t176 - t293) + (t51 * t212 + t52 * t213) * t154 + (-t75 - t52 - t83 + t309) * t51) * m(5) + (t56 * (t195 - t292) + t64 * (t188 - t245) + t57 * (t186 - t293) + (-t64 + t310) * t63) * m(4); t169 + (t18 * t170 + t19 * t177 + t161 - t165 + t306) * m(6) + (t57 * t186 + t56 * t195 + t310 * t63 + (-t245 + t90) * t64) * m(4) + (-(t88 * t105 + t89 * t208) * t154 - t72 * t105 - t73 * t208 - t88 * t93 + t89 * t92) * m(3) + (t48 * t176 + t47 * t179 - t51 * t83 + t52 * t82 + t307) * m(5); m(6) * t11; t169 + (t18 * t187 + t19 * t214 + t162 + t306) * m(6) + (t48 * t206 - t47 * t94 - (-t52 * t206 + t51 * t94) * t147 + t307) * m(5); t147 * ((t147 * t42 + t14) * t141 + (-t147 * t41 + t15) * t140) / 0.2e1 - t147 * ((-t244 * t156 + t243 * t158) * t147 + ((t277 * t140 + t278 * t141) * t158 + (-t279 * t140 - t280 * t141) * t156) * qJD(5)) / 0.2e1 + ((t76 * t240 + t258) * t141 + (t316 + (t298 * t140 + (-t297 - t259) * t141) * qJD(5)) * t140) * t222 + ((-t241 * t259 + t258) * t140 + (-t316 + (t297 * t141 + (-t298 + t76) * t140) * qJD(5)) * t141) * t224 + (t147 * t20 + ((t180 * t140 + t147 * t27) * t141 + (t181 * t140 - t147 * t26) * t140) * t314) * t140 / 0.2e1 + (t147 * t21 + ((t180 * t141 + t147 * t25) * t141 + (t181 * t141 - t147 * t24) * t140) * t314) * t141 / 0.2e1 - (t9 + t193) * t256 / 0.2e1 + (t10 + t192) * t253 / 0.2e1 + (t11 * t175 + t32 * t163 + t18 * t84 - t19 * t257 - (t257 * t29 - t28 * t84) * t147 - (t32 * (-t140 * t84 - t141 * t257) + t308 * t205) * qJD(5) + (t253 * t29 - t256 * t28) * t133 + t308 * t113) * m(6);];
tauc = t1(:);
