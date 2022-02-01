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
% Datum: 2022-01-20 10:34
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 10:33:59
% EndTime: 2022-01-20 10:34:06
% DurationCPUTime: 3.30s
% Computational Cost: add. (10079->303), mult. (5719->396), div. (0->0), fcn. (4238->10), ass. (0->199)
t154 = qJ(1) + qJ(2);
t147 = pkin(9) + t154;
t145 = qJ(4) + t147;
t140 = cos(t145);
t153 = qJD(1) + qJD(2);
t146 = qJD(4) + t153;
t242 = t140 * t146;
t101 = pkin(8) * t242;
t139 = sin(t145);
t155 = sin(qJ(5));
t244 = t139 * t155;
t114 = rSges(6,2) * t244;
t157 = cos(qJ(5));
t260 = rSges(6,2) * t157;
t225 = qJD(5) * t260;
t196 = rSges(6,3) * t242 + t146 * t114 - t140 * t225;
t232 = qJD(5) * t155;
t223 = t140 * t232;
t238 = t140 * rSges(6,3) + t114;
t243 = t139 * t157;
t69 = rSges(6,1) * t243 - t238;
t60 = t146 * t69;
t135 = t140 * pkin(8);
t93 = pkin(4) * t139 - t135;
t290 = -rSges(6,1) * t223 + t146 * t93 + t101 + t196 + t60;
t158 = cos(qJ(1));
t259 = pkin(1) * qJD(1);
t228 = t158 * t259;
t143 = cos(t147);
t149 = cos(t154);
t144 = pkin(2) * t149;
t282 = pkin(3) * t143 + t144;
t177 = t153 * t282 + t228;
t133 = t140 * rSges(5,1);
t91 = -rSges(5,2) * t139 + t133;
t256 = t146 * t91;
t52 = t177 + t256;
t90 = rSges(5,1) * t139 + rSges(5,2) * t140;
t81 = t146 * t90;
t289 = t52 * t81;
t137 = t143 * rSges(4,1);
t288 = -t137 - t144;
t156 = sin(qJ(1));
t227 = t156 * t259;
t148 = sin(t154);
t102 = rSges(3,1) * t148 + rSges(3,2) * t149;
t249 = t102 * t153;
t85 = -t227 - t249;
t287 = 0.2e1 * qJD(5);
t142 = sin(t147);
t261 = rSges(4,2) * t142;
t285 = -t261 - t288;
t271 = pkin(2) * t148;
t208 = -pkin(3) * t142 - t271;
t176 = t208 * t153 - t227;
t51 = t176 - t81;
t240 = t140 * t157;
t284 = rSges(6,1) * t240 + t139 * rSges(6,3);
t283 = -t140 * pkin(4) - t139 * pkin(8);
t150 = Icges(6,4) * t157;
t199 = -Icges(6,2) * t155 + t150;
t121 = Icges(6,1) * t155 + t150;
t262 = rSges(6,1) * t157;
t220 = -pkin(4) - t262;
t124 = rSges(6,1) * t155 + t260;
t233 = qJD(5) * t140;
t224 = t124 * t233;
t169 = t176 - t224;
t28 = (-t69 - t93) * t146 + t169;
t258 = t140 * t28;
t241 = t140 * t155;
t229 = rSges(6,2) * t241;
t70 = -t229 + t284;
t188 = t70 - t283;
t234 = qJD(5) * t139;
t92 = t124 * t234;
t280 = -t188 * t146 + t92;
t29 = t177 - t280;
t161 = (t220 * t258 + (t28 * (-rSges(6,3) - pkin(8)) + t29 * t220) * t139) * t146;
t226 = t146 * t229 + (rSges(6,1) * t232 + t225) * t139;
t281 = t161 + (t224 + t290) * t29 + (t226 - t280) * t28;
t118 = Icges(6,5) * t157 - Icges(6,6) * t155;
t117 = Icges(6,5) * t155 + Icges(6,6) * t157;
t180 = Icges(6,3) * t146 - t117 * qJD(5);
t191 = t199 * t140;
t66 = Icges(6,6) * t139 + t191;
t254 = t155 * t66;
t251 = Icges(6,4) * t155;
t122 = Icges(6,1) * t157 - t251;
t192 = t122 * t140;
t68 = Icges(6,5) * t139 + t192;
t201 = -t157 * t68 + t254;
t245 = t139 * t146;
t279 = -t118 * t245 + t180 * t140 + t201 * t146;
t190 = t118 * t140;
t65 = Icges(6,4) * t243 - Icges(6,2) * t244 - Icges(6,6) * t140;
t255 = t155 * t65;
t113 = Icges(6,4) * t244;
t67 = Icges(6,1) * t243 - Icges(6,5) * t140 - t113;
t202 = -t157 * t67 + t255;
t278 = t180 * t139 + (t190 + t202) * t146;
t119 = Icges(6,2) * t157 + t251;
t198 = t119 * t155 - t121 * t157;
t277 = t118 * qJD(5) + t198 * t146;
t63 = Icges(6,5) * t243 - Icges(6,6) * t244 - Icges(6,3) * t140;
t24 = -t202 * t139 - t140 * t63;
t264 = -Icges(6,2) * t243 - t113 + t67;
t266 = t121 * t139 + t65;
t276 = -t264 * t155 - t266 * t157;
t152 = t153 ^ 2;
t273 = t146 / 0.2e1;
t272 = pkin(1) * t156;
t270 = pkin(2) * t152;
t151 = t158 * pkin(1);
t269 = t52 * t90;
t268 = -t139 * t63 - t67 * t240;
t64 = Icges(6,3) * t139 + t190;
t267 = t139 * t64 + t68 * t240;
t265 = -t121 * t140 - t66;
t263 = -t119 * t140 + t68;
t141 = t149 * rSges(3,1);
t247 = t117 * t140;
t49 = -t198 * t139 - t247;
t253 = t49 * t146;
t248 = t117 * t139;
t246 = t118 * t146;
t239 = t148 * t153;
t237 = -t119 + t122;
t236 = t121 + t199;
t159 = qJD(1) ^ 2;
t231 = t159 * t272;
t230 = t159 * t151;
t97 = rSges(4,1) * t142 + rSges(4,2) * t143;
t189 = -t97 - t271;
t218 = -t234 / 0.2e1;
t215 = t233 / 0.2e1;
t53 = t68 * t243;
t213 = t140 * t64 - t53;
t212 = -t63 + t254;
t103 = -rSges(3,2) * t148 + t141;
t99 = rSges(5,2) * t245;
t74 = rSges(5,1) * t242 - t99;
t89 = -rSges(3,2) * t239 + t153 * t141;
t205 = -rSges(6,2) * t155 + t262;
t204 = -t139 * t29 - t258;
t203 = t139 * t69 + t140 * t70;
t41 = t155 * t67 + t157 * t65;
t42 = t155 * t68 + t157 * t66;
t197 = t91 + t282;
t25 = -t66 * t244 - t213;
t194 = (t139 * t25 - t140 * t24) * qJD(5);
t26 = -t65 * t241 - t268;
t27 = -t66 * t241 + t267;
t193 = (t139 * t27 - t140 * t26) * qJD(5);
t187 = -t263 * t155 + t265 * t157;
t186 = t220 * t139 + t135 + t238;
t185 = t208 - t90;
t184 = t188 + t282;
t183 = (-t236 * t155 + t237 * t157) * t146;
t182 = Icges(6,5) * t146 - qJD(5) * t121;
t181 = Icges(6,6) * t146 - t119 * qJD(5);
t179 = t208 * t152 - t231;
t178 = -t152 * t282 - t230;
t43 = (-t146 * t243 - t223) * rSges(6,1) + t196;
t44 = t284 * t146 - t226;
t175 = (t43 + t60) * t140 + (-t146 * t70 + t44) * t139;
t50 = -t198 * t140 + t248;
t46 = t50 * t146;
t10 = t46 + t193;
t107 = t199 * qJD(5);
t108 = t122 * qJD(5);
t14 = -t202 * qJD(5) + t155 * (t182 * t139 + t146 * t192) + t157 * (t181 * t139 + t146 * t191);
t15 = -t201 * qJD(5) + t155 * (-t122 * t245 + t182 * t140) + t157 * (t181 * t140 - t199 * t245);
t163 = -t107 * t155 + t108 * t157 + t117 * t146 + (-t119 * t157 - t121 * t155) * qJD(5);
t20 = t277 * t139 + t163 * t140;
t21 = t163 * t139 - t277 * t140;
t9 = t194 + t253;
t174 = (t46 + ((t25 - t53 + (t64 + t255) * t140 + t268) * t140 + t267 * t139) * qJD(5)) * t215 + (-t198 * qJD(5) + t107 * t157 + t108 * t155) * t146 + (t9 - t253 + ((t212 * t140 - t267 + t27) * t140 + (t212 * t139 + t213 + t26) * t139) * qJD(5)) * t218 + (t15 + t20) * t234 / 0.2e1 - (t10 + t14 + t21) * t233 / 0.2e1 + ((t41 + t49) * t139 + (t42 + t50) * t140) * qJD(5) * t273;
t170 = t186 + t208;
t61 = t189 * t153 - t227;
t62 = t153 * t285 + t228;
t162 = (t62 * t189 + t61 * t288) * t153;
t116 = t153 * t261;
t109 = t205 * qJD(5);
t87 = t153 * t97;
t86 = t103 * t153 + t228;
t83 = t124 * t140;
t82 = t124 * t139;
t72 = -t153 * t89 - t230;
t71 = -t153 * t249 - t231;
t57 = -t230 - t149 * t270 - t153 * (t153 * t137 - t116);
t56 = -t148 * t270 - t97 * t152 - t231;
t48 = -t146 * t74 + t178;
t47 = -t146 * t81 + t179;
t32 = t203 * qJD(5) + qJD(3);
t19 = -t109 * t233 + (t283 * t146 - t44 + t92) * t146 + t178;
t18 = -t109 * t234 + (-pkin(4) * t245 + t101 - t224 + t43) * t146 + t179;
t11 = t175 * qJD(5);
t1 = [m(3) * (t72 * (-t102 - t272) + t71 * (t103 + t151) + (-t89 - t228 + t86) * t85) + t174 + (t19 * (t170 - t272) + t28 * (t226 - t228) + t18 * (t151 + t184) + (t29 * t208 - t28 * t282) * t153 + t161 + (-t227 + t28 - t169 + t290) * t29) * m(6) + (t48 * (t185 - t272) + t51 * (-t74 - t228) + t47 * (t151 + t197) + (t52 * t208 - t282 * t51) * t153 + t52 * (-rSges(5,1) * t245 - rSges(5,2) * t242 - t227)) * m(5) + (t57 * (t189 - t272) + t61 * (t116 - t228) + t56 * (t151 + t285) + t162 - (-pkin(2) * t239 - t61 - t87) * t62) * m(4); t174 + (t61 * t116 + t57 * t189 + t56 * t285 + t162 + t62 * t87 - (-t62 * t271 - t285 * t61) * t153) * m(4) + (-t102 * t72 + t103 * t71 - t85 * t89 - t86 * t249 - (-t102 * t86 - t103 * t85) * t153) * m(3) + (t19 * t170 + t18 * t184 + t281) * m(6) + (-t269 * t146 + t48 * t185 + t47 * t197 + t289 + (-t133 * t146 + t256 + t99) * t51) * m(5); m(6) * t11; t174 + (t18 * t188 + t19 * t186 + t281) * m(6) + (t47 * t91 - t48 * t90 - t51 * t74 - t289 - (-t51 * t91 - t269) * t146) * m(5); ((t146 * t42 - t14) * t140 + (t146 * t41 + t15) * t139) * t273 + ((-t234 * t247 + t246) * t139 + (t183 + (-t276 * t140 + (t248 + t187) * t139) * qJD(5)) * t140) * t218 + ((-t233 * t248 - t246) * t140 + (t183 + (t187 * t139 + (-t276 + t247) * t140) * qJD(5)) * t139) * t215 - t146 * ((t237 * t155 + t236 * t157) * t146 + ((t263 * t139 - t264 * t140) * t157 + (t265 * t139 + t266 * t140) * t155) * qJD(5)) / 0.2e1 + (t146 * t20 + ((-t278 * t139 + t146 * t27) * t140 + (t279 * t139 + t146 * t26) * t139) * t287) * t139 / 0.2e1 - (t146 * t21 + ((t278 * t140 + t146 * t25) * t140 + (-t279 * t140 + t146 * t24) * t139) * t287) * t140 / 0.2e1 + (t9 + t194) * t245 / 0.2e1 + (t10 + t193) * t242 / 0.2e1 + (t11 * t203 + t32 * t175 + t204 * t109 + ((-t146 * t29 - t19) * t140 + (t146 * t28 - t18) * t139) * t124 - (t28 * t82 - t29 * t83) * t146 - (t32 * (-t139 * t82 - t140 * t83) + t204 * t205) * qJD(5)) * m(6);];
tauc = t1(:);
