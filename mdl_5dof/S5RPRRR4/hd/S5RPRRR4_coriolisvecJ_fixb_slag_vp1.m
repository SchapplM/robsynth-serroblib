% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRRR4
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRR4_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR4_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR4_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR4_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR4_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:52:02
% EndTime: 2020-01-03 11:52:10
% DurationCPUTime: 3.77s
% Computational Cost: add. (9883->295), mult. (5619->384), div. (0->0), fcn. (4190->10), ass. (0->187)
t158 = qJ(1) + pkin(9);
t151 = qJ(3) + t158;
t143 = sin(t151);
t157 = qJD(1) + qJD(3);
t251 = t143 * t157;
t236 = pkin(3) * t251;
t148 = sin(t158);
t160 = sin(qJ(1));
t154 = t160 * pkin(1);
t295 = pkin(2) * t148 + t154;
t291 = t295 * qJD(1);
t173 = t291 + t236;
t150 = qJD(4) + t157;
t145 = qJ(4) + t151;
t137 = sin(t145);
t138 = cos(t145);
t297 = t137 * rSges(5,1) + t138 * rSges(5,2);
t303 = t297 * t150;
t51 = t173 + t303;
t159 = sin(qJ(5));
t161 = cos(qJ(5));
t152 = Icges(6,4) * t161;
t202 = -Icges(6,2) * t159 + t152;
t191 = t202 * t150;
t263 = Icges(6,4) * t159;
t119 = Icges(6,1) * t161 - t263;
t192 = t119 * t150;
t65 = -Icges(6,5) * t138 + t119 * t137;
t267 = t161 * t65;
t63 = -Icges(6,6) * t138 + t137 * t202;
t269 = t159 * t63;
t206 = t267 - t269;
t294 = Icges(6,1) * t159 + t152;
t285 = -Icges(6,5) * t150 + qJD(5) * t294;
t116 = Icges(6,2) * t161 + t263;
t288 = -Icges(6,6) * t150 + t116 * qJD(5);
t252 = t138 * t161;
t253 = t138 * t159;
t64 = Icges(6,4) * t252 - Icges(6,2) * t253 + Icges(6,6) * t137;
t108 = Icges(6,4) * t253;
t66 = Icges(6,1) * t252 + Icges(6,5) * t137 - t108;
t307 = -t206 * qJD(5) - t159 * (-t137 * t285 + t138 * t192) - t161 * (-t137 * t288 + t138 * t191) + t150 * (t159 * t66 + t161 * t64);
t144 = cos(t151);
t296 = t143 * rSges(4,1) + t144 * rSges(4,2);
t81 = t296 * t157;
t306 = -t291 - t81;
t248 = t294 + t202;
t249 = t116 - t119;
t304 = (t248 * t159 + t249 * t161) * t150;
t128 = t137 * pkin(4);
t120 = rSges(6,1) * t159 + rSges(6,2) * t161;
t239 = qJD(5) * t138;
t227 = t120 * t239;
t256 = t137 * t159;
t234 = rSges(6,2) * t256;
t255 = t137 * t161;
t216 = rSges(6,1) * t255 - t234;
t67 = -rSges(6,3) * t138 + t216;
t186 = -(-pkin(8) * t138 + t128 + t67) * t150 - t227;
t28 = t173 - t186;
t235 = rSges(6,1) * t252;
t68 = -rSges(6,2) * t253 + rSges(6,3) * t137 + t235;
t89 = t138 * pkin(4) + t137 * pkin(8);
t178 = t68 + t89;
t240 = qJD(5) * t137;
t228 = t120 * t240;
t250 = t144 * t157;
t113 = pkin(3) * t250;
t162 = cos(qJ(1));
t149 = cos(t158);
t241 = qJD(1) * t149;
t271 = pkin(1) * qJD(1);
t246 = pkin(2) * t241 + t162 * t271;
t230 = t113 + t246;
t189 = -t228 + t230;
t29 = t178 * t150 + t189;
t78 = t120 * t137;
t79 = t120 * t138;
t301 = -t28 * t78 - t29 * t79;
t300 = 0.2e1 * qJD(5);
t58 = t150 * t68;
t298 = t150 * t89 + t58;
t87 = t138 * rSges(5,1) - rSges(5,2) * t137;
t77 = t150 * t87;
t52 = t230 + t77;
t114 = Icges(6,5) * t159 + Icges(6,6) * t161;
t289 = -Icges(6,3) * t150 + t114 * qJD(5);
t102 = t202 * qJD(5);
t103 = t119 * qJD(5);
t287 = t102 * t159 - t103 * t161 - t114 * t150 + (t116 * t161 + t159 * t294) * qJD(5);
t163 = qJD(1) ^ 2;
t281 = t150 / 0.2e1;
t142 = pkin(2) * t149;
t280 = pkin(3) * t157 ^ 2;
t155 = t162 * pkin(1);
t278 = t137 * t294 + t63;
t277 = t138 * t294 + t64;
t276 = -t116 * t137 + t65;
t275 = -Icges(6,2) * t252 - t108 + t66;
t254 = t138 * t150;
t273 = -rSges(6,3) * t254 - t150 * t234;
t257 = t137 * t150;
t272 = rSges(6,3) * t257 + t150 * t235;
t268 = t159 * t64;
t258 = t114 * t137;
t50 = -t116 * t253 + t252 * t294 + t258;
t266 = t50 * t150;
t265 = pkin(4) * t254 + pkin(8) * t257;
t72 = t114 * t138;
t115 = Icges(6,5) * t161 - Icges(6,6) * t159;
t190 = t115 * t150;
t146 = t163 * t155;
t247 = t163 * t142 + t146;
t243 = t142 + t155;
t242 = qJD(1) * t148;
t238 = qJD(5) * t159;
t237 = qJD(5) * t161;
t233 = t150 * t255;
t232 = t150 * t253;
t231 = t144 * t280 + t247;
t135 = pkin(3) * t143;
t229 = t135 + t297;
t223 = t240 / 0.2e1;
t221 = t239 / 0.2e1;
t219 = -pkin(4) * t257 + pkin(8) * t254;
t70 = rSges(5,1) * t254 - rSges(5,2) * t257;
t82 = rSges(4,1) * t250 - rSges(4,2) * t251;
t136 = pkin(3) * t144;
t215 = t136 + t87;
t93 = rSges(4,1) * t144 - rSges(4,2) * t143;
t60 = t157 * t93 + t246;
t94 = rSges(3,1) * t149 - rSges(3,2) * t148;
t209 = rSges(6,1) * t161 - rSges(6,2) * t159;
t208 = -t137 * t29 + t138 * t28;
t207 = t137 * t67 + t138 * t68;
t41 = t159 * t65 + t161 * t63;
t204 = -t161 * t66 + t268;
t200 = -t116 * t159 + t161 * t294;
t199 = t113 + t70;
t198 = -t228 + t298;
t197 = t295 * t163;
t53 = t65 * t255;
t61 = -Icges(6,3) * t138 + t115 * t137;
t24 = -t138 * t61 - t63 * t256 + t53;
t54 = t66 * t255;
t62 = Icges(6,5) * t252 - Icges(6,6) * t253 + Icges(6,3) * t137;
t25 = t138 * t62 + t64 * t256 - t54;
t194 = (-t137 * t25 - t138 * t24) * qJD(5);
t55 = t63 * t253;
t26 = -t137 * t61 - t65 * t252 + t55;
t27 = t137 * t62 - t138 * t204;
t193 = (-t137 * t27 - t138 * t26) * qJD(5);
t185 = t276 * t159 + t278 * t161;
t184 = t275 * t159 + t277 * t161;
t182 = -t137 * t190 - t138 * t289 + t204 * t150;
t181 = t137 * t289 - t138 * t190 + t206 * t150;
t180 = -rSges(6,2) * t232 + t265 + t272;
t179 = -t115 * qJD(5) + t200 * t150;
t177 = t113 + t180;
t176 = t128 + (-rSges(6,3) - pkin(8)) * t138 + t216;
t175 = t136 + t178;
t174 = -t143 * t280 - t197;
t172 = t135 + t176;
t170 = -rSges(6,1) * t233 + t219 - t273;
t43 = rSges(6,2) * t138 * t237 + (t138 * t238 + t233) * rSges(6,1) + t273;
t44 = -rSges(6,1) * t137 * t238 + (-t137 * t237 - t232) * rSges(6,2) + t272;
t169 = (t150 * t67 - t43) * t138 + (t44 - t58) * t137;
t10 = t193 - t266;
t15 = t204 * qJD(5) + t159 * (t137 * t192 + t138 * t285) + t161 * (t137 * t191 + t138 * t288);
t18 = t179 * t137 + t138 * t287;
t19 = -t137 * t287 + t179 * t138;
t49 = t137 * t200 - t72;
t46 = t49 * t150;
t9 = t46 + t194;
t168 = (t46 + ((t26 + t54 - t55 + (t61 - t268) * t137) * t137 + (-t53 - t27 + (-t204 + t61) * t138 + (t267 + t269) * t137) * t138) * qJD(5)) * t223 + (t200 * qJD(5) + t102 * t161 + t103 * t159) * t150 + (t10 + t266 + ((-t25 + t55 + (t62 - t267) * t138) * t138 + (t24 - t53 + (t62 + t269) * t137) * t137) * qJD(5)) * t221 - (t9 + t15 + t18) * t240 / 0.2e1 - (t19 - t307) * t239 / 0.2e1 + (t138 * t50 + (t41 + t49) * t137) * qJD(5) * t281;
t165 = t170 - t236;
t164 = t301 * qJD(5);
t104 = t209 * qJD(5);
t57 = t157 * t82 + t247;
t56 = -t157 * t81 - t197;
t48 = t150 * t70 + t231;
t47 = -t150 * t303 + t174;
t32 = t207 * qJD(5) + qJD(2);
t21 = t104 * t239 + (t44 - t228 + t265) * t150 + t231;
t20 = -t104 * t240 + (t219 - t43 - t227) * t150 + t174;
t11 = t169 * qJD(5);
t1 = [t168 + m(3) * (-t146 + t163 * (t155 + t94) + (-0.2e1 * rSges(3,1) * t241 + 0.2e1 * rSges(3,2) * t242 + qJD(1) * t94) * qJD(1)) * (-t148 * rSges(3,1) - t149 * rSges(3,2) - t154) + m(4) * (t56 * (t93 + t243) + t57 * (t295 + t296) + (t60 - t82 - t246) * t306) + m(5) * (t47 * (t215 + t243) + t48 * (t229 + t295) + (-t52 + t199 + t246) * t51) + (t20 * (t175 + t243) + t29 * (-pkin(2) * t242 - t160 * t271 + t165) + t21 * (t172 + t295) + t164 + (-t189 + t29 - t298 + t177 + t246) * t28) * m(6); m(6) * t11; t168 + (t172 * t21 + t175 * t20 + t164 + (-t186 + t236 + t165) * t29 + (-t113 - t198 + t177) * t28) * m(6) + (t47 * t215 + t48 * t229 + (-t113 - t77 + t199) * t51) * m(5) + (t56 * t93 + t57 * t296 - t306 * t82 - t60 * t81 - (-t296 * t60 - t306 * t93) * t157) * m(4); t168 + (t176 * t21 + t178 * t20 + t164 + (t170 - t186) * t29 + (t180 - t198) * t28) * m(6) + (-(-t297 * t52 + t51 * t87) * t150 + t47 * t87 + t48 * t297 + t51 * t70 - t52 * t303) * m(5); (t307 * t138 + (t150 * t41 - t15) * t137) * t281 - t150 * ((-t249 * t159 + t248 * t161) * t150 + ((t275 * t137 - t276 * t138) * t161 + (-t277 * t137 + t278 * t138) * t159) * qJD(5)) / 0.2e1 + ((-t239 * t258 - t190) * t138 + (-t304 + (-t184 * t137 + (t72 + t185) * t138) * qJD(5)) * t137) * t221 + ((t72 * t240 - t190) * t137 + (t304 + (-t185 * t138 + (-t258 + t184) * t137) * qJD(5)) * t138) * t223 - (t150 * t18 + ((-t181 * t137 - t150 * t27) * t138 + (-t182 * t137 + t150 * t26) * t137) * t300) * t137 / 0.2e1 - (t150 * t19 + ((-t181 * t138 - t150 * t25) * t138 + (-t182 * t138 + t150 * t24) * t137) * t300) * t138 / 0.2e1 + (t9 + t194) * t257 / 0.2e1 - (t10 + t193) * t254 / 0.2e1 + (t11 * t207 + t32 * t169 + t208 * t104 + ((-t150 * t29 + t21) * t138 + (-t150 * t28 - t20) * t137) * t120 - t301 * t150 - (t32 * (-t137 * t78 - t138 * t79) + t208 * t209) * qJD(5)) * m(6);];
tauc = t1(:);
