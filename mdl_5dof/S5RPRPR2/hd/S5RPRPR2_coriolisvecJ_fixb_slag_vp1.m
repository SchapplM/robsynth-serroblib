% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2022-01-23 09:19
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPR2_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR2_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR2_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR2_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR2_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:18:48
% EndTime: 2022-01-23 09:18:54
% DurationCPUTime: 3.51s
% Computational Cost: add. (8428->310), mult. (5174->398), div. (0->0), fcn. (3828->10), ass. (0->199)
t171 = qJD(1) ^ 2;
t305 = 2 * qJD(5);
t164 = qJD(1) + qJD(3);
t165 = qJ(1) + pkin(8);
t161 = qJ(3) + t165;
t155 = cos(t161);
t141 = t155 * qJ(4);
t154 = sin(t161);
t167 = cos(pkin(9));
t156 = pkin(4) * t167 + pkin(3);
t168 = -pkin(7) - qJ(4);
t257 = t155 * t168;
t65 = -t257 - t141 + (pkin(3) - t156) * t154;
t163 = pkin(9) + qJ(5);
t159 = cos(t163);
t263 = t154 * t159;
t242 = rSges(6,1) * t263;
t157 = sin(t163);
t264 = t154 * t157;
t125 = rSges(6,2) * t264;
t252 = t155 * rSges(6,3) + t125;
t78 = t242 - t252;
t67 = t164 * t78;
t304 = t164 * t65 - t67;
t124 = t155 * t156;
t261 = t154 * t168;
t259 = t155 * t159;
t207 = rSges(6,1) * t259 + t154 * rSges(6,3);
t260 = t155 * t157;
t240 = rSges(6,2) * t260;
t79 = -t240 + t207;
t303 = t124 - t261 + t79;
t216 = t155 * pkin(3) + t154 * qJ(4);
t279 = rSges(5,2) * sin(pkin(9));
t241 = t155 * t279;
t282 = rSges(5,1) * t167;
t301 = -t154 * rSges(5,3) - t155 * t282;
t194 = -t241 + t216 - t301;
t302 = t164 * t194;
t151 = Icges(6,4) * t159;
t210 = -Icges(6,2) * t157 + t151;
t111 = Icges(6,1) * t157 + t151;
t300 = pkin(2) * cos(t165) + cos(qJ(1)) * pkin(1);
t108 = Icges(6,5) * t159 - Icges(6,6) * t157;
t269 = Icges(6,4) * t157;
t109 = Icges(6,2) * t159 + t269;
t209 = t109 * t157 - t111 * t159;
t298 = t108 * qJD(5) + t164 * t209;
t107 = Icges(6,5) * t157 + Icges(6,6) * t159;
t183 = Icges(6,3) * t164 - qJD(5) * t107;
t197 = t210 * t155;
t74 = Icges(6,6) * t154 + t197;
t275 = t157 * t74;
t112 = Icges(6,1) * t159 - t269;
t198 = t112 * t155;
t76 = Icges(6,5) * t154 + t198;
t212 = -t159 * t76 + t275;
t262 = t154 * t164;
t297 = -t108 * t262 + t155 * t183 + t164 * t212;
t196 = t108 * t155;
t73 = Icges(6,4) * t263 - Icges(6,2) * t264 - Icges(6,6) * t155;
t276 = t157 * t73;
t123 = Icges(6,4) * t264;
t75 = Icges(6,1) * t263 - Icges(6,5) * t155 - t123;
t213 = -t159 * t75 + t276;
t296 = t154 * t183 + (t196 + t213) * t164;
t71 = Icges(6,5) * t263 - Icges(6,6) * t264 - Icges(6,3) * t155;
t28 = -t154 * t213 - t155 * t71;
t284 = -Icges(6,2) * t263 - t123 + t75;
t286 = t111 * t154 + t73;
t294 = -t157 * t284 - t159 * t286;
t293 = t154 / 0.2e1;
t292 = -t155 / 0.2e1;
t291 = t164 / 0.2e1;
t289 = pkin(3) * t154;
t288 = -t154 * t71 - t75 * t259;
t72 = Icges(6,3) * t154 + t196;
t287 = t154 * t72 + t76 * t259;
t285 = -t111 * t155 - t74;
t283 = -t109 * t155 + t76;
t281 = rSges(6,1) * t159;
t278 = rSges(6,2) * t159;
t149 = t155 * rSges(4,1);
t101 = -t141 + t289;
t138 = qJD(4) * t154;
t221 = -pkin(2) * sin(t165) - sin(qJ(1)) * pkin(1);
t204 = t221 * qJD(1);
t191 = t138 + t204;
t117 = rSges(6,1) * t157 + t278;
t246 = qJD(5) * t155;
t237 = t117 * t246;
t179 = t191 - t237;
t26 = (-t101 + t65 - t78) * t164 + t179;
t274 = t164 * t26;
t266 = t107 * t155;
t47 = -t154 * t209 - t266;
t273 = t47 * t164;
t102 = rSges(4,1) * t154 + rSges(4,2) * t155;
t104 = -rSges(4,2) * t154 + t149;
t203 = t300 * qJD(1);
t70 = t104 * t164 + t203;
t272 = t70 * t102;
t96 = t164 * t101;
t271 = t138 - t96;
t267 = t107 * t154;
t265 = t108 * t164;
t258 = t155 * t164;
t256 = t164 * t102;
t255 = -t109 + t112;
t254 = t111 + t210;
t134 = t154 * t279;
t253 = rSges(5,3) * t258 + t164 * t134;
t130 = qJ(4) * t258;
t251 = t130 + t138;
t250 = t155 * rSges(5,3) + t134;
t248 = qJD(4) * t164;
t247 = qJD(5) * t154;
t245 = qJD(5) * t157;
t243 = t154 * t282;
t239 = qJD(5) * t278;
t238 = -t164 * t240 + (-rSges(6,1) * t245 - t239) * t154;
t95 = t117 * t247;
t236 = t155 * t245;
t233 = -pkin(3) - t282;
t232 = -t247 / 0.2e1;
t229 = t246 / 0.2e1;
t57 = t76 * t263;
t227 = t155 * t72 - t57;
t226 = -t71 + t275;
t217 = -rSges(6,2) * t157 + t281;
t139 = qJD(4) * t155;
t190 = t203 - t139;
t27 = t164 * t303 + t190 - t95;
t215 = -t154 * t27 - t155 * t26;
t214 = t154 * t78 + t155 * t79;
t42 = t157 * t75 + t159 * t73;
t43 = t157 * t76 + t159 * t74;
t208 = -t154 * t156 - t257;
t206 = t221 * t171;
t205 = t300 * t171;
t29 = -t264 * t74 - t227;
t201 = (t154 * t29 - t155 * t28) * qJD(5);
t30 = -t260 * t73 - t288;
t31 = -t260 * t74 + t287;
t200 = (t154 * t31 - t155 * t30) * qJD(5);
t199 = rSges(6,3) * t258 + t164 * t125 - t155 * t239;
t193 = t155 * t248 - t205;
t192 = -t157 * t283 + t159 * t285;
t189 = t154 * t233 + t141 + t250;
t188 = t154 * t248 + t164 * (-pkin(3) * t262 + t251) + t206;
t69 = t204 - t256;
t186 = (-t157 * t254 + t159 * t255) * t164;
t185 = Icges(6,5) * t164 - qJD(5) * t111;
t184 = Icges(6,6) * t164 - qJD(5) * t109;
t182 = -t257 + (-t156 - t281) * t154 + t252;
t55 = (-t159 * t262 - t236) * rSges(6,1) + t199;
t56 = t164 * t207 + t238;
t181 = (t55 + t67) * t155 + (-t164 * t79 + t56) * t154;
t48 = -t155 * t209 + t267;
t44 = t48 * t164;
t10 = t44 + t200;
t52 = t154 * t184 + t164 * t197;
t54 = t154 * t185 + t164 * t198;
t16 = -qJD(5) * t213 + t157 * t54 + t159 * t52;
t51 = t155 * t184 - t210 * t262;
t53 = -t112 * t262 + t155 * t185;
t17 = -qJD(5) * t212 + t157 * t53 + t159 * t51;
t98 = t210 * qJD(5);
t99 = t112 * qJD(5);
t174 = t107 * t164 - t157 * t98 + t159 * t99 + (-t109 * t159 - t111 * t157) * qJD(5);
t20 = t154 * t298 + t174 * t155;
t21 = t174 * t154 - t155 * t298;
t9 = t201 + t273;
t180 = (t44 + ((t29 - t57 + (t72 + t276) * t155 + t288) * t155 + t287 * t154) * qJD(5)) * t229 + (-qJD(5) * t209 + t157 * t99 + t159 * t98) * t164 + (-t273 + ((t155 * t226 - t287 + t31) * t155 + (t154 * t226 + t227 + t30) * t154) * qJD(5) + t9) * t232 + (t17 + t20) * t247 / 0.2e1 - (t16 + t21 + t10) * t246 / 0.2e1 + ((t42 + t47) * t154 + (t43 + t48) * t155) * qJD(5) * t291;
t176 = -qJD(5) * t43 - t157 * t51 + t159 * t53 + t164 * t72;
t175 = -qJD(5) * t42 - t157 * t52 + t159 * t54 + t164 * t71;
t114 = t164 * t241;
t81 = t243 - t250;
t45 = (-t101 - t81) * t164 + t191;
t46 = t190 + t302;
t173 = t45 * (t114 + t139) + t46 * (t251 + t253) + (t45 * t233 * t155 + (t45 * (-rSges(5,3) - qJ(4)) + t46 * t233) * t154) * t164;
t120 = t164 * t261;
t172 = t26 * (t120 + t139 - t238) + t27 * (-rSges(6,1) * t236 + t138 + t199) + (t26 * (-t207 - t124) + t27 * (t208 - t242)) * t164;
t131 = rSges(4,2) * t262;
t100 = t217 * qJD(5);
t92 = rSges(4,1) * t258 - t131;
t90 = t117 * t155;
t89 = t117 * t154;
t80 = t164 * t216 - t139;
t77 = t164 * t81;
t62 = -t164 * t92 - t205;
t61 = -t164 * t256 + t206;
t41 = qJD(5) * t214 + qJD(2);
t33 = (t164 * t301 + t114 - t80) * t164 + t193;
t32 = t164 * (-t164 * t243 + t253) + t188;
t15 = -t100 * t246 + (t95 + t120 - t56 - t80 + (t216 - t124) * t164) * t164 + t193;
t14 = -t100 * t247 + (-t237 - t130 + t55 + (t208 + t289) * t164) * t164 + t188;
t11 = t181 * qJD(5);
t1 = [t180 + m(4) * (t62 * (-t102 + t221) + t69 * t131 + t61 * (t104 + t300) + (-t149 * t69 - t272) * t164 + (t221 * t70 - t300 * t69) * qJD(1)) + (-(-t26 - t96 + t179 + t304) * t27 + t15 * (t182 + t221) + t14 * (t303 + t300) + (t221 * t27 - t26 * t300) * qJD(1) + t172) * m(6) + (-(-t45 - t77 - t96 + t191) * t46 + t33 * (t189 + t221) + t32 * (t194 + t300) + (t221 * t46 - t300 * t45) * qJD(1) + t173) * m(5); m(6) * t11; t180 + (t15 * t182 + t172 - t26 * (t139 + t95) - t27 * (-t237 + t271 + t304) + (t14 + t274) * t303) * m(6) + (t33 * t189 + t194 * t32 + t173 - t45 * (t139 - t302) - t46 * (-t77 + t271)) * m(5) + (-(-t104 * t69 - t272) * t164 - t102 * t62 + t104 * t61 - t69 * t92 - t70 * t256) * m(4); 0.2e1 * (t14 * t292 + t15 * t293) * m(6) + 0.2e1 * (t292 * t32 + t293 * t33) * m(5); ((t164 * t43 - t16) * t155 + (t164 * t42 + t17) * t154) * t291 + ((-t247 * t266 + t265) * t154 + (t186 + (-t294 * t155 + (t267 + t192) * t154) * qJD(5)) * t155) * t232 + ((-t246 * t267 - t265) * t155 + (t186 + (t192 * t154 + (-t294 + t266) * t155) * qJD(5)) * t154) * t229 - t164 * ((t255 * t157 + t254 * t159) * t164 + ((t154 * t283 - t155 * t284) * t159 + (t154 * t285 + t155 * t286) * t157) * qJD(5)) / 0.2e1 + (t164 * t20 + ((-t154 * t296 - t175 * t155 + t164 * t31) * t155 + (t154 * t297 + t176 * t155 + t164 * t30) * t154) * t305) * t293 + (t164 * t21 + ((-t175 * t154 + t155 * t296 + t164 * t29) * t155 + (t176 * t154 - t155 * t297 + t164 * t28) * t154) * t305) * t292 + (t9 + t201) * t262 / 0.2e1 + (t10 + t200) * t258 / 0.2e1 + (t11 * t214 + t41 * t181 + t215 * t100 + ((-t164 * t27 - t15) * t155 + (-t14 + t274) * t154) * t117 - (t26 * t89 - t27 * t90) * t164 - (t41 * (-t154 * t89 - t155 * t90) + t215 * t217) * qJD(5)) * m(6);];
tauc = t1(:);
