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
% Datum: 2020-01-03 11:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:33:38
% EndTime: 2020-01-03 11:33:50
% DurationCPUTime: 4.37s
% Computational Cost: add. (8428->330), mult. (5174->415), div. (0->0), fcn. (3828->10), ass. (0->198)
t164 = qJD(1) + qJD(3);
t165 = qJ(1) + pkin(8);
t159 = qJ(3) + t165;
t151 = cos(t159);
t264 = t151 * t164;
t150 = sin(t159);
t267 = t150 * t164;
t253 = pkin(3) * t264 + qJ(4) * t267;
t167 = cos(pkin(9));
t287 = rSges(5,1) * t167;
t133 = t151 * t287;
t256 = rSges(5,3) * t267 + t164 * t133;
t251 = t150 * rSges(5,3) + t133;
t285 = rSges(5,2) * sin(pkin(9));
t78 = -t151 * t285 + t251;
t271 = qJ(4) * t150;
t293 = pkin(3) * t151;
t99 = t271 + t293;
t92 = t164 * t99;
t320 = -t164 * t78 + t253 + t256 - t92;
t152 = pkin(4) * t167 + pkin(3);
t101 = t152 * t264;
t163 = pkin(9) + qJ(5);
t155 = sin(t163);
t157 = cos(t163);
t112 = rSges(6,1) * t155 + rSges(6,2) * t157;
t241 = qJD(5) * t150;
t243 = qJD(4) * t151;
t194 = -t112 * t241 - t243;
t265 = t151 * t157;
t236 = rSges(6,1) * t265;
t260 = rSges(6,3) * t267 + t164 * t236;
t119 = t151 * t152;
t168 = -pkin(7) - qJ(4);
t62 = t293 - t119 + (qJ(4) + t168) * t150;
t266 = t151 * t155;
t121 = rSges(6,2) * t266;
t75 = rSges(6,3) * t150 - t121 + t236;
t63 = t164 * t75;
t319 = t164 * t62 + t101 - t194 + t260 - t63 - t92;
t145 = Icges(6,4) * t157;
t201 = -Icges(6,2) * t155 + t145;
t190 = t201 * t164;
t276 = Icges(6,4) * t155;
t109 = Icges(6,1) * t157 - t276;
t191 = t109 * t164;
t71 = -Icges(6,5) * t151 + t109 * t150;
t280 = t157 * t71;
t69 = -Icges(6,6) * t151 + t201 * t150;
t282 = t155 * t69;
t205 = t280 - t282;
t308 = Icges(6,1) * t155 + t145;
t298 = -Icges(6,5) * t164 + qJD(5) * t308;
t106 = Icges(6,2) * t157 + t276;
t301 = -Icges(6,6) * t164 + t106 * qJD(5);
t70 = Icges(6,4) * t265 - Icges(6,2) * t266 + Icges(6,6) * t150;
t117 = Icges(6,4) * t266;
t72 = Icges(6,1) * t265 + Icges(6,5) * t150 - t117;
t318 = -t205 * qJD(5) - t155 * (-t150 * t298 + t151 * t191) - t157 * (-t150 * t301 + t151 * t190) + t164 * (t155 * t72 + t157 * t70);
t317 = -rSges(5,3) - qJ(4);
t156 = sin(t165);
t169 = sin(qJ(1));
t161 = t169 * pkin(1);
t307 = pkin(2) * t156 + t161;
t304 = t307 * qJD(1);
t309 = t150 * rSges(4,1) + t151 * rSges(4,2);
t88 = t309 * t164;
t316 = -t304 - t88;
t258 = t308 + t201;
t259 = t106 - t109;
t314 = (t258 * t155 + t259 * t157) * t164;
t312 = 0.2e1 * qJD(5);
t132 = t150 * t287;
t144 = t150 * pkin(3);
t198 = -t150 * t285 + t132 + t144;
t310 = (t317 * t151 + t198) * t164;
t255 = t150 * t152 + t151 * t168;
t104 = Icges(6,5) * t155 + Icges(6,6) * t157;
t302 = -Icges(6,3) * t164 + t104 * qJD(5);
t94 = t201 * qJD(5);
t95 = t109 * qJD(5);
t300 = (t106 * t157 + t155 * t308) * qJD(5) - t104 * t164 + t155 * t94 - t157 * t95;
t171 = qJD(1) ^ 2;
t294 = t164 / 0.2e1;
t158 = cos(t165);
t149 = pkin(2) * t158;
t170 = cos(qJ(1));
t162 = t170 * pkin(1);
t292 = t150 * t308 + t69;
t291 = t151 * t308 + t70;
t290 = -t106 * t150 + t71;
t289 = -Icges(6,2) * t265 - t117 + t72;
t286 = rSges(6,1) * t157;
t284 = rSges(6,2) * t155;
t283 = pkin(1) * qJD(1);
t281 = t155 * t70;
t244 = qJD(1) * t158;
t249 = pkin(2) * t244 + t170 * t283;
t26 = (-t62 + t75 + t99) * t164 + t194 + t249;
t279 = t164 * t26;
t270 = t104 * t150;
t46 = -t106 * t266 + t265 * t308 + t270;
t278 = t46 * t164;
t80 = t104 * t151;
t105 = Icges(6,5) * t157 - Icges(6,6) * t155;
t189 = t105 * t164;
t269 = t150 * t155;
t268 = t150 * t157;
t263 = t155 * t164;
t262 = t164 * t168;
t234 = rSges(6,2) * t263;
t261 = rSges(6,3) * t264 + t150 * t234;
t235 = t164 * t285;
t257 = rSges(5,3) * t264 + t150 * t235;
t254 = t119 - t121;
t125 = qJ(4) * t264;
t137 = qJD(4) * t150;
t252 = t125 + t137;
t153 = t171 * t162;
t250 = t171 * t149 + t153;
t246 = t149 + t162;
t245 = qJD(1) * t156;
t242 = qJD(5) * t112;
t240 = qJD(5) * t151;
t239 = qJD(5) * t155;
t238 = qJD(5) * t157;
t120 = rSges(6,1) * t268;
t74 = -rSges(6,2) * t269 - rSges(6,3) * t151 + t120;
t237 = t255 + t74;
t233 = t164 * (-t243 + t253) + t250;
t232 = rSges(6,1) * t239;
t231 = rSges(6,2) * t238;
t229 = t137 + t261;
t228 = t120 + t255;
t227 = t112 * t240;
t223 = t241 / 0.2e1;
t221 = t240 / 0.2e1;
t53 = t151 * t231 + (t151 * t239 + t157 * t267) * rSges(6,1) - t261;
t219 = t164 * t74 - t53;
t54 = -t150 * t232 + (-t150 * t238 - t151 * t263) * rSges(6,2) + t260;
t218 = t54 - t63;
t89 = rSges(4,1) * t264 - rSges(4,2) * t267;
t212 = -qJD(4) - t235;
t210 = t251 + t271;
t100 = rSges(4,1) * t151 - rSges(4,2) * t150;
t66 = t100 * t164 + t249;
t209 = -t243 + t249;
t114 = rSges(3,1) * t158 - rSges(3,2) * t156;
t206 = -t284 + t286;
t40 = t155 * t71 + t157 * t69;
t203 = -t157 * t72 + t281;
t199 = -t106 * t155 + t157 * t308;
t197 = t307 * t171;
t55 = t71 * t268;
t67 = -Icges(6,3) * t151 + t105 * t150;
t27 = -t151 * t67 - t69 * t269 + t55;
t56 = t72 * t268;
t68 = Icges(6,5) * t265 - Icges(6,6) * t266 + Icges(6,3) * t150;
t28 = t151 * t68 + t70 * t269 - t56;
t193 = (-t150 * t28 - t151 * t27) * qJD(5);
t57 = t69 * t266;
t29 = -t150 * t67 - t71 * t265 + t57;
t30 = t150 * t68 - t203 * t151;
t192 = (-t150 * t30 - t151 * t29) * qJD(5);
t185 = t164 * t137 - t197;
t184 = t290 * t155 + t292 * t157;
t183 = t289 * t155 + t291 * t157;
t182 = -t150 * t189 - t151 * t302 + t203 * t164;
t181 = t150 * t302 - t151 * t189 + t205 * t164;
t180 = t304 - t137;
t179 = -t105 * qJD(5) + t199 * t164;
t177 = (-pkin(3) - t287) * t267 + t252 + t257;
t10 = t192 - t278;
t17 = t203 * qJD(5) + t155 * (t150 * t191 + t151 * t298) + t157 * (t150 * t190 + t151 * t301);
t20 = t179 * t150 + t151 * t300;
t21 = -t150 * t300 + t179 * t151;
t45 = t199 * t150 - t80;
t42 = t45 * t164;
t9 = t42 + t193;
t176 = (t42 + ((t29 + t56 - t57 + (t67 - t281) * t150) * t150 + (-t55 - t30 + (-t203 + t67) * t151 + (t280 + t282) * t150) * t151) * qJD(5)) * t223 + (t199 * qJD(5) + t155 * t95 + t157 * t94) * t164 + (t278 + ((-t28 + t57 + (t68 - t280) * t151) * t151 + (t27 - t55 + (t68 + t282) * t150) * t150) * qJD(5) + t10) * t221 - (t17 + t20 + t9) * t241 / 0.2e1 - (t21 - t318) * t240 / 0.2e1 + (t151 * t46 + (t40 + t45) * t150) * qJD(5) * t294;
t76 = pkin(3) * t267 - t252;
t31 = (-t164 * t132 + t257 - t76) * t164 + t185;
t32 = (t212 * t151 + t256) * t164 + t233;
t43 = t180 + t310;
t173 = (t31 * (pkin(3) - t285) + t32 * t317 + t43 * t212) * t151;
t96 = t206 * qJD(5);
t14 = -t96 * t241 + (-t227 - t125 - t53 - t76 + (t144 - t255) * t164) * t164 + t185;
t15 = t96 * t240 + (-t243 + t101 + t54 + (-t242 - t262) * t150 - t253) * t164 + t233;
t25 = t237 * t164 + t180 + t227;
t172 = (t14 * (rSges(6,3) - t168) - t15 * t284 - t25 * t242 + (t26 * (-t152 - t286) - t25 * t168) * t164) * t150 + (t14 * t286 + t26 * (-t231 - t232 - t262) - t15 * rSges(6,3) + t25 * (-qJD(4) - t234)) * t151;
t86 = t112 * t151;
t85 = t112 * t150;
t60 = t164 * t89 + t250;
t59 = -t164 * t88 - t197;
t44 = (t78 + t99) * t164 + t209;
t39 = qJD(2) + (t150 * t74 + t151 * t75) * qJD(5);
t11 = (t218 * t150 + t219 * t151) * qJD(5);
t1 = [m(3) * (-t153 + t171 * (t114 + t162) + (-0.2e1 * rSges(3,1) * t244 + 0.2e1 * rSges(3,2) * t245 + qJD(1) * t114) * qJD(1)) * (-t156 * rSges(3,1) - t158 * rSges(3,2) - t161) + m(4) * (t59 * (t100 + t246) + t60 * (t307 + t309) + (t66 - t89 - t249) * t316) + t176 + (t14 * (t246 + t254) + t26 * (-pkin(2) * t245 - t169 * t283 + t229) + t15 * (t228 + t307) + t172 + (t26 + t319) * t25) * m(6) + (t31 * (t210 + t246) + t44 * (-t304 + t177) + t32 * (t198 + t307) + t173 + (t249 - t209 + t44 + t320) * t43) * m(5); m(6) * t11; t176 + (t14 * t254 + t15 * t228 + t237 * t279 + t172 + (t229 - t137 + t227) * t26 + t319 * t25) * m(6) + (t32 * t198 + t31 * t210 + t173 + (-t137 + t177 + t310) * t44 + (t243 + t320) * t43) * m(5) + (t100 * t59 + t60 * t309 - t316 * t89 - t66 * t88 - (-t100 * t316 - t309 * t66) * t164) * m(4); m(5) * (-t150 * t32 - t31 * t151) + m(6) * (-t14 * t151 - t15 * t150); (t318 * t151 + (t164 * t40 - t17) * t150) * t294 - t164 * ((-t259 * t155 + t258 * t157) * t164 + ((t289 * t150 - t290 * t151) * t157 + (-t291 * t150 + t292 * t151) * t155) * qJD(5)) / 0.2e1 + ((-t240 * t270 - t189) * t151 + (-t314 + (-t183 * t150 + (t80 + t184) * t151) * qJD(5)) * t150) * t221 + ((t80 * t241 - t189) * t150 + (t314 + (-t184 * t151 + (-t270 + t183) * t150) * qJD(5)) * t151) * t223 - (t164 * t20 + ((-t181 * t150 - t164 * t30) * t151 + (-t182 * t150 + t164 * t29) * t150) * t312) * t150 / 0.2e1 - (t164 * t21 + ((-t181 * t151 - t164 * t28) * t151 + (-t182 * t151 + t164 * t27) * t150) * t312) * t151 / 0.2e1 + (t9 + t193) * t267 / 0.2e1 - (t10 + t192) * t264 / 0.2e1 + ((t11 * t75 + t39 * t219 + t25 * t96) * t151 + (t11 * t74 + t39 * t218 - t26 * t96) * t150 + ((t15 - t279) * t151 + (-t164 * t25 - t14) * t150) * t112 - (-t25 * t85 - t26 * t86) * t164 - (t39 * (-t150 * t85 - t151 * t86) + (-t150 * t26 + t151 * t25) * t206) * qJD(5)) * m(6);];
tauc = t1(:);
