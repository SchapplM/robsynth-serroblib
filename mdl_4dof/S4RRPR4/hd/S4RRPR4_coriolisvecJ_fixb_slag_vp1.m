% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRPR4_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR4_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR4_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR4_coriolisvecJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR4_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR4_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR4_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:02:26
% EndTime: 2019-12-31 17:02:30
% DurationCPUTime: 2.72s
% Computational Cost: add. (6057->295), mult. (4880->375), div. (0->0), fcn. (3686->8), ass. (0->192)
t288 = 2 * qJD(4);
t160 = sin(qJ(1));
t263 = pkin(1) * qJD(1);
t226 = t160 * t263;
t156 = qJ(1) + qJ(2);
t151 = sin(t156);
t152 = cos(t156);
t109 = rSges(3,1) * t151 + rSges(3,2) * t152;
t155 = qJD(1) + qJD(2);
t251 = t109 * t155;
t80 = -t226 - t251;
t158 = cos(pkin(7));
t147 = pkin(3) * t158 + pkin(2);
t121 = t152 * t147;
t159 = -pkin(6) - qJ(3);
t245 = t151 * t159;
t154 = pkin(7) + qJ(4);
t150 = cos(t154);
t247 = t150 * t152;
t190 = rSges(5,1) * t247 + t151 * rSges(5,3);
t149 = sin(t154);
t249 = t149 * t152;
t227 = rSges(5,2) * t249;
t76 = -t227 + t190;
t287 = t121 - t245 + t76;
t198 = t152 * pkin(2) + t151 * qJ(3);
t265 = rSges(4,2) * sin(pkin(7));
t228 = t152 * t265;
t268 = rSges(4,1) * t158;
t285 = -t151 * rSges(4,3) - t152 * t268;
t182 = -t228 + t198 - t285;
t286 = t182 * t155;
t134 = Icges(5,4) * t150;
t194 = -Icges(5,2) * t149 + t134;
t104 = Icges(5,1) * t149 + t134;
t132 = t151 * t265;
t244 = t152 * t155;
t240 = rSges(4,3) * t244 + t155 * t132;
t230 = t151 * t268;
t237 = t152 * rSges(4,3) + t132;
t78 = t230 - t237;
t138 = t152 * qJ(3);
t273 = pkin(2) * t151;
t108 = -t138 + t273;
t93 = t155 * t108;
t284 = t155 * t78 + t240 + t93;
t101 = Icges(5,5) * t150 - Icges(5,6) * t149;
t256 = Icges(5,4) * t149;
t102 = Icges(5,2) * t150 + t256;
t193 = t149 * t102 - t150 * t104;
t283 = t101 * qJD(4) + t193 * t155;
t100 = Icges(5,5) * t149 + Icges(5,6) * t150;
t174 = Icges(5,3) * t155 - t100 * qJD(4);
t184 = t194 * t152;
t71 = Icges(5,6) * t151 + t184;
t261 = t149 * t71;
t105 = Icges(5,1) * t150 - t256;
t185 = t105 * t152;
t73 = Icges(5,5) * t151 + t185;
t196 = -t150 * t73 + t261;
t246 = t151 * t155;
t282 = -t101 * t246 + t174 * t152 + t196 * t155;
t183 = t101 * t152;
t248 = t150 * t151;
t250 = t149 * t151;
t70 = Icges(5,4) * t248 - Icges(5,2) * t250 - Icges(5,6) * t152;
t262 = t149 * t70;
t120 = Icges(5,4) * t250;
t72 = Icges(5,1) * t248 - Icges(5,5) * t152 - t120;
t197 = -t150 * t72 + t262;
t281 = t174 * t151 + (t183 + t197) * t155;
t135 = qJD(3) * t151;
t122 = rSges(5,2) * t250;
t264 = rSges(5,2) * t150;
t223 = qJD(4) * t264;
t189 = rSges(5,3) * t244 + t155 * t122 - t152 * t223;
t234 = qJD(4) * t152;
t222 = t149 * t234;
t229 = rSges(5,1) * t248;
t239 = t152 * rSges(5,3) + t122;
t75 = t229 - t239;
t62 = t155 * t75;
t243 = t152 * t159;
t63 = -t243 - t138 + (pkin(2) - t147) * t151;
t280 = -rSges(5,1) * t222 - t155 * t63 + t135 + t189 + t62 + t93;
t68 = Icges(5,5) * t248 - Icges(5,6) * t250 - Icges(5,3) * t152;
t27 = -t197 * t151 - t152 * t68;
t267 = rSges(5,1) * t149;
t106 = t264 + t267;
t221 = t106 * t234;
t202 = t135 - t221;
t180 = t202 - t226;
t25 = (-t108 + t63 - t75) * t155 + t180;
t136 = qJD(3) * t152;
t161 = cos(qJ(1));
t225 = t161 * t263;
t203 = -t136 + t225;
t235 = qJD(4) * t151;
t90 = t106 * t235;
t26 = t155 * t287 + t203 - t90;
t279 = -t151 * t26 - t152 * t25;
t278 = (-t102 * t152 + t73) * t151 - (-Icges(5,2) * t248 - t120 + t72) * t152;
t277 = t151 / 0.2e1;
t276 = -t152 / 0.2e1;
t275 = t155 / 0.2e1;
t274 = pkin(1) * t160;
t153 = t161 * pkin(1);
t272 = -t151 * t68 - t72 * t247;
t69 = Icges(5,3) * t151 + t183;
t271 = t151 * t69 + t73 * t247;
t266 = rSges(5,1) * t150;
t259 = t155 * t25;
t253 = t100 * t152;
t42 = -t193 * t151 - t253;
t258 = t42 * t155;
t254 = t100 * t151;
t252 = t101 * t155;
t242 = -t102 + t105;
t241 = t104 + t194;
t126 = qJ(3) * t244;
t238 = t126 + t135;
t236 = qJD(3) * t155;
t162 = qJD(1) ^ 2;
t233 = t162 * t274;
t232 = t162 * t153;
t224 = -t151 * t223 - t155 * t227 - t235 * t267;
t218 = -pkin(2) - t268;
t217 = -t235 / 0.2e1;
t214 = t234 / 0.2e1;
t212 = -t68 + t261;
t56 = t73 * t248;
t211 = t152 * t69 - t56;
t111 = t152 * rSges(3,1) - rSges(3,2) * t151;
t205 = t152 * t236 - t232;
t92 = rSges(3,1) * t244 - rSges(3,2) * t246;
t117 = t155 * t245;
t201 = t117 + t136 - t224;
t199 = -rSges(5,2) * t149 + t266;
t40 = t149 * t72 + t150 * t70;
t41 = t149 * t73 + t150 * t71;
t192 = -t147 * t151 - t243;
t191 = t151 * t236 + t155 * (-pkin(2) * t246 + t238) - t233;
t28 = -t71 * t250 - t211;
t187 = (t151 * t28 - t152 * t27) * qJD(4);
t29 = -t70 * t249 - t272;
t30 = -t71 * t249 + t271;
t186 = (t151 * t30 - t152 * t29) * qJD(4);
t39 = (t151 * t75 + t152 * t76) * qJD(4);
t181 = -t71 * t151 + t70 * t152;
t179 = t218 * t151 + t138 + t237;
t177 = (-t241 * t149 + t242 * t150) * t155;
t176 = Icges(5,5) * t155 - qJD(4) * t104;
t175 = Icges(5,6) * t155 - t102 * qJD(4);
t173 = -t243 + (-t147 - t266) * t151 + t239;
t43 = -t193 * t152 + t254;
t38 = t43 * t155;
t10 = t38 + t186;
t49 = t175 * t151 + t155 * t184;
t51 = t176 * t151 + t155 * t185;
t13 = -t197 * qJD(4) + t149 * t51 + t150 * t49;
t48 = t175 * t152 - t194 * t246;
t50 = -t105 * t246 + t176 * t152;
t14 = -t196 * qJD(4) + t149 * t50 + t150 * t48;
t95 = t194 * qJD(4);
t96 = t105 * qJD(4);
t166 = t100 * t155 - t149 * t95 + t150 * t96 + (-t102 * t150 - t104 * t149) * qJD(4);
t19 = t283 * t151 + t166 * t152;
t20 = t166 * t151 - t283 * t152;
t9 = t187 + t258;
t172 = (t38 + ((t28 - t56 + (t69 + t262) * t152 + t272) * t152 + t271 * t151) * qJD(4)) * t214 + (-t193 * qJD(4) + t149 * t96 + t150 * t95) * t155 + (-t258 + ((t212 * t152 - t271 + t30) * t152 + (t212 * t151 + t211 + t29) * t151) * qJD(4) + t9) * t217 + (t14 + t19) * t235 / 0.2e1 - (t13 + t20 + t10) * t234 / 0.2e1 + ((t40 + t42) * t151 + (t41 + t43) * t152) * qJD(4) * t275;
t168 = -t41 * qJD(4) - t149 * t48 + t150 * t50 + t155 * t69;
t167 = -t40 * qJD(4) - t149 * t49 + t150 * t51 + t155 * t68;
t54 = (-t108 - t78) * t155 + t135 - t226;
t55 = t203 + t286;
t165 = (t54 * t218 * t152 + (t54 * (-rSges(4,3) - qJ(3)) + t55 * t218) * t151) * t155;
t164 = -t278 * t149 + t181 * t150;
t163 = (t25 * (-t190 - t121) + t26 * (t192 - t229)) * t155;
t113 = t155 * t228;
t97 = t199 * qJD(4);
t89 = t106 * t152;
t88 = t106 * t151;
t81 = t111 * t155 + t225;
t77 = t198 * t155 - t136;
t67 = -t155 * t92 - t232;
t66 = -t155 * t251 - t233;
t53 = t190 * t155 + t224;
t52 = (-t150 * t246 - t222) * rSges(5,1) + t189;
t33 = (t285 * t155 + t113 - t77) * t155 + t205;
t32 = t155 * (-t155 * t230 + t240) + t191;
t18 = -t97 * t234 + (t90 + t117 - t53 - t77 + (t198 - t121) * t155) * t155 + t205;
t17 = -t97 * t235 + (-t221 - t126 + t52 + (t192 + t273) * t155) * t155 + t191;
t1 = [m(3) * (t67 * (-t109 - t274) + t66 * (t111 + t153) + (-t92 - t225 + t81) * t80) + t172 + (t18 * (t173 - t274) + t25 * (t201 - t225) + t17 * (t153 + t287) + t163 + (-t180 + t25 - t226 + t280) * t26) * m(5) + (t33 * (t179 - t274) + t54 * (t113 - t203) + t32 * (t153 + t182) + t165 + (t54 + t126 + t284) * t55) * m(4); t172 + (t18 * t173 + t163 + (-t202 + t280) * t26 + (-t136 - t90 + t201) * t25 + (t17 + t259) * t287) * m(5) + (t33 * t179 + t32 * t182 + t165 + (-t135 + t238 + t284) * t55 + (t113 + t286) * t54) * m(4) + (-(-t109 * t81 - t111 * t80) * t155 - t109 * t67 + t111 * t66 - t80 * t92 - t81 * t251) * m(3); 0.2e1 * (t17 * t276 + t18 * t277) * m(5) + 0.2e1 * (t32 * t276 + t33 * t277) * m(4); ((t155 * t41 - t13) * t152 + (t155 * t40 + t14) * t151) * t275 + ((-t235 * t253 + t252) * t151 + (t177 + (t151 * t254 + t164) * qJD(4)) * t152) * t217 + ((-t234 * t254 - t252) * t152 + (t177 + (t152 * t253 + t164) * qJD(4)) * t151) * t214 - t155 * ((t242 * t149 + t241 * t150) * t155 + (t181 * t149 + t278 * t150) * qJD(4)) / 0.2e1 + (t155 * t19 + ((-t281 * t151 - t167 * t152 + t155 * t30) * t152 + (t282 * t151 + t168 * t152 + t155 * t29) * t151) * t288) * t277 + (t155 * t20 + ((-t167 * t151 + t281 * t152 + t155 * t28) * t152 + (t168 * t151 - t282 * t152 + t155 * t27) * t151) * t288) * t276 + (t9 + t187) * t246 / 0.2e1 + (t10 + t186) * t244 / 0.2e1 + (((-t155 * t26 - t18) * t152 + (-t17 + t259) * t151) * t106 + t279 * t97 + 0.2e1 * ((-t155 * t76 + t53) * t151 + (t52 + t62) * t152) * t39 - (t25 * t88 - t26 * t89) * t155 - (t39 * (-t151 * t88 - t152 * t89) + t279 * t199) * qJD(4)) * m(5);];
tauc = t1(:);
