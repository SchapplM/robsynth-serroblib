% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:16
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRPR1_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR1_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR1_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR1_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR1_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPR1_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPR1_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:15:34
% EndTime: 2019-12-05 16:15:40
% DurationCPUTime: 2.77s
% Computational Cost: add. (8320->306), mult. (4934->389), div. (0->0), fcn. (3718->8), ass. (0->198)
t159 = pkin(8) + qJ(2);
t154 = sin(t159);
t265 = pkin(2) * qJD(2);
t227 = t154 * t265;
t157 = qJ(3) + t159;
t150 = sin(t157);
t151 = cos(t157);
t101 = rSges(4,1) * t150 + rSges(4,2) * t151;
t160 = qJD(2) + qJD(3);
t256 = t101 * t160;
t81 = -t227 - t256;
t291 = 2 * qJD(5);
t162 = cos(pkin(9));
t152 = pkin(4) * t162 + pkin(3);
t124 = t151 * t152;
t163 = -pkin(7) - qJ(4);
t249 = t150 * t163;
t158 = pkin(9) + qJ(5);
t155 = cos(t158);
t247 = t151 * t155;
t191 = rSges(6,1) * t247 + t150 * rSges(6,3);
t153 = sin(t158);
t248 = t151 * t153;
t228 = rSges(6,2) * t248;
t77 = -t228 + t191;
t290 = t124 - t249 + t77;
t200 = t151 * pkin(3) + t150 * qJ(4);
t267 = rSges(5,2) * sin(pkin(9));
t229 = t151 * t267;
t269 = rSges(5,1) * t162;
t288 = -t150 * rSges(5,3) - t151 * t269;
t183 = -t229 + t200 - t288;
t289 = t160 * t183;
t148 = Icges(6,4) * t155;
t195 = -Icges(6,2) * t153 + t148;
t110 = Icges(6,1) * t153 + t148;
t133 = t150 * t267;
t246 = t151 * t160;
t242 = rSges(5,3) * t246 + t160 * t133;
t231 = t150 * t269;
t239 = t151 * rSges(5,3) + t133;
t79 = t231 - t239;
t138 = t151 * qJ(4);
t276 = pkin(3) * t150;
t100 = -t138 + t276;
t94 = t160 * t100;
t287 = t160 * t79 + t242 + t94;
t107 = Icges(6,5) * t155 - Icges(6,6) * t153;
t258 = Icges(6,4) * t153;
t108 = Icges(6,2) * t155 + t258;
t194 = t153 * t108 - t155 * t110;
t286 = t107 * qJD(5) + t160 * t194;
t106 = Icges(6,5) * t153 + Icges(6,6) * t155;
t175 = Icges(6,3) * t160 - qJD(5) * t106;
t185 = t195 * t151;
t72 = Icges(6,6) * t150 + t185;
t262 = t153 * t72;
t111 = Icges(6,1) * t155 - t258;
t186 = t111 * t151;
t74 = Icges(6,5) * t150 + t186;
t198 = -t155 * t74 + t262;
t250 = t150 * t160;
t285 = -t107 * t250 + t151 * t175 + t160 * t198;
t184 = t107 * t151;
t251 = t150 * t155;
t252 = t150 * t153;
t71 = Icges(6,4) * t251 - Icges(6,2) * t252 - Icges(6,6) * t151;
t263 = t153 * t71;
t123 = Icges(6,4) * t252;
t73 = Icges(6,1) * t251 - Icges(6,5) * t151 - t123;
t199 = -t155 * t73 + t263;
t284 = t150 * t175 + (t184 + t199) * t160;
t135 = qJD(4) * t150;
t125 = rSges(6,2) * t252;
t266 = rSges(6,2) * t155;
t225 = qJD(5) * t266;
t187 = rSges(6,3) * t246 + t160 * t125 - t151 * t225;
t235 = qJD(5) * t153;
t222 = t151 * t235;
t245 = t151 * t163;
t63 = -t245 - t138 + (pkin(3) - t152) * t150;
t230 = rSges(6,1) * t251;
t241 = t151 * rSges(6,3) + t125;
t76 = t230 - t241;
t65 = t160 * t76;
t283 = -rSges(6,1) * t222 - t160 * t63 + t135 + t187 + t65 + t94;
t69 = Icges(6,5) * t251 - Icges(6,6) * t252 - Icges(6,3) * t151;
t28 = -t150 * t199 - t151 * t69;
t271 = -Icges(6,2) * t251 - t123 + t73;
t273 = t110 * t150 + t71;
t282 = -t153 * t271 - t155 * t273;
t281 = t150 / 0.2e1;
t280 = -t151 / 0.2e1;
t279 = t160 / 0.2e1;
t278 = pkin(2) * t154;
t277 = pkin(2) * qJD(2) ^ 2;
t275 = -t150 * t69 - t73 * t247;
t70 = Icges(6,3) * t150 + t184;
t274 = t150 * t70 + t74 * t247;
t272 = -t110 * t151 - t72;
t270 = -t108 * t151 + t74;
t268 = rSges(6,1) * t155;
t116 = rSges(6,1) * t153 + t266;
t236 = qJD(5) * t151;
t223 = t116 * t236;
t203 = t135 - t223;
t181 = t203 - t227;
t26 = (-t100 + t63 - t76) * t160 + t181;
t261 = t160 * t26;
t254 = t106 * t151;
t45 = -t150 * t194 - t254;
t260 = t45 * t160;
t255 = t106 * t150;
t253 = t107 * t160;
t244 = -t108 + t111;
t243 = t110 + t195;
t129 = qJ(4) * t246;
t240 = t129 + t135;
t238 = qJD(4) * t160;
t237 = qJD(5) * t150;
t234 = t154 * t277;
t156 = cos(t159);
t233 = t156 * t277;
t226 = t156 * t265;
t224 = -t160 * t228 + (-rSges(6,1) * t235 - t225) * t150;
t93 = t116 * t237;
t219 = -pkin(3) - t269;
t218 = -t237 / 0.2e1;
t215 = t236 / 0.2e1;
t57 = t74 * t251;
t213 = t151 * t70 - t57;
t212 = -t69 + t262;
t53 = (-t155 * t250 - t222) * rSges(6,1) + t187;
t211 = t53 + t65;
t54 = t160 * t191 + t224;
t210 = -t160 * t77 + t54;
t103 = t151 * rSges(4,1) - rSges(4,2) * t150;
t206 = t151 * t238 - t233;
t92 = rSges(4,1) * t246 - rSges(4,2) * t250;
t136 = qJD(4) * t151;
t204 = -t136 + t226;
t201 = -rSges(6,2) * t153 + t268;
t42 = t153 * t73 + t155 * t71;
t43 = t153 * t74 + t155 * t72;
t120 = t160 * t249;
t197 = t120 + t136 - t224;
t193 = -t150 * t152 - t245;
t192 = t150 * t238 + t160 * (-pkin(3) * t250 + t240) - t234;
t29 = -t252 * t72 - t213;
t189 = (t150 * t29 - t151 * t28) * qJD(5);
t30 = -t248 * t71 - t275;
t31 = -t248 * t72 + t274;
t188 = (t150 * t31 - t151 * t30) * qJD(5);
t182 = -t153 * t270 + t155 * t272;
t180 = t150 * t219 + t138 + t239;
t178 = (-t153 * t243 + t155 * t244) * t160;
t177 = Icges(6,5) * t160 - qJD(5) * t110;
t176 = Icges(6,6) * t160 - qJD(5) * t108;
t174 = -t245 + (-t152 - t268) * t150 + t241;
t46 = -t151 * t194 + t255;
t44 = t46 * t160;
t10 = t44 + t188;
t50 = t150 * t176 + t160 * t185;
t52 = t150 * t177 + t160 * t186;
t14 = -qJD(5) * t199 + t153 * t52 + t155 * t50;
t49 = t151 * t176 - t195 * t250;
t51 = -t111 * t250 + t151 * t177;
t15 = -qJD(5) * t198 + t153 * t51 + t155 * t49;
t96 = t195 * qJD(5);
t97 = t111 * qJD(5);
t167 = t106 * t160 - t153 * t96 + t155 * t97 + (-t108 * t155 - t110 * t153) * qJD(5);
t20 = t150 * t286 + t167 * t151;
t21 = t167 * t150 - t151 * t286;
t9 = t189 + t260;
t173 = (-qJD(5) * t194 + t153 * t97 + t155 * t96) * t160 + (t44 + ((t29 - t57 + (t70 + t263) * t151 + t275) * t151 + t274 * t150) * qJD(5)) * t215 + (t9 - t260 + ((t151 * t212 - t274 + t31) * t151 + (t150 * t212 + t213 + t30) * t150) * qJD(5)) * t218 + (t15 + t20) * t237 / 0.2e1 - (t14 + t21 + t10) * t236 / 0.2e1 + ((t42 + t45) * t150 + (t43 + t46) * t151) * qJD(5) * t279;
t169 = -qJD(5) * t43 - t153 * t49 + t155 * t51 + t160 * t70;
t168 = -qJD(5) * t42 - t153 * t50 + t155 * t52 + t160 * t69;
t55 = (-t100 - t79) * t160 + t135 - t227;
t56 = t204 + t289;
t166 = (t55 * t219 * t151 + (t55 * (-rSges(5,3) - qJ(4)) + t56 * t219) * t150) * t160;
t27 = t160 * t290 + t204 - t93;
t165 = (t26 * (-t191 - t124) + t27 * (t193 - t230)) * t160;
t149 = pkin(2) * t156;
t113 = t160 * t229;
t98 = t201 * qJD(5);
t90 = t116 * t151;
t89 = t116 * t150;
t82 = t103 * t160 + t226;
t78 = t160 * t200 - t136;
t67 = -t160 * t92 - t233;
t66 = -t160 * t256 - t234;
t39 = qJD(1) + (t150 * t76 + t151 * t77) * qJD(5);
t33 = (t160 * t288 + t113 - t78) * t160 + t206;
t32 = t160 * (-t160 * t231 + t242) + t192;
t17 = -t98 * t236 + (t93 + t120 - t54 - t78 + (t200 - t124) * t160) * t160 + t206;
t16 = -t98 * t237 + (-t223 - t129 + t53 + (t193 + t276) * t160) * t160 + t192;
t11 = (t210 * t150 + t211 * t151) * qJD(5);
t1 = [m(6) * t11; m(4) * (t67 * (-t101 - t278) + t66 * (t103 + t149) + (-t92 - t226 + t82) * t81) + t173 + (t17 * (t174 - t278) + t26 * (t197 - t226) + t16 * (t149 + t290) + t165 + (-t181 + t26 - t227 + t283) * t27) * m(6) + (t33 * (t180 - t278) + t55 * (t113 - t204) + t32 * (t149 + t183) + t166 + (t55 + t129 + t287) * t56) * m(5); t173 + (t17 * t174 + t165 + (-t203 + t283) * t27 + (-t136 - t93 + t197) * t26 + (t16 + t261) * t290) * m(6) + (t33 * t180 + t32 * t183 + t166 + (-t135 + t240 + t287) * t56 + (t113 + t289) * t55) * m(5) + (-(-t82 * t101 - t103 * t81) * t160 - t101 * t67 + t103 * t66 - t81 * t92 - t82 * t256) * m(4); 0.2e1 * (t16 * t280 + t17 * t281) * m(6) + 0.2e1 * (t32 * t280 + t281 * t33) * m(5); ((t160 * t43 - t14) * t151 + (t160 * t42 + t15) * t150) * t279 + ((-t237 * t254 + t253) * t150 + (t178 + (-t282 * t151 + (t255 + t182) * t150) * qJD(5)) * t151) * t218 + ((-t236 * t255 - t253) * t151 + (t178 + (t182 * t150 + (-t282 + t254) * t151) * qJD(5)) * t150) * t215 - t160 * ((t244 * t153 + t243 * t155) * t160 + ((t150 * t270 - t151 * t271) * t155 + (t150 * t272 + t151 * t273) * t153) * qJD(5)) / 0.2e1 + (t160 * t20 + ((-t150 * t284 - t168 * t151 + t160 * t31) * t151 + (t150 * t285 + t169 * t151 + t160 * t30) * t150) * t291) * t281 + (t160 * t21 + ((-t168 * t150 + t151 * t284 + t160 * t29) * t151 + (t169 * t150 - t151 * t285 + t160 * t28) * t150) * t291) * t280 + (t9 + t189) * t250 / 0.2e1 + (t10 + t188) * t246 / 0.2e1 + ((t11 * t77 + t39 * t211 - t26 * t98) * t151 + (t11 * t76 + t39 * t210 - t27 * t98) * t150 + ((-t160 * t27 - t17) * t151 + (-t16 + t261) * t150) * t116 - (t26 * t89 - t27 * t90) * t160 - (t39 * (-t150 * t89 - t151 * t90) + (-t150 * t27 - t151 * t26) * t201) * qJD(5)) * m(6);];
tauc = t1(:);
