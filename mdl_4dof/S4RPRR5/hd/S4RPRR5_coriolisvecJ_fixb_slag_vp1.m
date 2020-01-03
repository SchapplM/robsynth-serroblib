% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% Datum: 2019-12-31 16:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPRR5_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR5_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR5_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR5_coriolisvecJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR5_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR5_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR5_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:51:32
% EndTime: 2019-12-31 16:51:38
% DurationCPUTime: 4.49s
% Computational Cost: add. (4330->311), mult. (9695->413), div. (0->0), fcn. (10012->6), ass. (0->171)
t149 = cos(qJ(4));
t148 = sin(qJ(4));
t271 = sin(qJ(3));
t272 = sin(qJ(1));
t273 = cos(qJ(3));
t274 = cos(qJ(1));
t109 = -t272 * t271 - t274 * t273;
t110 = t274 * t271 - t272 * t273;
t251 = Icges(5,4) * t149;
t190 = -Icges(5,2) * t148 + t251;
t60 = -Icges(5,6) * t110 + t109 * t190;
t257 = t148 * t60;
t252 = Icges(5,4) * t148;
t192 = Icges(5,1) * t149 - t252;
t64 = -Icges(5,5) * t110 + t109 * t192;
t193 = t149 * t64 - t257;
t61 = Icges(5,6) * t109 + t110 * t190;
t258 = t148 * t61;
t65 = Icges(5,5) * t109 + t110 * t192;
t301 = -t149 * t65 + t258;
t239 = t110 * t149;
t188 = Icges(5,5) * t149 - Icges(5,6) * t148;
t57 = Icges(5,3) * t109 + t110 * t188;
t267 = -t109 * t57 - t239 * t65;
t56 = -Icges(5,3) * t110 + t109 * t188;
t300 = -(t56 - t258) * t110 + t267;
t230 = qJD(4) * t110;
t231 = qJD(4) * t109;
t147 = qJD(1) - qJD(3);
t275 = -t147 / 0.2e1;
t33 = -t148 * t65 - t149 * t61;
t34 = -t148 * t64 - t149 * t60;
t299 = ((-t34 * t109 + t33 * t110) * qJD(4) - t33 * t230 + t34 * t231) * t275;
t297 = t110 * t56;
t292 = t57 * t110;
t290 = t57 + t257;
t200 = rSges(5,1) * t149 - rSges(5,2) * t148;
t111 = t200 * qJD(4);
t150 = qJD(1) ^ 2;
t223 = t272 * pkin(2);
t138 = qJD(2) * t272;
t211 = qJD(1) * t272;
t212 = qJD(1) * t274;
t234 = qJ(2) * t212 + t138;
t253 = (-pkin(1) * t211 + t138 + t234) * qJD(1);
t179 = -t150 * t223 + t253;
t199 = rSges(5,1) * t148 + rSges(5,2) * t149;
t90 = t147 * t109;
t91 = t147 * t110;
t219 = t91 * pkin(3) + t90 * pkin(6);
t228 = qJD(4) * t148;
t182 = t109 * t228 + t149 * t91;
t178 = rSges(5,1) * t182 + t90 * rSges(5,3);
t227 = qJD(4) * t149;
t183 = t109 * t227 - t148 * t91;
t25 = rSges(5,2) * t183 + t178;
t10 = (t219 + t25) * t147 + (t110 * t111 + t199 * t90) * qJD(4) + t179;
t145 = t274 * pkin(2);
t139 = qJD(2) * t274;
t232 = t274 * pkin(1) + t272 * qJ(2);
t207 = (-qJD(1) * t232 + 0.2e1 * t139) * qJD(1);
t164 = -t145 * t150 + t207;
t220 = t90 * pkin(3) - t91 * pkin(6);
t180 = t110 * t228 - t149 * t90;
t177 = rSges(5,1) * t180 + t91 * rSges(5,3);
t181 = t110 * t227 + t148 * t90;
t24 = rSges(5,2) * t181 + t177;
t11 = (t220 - t24) * t147 + (t109 * t111 - t199 * t91) * qJD(4) + t164;
t141 = t274 * qJ(2);
t224 = t272 * pkin(1);
t121 = t224 - t141;
t168 = t138 + (-t223 - t121) * qJD(1);
t229 = qJD(4) * t199;
t240 = t110 * t148;
t103 = t109 * rSges(5,3);
t254 = rSges(5,1) * t239 + t103;
t70 = rSges(5,2) * t240 - t254;
t94 = -t110 * pkin(3) - t109 * pkin(6);
t284 = t109 * t229 - (t70 + t94) * t147;
t28 = t168 + t284;
t221 = t145 + t232;
t166 = qJD(1) * t221 - t139;
t208 = -t109 * pkin(3) + pkin(6) * t110;
t216 = t110 * t229;
t102 = t110 * rSges(5,3);
t242 = t109 * t149;
t237 = -rSges(5,1) * t242 + t102;
t243 = t109 * t148;
t71 = rSges(5,2) * t243 + t237;
t29 = t216 + (t208 + t71) * t147 + t166;
t289 = (t148 * (t10 * t109 - t11 * t110 - t28 * t90 - t29 * t91) + (t109 * t29 - t110 * t28) * t227) * rSges(5,2);
t202 = t110 * rSges(4,1) - t109 * rSges(4,2);
t288 = t147 * t202;
t187 = Icges(5,5) * t148 + Icges(5,6) * t149;
t75 = t187 * t109;
t74 = t187 * t110;
t287 = qJD(1) * t121 - t138 + t234;
t286 = 0.2e1 * qJD(4);
t189 = Icges(5,2) * t149 + t252;
t191 = Icges(5,1) * t148 + t251;
t186 = -t148 * t189 + t149 * t191;
t39 = t110 * t186 + t75;
t37 = t39 * t147;
t40 = t109 * t186 - t74;
t38 = t40 * t147;
t171 = t274 * rSges(3,1) + t272 * rSges(3,3);
t285 = t232 + t171;
t173 = -t224 - t223;
t282 = pkin(2) * t211 + t173 * qJD(1) + t287;
t261 = t189 * t110 - t65;
t263 = -t191 * t110 - t61;
t281 = t148 * t261 + t149 * t263;
t266 = t109 * t56 + t239 * t64;
t265 = t242 * t65 - t292;
t264 = t242 * t64 - t297;
t262 = -t191 * t109 - t60;
t260 = t189 * t109 - t64;
t238 = t188 * t147;
t236 = -t189 + t192;
t235 = -t190 - t191;
t226 = -t274 / 0.2e1;
t225 = t272 / 0.2e1;
t222 = t272 * rSges(3,1);
t210 = t231 / 0.2e1;
t209 = -t230 / 0.2e1;
t46 = t91 * rSges(4,1) - t90 * rSges(4,2);
t45 = -t90 * rSges(4,1) - t91 * rSges(4,2);
t201 = rSges(4,1) * t109 + rSges(4,2) * t110;
t197 = -t109 * t28 - t110 * t29;
t185 = -t94 + t254;
t184 = -t208 - t237;
t14 = -t240 * t61 - t267;
t15 = -t240 * t60 + t266;
t176 = (-t109 * t14 + t110 * t15) * qJD(4);
t16 = -t243 * t61 + t265;
t17 = -t243 * t60 + t264;
t175 = (-t109 * t16 + t110 * t17) * qJD(4);
t30 = (t109 * t71 + t110 * t70) * qJD(4);
t172 = -t222 - t224;
t167 = t148 * t260 + t149 * t262;
t165 = t141 + t173;
t162 = -t178 - t219;
t161 = t177 - t220;
t158 = (t148 * t235 + t149 * t236) * t147;
t107 = t190 * qJD(4);
t108 = t192 * qJD(4);
t152 = -t107 * t148 + t108 * t149 + (-t148 * t191 - t149 * t189) * qJD(4);
t106 = t188 * qJD(4);
t12 = t106 * t109 + t110 * t152 + t186 * t90 - t187 * t91;
t13 = -t106 * t110 + t109 * t152 - t186 * t91 - t187 * t90;
t5 = t37 + t176;
t6 = t38 + t175;
t7 = -t301 * qJD(4) - t148 * (Icges(5,1) * t180 + Icges(5,4) * t181 + Icges(5,5) * t91) - t149 * (Icges(5,4) * t180 + Icges(5,2) * t181 + Icges(5,6) * t91);
t8 = qJD(4) * t193 - t148 * (Icges(5,1) * t182 + Icges(5,4) * t183 + Icges(5,5) * t90) - t149 * (Icges(5,4) * t182 + Icges(5,2) * t183 + Icges(5,6) * t90);
t151 = t147 * (-t108 * t148 + t189 * t228 + (-qJD(4) * t191 - t107) * t149) + t5 * t230 / 0.2e1 + (t13 + t8) * t209 + (t12 + t6 + t7) * t210 - ((t40 - t34) * t90 + (t39 - t33) * t91) * qJD(4) / 0.2e1;
t143 = t274 * rSges(3,3);
t137 = rSges(3,3) * t212;
t81 = t199 * t109;
t80 = t199 * t110;
t73 = -t150 * t171 + t207;
t72 = qJD(1) * (-rSges(3,1) * t211 + t137) + t253;
t69 = t110 * t200 + t103;
t68 = t109 * t200 - t102;
t55 = -t201 * t147 + t166;
t54 = t168 + t288;
t36 = -t147 * t45 + t164;
t35 = t147 * t46 + t179;
t19 = Icges(5,5) * t182 + Icges(5,6) * t183 + Icges(5,3) * t90;
t18 = Icges(5,5) * t180 + Icges(5,6) * t181 + Icges(5,3) * t91;
t9 = (-t109 * t61 - t110 * t60) * t148 + t265 + t266;
t1 = [-t151 + t299 + (-t37 + ((t16 + (-t193 + t57) * t110) * t110 + (t290 * t109 + t17 - t264 - t297) * t109) * qJD(4)) * t209 + (t38 + ((t301 * t110 + t14 + t264) * t110 + (-t110 * t290 + t15 - t9) * t109) * qJD(4)) * t210 + (t11 * (t165 + t185) + t10 * (-t184 + t221) + t289 - t30 * (t69 + t70) * t231 + (-t161 - t166) * t28 + (-t162 + t28 - (t69 - t94) * t147 - t199 * t231 + t282) * t29) * m(5) + (t36 * (t165 + t202) + t35 * (-t201 + t221) + (-t166 - t45) * t54 + (t282 + t46 + t54 - t288) * t55) * m(4) + (t73 * (t141 + t143 + t172) + t72 * t285 + (t137 + (t222 - t143 + t172) * qJD(1) + t287) * (qJD(1) * t285 - t139)) * m(3); 0.2e1 * (t10 * t226 + t11 * t225) * m(5) + 0.2e1 * (t225 * t36 + t226 * t35) * m(4) + 0.2e1 * (t225 * t73 + t226 * t72) * m(3); t151 + t299 + (t37 + ((-t16 + t9) * t110 + (t109 * t193 - t17 + t300) * t109) * qJD(4)) * t209 + (-t38 + ((-t14 - t300) * t110 + (-t15 + (-t301 + t56) * t109 - t292) * t109) * qJD(4)) * t210 + (-t11 * t185 + t10 * t184 - t289 - t30 * (t68 + t71) * t230 + (t162 + t284) * t29 + (t161 - t216 - (-t68 + t208) * t147) * t28) * m(5) + (t35 * t201 - t36 * t202 + t45 * t54 - t46 * t55 - (-t201 * t54 - t202 * t55) * t147) * m(4); t147 * (-t7 * t109 + t8 * t110 - t33 * t91 - t34 * t90) / 0.2e1 + ((t75 * t230 - t238) * t110 + (t158 + (-t281 * t109 + (-t74 + t167) * t110) * qJD(4)) * t109) * t209 + ((t74 * t231 + t238) * t109 + (t158 + (t167 * t110 + (-t281 - t75) * t109) * qJD(4)) * t110) * t210 + ((t148 * t236 - t149 * t235) * t147 + ((t109 * t261 - t110 * t260) * t149 + (-t109 * t263 + t110 * t262) * t148) * qJD(4)) * t275 - (t12 * t147 + (-(-t109 * t18 - t301 * t90 - t57 * t91) * t109 + t110 * (-t109 * t19 + t193 * t90 - t56 * t91) + t14 * t91 + t15 * t90) * t286) * t109 / 0.2e1 + (t13 * t147 + (-t109 * (t110 * t18 + t301 * t91 - t57 * t90) + t110 * (t110 * t19 - t193 * t91 - t56 * t90) + t16 * t91 + t17 * t90) * t286) * t110 / 0.2e1 + (t6 + t175) * t90 / 0.2e1 + (t5 + t176) * t91 / 0.2e1 + (0.2e1 * t30 * (t109 * t25 + t110 * t24 + t70 * t90 - t71 * t91) - t197 * t111 - (-t10 * t110 - t109 * t11 + t28 * t91 - t29 * t90) * t199 - (-t28 * t80 + t29 * t81) * t147 - (t30 * (t109 * t81 + t110 * t80) - t197 * t200) * qJD(4)) * m(5);];
tauc = t1(:);
