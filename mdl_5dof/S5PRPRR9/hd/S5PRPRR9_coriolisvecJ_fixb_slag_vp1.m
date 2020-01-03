% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
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
% Datum: 2019-12-31 17:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRPRR9_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR9_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR9_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR9_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR9_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR9_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRR9_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:39:37
% EndTime: 2019-12-31 17:39:43
% DurationCPUTime: 4.55s
% Computational Cost: add. (8909->314), mult. (9829->415), div. (0->0), fcn. (10136->6), ass. (0->175)
t150 = cos(qJ(5));
t149 = sin(qJ(5));
t231 = pkin(8) + qJ(2);
t217 = sin(t231);
t218 = cos(t231);
t276 = sin(qJ(4));
t277 = cos(qJ(4));
t106 = -t217 * t276 - t218 * t277;
t107 = -t217 * t277 + t218 * t276;
t257 = Icges(6,4) * t150;
t195 = -Icges(6,2) * t149 + t257;
t60 = -Icges(6,6) * t107 + t106 * t195;
t262 = t149 * t60;
t258 = Icges(6,4) * t149;
t197 = Icges(6,1) * t150 - t258;
t64 = -Icges(6,5) * t107 + t106 * t197;
t200 = t150 * t64 - t262;
t61 = Icges(6,6) * t106 + t107 * t195;
t263 = t149 * t61;
t65 = Icges(6,5) * t106 + t107 * t197;
t304 = -t150 * t65 + t263;
t245 = t107 * t150;
t193 = Icges(6,5) * t150 - Icges(6,6) * t149;
t57 = Icges(6,3) * t106 + t107 * t193;
t272 = -t106 * t57 - t245 * t65;
t56 = -Icges(6,3) * t107 + t106 * t193;
t303 = -(t56 - t263) * t107 + t272;
t235 = qJD(5) * t107;
t236 = qJD(5) * t106;
t148 = qJD(2) - qJD(4);
t278 = -t148 / 0.2e1;
t34 = -t149 * t65 - t150 * t61;
t35 = -t149 * t64 - t150 * t60;
t302 = ((-t35 * t106 + t34 * t107) * qJD(5) - t34 * t235 + t35 * t236) * t278;
t300 = t107 * t56;
t295 = t57 * t107;
t293 = t57 + t262;
t208 = rSges(6,1) * t150 - rSges(6,2) * t149;
t119 = t208 * qJD(5);
t151 = qJD(2) ^ 2;
t212 = t217 * pkin(3);
t139 = qJD(3) * t217;
t198 = qJD(2) * t217;
t199 = qJD(2) * t218;
t239 = qJ(3) * t199 + t139;
t259 = (-pkin(2) * t198 + t139 + t239) * qJD(2);
t175 = -t151 * t212 + t259;
t207 = rSges(6,1) * t149 + rSges(6,2) * t150;
t93 = t148 * t106;
t94 = t148 * t107;
t228 = t94 * pkin(4) + t93 * pkin(7);
t233 = qJD(5) * t149;
t183 = t106 * t233 + t150 * t94;
t180 = rSges(6,1) * t183 + t93 * rSges(6,3);
t232 = qJD(5) * t150;
t184 = t106 * t232 - t149 * t94;
t28 = rSges(6,2) * t184 + t180;
t11 = (t228 + t28) * t148 + (t107 * t119 + t207 * t93) * qJD(5) + t175;
t146 = t218 * pkin(3);
t140 = qJD(3) * t218;
t237 = t218 * pkin(2) + t217 * qJ(3);
t216 = (-qJD(2) * t237 + 0.2e1 * t140) * qJD(2);
t166 = -t146 * t151 + t216;
t229 = t93 * pkin(4) - t94 * pkin(7);
t181 = t107 * t233 - t150 * t93;
t179 = rSges(6,1) * t181 + t94 * rSges(6,3);
t182 = t107 * t232 + t149 * t93;
t27 = rSges(6,2) * t182 + t179;
t12 = (t229 - t27) * t148 + (t106 * t119 - t207 * t94) * qJD(5) + t166;
t142 = t218 * qJ(3);
t213 = t217 * pkin(2);
t111 = t213 - t142;
t169 = t139 + (-t212 - t111) * qJD(2);
t234 = qJD(5) * t207;
t103 = t106 * rSges(6,3);
t243 = rSges(6,1) * t245 + t103;
t246 = t107 * t149;
t70 = rSges(6,2) * t246 - t243;
t97 = -t107 * pkin(4) - t106 * pkin(7);
t287 = t106 * t234 - t148 * (t70 + t97);
t29 = t169 + t287;
t230 = t146 + t237;
t168 = qJD(2) * t230 - t140;
t219 = -t106 * pkin(4) + pkin(7) * t107;
t225 = t107 * t234;
t102 = t107 * rSges(6,3);
t248 = t106 * t150;
t242 = -rSges(6,1) * t248 + t102;
t249 = t106 * t149;
t71 = rSges(6,2) * t249 + t242;
t30 = t225 + (t219 + t71) * t148 + t168;
t292 = (t149 * (t106 * t11 - t107 * t12 - t29 * t93 - t30 * t94) + (t106 * t30 - t107 * t29) * t232) * rSges(6,2);
t210 = t107 * rSges(5,1) - t106 * rSges(5,2);
t291 = t148 * t210;
t192 = Icges(6,5) * t149 + Icges(6,6) * t150;
t78 = t192 * t106;
t77 = t192 * t107;
t290 = qJD(2) * t111 - t139 + t239;
t289 = 0.2e1 * qJD(5);
t194 = Icges(6,2) * t150 + t258;
t196 = Icges(6,1) * t149 + t257;
t191 = -t149 * t194 + t150 * t196;
t40 = t107 * t191 + t78;
t38 = t40 * t148;
t41 = t106 * t191 - t77;
t39 = t41 * t148;
t161 = t218 * rSges(4,1) + t217 * rSges(4,3);
t288 = t237 + t161;
t165 = -t213 - t212;
t285 = pkin(3) * t198 + t165 * qJD(2) + t290;
t266 = t194 * t107 - t65;
t268 = -t196 * t107 - t61;
t284 = t149 * t266 + t268 * t150;
t271 = t106 * t56 + t245 * t64;
t270 = t248 * t65 - t295;
t269 = t248 * t64 - t300;
t267 = -t196 * t106 - t60;
t265 = t194 * t106 - t64;
t244 = t193 * t148;
t241 = -t194 + t197;
t240 = -t195 - t196;
t221 = t236 / 0.2e1;
t220 = -t235 / 0.2e1;
t215 = -t218 / 0.2e1;
t214 = t217 / 0.2e1;
t47 = t94 * rSges(5,1) - t93 * rSges(5,2);
t46 = -t93 * rSges(5,1) - t94 * rSges(5,2);
t211 = t217 * rSges(4,1);
t209 = rSges(5,1) * t106 + rSges(5,2) * t107;
t205 = -t106 * t29 - t107 * t30;
t204 = t106 * t71 + t107 * t70;
t190 = -t97 + t243;
t189 = -t219 - t242;
t15 = -t246 * t61 - t272;
t16 = -t246 * t60 + t271;
t178 = qJD(5) * (-t106 * t15 + t107 * t16);
t17 = -t249 * t61 + t270;
t18 = -t249 * t60 + t269;
t177 = (-t106 * t17 + t107 * t18) * qJD(5);
t176 = t149 * t265 + t150 * t267;
t173 = -t180 - t228;
t172 = t179 - t229;
t170 = t106 * t28 + t107 * t27 + t70 * t93 - t71 * t94;
t167 = (t149 * t240 + t150 * t241) * t148;
t164 = -t211 - t213;
t158 = t142 + t165;
t117 = t195 * qJD(5);
t118 = t197 * qJD(5);
t154 = -t117 * t149 + t118 * t150 + (-t149 * t196 - t150 * t194) * qJD(5);
t10 = qJD(5) * t200 - t149 * (Icges(6,1) * t183 + Icges(6,4) * t184 + Icges(6,5) * t93) - t150 * (Icges(6,4) * t183 + Icges(6,2) * t184 + Icges(6,6) * t93);
t116 = t193 * qJD(5);
t13 = t106 * t116 + t107 * t154 + t191 * t93 - t192 * t94;
t14 = t106 * t154 - t107 * t116 - t191 * t94 - t192 * t93;
t5 = t38 + t178;
t6 = t39 + t177;
t9 = -qJD(5) * t304 - t149 * (Icges(6,1) * t181 + Icges(6,4) * t182 + Icges(6,5) * t94) - t150 * (Icges(6,4) * t181 + Icges(6,2) * t182 + Icges(6,6) * t94);
t152 = t148 * (-t118 * t149 + t194 * t233 + (-qJD(5) * t196 - t117) * t150) + t5 * t235 / 0.2e1 + (t14 + t10) * t220 + (t13 + t6 + t9) * t221 - ((t41 - t35) * t93 + (t40 - t34) * t94) * qJD(5) / 0.2e1;
t144 = t218 * rSges(4,3);
t138 = rSges(4,3) * t199;
t84 = t207 * t106;
t83 = t207 * t107;
t74 = -t151 * t161 + t216;
t73 = qJD(2) * (-rSges(4,1) * t198 + t138) + t259;
t69 = t107 * t208 + t103;
t68 = t106 * t208 - t102;
t55 = -t209 * t148 + t168;
t54 = t169 + t291;
t37 = -t148 * t46 + t166;
t36 = t148 * t47 + t175;
t31 = qJD(5) * t204 + qJD(1);
t22 = Icges(6,5) * t183 + Icges(6,6) * t184 + Icges(6,3) * t93;
t21 = Icges(6,5) * t181 + Icges(6,6) * t182 + Icges(6,3) * t94;
t8 = (-t106 * t61 - t107 * t60) * t149 + t270 + t271;
t7 = t170 * qJD(5);
t1 = [m(6) * t7; (-t38 + ((t17 + (-t200 + t57) * t107) * t107 + (t293 * t106 + t18 - t269 - t300) * t106) * qJD(5)) * t220 + (t39 + ((t304 * t107 + t15 + t269) * t107 + (-t293 * t107 + t16 - t8) * t106) * qJD(5)) * t221 - t152 + t302 + (t12 * (t158 + t190) + t11 * (-t189 + t230) + t292 - t31 * (t69 + t70) * t236 + (-t172 - t168) * t29 + (-t173 + t29 - (t69 - t97) * t148 - t207 * t236 + t285) * t30) * m(6) + (t37 * (t158 + t210) + t36 * (-t209 + t230) + (-t168 - t46) * t54 + (t285 + t47 + t54 - t291) * t55) * m(5) + (t74 * (t142 + t144 + t164) + t73 * t288 + (t138 + (t211 - t144 + t164) * qJD(2) + t290) * (qJD(2) * t288 - t140)) * m(4); 0.2e1 * (t11 * t215 + t12 * t214) * m(6) + 0.2e1 * (t214 * t37 + t215 * t36) * m(5) + 0.2e1 * (t214 * t74 + t215 * t73) * m(4); (t38 + ((-t17 + t8) * t107 + (t200 * t106 - t18 + t303) * t106) * qJD(5)) * t220 + (-t39 + ((-t15 - t303) * t107 + (-t16 + (-t304 + t56) * t106 - t295) * t106) * qJD(5)) * t221 + t152 + t302 + (-t12 * t190 + t11 * t189 - t292 - t31 * (t68 + t71) * t235 + (t173 + t287) * t30 + (t172 - t225 - (-t68 + t219) * t148) * t29) * m(6) + (-(-t209 * t54 - t210 * t55) * t148 + t36 * t209 - t37 * t210 + t46 * t54 - t47 * t55) * m(5); t148 * (t10 * t107 - t106 * t9 - t34 * t94 - t35 * t93) / 0.2e1 + ((t235 * t78 - t244) * t107 + (t167 + (-t284 * t106 + (-t77 + t176) * t107) * qJD(5)) * t106) * t220 + ((t236 * t77 + t244) * t106 + (t167 + (t176 * t107 + (-t284 - t78) * t106) * qJD(5)) * t107) * t221 + ((t149 * t241 - t240 * t150) * t148 + ((t106 * t266 - t107 * t265) * t150 + (-t106 * t268 + t107 * t267) * t149) * qJD(5)) * t278 - (t13 * t148 + (-(-t106 * t21 - t304 * t93 - t57 * t94) * t106 + t107 * (-t106 * t22 + t200 * t93 - t56 * t94) + t15 * t94 + t16 * t93) * t289) * t106 / 0.2e1 + (t14 * t148 + (-t106 * (t107 * t21 + t304 * t94 - t57 * t93) + t107 * (t107 * t22 - t200 * t94 - t56 * t93) + t17 * t94 + t18 * t93) * t289) * t107 / 0.2e1 + (t6 + t177) * t93 / 0.2e1 + (t5 + t178) * t94 / 0.2e1 + (t7 * t204 + t31 * t170 - t205 * t119 - (-t106 * t12 - t107 * t11 + t29 * t94 - t30 * t93) * t207 - (-t29 * t83 + t30 * t84) * t148 - (t31 * (t106 * t84 + t107 * t83) - t205 * t208) * qJD(5)) * m(6);];
tauc = t1(:);
