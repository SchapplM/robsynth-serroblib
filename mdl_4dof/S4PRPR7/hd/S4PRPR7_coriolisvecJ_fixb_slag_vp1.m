% Calculate vector of centrifugal and Coriolis load on the joints for
% S4PRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
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
% Datum: 2019-12-31 16:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PRPR7_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR7_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR7_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR7_coriolisvecJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR7_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPR7_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPR7_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:25:37
% EndTime: 2019-12-31 16:25:50
% DurationCPUTime: 11.52s
% Computational Cost: add. (3931->430), mult. (11144->689), div. (0->0), fcn. (10513->6), ass. (0->223)
t189 = sin(pkin(6));
t192 = sin(qJ(2));
t194 = cos(qJ(2));
t371 = (Icges(4,5) - Icges(3,6)) * t194 + (Icges(4,4) - Icges(3,5)) * t192;
t357 = t371 * t189;
t190 = cos(pkin(6));
t356 = t371 * t190;
t370 = -Icges(3,2) - Icges(4,3);
t369 = (Icges(3,4) + Icges(4,6)) * t194;
t317 = Icges(3,4) * t192;
t365 = t369 * t194 + (-t317 + (Icges(3,1) + t370) * t194) * t192;
t283 = qJD(4) * t194;
t290 = qJD(2) * t189;
t164 = t190 * t283 + t290;
t191 = sin(qJ(4));
t193 = cos(qJ(4));
t258 = rSges(5,1) * t191 + rSges(5,2) * t193;
t211 = -rSges(5,3) * t192 + t194 * t258;
t364 = t164 * t211;
t288 = qJD(2) * t190;
t165 = t189 * t283 - t288;
t363 = t165 * t211;
t187 = t189 ^ 2;
t188 = t190 ^ 2;
t291 = t187 + t188;
t362 = t371 * qJD(2);
t286 = qJD(2) * t194;
t287 = qJD(2) * t192;
t360 = (-Icges(4,6) * t192 + t370 * t194 - t317) * t287 + ((Icges(3,1) + Icges(4,2)) * t192 + t369) * t286;
t359 = t362 * t189;
t358 = t362 * t190;
t303 = t190 * t194;
t304 = t190 * t192;
t355 = t360 * t190 + ((Icges(4,2) * t303 - Icges(4,6) * t304) * t192 + t365 * t190 - t357) * qJD(2);
t305 = t189 * t194;
t306 = t189 * t192;
t354 = t360 * t189 + ((Icges(4,2) * t305 - Icges(4,6) * t306) * t192 + t365 * t189 + t356) * qJD(2);
t301 = t192 * t193;
t142 = -t189 * t191 + t190 * t301;
t302 = t191 * t192;
t143 = t189 * t193 + t190 * t302;
t56 = Icges(5,5) * t143 + Icges(5,6) * t142 + Icges(5,3) * t303;
t315 = Icges(5,4) * t143;
t58 = Icges(5,2) * t142 + Icges(5,6) * t303 + t315;
t138 = Icges(5,4) * t142;
t60 = Icges(5,1) * t143 + Icges(5,5) * t303 + t138;
t19 = t142 * t58 + t143 * t60 + t303 * t56;
t144 = t189 * t301 + t190 * t191;
t145 = -t189 * t302 + t190 * t193;
t57 = -Icges(5,5) * t145 + Icges(5,6) * t144 + Icges(5,3) * t305;
t314 = Icges(5,4) * t145;
t59 = Icges(5,2) * t144 + Icges(5,6) * t305 - t314;
t139 = Icges(5,4) * t144;
t61 = -Icges(5,1) * t145 + Icges(5,5) * t305 + t139;
t20 = t142 * t59 + t143 * t61 + t303 * t57;
t21 = t144 * t58 - t145 * t60 + t305 * t56;
t22 = t144 * t59 - t145 * t61 + t305 * t57;
t284 = qJD(4) * t192;
t313 = Icges(5,4) * t191;
t241 = Icges(5,2) * t193 + t313;
t209 = -Icges(5,6) * t192 + t194 * t241;
t312 = Icges(5,4) * t193;
t245 = Icges(5,1) * t191 + t312;
t210 = -Icges(5,5) * t192 + t194 * t245;
t239 = Icges(5,5) * t191 + Icges(5,6) * t193;
t208 = -Icges(5,3) * t192 + t194 * t239;
t300 = t194 * t208;
t38 = -t142 * t209 - t143 * t210 - t190 * t300;
t39 = -t144 * t209 + t145 * t210 - t189 * t300;
t349 = t190 * (t164 * t19 + t165 * t20 + t284 * t38) + t189 * (t164 * t21 + t165 * t22 + t284 * t39);
t257 = rSges(4,2) * t192 + rSges(4,3) * t194;
t348 = t291 * qJD(2) * t257;
t334 = t164 / 0.2e1;
t332 = t165 / 0.2e1;
t347 = -t287 / 0.2e1;
t172 = rSges(3,1) * t192 + rSges(3,2) * t194;
t224 = qJD(2) * t172;
t161 = (Icges(5,2) * t191 - t312) * t194;
t201 = t164 * (-Icges(5,2) * t143 + t138 + t60) + t165 * (Icges(5,2) * t145 + t139 + t61) + t284 * (-t210 + t161);
t162 = (-Icges(5,1) * t193 + t313) * t194;
t202 = t164 * (-Icges(5,1) * t142 + t315 + t58) + t165 * (-Icges(5,1) * t144 - t314 + t59) + t284 * (-t209 - t162);
t195 = qJD(2) ^ 2;
t112 = Icges(5,3) * t194 + t192 * t239;
t160 = (-Icges(5,5) * t193 + Icges(5,6) * t191) * t194;
t64 = qJD(2) * t112 + qJD(4) * t160;
t227 = t194 * t64 + t208 * t287;
t323 = t189 * t20;
t255 = t190 * t19 + t323;
t321 = t38 * t194;
t114 = Icges(5,6) * t194 + t192 * t241;
t65 = qJD(2) * t114 + qJD(4) * t161;
t116 = Icges(5,5) * t194 + t192 * t245;
t66 = qJD(2) * t116 + qJD(4) * t162;
t273 = t191 * t286;
t82 = qJD(4) * t142 + t190 * t273;
t272 = t193 * t286;
t83 = -qJD(4) * t143 + t190 * t272;
t200 = (t142 * t65 + t143 * t66 + t190 * t227 - t209 * t83 - t210 * t82) * t192 + (-t192 * t255 + t321) * qJD(2);
t274 = t190 * t287;
t42 = Icges(5,5) * t82 + Icges(5,6) * t83 - Icges(5,3) * t274;
t229 = t194 * t42 - t287 * t56;
t44 = Icges(5,4) * t82 + Icges(5,2) * t83 - Icges(5,6) * t274;
t46 = Icges(5,1) * t82 + Icges(5,4) * t83 - Icges(5,5) * t274;
t7 = t142 * t44 + t143 * t46 + t190 * t229 + t58 * t83 + t60 * t82;
t275 = t189 * t287;
t84 = qJD(4) * t144 + t189 * t273;
t85 = qJD(4) * t145 + t189 * t272;
t43 = Icges(5,5) * t84 + Icges(5,6) * t85 - Icges(5,3) * t275;
t228 = t194 * t43 - t287 * t57;
t45 = Icges(5,4) * t84 + Icges(5,2) * t85 - Icges(5,6) * t275;
t47 = Icges(5,1) * t84 + Icges(5,4) * t85 - Icges(5,5) * t275;
t8 = t142 * t45 + t143 * t47 + t190 * t228 + t59 * t83 + t61 * t82;
t337 = qJD(4) * t200 / 0.2e1 + t7 * t334 + t8 * t332;
t336 = t291 * t347;
t335 = -t164 / 0.2e1;
t333 = -t165 / 0.2e1;
t331 = -t194 / 0.2e1;
t330 = pkin(5) * t195;
t186 = qJD(3) * t192;
t179 = t189 * t186;
t170 = pkin(2) * t192 - qJ(3) * t194;
t222 = qJD(2) * t170;
t94 = -t189 * t222 + t179;
t181 = t190 * t186;
t95 = -t190 * t222 + t181;
t325 = t189 * t94 + t190 * t95;
t322 = t190 * t21;
t320 = t39 * t194;
t307 = t208 * t192;
t173 = pkin(2) * t194 + qJ(3) * t192;
t158 = t173 * t189;
t159 = t173 * t190;
t295 = t189 * t158 + t190 * t159;
t285 = qJD(3) * t194;
t140 = qJD(2) * t173 - t285;
t174 = -rSges(4,2) * t194 + rSges(4,3) * t192;
t294 = -t174 * qJD(2) - t140;
t293 = -t170 + t257;
t292 = -t173 - t174;
t281 = qJD(2) * qJD(3);
t280 = t194 * t330;
t279 = t192 * t281 + t95 * t288 + t94 * t290;
t276 = -t222 * t291 + t186;
t271 = t194 * t281;
t267 = t286 / 0.2e1;
t266 = -t284 / 0.2e1;
t265 = t284 / 0.2e1;
t264 = -pkin(5) * t192 - t170;
t263 = qJD(2) * t294;
t262 = qJD(2) * t293;
t261 = t211 + t264;
t260 = qJD(2) * t266;
t259 = qJD(4) * t267;
t175 = rSges(3,1) * t194 - rSges(3,2) * t192;
t256 = -t56 * t164 - t57 * t165;
t254 = t189 * t22 + t322;
t250 = t191 * t60 + t193 * t58;
t23 = t192 * t56 - t194 * t250;
t249 = t191 * t61 + t193 * t59;
t24 = t192 * t57 - t194 * t249;
t253 = t24 * t189 + t23 * t190;
t248 = qJD(2) * t264;
t62 = rSges(5,1) * t143 + rSges(5,2) * t142 + rSges(5,3) * t303;
t36 = t189 * t248 + t284 * t62 + t179 + t364;
t63 = -rSges(5,1) * t145 + rSges(5,2) * t144 + rSges(5,3) * t305;
t37 = t190 * t248 - t284 * t63 + t181 - t363;
t252 = -t189 * t36 - t190 * t37;
t251 = t189 * t62 - t190 * t63;
t236 = t291 * t175;
t235 = -t191 * t210 - t193 * t209;
t234 = t291 * t224;
t118 = rSges(5,3) * t194 + t192 * t258;
t163 = (-rSges(5,1) * t193 + rSges(5,2) * t191) * t194;
t67 = qJD(2) * t118 + qJD(4) * t163;
t233 = -pkin(5) * t286 - t140 - t67;
t232 = t189 * t260;
t231 = t190 * t260;
t230 = t158 * t290 + t159 * t288 + qJD(1) - t285;
t166 = pkin(3) * t189 + pkin(5) * t303;
t167 = -pkin(3) * t190 + pkin(5) * t305;
t18 = t164 * t63 - t165 * t62 + (t166 * t190 + t167 * t189) * qJD(2) + t230;
t226 = t18 * t251;
t110 = rSges(4,1) * t189 + t174 * t190;
t111 = -rSges(4,1) * t190 + t174 * t189;
t40 = (t110 * t190 + t111 * t189) * qJD(2) + t230;
t225 = t40 * t257;
t215 = (t112 + t235) * t192;
t212 = t160 * t284 + t164 * (Icges(5,5) * t142 - Icges(5,6) * t143) + t165 * (Icges(5,5) * t144 + Icges(5,6) * t145);
t207 = t194 * t212;
t199 = (t144 * t65 - t145 * t66 + t189 * t227 - t209 * t85 - t210 * t84) * t192 + (-t192 * t254 + t320) * qJD(2);
t41 = -t194 * t235 - t307;
t198 = ((qJD(2) * t235 + t64) * t192 + (-qJD(2) * t208 - t191 * t66 - t193 * t65 + (-t191 * t209 + t193 * t210) * qJD(4)) * t194) * t192 + (-t192 * t253 + t194 * t41) * qJD(2);
t197 = (t208 * t190 + t250) * t164 + (t208 * t189 + t249) * t165;
t196 = (qJD(4) * t215 + t197) * t194;
t182 = t190 * t285;
t180 = t189 * t285;
t177 = t190 * t271;
t176 = t189 * t271;
t91 = t210 * t190;
t90 = t210 * t189;
t89 = t209 * t190;
t88 = t209 * t189;
t77 = rSges(5,1) * t144 + rSges(5,2) * t145;
t76 = rSges(5,1) * t142 - rSges(5,2) * t143;
t69 = t190 * t262 + t181;
t68 = t189 * t262 + t179;
t53 = t234 * qJD(2);
t52 = t190 * t263 + t177;
t51 = t189 * t263 + t176;
t50 = qJD(2) * t236 + qJD(1);
t49 = rSges(5,1) * t84 + rSges(5,2) * t85 - rSges(5,3) * t275;
t48 = rSges(5,1) * t82 + rSges(5,2) * t83 - rSges(5,3) * t274;
t35 = qJD(2) * t348 + t279;
t17 = -t190 * t280 - t49 * t284 + t165 * t67 + t177 + (-t140 * t190 + (-t194 * t63 + t211 * t306) * qJD(4)) * qJD(2);
t16 = -t189 * t280 + t48 * t284 - t164 * t67 + t176 + (-t140 * t189 + (t194 * t62 - t211 * t304) * qJD(4)) * qJD(2);
t12 = t164 * t49 - t165 * t48 + (qJD(2) * qJD(4) * t251 - t291 * t330) * t192 + t279;
t11 = t164 * t23 + t165 * t24 + t284 * t41;
t10 = t144 * t45 - t145 * t47 + t189 * t228 + t59 * t85 + t61 * t84;
t9 = t144 * t44 - t145 * t46 + t189 * t229 + t58 * t85 + t60 * t84;
t4 = (qJD(2) * t249 + t43) * t192 + (qJD(2) * t57 - t191 * t47 - t193 * t45 + (t191 * t59 - t193 * t61) * qJD(4)) * t194;
t3 = (qJD(2) * t250 + t42) * t192 + (qJD(2) * t56 - t191 * t46 - t193 * t44 + (t191 * t58 - t193 * t60) * qJD(4)) * t194;
t2 = qJD(4) * t199 + t10 * t165 + t164 * t9;
t1 = [-m(3) * t53 + m(4) * t35 + m(5) * t12; -t11 * t283 / 0.2e1 + (t189 * t23 - t190 * t24) * t259 + (t189 * t19 - t190 * t20) * t231 + (t189 * t21 - t190 * t22) * t232 + (((-t191 * t91 - t193 * t89 + t56) * t164 + (-t191 * t90 - t193 * t88 + t57) * t165 + t41 * qJD(4)) * t194 + ((t215 + (-t114 * t193 - t116 * t191 - t208) * t194 - t253) * qJD(4) + t197) * t192) * t266 + (-t10 * t190 + t189 * t9) * t332 + ((t144 * t89 - t145 * t91) * t164 + (t144 * t88 - t145 * t90) * t165 + (t320 + (t114 * t144 - t116 * t145 - t322) * t192) * qJD(4) + (((-t22 + t307) * qJD(4) + t256) * t192 + t196) * t189) * t333 + (t189 * t7 - t190 * t8) * t334 + ((t142 * t89 + t143 * t91) * t164 + (t142 * t88 + t143 * t90) * t165 + (t321 + (t114 * t142 + t116 * t143 - t323) * t192) * qJD(4) + (((-t19 + t307) * qJD(4) + t256) * t192 + t196) * t190) * t335 + t189 * t337 - t190 * t2 / 0.2e1 + (t354 * t188 + (t358 * t189 + (-t355 - t359) * t190) * t189) * t290 + (-t359 * t188 + (t355 * t189 + (-t354 + t358) * t190) * t189) * t288 - (t356 * qJD(2) * t187 - t357 * t189 * t288) * t290 / 0.2e1 + (t357 * t188 * qJD(2) - t190 * t356 * t290) * t288 / 0.2e1 + (-t37 * (t118 * t165 + t182) - t36 * (-t118 * t164 + t180) - t18 * (t189 * t364 - t190 * t363 + t276) - (t252 * t173 + (-t18 * t192 * t291 + t194 * t252) * pkin(5)) * qJD(2) - ((t36 * t62 - t37 * t63) * t194 + t226 * t192) * qJD(4) + t12 * t295 + t18 * (-pkin(5) * t287 * t291 + t325) + (t17 * t261 + t37 * t233 + t12 * (t166 + t62) + t18 * t48) * t190 + (t16 * t261 + t36 * t233 + t12 * (t167 + t63) + t18 * t49) * t189) * m(5) + (-t69 * t182 - t68 * t180 - ((t190 * t225 + t292 * t69) * t190 + (t189 * t225 + t292 * t68) * t189) * qJD(2) + t35 * t295 + (t35 * t110 + t293 * t52 + t294 * t69) * t190 + (t35 * t111 + t293 * t51 + t294 * t68) * t189 + (-t276 + t325 + t348) * t40) * m(4) + (t189 * t3 - t190 * t4 + t349) * t265 + (-t234 * t50 - t236 * t53 + (t172 * t175 * t195 + t224 * t50) * t291) * m(3); 0.2e1 * (t12 * t331 + t18 * t336) * m(5) + 0.2e1 * (t331 * t35 + t336 * t40) * m(4) + 0.2e1 * (m(4) * (qJD(2) * t40 + t189 * t51 + t190 * t52) / 0.2e1 + m(5) * (qJD(2) * t18 + t16 * t189 + t17 * t190) / 0.2e1) * t192; t303 * t337 + (t192 * t38 + t194 * t255) * t231 + ((t189 * t8 + t190 * t7) * t194 + t200) * t334 + t2 * t305 / 0.2e1 + (t192 * t39 + t194 * t254) * t232 + ((t10 * t189 + t190 * t9) * t194 + t199) * t332 + t11 * t267 + t192 * (qJD(4) * t198 + t164 * t3 + t165 * t4) / 0.2e1 + (t192 * t41 + t194 * t253) * t259 + ((t189 * t4 + t190 * t3) * t194 + t198) * t265 + (t201 * t142 - t143 * t202 + t190 * t207) * t335 + (t144 * t201 + t145 * t202 + t189 * t207) * t333 + (t212 * t192 + (t202 * t191 - t193 * t201) * t194) * t266 + t349 * t347 + ((t16 * t62 - t17 * t63 + t36 * t48 - t37 * t49 + (t226 - (-t189 * t37 + t190 * t36) * t211) * qJD(2)) * t192 + (t37 * (-qJD(2) * t63 + t189 * t67) + t36 * (qJD(2) * t62 - t190 * t67) - t12 * t251 + t18 * (-t189 * t48 + t190 * t49) - (-t16 * t190 + t17 * t189) * t211) * t194 - t37 * (t163 * t165 - t284 * t77) - t36 * (-t163 * t164 + t284 * t76) - t18 * (t164 * t77 - t165 * t76)) * m(5);];
tauc = t1(:);
