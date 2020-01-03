% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-31 18:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRR8_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR8_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR8_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR8_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR8_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR8_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR8_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:01:04
% EndTime: 2019-12-31 18:01:11
% DurationCPUTime: 5.34s
% Computational Cost: add. (8893->358), mult. (10493->455), div. (0->0), fcn. (10624->8), ass. (0->202)
t178 = cos(qJ(5));
t177 = sin(qJ(5));
t271 = pkin(8) + qJ(4);
t247 = sin(t271);
t248 = cos(t271);
t318 = sin(qJ(1));
t319 = cos(qJ(1));
t118 = -t318 * t247 - t319 * t248;
t119 = t319 * t247 - t318 * t248;
t299 = Icges(6,4) * t178;
t229 = -Icges(6,2) * t177 + t299;
t62 = -Icges(6,6) * t119 + t118 * t229;
t304 = t177 * t62;
t300 = Icges(6,4) * t177;
t231 = Icges(6,1) * t178 - t300;
t66 = -Icges(6,5) * t119 + t118 * t231;
t232 = t178 * t66 - t304;
t63 = Icges(6,6) * t118 + t119 * t229;
t305 = t177 * t63;
t67 = Icges(6,5) * t118 + t119 * t231;
t349 = -t178 * t67 + t305;
t287 = t119 * t178;
t227 = Icges(6,5) * t178 - Icges(6,6) * t177;
t59 = Icges(6,3) * t118 + t119 * t227;
t314 = -t118 * t59 - t287 * t67;
t58 = -Icges(6,3) * t119 + t118 * t227;
t348 = -(t58 - t305) * t119 + t314;
t275 = qJD(5) * t119;
t276 = qJD(5) * t118;
t175 = qJD(1) - qJD(4);
t320 = -t175 / 0.2e1;
t36 = -t177 * t67 - t178 * t63;
t37 = -t177 * t66 - t178 * t62;
t347 = ((-t37 * t118 + t36 * t119) * qJD(5) - t36 * t275 + t37 * t276) * t320;
t345 = t119 * t58;
t340 = t59 * t119;
t338 = t59 + t304;
t240 = rSges(6,1) * t178 - rSges(6,2) * t177;
t131 = t240 * qJD(5);
t301 = cos(pkin(8));
t162 = pkin(3) * t301 + pkin(2);
t151 = t318 * t162;
t176 = sin(pkin(8));
t253 = qJD(1) * t319;
t245 = t176 * t253;
t153 = pkin(3) * t245;
t179 = qJD(1) ^ 2;
t266 = t318 * pkin(2);
t167 = qJD(2) * t318;
t252 = qJD(1) * t318;
t279 = qJ(2) * t253 + t167;
t284 = (-pkin(1) * t252 + t167 + t279) * qJD(1);
t209 = -t179 * t266 + t284;
t201 = qJD(1) * (t153 + (t266 - t151) * qJD(1)) + t209;
t239 = rSges(6,1) * t177 + rSges(6,2) * t178;
t95 = t175 * t118;
t96 = t175 * t119;
t262 = t96 * pkin(4) + t95 * pkin(7);
t273 = qJD(5) * t177;
t214 = t118 * t273 + t178 * t96;
t211 = t214 * rSges(6,1) + t95 * rSges(6,3);
t272 = qJD(5) * t178;
t215 = t118 * t272 - t177 * t96;
t30 = rSges(6,2) * t215 + t211;
t11 = (t262 + t30) * t175 + (t119 * t131 + t239 * t95) * qJD(5) + t201;
t168 = qJD(2) * t319;
t277 = t319 * pkin(1) + t318 * qJ(2);
t120 = qJD(1) * t277 - t168;
t260 = t318 * t176;
t198 = pkin(3) * t260 + t319 * t162;
t164 = qJD(1) * t168;
t268 = t319 * pkin(2);
t222 = -t179 * t268 + t164;
t189 = (-(-t268 + t198) * qJD(1) - t120) * qJD(1) + t222;
t263 = t95 * pkin(4) - t96 * pkin(7);
t212 = t119 * t273 - t178 * t95;
t210 = t212 * rSges(6,1) + t96 * rSges(6,3);
t213 = t119 * t272 + t177 * t95;
t29 = rSges(6,2) * t213 + t210;
t12 = (t263 - t29) * t175 + (t118 * t131 - t239 * t96) * qJD(5) + t189;
t261 = t319 * t176;
t221 = pkin(3) * t261 - t151;
t107 = t266 + t221;
t170 = t319 * qJ(2);
t267 = t318 * pkin(1);
t145 = t267 - t170;
t196 = t167 + (-t266 + t107 - t145) * qJD(1);
t274 = qJD(5) * t239;
t113 = t118 * rSges(6,3);
t285 = rSges(6,1) * t287 + t113;
t288 = t119 * t177;
t72 = rSges(6,2) * t288 - t285;
t99 = -t119 * pkin(4) - t118 * pkin(7);
t329 = t118 * t274 - t175 * (t72 + t99);
t21 = t196 + t329;
t246 = t277 + t198;
t195 = qJD(1) * t246 - t168;
t249 = -t118 * pkin(4) + pkin(7) * t119;
t257 = t119 * t274;
t112 = t119 * rSges(6,3);
t290 = t118 * t178;
t283 = -rSges(6,1) * t290 + t112;
t291 = t118 * t177;
t73 = rSges(6,2) * t291 + t283;
t22 = t257 + (t249 + t73) * t175 + t195;
t337 = (t177 * (t11 * t118 - t119 * t12 - t21 * t95 - t22 * t96) + (t118 * t22 - t119 * t21) * t272) * rSges(6,2);
t242 = t119 * rSges(5,1) - t118 * rSges(5,2);
t336 = t175 * t242;
t226 = Icges(6,5) * t177 + Icges(6,6) * t178;
t80 = t226 * t118;
t79 = t226 * t119;
t335 = -t268 - t277;
t330 = qJD(1) * t145 - t167 + t279;
t334 = 0.2e1 * qJD(5);
t228 = Icges(6,2) * t178 + t300;
t230 = Icges(6,1) * t177 + t299;
t225 = -t177 * t228 + t178 * t230;
t44 = t225 * t119 + t80;
t40 = t44 * t175;
t45 = t118 * t225 - t79;
t41 = t45 * t175;
t126 = -t301 * t319 - t260;
t244 = t318 * t301;
t127 = -t244 + t261;
t197 = t127 * rSges(4,1) - t126 * rSges(4,2) - t266;
t333 = -t145 + t197;
t332 = -t126 * rSges(4,1) - t127 * rSges(4,2) - t335;
t205 = t319 * rSges(3,1) + t318 * rSges(3,3);
t331 = t277 + t205;
t327 = pkin(2) * t252 + t153 + t330 + (-t107 - t267 - t151) * qJD(1);
t308 = t228 * t119 - t67;
t310 = -t230 * t119 - t63;
t326 = t177 * t308 + t178 * t310;
t313 = t118 * t58 + t287 * t66;
t312 = t290 * t67 - t340;
t311 = t290 * t66 - t345;
t309 = -t230 * t118 - t62;
t307 = t228 * t118 - t66;
t286 = t227 * t175;
t121 = t126 * qJD(1);
t122 = -qJD(1) * t244 + t245;
t282 = t122 * rSges(4,1) - t121 * rSges(4,2);
t281 = -t228 + t231;
t280 = -t229 - t230;
t270 = -t319 / 0.2e1;
t269 = t318 / 0.2e1;
t265 = t318 * rSges(3,1);
t251 = t276 / 0.2e1;
t250 = -t275 / 0.2e1;
t51 = t96 * rSges(5,1) - t95 * rSges(5,2);
t50 = -t95 * rSges(5,1) - t96 * rSges(5,2);
t243 = t121 * rSges(4,1) + t122 * rSges(4,2);
t241 = rSges(5,1) * t118 + rSges(5,2) * t119;
t237 = -t118 * t21 - t119 * t22;
t236 = t118 * t73 + t119 * t72;
t224 = -t99 + t285;
t223 = -t249 - t283;
t15 = -t288 * t63 - t314;
t16 = -t288 * t62 + t313;
t208 = (-t118 * t15 + t119 * t16) * qJD(5);
t17 = -t291 * t63 + t312;
t18 = -t291 * t62 + t311;
t207 = (-t118 * t17 + t119 * t18) * qJD(5);
t206 = -t265 - t267;
t199 = t177 * t307 + t178 * t309;
t193 = -t211 - t262;
t192 = t210 - t263;
t190 = t118 * t30 + t119 * t29 + t72 * t95 - t73 * t96;
t188 = (t177 * t280 + t178 * t281) * t175;
t187 = t221 - t145;
t129 = t229 * qJD(5);
t130 = t231 * qJD(5);
t182 = -t129 * t177 + t130 * t178 + (-t177 * t230 - t178 * t228) * qJD(5);
t10 = qJD(5) * t232 - t177 * (Icges(6,1) * t214 + Icges(6,4) * t215 + Icges(6,5) * t95) - t178 * (Icges(6,4) * t214 + Icges(6,2) * t215 + Icges(6,6) * t95);
t128 = t227 * qJD(5);
t13 = t118 * t128 + t119 * t182 + t225 * t95 - t226 * t96;
t14 = t118 * t182 - t119 * t128 - t225 * t96 - t226 * t95;
t5 = t40 + t208;
t6 = t41 + t207;
t9 = -qJD(5) * t349 - t177 * (Icges(6,1) * t212 + Icges(6,4) * t213 + Icges(6,5) * t96) - t178 * (Icges(6,4) * t212 + Icges(6,2) * t213 + Icges(6,6) * t96);
t180 = t175 * (-t130 * t177 + t228 * t273 + (-qJD(5) * t230 - t129) * t178) + t5 * t275 / 0.2e1 + (t14 + t10) * t250 + (t13 + t6 + t9) * t251 - ((t45 - t37) * t95 + (t44 - t36) * t96) * qJD(5) / 0.2e1;
t172 = t319 * rSges(3,3);
t166 = rSges(3,3) * t253;
t94 = -qJD(1) * t120 - t179 * t205 + t164;
t93 = qJD(1) * (-rSges(3,1) * t252 + t166) + t284;
t86 = t239 * t118;
t85 = t239 * t119;
t77 = qJD(1) * t333 + t167;
t71 = t119 * t240 + t113;
t70 = t118 * t240 - t112;
t57 = (-t120 + t243) * qJD(1) + t222;
t56 = qJD(1) * t282 + t209;
t43 = -t175 * t241 + t195;
t42 = t196 + t336;
t35 = -t175 * t50 + t189;
t34 = t175 * t51 + t201;
t31 = qJD(5) * t236 - qJD(3);
t24 = Icges(6,5) * t214 + Icges(6,6) * t215 + Icges(6,3) * t95;
t23 = Icges(6,5) * t212 + Icges(6,6) * t213 + Icges(6,3) * t96;
t8 = (-t118 * t63 - t119 * t62) * t177 + t312 + t313;
t7 = t190 * qJD(5);
t1 = [(-t40 + ((t17 + (-t232 + t59) * t119) * t119 + (t338 * t118 + t18 - t311 - t345) * t118) * qJD(5)) * t250 + (t41 + ((t119 * t349 + t15 + t311) * t119 + (-t119 * t338 + t16 - t8) * t118) * qJD(5)) * t251 - t180 + t347 + (t12 * (t187 + t224) + t11 * (-t223 + t246) + t337 - t31 * (t71 + t72) * t276 + (-t192 - t195) * t21 + (-t193 + t21 - (t71 - t99) * t175 - t239 * t276 + t327) * t22) * m(6) + (t35 * (t187 + t242) + t34 * (-t241 + t246) + (-t195 - t50) * t42 + (t327 + t42 + t51 - t336) * t43) * m(5) + (t57 * t333 + t56 * t332 + (t335 * qJD(1) + t168 + t243) * t77 + (t77 + t282 + (-t197 - t267 - t266) * qJD(1) + t330) * (t332 * qJD(1) - t168)) * m(4) + (t94 * (t170 + t172 + t206) + t93 * t331 + (t166 + (t265 - t172 + t206) * qJD(1) + t330) * (qJD(1) * t331 - t168)) * m(3); 0.2e1 * (t11 * t270 + t12 * t269) * m(6) + 0.2e1 * (t269 * t35 + t270 * t34) * m(5) + 0.2e1 * (t269 * t57 + t270 * t56) * m(4) + 0.2e1 * (t269 * t94 + t270 * t93) * m(3); -m(6) * t7; (t40 + ((-t17 + t8) * t119 + (t118 * t232 - t18 + t348) * t118) * qJD(5)) * t250 + (-t41 + ((-t15 - t348) * t119 + (-t16 + (-t349 + t58) * t118 - t340) * t118) * qJD(5)) * t251 + t180 + t347 + (-t12 * t224 + t11 * t223 - t337 - t31 * (t70 + t73) * t275 + (t193 + t329) * t22 + (t192 - t257 - (-t70 + t249) * t175) * t21) * m(6) + (-(-t241 * t42 - t242 * t43) * t175 + t34 * t241 - t35 * t242 + t42 * t50 - t43 * t51) * m(5); t175 * (t10 * t119 - t118 * t9 - t36 * t96 - t37 * t95) / 0.2e1 + ((t80 * t275 - t286) * t119 + (t188 + (-t326 * t118 + (-t79 + t199) * t119) * qJD(5)) * t118) * t250 + ((t79 * t276 + t286) * t118 + (t188 + (t199 * t119 + (-t326 - t80) * t118) * qJD(5)) * t119) * t251 + ((t177 * t281 - t178 * t280) * t175 + ((t118 * t308 - t119 * t307) * t178 + (-t118 * t310 + t119 * t309) * t177) * qJD(5)) * t320 - (t13 * t175 + (-(-t118 * t23 - t349 * t95 - t59 * t96) * t118 + t119 * (-t118 * t24 + t232 * t95 - t58 * t96) + t15 * t96 + t16 * t95) * t334) * t118 / 0.2e1 + (t14 * t175 + (-t118 * (t119 * t23 + t349 * t96 - t59 * t95) + t119 * (t119 * t24 - t232 * t96 - t58 * t95) + t17 * t96 + t18 * t95) * t334) * t119 / 0.2e1 + (t6 + t207) * t95 / 0.2e1 + (t5 + t208) * t96 / 0.2e1 + (t7 * t236 + t31 * t190 - t237 * t131 - (-t11 * t119 - t12 * t118 + t21 * t96 - t22 * t95) * t239 - (-t21 * t85 + t22 * t86) * t175 - (t31 * (t118 * t86 + t119 * t85) - t237 * t240) * qJD(5)) * m(6);];
tauc = t1(:);
