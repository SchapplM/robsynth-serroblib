% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPPR2_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR2_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR2_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR2_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR2_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:57:21
% EndTime: 2020-01-03 11:57:27
% DurationCPUTime: 3.27s
% Computational Cost: add. (23412->249), mult. (17995->348), div. (0->0), fcn. (18859->10), ass. (0->156)
t236 = qJ(1) + qJ(2);
t232 = pkin(8) + t236;
t229 = sin(t232);
t237 = sin(pkin(9));
t238 = cos(pkin(9));
t230 = cos(t232);
t241 = cos(qJ(5));
t239 = sin(qJ(5));
t274 = t238 * t239;
t195 = -t229 * t274 - t230 * t241;
t273 = t238 * t241;
t196 = t229 * t273 - t230 * t239;
t251 = t196 * rSges(6,1) + t195 * rSges(6,2);
t233 = sin(t236);
t304 = pkin(2) * t233;
t252 = -t230 * qJ(4) + t304;
t124 = (pkin(4) * t238 + pkin(3) + (rSges(6,3) + pkin(7)) * t237) * t229 + t251 + t252;
t303 = sin(qJ(1)) * pkin(1);
t119 = t124 + t303;
t360 = t119 - t124;
t170 = t195 * rSges(6,1) - t196 * rSges(6,2);
t197 = -t229 * t241 + t230 * t274;
t198 = t229 * t239 + t230 * t273;
t171 = t197 * rSges(6,1) + t198 * rSges(6,2);
t279 = t230 * t237;
t154 = t198 * rSges(6,1) - t197 * rSges(6,2) + rSges(6,3) * t279;
t234 = cos(t236);
t231 = pkin(2) * t234;
t260 = t230 * pkin(3) + t229 * qJ(4) + t231;
t278 = t230 * t238;
t336 = pkin(4) * t278 + pkin(7) * t279 + t154 + t260;
t54 = t124 * t170 - t171 * t336;
t359 = m(6) * t54;
t235 = cos(qJ(1)) * pkin(1);
t341 = t235 + t336;
t113 = t341 * t229;
t161 = rSges(5,1) * t278 - rSges(5,2) * t279 + t229 * rSges(5,3) + t260;
t158 = t161 + t235;
t143 = t158 * t229;
t156 = t161 * t229;
t280 = t229 * t237;
t160 = -rSges(5,2) * t280 - t230 * rSges(5,3) + (rSges(5,1) * t238 + pkin(3)) * t229 + t252;
t157 = t160 + t303;
t330 = m(6) / 0.2e1;
t343 = t229 * t336;
t350 = m(5) / 0.2e1;
t300 = (t113 + t343 + (-t119 - t124) * t230) * t330 + (t143 + t156 + (-t157 - t160) * t230) * t350;
t301 = (-t360 * t230 + t113 - t343) * t330 + (t143 - t156 + (-t157 + t160) * t230) * t350;
t340 = t300 - t301;
t358 = t340 * qJD(1);
t51 = t119 * t170 - t171 * t341;
t357 = m(6) * t51;
t145 = Icges(6,5) * t198 - Icges(6,6) * t197 + Icges(6,3) * t279;
t186 = Icges(6,4) * t198;
t149 = Icges(6,2) * t197 - Icges(6,6) * t279 - t186;
t185 = Icges(6,4) * t197;
t151 = Icges(6,1) * t198 + Icges(6,5) * t279 - t185;
t249 = t197 * t149 + t198 * t151;
t64 = t145 * t279 + t249;
t356 = t230 * t64;
t290 = Icges(6,4) * t239;
t203 = -Icges(6,5) * t238 + (Icges(6,1) * t241 - t290) * t237;
t266 = t203 + (-Icges(6,2) * t241 - t290) * t237;
t354 = t266 * t239;
t43 = t119 * t336 - t124 * t341;
t204 = -t238 * rSges(6,3) + (rSges(6,1) * t241 - rSges(6,2) * t239) * t237;
t353 = -t238 * t154 - t204 * t279;
t352 = -t195 * t149 + t196 * t151;
t327 = m(5) * (t157 * t161 - t160 * t158);
t312 = m(3) * (-t235 * (t233 * rSges(3,1) + t234 * rSges(3,2)) + t303 * (t234 * rSges(3,1) - t233 * rSges(3,2)));
t310 = m(4) * (-t235 * (t229 * rSges(4,1) + t230 * rSges(4,2) + t304) + t303 * (t230 * rSges(4,1) - t229 * rSges(4,2) + t231));
t201 = -Icges(6,3) * t238 + (Icges(6,5) * t241 - Icges(6,6) * t239) * t237;
t289 = Icges(6,4) * t241;
t202 = -Icges(6,6) * t238 + (-Icges(6,2) * t239 + t289) * t237;
t347 = (-t197 * t202 + t198 * t203 + t201 * t279) * t238;
t265 = qJD(1) + qJD(2);
t207 = (-Icges(6,5) * t239 - Icges(6,6) * t241) * t237;
t275 = t238 * t207;
t339 = -t275 / 0.2e1 - t237 * t354 / 0.2e1;
t291 = Icges(6,4) * t196;
t147 = Icges(6,2) * t195 + Icges(6,6) * t280 + t291;
t184 = Icges(6,4) * t195;
t150 = Icges(6,1) * t196 + Icges(6,5) * t280 + t184;
t272 = -t197 * t147 + t198 * t150;
t335 = 0.4e1 * qJD(1);
t153 = rSges(6,3) * t280 + t251;
t128 = t238 * t153 + t204 * t280;
t325 = m(6) * (t360 * t353 + (-t336 + t341) * t128);
t322 = m(6) * (t54 + t51);
t318 = m(6) * t43;
t315 = m(6) * (-t119 * t230 + t113);
t314 = m(6) * (-t124 * t230 + t343);
t313 = m(6) * (t128 * t230 + t353 * t229);
t308 = m(5) * (-t157 * t230 + t143);
t307 = m(5) * (-t160 * t230 + t156);
t305 = m(6) * (-t229 * t170 + t230 * t171);
t112 = t305 / 0.2e1;
t73 = t313 / 0.2e1;
t25 = t112 + t73;
t302 = t25 * qJD(4);
t296 = -t128 * t170 - t171 * t353;
t294 = m(6) * qJD(5);
t288 = (t195 * t202 + t196 * t203 + t201 * t280) * t238;
t209 = (-Icges(6,1) * t239 - t289) * t237;
t267 = t202 - t209;
t287 = (-t275 + (-t267 * t241 - t354) * t237) * t238;
t144 = Icges(6,5) * t196 + Icges(6,6) * t195 + Icges(6,3) * t280;
t286 = t144 * t230;
t285 = t145 * t229;
t276 = t237 * t241;
t271 = -Icges(6,1) * t195 + t147 + t291;
t270 = -Icges(6,1) * t197 + t149 - t186;
t269 = -Icges(6,2) * t196 + t150 + t184;
t268 = Icges(6,2) * t198 - t151 + t185;
t61 = t144 * t280 + t195 * t147 + t196 * t150;
t62 = -t145 * t280 - t352;
t63 = -t144 * t279 - t272;
t11 = -t288 + ((-t249 + t61 + t64) * t229 + (t63 + (-t285 + t286) * t237 - t62 + t272) * t230) * t237;
t12 = -t347 + (t356 + (t272 + t62 + (t285 + t286) * t237 + t352) * t229) * t237;
t33 = -t288 + (t229 * t61 - t230 * t62) * t237;
t34 = t347 + (t229 * t63 - t356) * t237;
t2 = ((t33 / 0.2e1 - t11 / 0.2e1) * t230 + (t12 / 0.2e1 + t34 / 0.2e1) * t229) * t237;
t26 = t73 - t305 / 0.2e1;
t264 = t26 * qJD(4) + t2 * qJD(5);
t258 = t280 / 0.2e1;
t257 = -t279 / 0.2e1;
t254 = -t276 / 0.2e1;
t253 = t276 / 0.2e1;
t250 = t202 * t254 + t209 * t253 + t339;
t248 = t322 / 0.2e1 + t250;
t244 = t202 * t253 + t209 * t254 - t339;
t164 = Icges(6,5) * t195 - Icges(6,6) * t196;
t55 = -t238 * t164 + (-t269 * t239 - t271 * t241) * t237;
t165 = Icges(6,5) * t197 + Icges(6,6) * t198;
t56 = -t238 * t165 + (-t268 * t239 - t270 * t241) * t237;
t85 = t266 * t195 - t267 * t196 + t207 * t280;
t86 = t266 * t197 + t267 * t198 - t207 * t279;
t243 = t11 * t279 / 0.2e1 - t287 - (t12 + t34) * t280 / 0.2e1 + (t55 + t85) * t258 + (t33 + t56 + t86) * t257;
t210 = (-rSges(6,1) * t239 - rSges(6,2) * t241) * t237;
t139 = t238 * t171 - t210 * t279;
t138 = -t238 * t170 - t210 * t280;
t108 = (t170 * t230 + t171 * t229) * t237;
t50 = t307 + t314;
t49 = t308 + t315;
t42 = t250 + t359;
t40 = t250 + t357;
t24 = t112 - t313 / 0.2e1;
t21 = t25 * qJD(5);
t20 = t24 * qJD(5);
t14 = t325 / 0.2e1;
t13 = t310 + t312 + t318 + t327;
t8 = t300 + t301;
t5 = -t325 / 0.2e1 + t248;
t4 = t14 + t248;
t3 = t14 - t322 / 0.2e1 + t244;
t1 = [t13 * qJD(2) + t49 * qJD(4) + t40 * qJD(5), t13 * qJD(1) + t8 * qJD(4) + t4 * qJD(5) + 0.2e1 * (t312 / 0.2e1 + t310 / 0.2e1 + t327 / 0.2e1 + t43 * t330) * qJD(2), 0, t49 * qJD(1) + t8 * qJD(2) + t21, t40 * qJD(1) + t4 * qJD(2) + (m(6) * (t138 * t119 + t139 * t341 + t296) + t243) * qJD(5) + t302; t340 * qJD(4) + t5 * qJD(5) + (-t312 / 0.4e1 - t310 / 0.4e1 - t327 / 0.4e1 - t318 / 0.4e1) * t335, t50 * qJD(4) + t42 * qJD(5), 0, t50 * qJD(2) + t21 + t358, t5 * qJD(1) + t42 * qJD(2) + (m(6) * (t138 * t124 + t139 * t336 + t296) + t243) * qJD(5) + t302; 0, 0, 0, 0, t108 * t294; -t340 * qJD(2) + t20 + (-t308 / 0.4e1 - t315 / 0.4e1) * t335, -t358 + t20 + 0.4e1 * (-t314 / 0.4e1 - t307 / 0.4e1) * qJD(2), 0, 0, (-t138 * t229 - t139 * t230) * t294 + t265 * t24; (t244 - t357) * qJD(1) + t3 * qJD(2) + t264, t3 * qJD(1) + (t244 - t359) * qJD(2) + t264, 0, t265 * t26, (m(6) * ((t153 * t230 - t154 * t229) * t237 * t108 - t128 * t138 + t139 * t353) - t238 * (-t287 + (t229 * t55 - t230 * t56) * t237) / 0.2e1 + (-t85 * t238 + (t164 * t280 + t195 * t269 - t196 * t271) * t280 - (t165 * t280 + t195 * t268 - t196 * t270) * t279) * t258 + (-t86 * t238 + (-t164 * t279 + t197 * t269 + t198 * t271) * t280 - (-t165 * t279 + t197 * t268 + t198 * t270) * t279) * t257) * qJD(5) + t265 * t2;];
Cq = t1;
