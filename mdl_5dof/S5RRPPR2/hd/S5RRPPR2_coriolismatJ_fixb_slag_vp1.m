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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:06
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 10:05:27
% EndTime: 2022-01-20 10:05:36
% DurationCPUTime: 3.75s
% Computational Cost: add. (23412->236), mult. (17995->337), div. (0->0), fcn. (18859->10), ass. (0->155)
t236 = qJ(1) + qJ(2);
t233 = pkin(8) + t236;
t231 = sin(t233);
t232 = cos(t233);
t234 = sin(t236);
t308 = pkin(2) * t234;
t254 = t232 * qJ(4) - t308;
t238 = cos(pkin(9));
t255 = rSges(5,1) * t238 + pkin(3);
t237 = sin(pkin(9));
t283 = t231 * t237;
t168 = rSges(5,2) * t283 + t232 * rSges(5,3) - t231 * t255 + t254;
t306 = sin(qJ(1)) * pkin(1);
t165 = t168 - t306;
t282 = t232 * t237;
t235 = cos(t236);
t307 = pkin(2) * t235;
t169 = t307 - rSges(5,2) * t282 + t255 * t232 + (rSges(5,3) + qJ(4)) * t231;
t305 = cos(qJ(1)) * pkin(1);
t166 = t169 + t305;
t103 = t165 * t232 + t166 * t231;
t110 = t168 * t232 + t169 * t231;
t334 = m(6) / 0.2e1;
t354 = m(5) / 0.2e1;
t241 = cos(qJ(5));
t239 = sin(qJ(5));
t278 = t238 * t239;
t206 = t231 * t241 - t232 * t278;
t277 = t238 * t241;
t207 = t231 * t239 + t232 * t277;
t251 = t207 * rSges(6,1) + t206 * rSges(6,2);
t341 = pkin(4) * t238 + pkin(3) + (rSges(6,3) + pkin(7)) * t237;
t128 = t231 * qJ(4) + t341 * t232 + t251 + t307;
t123 = t128 + t305;
t204 = t231 * t278 + t232 * t241;
t205 = t231 * t277 - t232 * t239;
t344 = -t205 * rSges(6,1) + t204 * rSges(6,2);
t348 = -t341 * t231 + t254 + t344;
t358 = t348 - t306;
t65 = t123 * t231 + t232 * t358;
t66 = t128 * t231 + t232 * t348;
t302 = (t66 + t65) * t334 + (t110 + t103) * t354;
t303 = (-t66 + t65) * t334 + (t103 - t110) * t354;
t6 = t303 - t302;
t366 = t6 * qJD(1);
t178 = -t204 * rSges(6,1) - t205 * rSges(6,2);
t179 = t206 * rSges(6,1) - t207 * rSges(6,2);
t54 = t128 * t179 - t178 * t348;
t365 = m(6) * t54;
t51 = t123 * t179 - t178 * t358;
t364 = m(6) * t51;
t43 = t123 * t348 - t128 * t358;
t159 = rSges(6,3) * t283 - t344;
t213 = -t238 * rSges(6,3) + (rSges(6,1) * t241 - rSges(6,2) * t239) * t237;
t362 = t238 * t159 + t213 * t283;
t292 = Icges(6,4) * t239;
t212 = -Icges(6,5) * t238 + (Icges(6,1) * t241 - t292) * t237;
t270 = t212 + (-Icges(6,2) * t241 - t292) * t237;
t359 = t270 * t239;
t194 = Icges(6,4) * t205;
t153 = -Icges(6,2) * t204 + Icges(6,6) * t283 + t194;
t193 = Icges(6,4) * t204;
t157 = -Icges(6,1) * t205 - Icges(6,5) * t283 + t193;
t357 = -t204 * t153 - t205 * t157;
t356 = t206 * t153 - t207 * t157;
t331 = m(5) * (-t169 * t165 + t166 * t168);
t315 = m(3) * (t305 * (-t234 * rSges(3,1) - t235 * rSges(3,2)) + (t235 * rSges(3,1) - t234 * rSges(3,2)) * t306);
t313 = m(4) * (t305 * (-t231 * rSges(4,1) - t232 * rSges(4,2) - t308) + (t232 * rSges(4,1) - t231 * rSges(4,2) + t307) * t306);
t210 = -Icges(6,3) * t238 + (Icges(6,5) * t241 - Icges(6,6) * t239) * t237;
t291 = Icges(6,4) * t241;
t211 = -Icges(6,6) * t238 + (-Icges(6,2) * t239 + t291) * t237;
t351 = (-t204 * t211 + t205 * t212 + t210 * t283) * t238;
t150 = Icges(6,5) * t205 - Icges(6,6) * t204 + Icges(6,3) * t283;
t350 = t150 * t282;
t269 = qJD(1) + qJD(2);
t216 = (-Icges(6,5) * t239 - Icges(6,6) * t241) * t237;
t279 = t238 * t216;
t346 = -t279 / 0.2e1 - t237 * t359 / 0.2e1;
t63 = t350 + t356;
t345 = t63 - t350;
t152 = Icges(6,5) * t207 + Icges(6,6) * t206 + Icges(6,3) * t282;
t293 = Icges(6,4) * t207;
t155 = Icges(6,2) * t206 + Icges(6,6) * t282 + t293;
t195 = Icges(6,4) * t206;
t158 = Icges(6,1) * t207 + Icges(6,5) * t282 + t195;
t342 = t152 * t283 - t204 * t155 + t205 * t158;
t340 = 0.4e1 * qJD(1);
t161 = rSges(6,3) * t282 + t251;
t134 = t238 * t161 + t213 * t282;
t329 = m(6) * ((t123 - t128) * t362 + (-t348 + t358) * t134);
t326 = m(6) * (t54 + t51);
t323 = m(6) * t43;
t320 = m(6) * t65;
t319 = m(6) * t66;
t318 = m(6) * (-t134 * t231 + t232 * t362);
t311 = m(5) * t103;
t310 = m(5) * t110;
t309 = m(6) * (-t231 * t178 - t232 * t179);
t112 = t309 / 0.2e1;
t73 = t318 / 0.2e1;
t25 = t112 + t73;
t304 = t25 * qJD(4);
t298 = -t134 * t179 - t178 * t362;
t296 = m(6) * qJD(5);
t61 = t150 * t283 + t357;
t295 = t231 * t61;
t290 = (t206 * t211 + t207 * t212 + t210 * t282) * t238;
t218 = (-Icges(6,1) * t239 - t291) * t237;
t271 = -t211 + t218;
t289 = (-t279 + (t271 * t241 - t359) * t237) * t238;
t280 = t237 * t241;
t276 = -Icges(6,1) * t204 - t153 - t194;
t275 = Icges(6,1) * t206 - t155 - t293;
t274 = -Icges(6,2) * t205 - t157 - t193;
t273 = -Icges(6,2) * t207 + t158 + t195;
t11 = t351 + ((-t342 + t345 - t356) * t232 - t295) * t237;
t64 = t152 * t282 + t206 * t155 + t207 * t158;
t12 = -t290 + ((-t357 + t61 + t64) * t232 + t345 * t231) * t237;
t33 = -t351 + (t232 * t342 + t295) * t237;
t34 = -t290 + (t231 * t63 + t232 * t64) * t237;
t2 = ((t11 / 0.2e1 + t33 / 0.2e1) * t232 + (-t34 / 0.2e1 + t12 / 0.2e1) * t231) * t237;
t26 = t73 - t309 / 0.2e1;
t268 = t26 * qJD(4) + t2 * qJD(5);
t262 = t283 / 0.2e1;
t260 = t282 / 0.2e1;
t257 = -t280 / 0.2e1;
t256 = t280 / 0.2e1;
t250 = t211 * t257 + t218 * t256 + t346;
t248 = t326 / 0.2e1 + t250;
t245 = t211 * t256 + t218 * t257 - t346;
t172 = -Icges(6,5) * t204 - Icges(6,6) * t205;
t55 = -t238 * t172 + (-t274 * t239 + t276 * t241) * t237;
t173 = Icges(6,5) * t206 - Icges(6,6) * t207;
t56 = -t238 * t173 + (-t273 * t239 + t275 * t241) * t237;
t85 = -t270 * t204 + t271 * t205 + t216 * t283;
t86 = t270 * t206 + t271 * t207 + t216 * t282;
t244 = -t12 * t283 / 0.2e1 - t289 - (t11 + t33) * t282 / 0.2e1 + (t56 + t86) * t260 + (t34 + t55 + t85) * t262;
t219 = (-rSges(6,1) * t239 - rSges(6,2) * t241) * t237;
t143 = -t238 * t179 - t219 * t282;
t142 = t238 * t178 + t219 * t283;
t108 = (t178 * t232 - t179 * t231) * t237;
t50 = t310 + t319;
t49 = t311 + t320;
t42 = t250 + t365;
t40 = t250 + t364;
t24 = t112 - t318 / 0.2e1;
t21 = t25 * qJD(5);
t20 = t24 * qJD(5);
t14 = t329 / 0.2e1;
t13 = t313 + t315 + t323 + t331;
t7 = t302 + t303;
t5 = -t329 / 0.2e1 + t248;
t4 = t14 + t248;
t3 = t14 - t326 / 0.2e1 + t245;
t1 = [t13 * qJD(2) + t49 * qJD(4) + t40 * qJD(5), t13 * qJD(1) + t7 * qJD(4) + t4 * qJD(5) + 0.2e1 * (t43 * t334 + t313 / 0.2e1 + t331 / 0.2e1 + t315 / 0.2e1) * qJD(2), 0, t49 * qJD(1) + t7 * qJD(2) + t21, t40 * qJD(1) + t4 * qJD(2) + (m(6) * (t143 * t123 + t142 * t358 + t298) + t244) * qJD(5) + t304; -t6 * qJD(4) + t5 * qJD(5) + (-t323 / 0.4e1 - t313 / 0.4e1 - t331 / 0.4e1 - t315 / 0.4e1) * t340, t50 * qJD(4) + t42 * qJD(5), 0, t50 * qJD(2) + t21 - t366, t5 * qJD(1) + t42 * qJD(2) + (m(6) * (t143 * t128 + t142 * t348 + t298) + t244) * qJD(5) + t304; 0, 0, 0, 0, t108 * t296; t6 * qJD(2) + t20 + (-t320 / 0.4e1 - t311 / 0.4e1) * t340, t366 + t20 + 0.4e1 * (-t319 / 0.4e1 - t310 / 0.4e1) * qJD(2), 0, 0, (t142 * t231 - t143 * t232) * t296 + t269 * t24; (t245 - t364) * qJD(1) + t3 * qJD(2) + t268, t3 * qJD(1) + (t245 - t365) * qJD(2) + t268, 0, t269 * t26, (m(6) * ((t159 * t232 - t161 * t231) * t237 * t108 - t134 * t143 + t142 * t362) + ((t173 * t282 + t273 * t206 + t275 * t207) * t282 + (t172 * t282 + t274 * t206 + t276 * t207) * t283 - t86 * t238) * t260 + ((t173 * t283 - t273 * t204 + t275 * t205) * t282 + (t172 * t283 - t274 * t204 + t276 * t205) * t283 - t85 * t238) * t262 - t238 * (-t289 + (t231 * t55 + t232 * t56) * t237) / 0.2e1) * qJD(5) + t269 * t2;];
Cq = t1;
