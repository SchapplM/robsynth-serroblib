% Calculate matrix of centrifugal and coriolis load on the joints for
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
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RPRR5_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR5_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR5_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR5_coriolismatJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR5_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR5_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR5_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:51:33
% EndTime: 2019-12-31 16:51:37
% DurationCPUTime: 2.89s
% Computational Cost: add. (6986->242), mult. (14628->325), div. (0->0), fcn. (17082->6), ass. (0->171)
t287 = sin(qJ(3));
t288 = sin(qJ(1));
t289 = cos(qJ(3));
t290 = cos(qJ(1));
t173 = t290 * t287 - t288 * t289;
t172 = -t288 * t287 - t290 * t289;
t189 = sin(qJ(4));
t190 = cos(qJ(4));
t265 = Icges(5,4) * t190;
t209 = -Icges(5,2) * t189 + t265;
t119 = Icges(5,6) * t172 + t209 * t173;
t266 = Icges(5,4) * t189;
t211 = Icges(5,1) * t190 - t266;
t123 = Icges(5,5) * t172 + t211 * t173;
t345 = t119 * t190 + t189 * t123;
t268 = t345 * t173;
t118 = -Icges(5,6) * t173 + t209 * t172;
t122 = -Icges(5,5) * t173 + t211 * t172;
t346 = t118 * t190 + t189 * t122;
t359 = t346 * t172;
t360 = t268 / 0.2e1 - t359 / 0.2e1;
t358 = -t173 / 0.4e1;
t181 = -t189 * rSges(5,1) - rSges(5,2) * t190;
t141 = t173 * t181;
t142 = t181 * t172;
t354 = m(5) * (t288 * t141 + t142 * t290);
t230 = -t354 / 0.2e1;
t84 = t354 / 0.2e1;
t312 = m(5) / 0.2e1;
t314 = m(4) / 0.2e1;
t145 = t172 * rSges(4,1) + t173 * rSges(4,2);
t212 = -t173 * rSges(4,1) + t172 * rSges(4,2);
t324 = t145 * t288 + t290 * t212;
t169 = t172 * pkin(6);
t269 = t189 * rSges(5,2);
t158 = t173 * t269;
t168 = t172 * rSges(5,3);
t223 = -t158 + t168;
t273 = rSges(5,1) * t190;
t225 = pkin(3) + t273;
t90 = -t225 * t173 - t169 - t223;
t86 = t290 * t90;
t341 = t173 * rSges(5,3) + t172 * t269;
t88 = -t173 * pkin(6) + t225 * t172 - t341;
t331 = t88 * t288 + t86;
t278 = t331 * t312 + t324 * t314;
t353 = t84 + t230;
t356 = qJD(4) * t353;
t157 = t173 * t273;
t234 = t157 + t168;
t89 = (-pkin(3) + t269) * t173 - t169 - t234;
t275 = t89 - t90;
t332 = m(5) * t275;
t299 = t88 * t332;
t355 = t299 * qJD(3);
t246 = t353 * qJD(2);
t252 = t118 * t189;
t253 = t119 * t189;
t249 = t122 * t190;
t250 = t123 * t190;
t207 = Icges(5,5) * t190 - Icges(5,6) * t189;
t116 = -Icges(5,3) * t172 - t207 * t173;
t352 = t173 * t116;
t117 = Icges(5,3) * t173 - t207 * t172;
t351 = t173 * t117;
t215 = -t288 * pkin(1) + t290 * qJ(2);
t198 = -t288 * pkin(2) + t215;
t73 = t198 - t90;
t199 = t290 * pkin(1) + t288 * qJ(2);
t193 = t290 * pkin(2) + t199;
t74 = t193 - t88;
t348 = -t73 * t88 + t74 * t90;
t245 = t172 * t116 - t173 * t250;
t326 = -t117 - t253;
t347 = -t326 * t173 + t245;
t344 = -t189 / 0.2e1;
t333 = -t190 / 0.2e1;
t206 = Icges(5,5) * t189 + Icges(5,6) * t190;
t136 = t206 * t172;
t343 = t206 * t173;
t342 = 0.2e1 * t312 * qJD(1);
t279 = -m(5) * t331 / 0.2e1 - m(4) * t324 / 0.2e1;
t280 = (-t89 * t290 + t86) * t312;
t340 = t279 - t280;
t338 = 0.2e1 * t230;
t133 = -t145 + t193;
t191 = t198 - t212;
t337 = t133 * t212 - t191 * t145;
t262 = Icges(5,2) * t190;
t208 = t262 + t266;
t200 = t189 * t208;
t210 = Icges(5,1) * t189 + t265;
t247 = t190 * t210;
t325 = t200 - t247;
t93 = t325 * t172 + t343;
t336 = -t172 / 0.2e1;
t335 = t173 / 0.2e1;
t334 = -t207 / 0.2e1;
t282 = m(5) * t181;
t329 = t172 * t173;
t327 = -t116 + t252;
t94 = t325 * t173 - t136;
t270 = t173 * t94;
t272 = t172 * t93;
t323 = t270 / 0.4e1 + t272 / 0.4e1;
t322 = t359 / 0.4e1 + t141 * t332 / 0.2e1;
t321 = t268 / 0.4e1 - ((-t73 - t90) * t173 + (t74 + t88) * t172) * t282 / 0.2e1;
t219 = t200 / 0.2e1 + t211 * t344 + t209 * t333 - t247 / 0.2e1;
t320 = t172 ^ 2;
t319 = t173 ^ 2;
t317 = 0.4e1 * qJD(1);
t316 = 0.2e1 * qJD(3);
t35 = -t173 * t253 - t245;
t244 = -t172 * t117 + t173 * t249;
t36 = -t173 * t252 + t244;
t19 = -t172 * t35 + t173 * t36;
t311 = -t19 / 0.2e1;
t243 = t172 * t250 + t352;
t37 = -t172 * t253 + t243;
t242 = t172 * t249 + t351;
t38 = -t172 * t252 + t242;
t20 = -t172 * t37 + t173 * t38;
t310 = -t20 / 0.2e1;
t309 = m(4) * t337;
t306 = m(4) * (t133 * t288 + t191 * t290);
t304 = m(5) * ((t74 - t88) * t142 + (-t73 + t90) * t141);
t301 = m(5) * t348;
t298 = m(5) * (t141 * t73 - t142 * t74);
t297 = m(5) * (-t141 * t90 + t142 * t88);
t70 = t73 * t290;
t296 = m(5) * (t74 * t288 + t70);
t286 = m(3) * ((t290 * rSges(3,3) + t215) * t290 + (t288 * rSges(3,3) + t199) * t288);
t75 = t198 - t89;
t277 = t73 - t75;
t274 = m(5) * qJD(4);
t256 = t141 * t181;
t255 = t142 * t181;
t233 = qJD(4) * t172;
t232 = qJD(4) * t173;
t64 = t359 / 0.2e1;
t192 = (-t211 + t262) * t190 + (t209 + t210 + t265) * t189;
t218 = t311 + (-t210 * t172 - t118) * t189 / 0.2e1 + (t208 * t172 - t122) * t333 + t192 * t336 + t173 * t334;
t217 = t310 + (-t210 * t173 - t119) * t344 + (-t208 * t173 + t123) * t333 + t172 * t334 + t192 * t335;
t205 = t219 + t360;
t204 = -t268 / 0.2e1 + t64 + t219;
t202 = -t270 / 0.4e1 - t321 + t323 - (-t346 + t93) * t172 / 0.4e1;
t201 = -t272 / 0.4e1 - t322 + t323 + (-t345 + t94) * t358;
t183 = t269 - t273;
t126 = -t223 - t157;
t56 = t126 * t173 + (-t172 * t273 + t341) * t172;
t40 = 0.2e1 * t84;
t32 = -t219 + t297;
t30 = -t219 + t298;
t28 = t56 * (t126 - t158 + t234) * t172;
t25 = t286 + t296 + t306;
t23 = (-t118 * t173 - t119 * t172) * t189 + t243 + t244;
t14 = t304 / 0.2e1;
t13 = t301 + t309;
t12 = t278 - t340;
t11 = t278 + t340;
t10 = t279 + t280 - t278;
t9 = (t35 + (-t250 + t253) * t173 + t242) * t173 + (-t173 * t327 - t23 + t36) * t172;
t8 = (t23 - t37) * t173 + (-t38 + (t249 - t252) * t172 + t347) * t172;
t7 = (t37 + (t327 - t249) * t173) * t173 + (t327 * t172 - t242 + t351 + t38) * t172;
t6 = (-t347 - t35) * t173 + (-t36 + (t326 + t250) * t172 + t352) * t172;
t5 = t201 - t304 / 0.2e1 + t202 - t219;
t4 = t345 * t358 + t14 + t202 + t205 + t322;
t3 = t201 + t14 - t359 / 0.4e1 + t204 + t321;
t2 = (t8 / 0.2e1 + t311) * t173 + (t310 - t6 / 0.2e1) * t172;
t1 = m(5) * t28 + (t7 / 0.2e1 + t19 / 0.2e1) * t173 + (t20 / 0.2e1 - t9 / 0.2e1) * t172;
t15 = [-m(5) * t277 * t74 * t317 / 0.4e1 + t25 * qJD(2) + t13 * qJD(3) + t30 * qJD(4), qJD(1) * t25 + qJD(3) * t11 + t356, t13 * qJD(1) + t11 * qJD(2) + t4 * qJD(4) - t355 + (-t348 * t312 - t337 * t314) * t316, t30 * qJD(1) + t246 + t4 * qJD(3) - t28 * t274 + (m(5) * (-t183 * t74 + t255) - t7 / 0.2e1 + t218) * t232 + (m(5) * (-t183 * t73 - t256) + t9 / 0.2e1 + t217) * t233; t10 * qJD(3) + t40 * qJD(4) + (-t286 / 0.4e1 - t306 / 0.4e1 - t296 / 0.4e1) * t317 + (-t290 * t75 + t70) * t342, 0, t10 * qJD(1) + t338 * qJD(4) + t278 * t316, t40 * qJD(1) + t338 * qJD(3) + (-t288 * t172 + t290 * t173) * t183 * t274; t12 * qJD(2) + t355 + t3 * qJD(4) + (-t309 / 0.4e1 - t301 / 0.4e1) * t317 + (t275 * t74 - t277 * t88) * t342, qJD(1) * t12 - t356, qJD(1) * t299 + t32 * qJD(4), t3 * qJD(1) - t246 + t32 * qJD(3) + (m(5) * (-t183 * t88 - t255) - t8 / 0.2e1 - t218) * t232 + (m(5) * (-t183 * t90 + t256) + t6 / 0.2e1 - t217) * t233; t338 * qJD(2) + t5 * qJD(3) + t1 * qJD(4) + (t205 + t64 - t298 + (-t345 / 0.2e1 + t277 * t282) * t173) * qJD(1), qJD(1) * t338 + qJD(3) * t353, t5 * qJD(1) + t2 * qJD(4) + t246 + (t204 - t297 + t360) * qJD(3), t1 * qJD(1) + t2 * qJD(3) + (m(5) * (t56 * (-t141 * t173 - t142 * t172) + (t319 + t320) * t183 * t181) + (t319 * t136 - t329 * t343) * t335 + (-t136 * t329 + t320 * t343) * t336) * qJD(4);];
Cq = t15;
