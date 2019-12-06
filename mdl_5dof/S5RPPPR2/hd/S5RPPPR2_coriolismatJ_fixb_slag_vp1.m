% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% Datum: 2019-12-05 17:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPPPR2_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR2_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR2_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPPR2_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:31:00
% EndTime: 2019-12-05 17:31:11
% DurationCPUTime: 4.66s
% Computational Cost: add. (11592->269), mult. (29923->407), div. (0->0), fcn. (37402->10), ass. (0->179)
t215 = sin(qJ(1));
t217 = cos(qJ(1));
t271 = cos(pkin(8));
t272 = cos(pkin(7));
t235 = t272 * t271;
t269 = sin(pkin(8));
t204 = t215 * t269 + t217 * t235;
t213 = sin(pkin(7));
t270 = cos(pkin(9));
t240 = t213 * t270;
t268 = sin(pkin(9));
t179 = t204 * t268 - t217 * t240;
t238 = t213 * t268;
t180 = t204 * t270 + t217 * t238;
t234 = t272 * t269;
t203 = -t215 * t271 + t217 * t234;
t214 = sin(qJ(5));
t216 = cos(qJ(5));
t148 = t180 * t216 + t203 * t214;
t149 = t180 * t214 - t203 * t216;
t233 = -t148 * rSges(6,1) + t149 * rSges(6,2);
t237 = t272 * pkin(2) + pkin(1);
t229 = qJ(3) * t213 + t237;
t259 = t215 * qJ(2);
t222 = -t204 * pkin(3) - t229 * t217 - t259;
t220 = -t203 * qJ(4) + t222;
t339 = -t180 * pkin(4) + t220;
t72 = (-rSges(6,3) - pkin(6)) * t179 + t233 + t339;
t109 = t179 * rSges(6,3) - t233;
t73 = -t179 * pkin(6) - t109 + t339;
t334 = (t72 - t73) * m(6);
t266 = Icges(6,4) * t148;
t103 = Icges(6,2) * t149 - Icges(6,6) * t179 - t266;
t140 = Icges(6,4) * t149;
t106 = -Icges(6,1) * t148 - Icges(6,5) * t179 + t140;
t202 = t215 * t234 + t217 * t271;
t224 = -t215 * t235 + t217 * t269;
t218 = -t215 * t238 + t224 * t270;
t145 = -t202 * t216 - t218 * t214;
t146 = -t202 * t214 + t218 * t216;
t350 = t103 * t145 + t106 * t146;
t349 = -t103 * t149 + t106 * t148;
t178 = t215 * t240 + t224 * t268;
t99 = Icges(6,5) * t148 - Icges(6,6) * t149 + Icges(6,3) * t179;
t348 = t178 * t99;
t200 = t271 * t238 + t272 * t270;
t201 = t271 * t240 - t272 * t268;
t239 = t213 * t269;
t176 = -t201 * t214 + t216 * t239;
t177 = t201 * t216 + t214 * t239;
t228 = rSges(6,1) * t177 + rSges(6,2) * t176 + rSges(6,3) * t200;
t346 = -t200 * t109 + t179 * t228;
t35 = t348 - t350;
t345 = t35 - t348;
t342 = -t179 * t99 + t349;
t166 = t215 * t202 + t203 * t217;
t242 = -m(5) / 0.4e1 - m(6) / 0.4e1;
t320 = 0.2e1 * t242;
t341 = t166 * t320;
t336 = m(5) + m(6);
t340 = t336 * t166;
t337 = t215 ^ 2;
t282 = t179 / 0.2e1;
t267 = Icges(6,4) * t146;
t101 = Icges(6,2) * t145 + Icges(6,6) * t178 + t267;
t139 = Icges(6,4) * t145;
t104 = Icges(6,1) * t146 + Icges(6,5) * t178 + t139;
t98 = Icges(6,5) * t146 + Icges(6,6) * t145 + Icges(6,3) * t178;
t34 = t145 * t101 + t146 * t104 + t178 * t98;
t333 = t34 - t342;
t313 = (t202 * t217 - t203 * t215) * t213;
t323 = t336 * t313;
t122 = t145 * rSges(6,1) - t146 * rSges(6,2);
t123 = -rSges(6,1) * t149 - t148 * rSges(6,2);
t265 = Icges(6,4) * t177;
t127 = Icges(6,2) * t176 + Icges(6,6) * t200 + t265;
t136 = Icges(6,1) * t176 - t265;
t134 = Icges(6,5) * t176 - Icges(6,6) * t177;
t264 = t200 * t134;
t108 = t146 * rSges(6,1) + t145 * rSges(6,2) + t178 * rSges(6,3);
t211 = t217 * qJ(2);
t221 = -t224 * pkin(3) + t202 * qJ(4) + t229 * t215 - t211;
t71 = -t218 * pkin(4) - t178 * pkin(6) - t108 + t221;
t322 = (t136 / 0.2e1 - t127 / 0.2e1) * t177 + m(6) * (-t122 * t71 - t123 * t72) + t264 / 0.2e1;
t321 = -t180 * rSges(5,1) + t179 * rSges(5,2);
t319 = t178 / 0.2e1;
t318 = -t179 / 0.2e1;
t113 = -t203 * rSges(5,3) + t220 + t321;
t114 = (-rSges(5,3) - qJ(4)) * t203 + t222 + t321;
t316 = m(5) * (t113 - t114);
t315 = m(6) * t213;
t169 = Icges(6,4) * t176;
t128 = Icges(6,1) * t177 + Icges(6,5) * t200 + t169;
t135 = -Icges(6,2) * t177 + t169;
t251 = t128 + t135;
t312 = t251 * t176;
t309 = (rSges(4,3) + qJ(3)) * t213 + t237;
t308 = 0.2e1 * qJD(1) * (-t316 / 0.2e1 - t334 / 0.2e1);
t205 = (t217 ^ 2 + t337) * t213;
t307 = 0.2e1 * t205;
t305 = 0.4e1 * qJD(1);
t131 = -t224 * rSges(4,1) - t202 * rSges(4,2) + t309 * t215 - t211;
t219 = -t204 * rSges(4,1) + t203 * rSges(4,2) - t309 * t217 - t259;
t302 = m(4) * (t131 * t215 - t217 * t219) * t213;
t112 = -t218 * rSges(5,1) + t178 * rSges(5,2) + t202 * rSges(5,3) + t221;
t301 = m(5) * (t112 * t215 - t113 * t217) * t213;
t300 = m(5) * (-t217 * t112 - t113 * t215);
t67 = t200 * t108 - t178 * t228;
t296 = m(6) * (-t67 * t202 - t346 * t203);
t295 = (-t215 * t67 - t217 * t346) * t315;
t294 = m(6) * (-t215 * t346 + t67 * t217);
t293 = (t215 * t71 - t217 * t72) * t315;
t292 = m(6) * (-t72 * t215 - t217 * t71);
t291 = m(6) * (t122 * t203 + t123 * t202);
t290 = (t122 * t217 + t123 * t215) * t315;
t289 = m(6) * (t215 * t122 - t123 * t217);
t279 = m(3) * (-t217 * (-t217 * rSges(3,3) - t211) + (rSges(3,3) + qJ(2)) * t337);
t278 = m(4) * (-t217 * t131 - t215 * t219);
t275 = m(6) * qJD(5);
t257 = Icges(6,1) * t145 - t101 - t267;
t256 = -Icges(6,1) * t149 + t103 - t266;
t255 = -Icges(6,2) * t146 + t104 + t139;
t254 = -Icges(6,2) * t148 - t106 - t140;
t252 = -t127 + t136;
t249 = -t323 / 0.2e1;
t247 = t323 / 0.2e1;
t246 = t340 / 0.2e1;
t245 = -t340 / 0.2e1;
t124 = (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * t307;
t243 = t124 * qJD(1);
t241 = t34 * t318 + (t345 + t350) * t319 + (t333 + t342) * t282;
t65 = t72 * t203;
t95 = t113 * t203;
t26 = m(5) * (t202 * t112 - t95) + m(6) * (t202 * t71 - t65);
t116 = Icges(6,5) * t145 - Icges(6,6) * t146;
t117 = -Icges(6,5) * t149 - Icges(6,6) * t148;
t232 = t178 * t116 + t179 * t117;
t50 = (Icges(6,5) * t177 + Icges(6,6) * t176 + Icges(6,3) * t200) * t178 + t127 * t145 + t128 * t146;
t49 = t50 * t200;
t227 = (t178 * t34 + t179 * t35 + t49) * t282;
t226 = t257 * t178 + t256 * t179;
t225 = t255 * t178 + t254 * t179;
t137 = t176 * rSges(6,1) - t177 * rSges(6,2);
t125 = (m(4) / 0.4e1 - t242) * t307 - (m(4) + t336) * t205 / 0.2e1;
t84 = -t123 * t200 + t137 * t179;
t83 = t122 * t200 - t137 * t178;
t81 = t289 / 0.2e1;
t79 = t290 / 0.2e1;
t77 = t291 / 0.2e1;
t76 = t245 + t341;
t75 = t245 + t246;
t74 = t246 - t341;
t63 = -t122 * t179 + t123 * t178;
t58 = -0.2e1 * t242 * t313 + t247;
t57 = t247 + t249;
t56 = t313 * t320 + t249;
t45 = (t252 * t177 + t264 + t312) * t200;
t41 = t294 / 0.2e1;
t38 = t295 / 0.2e1;
t32 = t296 / 0.2e1;
t31 = t179 * t134 + t252 * t148 - t149 * t251;
t30 = t178 * t134 + t251 * t145 + t252 * t146;
t28 = t200 * t117 + t254 * t176 + t256 * t177;
t27 = t200 * t116 + t255 * t176 + t257 * t177;
t25 = t293 + t301 + t302;
t24 = t278 + t279 + t292 + t300;
t21 = (t128 / 0.2e1 + t135 / 0.2e1) * t176 + t322;
t13 = t41 - t289 / 0.2e1;
t12 = t81 + t41;
t11 = t81 - t294 / 0.2e1;
t9 = t38 - t290 / 0.2e1;
t8 = t79 + t38;
t7 = t79 - t295 / 0.2e1;
t6 = t32 - t291 / 0.2e1;
t5 = t77 + t32;
t4 = t77 - t296 / 0.2e1;
t2 = t49 + (t333 + t349) * t178 + t345 * t179;
t1 = t241 * t178 + t2 * t282 - t227;
t3 = [(t71 * t334 / 0.4e1 + t112 * t316 / 0.4e1) * t305 + t24 * qJD(2) + t25 * qJD(3) + t26 * qJD(4) + t21 * qJD(5), qJD(1) * t24 + qJD(3) * t125 + qJD(4) * t75 + qJD(5) * t12, qJD(1) * t25 + qJD(2) * t125 + qJD(4) * t57 + qJD(5) * t8, t26 * qJD(1) + t75 * qJD(2) + t57 * qJD(3) + t5 * qJD(5), t21 * qJD(1) + t12 * qJD(2) + t8 * qJD(3) + t5 * qJD(4) + (t45 + m(6) * (t122 * t67 - t123 * t346 - t71 * t83 + t72 * t84) + t227 + (t28 / 0.2e1 + t31 / 0.2e1 - t2 / 0.2e1) * t179 + (t27 / 0.2e1 + t30 / 0.2e1 - t241) * t178) * qJD(5); -t124 * qJD(3) + t76 * qJD(4) + t11 * qJD(5) + (-t279 / 0.4e1 - t278 / 0.4e1 - t300 / 0.4e1 - t292 / 0.4e1) * t305 + t215 * t308, 0, -t243, t76 * qJD(1), t11 * qJD(1) + (t83 * t215 + t217 * t84) * t275; t124 * qJD(2) + t56 * qJD(4) + t7 * qJD(5) + (-t302 / 0.4e1 - t301 / 0.4e1 - t293 / 0.4e1) * t305 + t213 * t217 * t308, t243, 0, t56 * qJD(1), t7 * qJD(1) + (-t63 * t272 + (-t215 * t84 + t217 * t83) * t213) * t275; (m(5) * (t114 * t203 - t95) + m(6) * (t203 * t73 - t65) - t26) * qJD(1) + t74 * qJD(2) + t58 * qJD(3) + t4 * qJD(5), t74 * qJD(1), t58 * qJD(1), 0, t4 * qJD(1) + (-t84 * t202 + t83 * t203 + t63 * t239) * t275; t13 * qJD(2) + t9 * qJD(3) + t6 * qJD(4) + t1 * qJD(5) + (-t67 * t334 - t312 / 0.2e1 + (t282 + t318) * (t101 * t176 + t104 * t177 + t200 * t98 + t50) - t322) * qJD(1), t13 * qJD(1), t9 * qJD(1), t6 * qJD(1), t1 * qJD(1) + (m(6) * (t346 * t84 + (-t108 * t179 + t109 * t178) * t63 + t67 * t83) + t200 * (t27 * t178 + t28 * t179 + t45) / 0.2e1 + (t225 * t145 + t226 * t146 + t232 * t178 + t30 * t200) * t319 + (t226 * t148 - t149 * t225 + t232 * t179 + t31 * t200) * t282) * qJD(5);];
Cq = t3;
