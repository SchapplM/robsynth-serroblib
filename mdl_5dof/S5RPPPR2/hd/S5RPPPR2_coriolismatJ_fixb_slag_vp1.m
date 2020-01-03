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
% Datum: 2020-01-03 11:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:22:23
% EndTime: 2020-01-03 11:22:34
% DurationCPUTime: 4.61s
% Computational Cost: add. (11592->289), mult. (29923->427), div. (0->0), fcn. (37402->10), ass. (0->190)
t224 = cos(pkin(7));
t228 = cos(qJ(1));
t274 = cos(pkin(8));
t241 = t228 * t274;
t221 = sin(pkin(8));
t226 = sin(qJ(1));
t264 = t226 * t221;
t205 = t224 * t241 + t264;
t220 = sin(pkin(9));
t223 = cos(pkin(9));
t222 = sin(pkin(7));
t265 = t222 * t228;
t184 = -t205 * t220 + t223 * t265;
t183 = t205 * t223 + t220 * t265;
t242 = t226 * t274;
t263 = t228 * t221;
t204 = t224 * t263 - t242;
t225 = sin(qJ(5));
t227 = cos(qJ(5));
t154 = t183 * t227 + t204 * t225;
t155 = t183 * t225 - t204 * t227;
t239 = t154 * rSges(6,1) - t155 * rSges(6,2);
t250 = t228 * pkin(1) + t226 * qJ(2);
t279 = pkin(2) * t224;
t240 = qJ(3) * t265 + t228 * t279 + t250;
t238 = t205 * pkin(3) + t240;
t234 = t204 * qJ(4) + t238;
t335 = t183 * pkin(4) + t234;
t76 = (-rSges(6,3) - pkin(6)) * t184 + t239 + t335;
t116 = t184 * rSges(6,3) - t239;
t77 = -pkin(6) * t184 - t116 + t335;
t347 = m(6) * (-t76 + t77);
t272 = Icges(6,4) * t154;
t108 = -Icges(6,2) * t155 - Icges(6,6) * t184 + t272;
t148 = Icges(6,4) * t155;
t111 = Icges(6,1) * t154 - Icges(6,5) * t184 - t148;
t203 = t224 * t242 - t263;
t266 = t222 * t226;
t181 = t203 * t223 + t220 * t266;
t202 = t224 * t264 + t241;
t151 = -t181 * t225 + t202 * t227;
t152 = t181 * t227 + t202 * t225;
t346 = t151 * t108 + t152 * t111;
t345 = t155 * t108 - t111 * t154;
t105 = Icges(6,5) * t154 - Icges(6,6) * t155 - Icges(6,3) * t184;
t180 = t203 * t220 - t223 * t266;
t344 = t180 * t105;
t243 = t222 * t274;
t200 = t220 * t243 + t224 * t223;
t201 = -t224 * t220 + t223 * t243;
t267 = t221 * t222;
t178 = -t201 * t225 + t227 * t267;
t179 = t201 * t227 + t225 * t267;
t235 = t179 * rSges(6,1) + t178 * rSges(6,2) + t200 * rSges(6,3);
t342 = -t200 * t116 + t184 * t235;
t35 = -t344 - t346;
t341 = t35 + t344;
t338 = t105 * t184 + t345;
t163 = -t202 * t265 + t204 * t266;
t248 = -m(5) / 0.4e1 - m(6) / 0.4e1;
t318 = 0.2e1 * t248;
t337 = t163 * t318;
t333 = m(5) + m(6);
t336 = t333 * t163;
t284 = t184 / 0.2e1;
t104 = Icges(6,5) * t152 + Icges(6,6) * t151 + Icges(6,3) * t180;
t273 = Icges(6,4) * t152;
t107 = Icges(6,2) * t151 + Icges(6,6) * t180 + t273;
t147 = Icges(6,4) * t151;
t110 = Icges(6,1) * t152 + Icges(6,5) * t180 + t147;
t34 = t180 * t104 + t151 * t107 + t152 * t110;
t331 = -t338 + t34;
t172 = t226 * t202 + t228 * t204;
t321 = t333 * t172;
t128 = t151 * rSges(6,1) - t152 * rSges(6,2);
t129 = t155 * rSges(6,1) + rSges(6,2) * t154;
t271 = Icges(6,4) * t179;
t133 = Icges(6,2) * t178 + Icges(6,6) * t200 + t271;
t144 = Icges(6,1) * t178 - t271;
t142 = Icges(6,5) * t178 - Icges(6,6) * t179;
t270 = t200 * t142;
t114 = t152 * rSges(6,1) + t151 * rSges(6,2) + t180 * rSges(6,3);
t217 = t228 * qJ(2);
t244 = pkin(1) + t279;
t229 = t202 * qJ(4) + t203 * pkin(3) + (qJ(3) * t222 + t244) * t226 - t217;
t75 = t181 * pkin(4) + t180 * pkin(6) + t114 + t229;
t320 = (t144 / 0.2e1 - t133 / 0.2e1) * t179 + m(6) * (t75 * t128 - t76 * t129) + t270 / 0.2e1;
t319 = t183 * rSges(5,1) + t184 * rSges(5,2);
t317 = t180 / 0.2e1;
t316 = -t184 / 0.2e1;
t175 = Icges(6,4) * t178;
t134 = Icges(6,1) * t179 + Icges(6,5) * t200 + t175;
t143 = -Icges(6,2) * t179 + t175;
t257 = t134 + t143;
t313 = t257 * t178;
t206 = (t226 ^ 2 + t228 ^ 2) * t222;
t311 = 0.2e1 * t206;
t310 = 2 * qJD(1);
t309 = 4 * qJD(1);
t308 = m(4) / 0.2e1;
t306 = m(5) / 0.2e1;
t305 = m(6) / 0.2e1;
t118 = t181 * rSges(5,1) - t180 * rSges(5,2) + t202 * rSges(5,3) + t229;
t119 = t204 * rSges(5,3) + t234 + t319;
t59 = t118 * t266 + t119 * t265;
t304 = m(5) * t59;
t103 = t119 * t226;
t303 = m(5) * (-t118 * t228 + t103);
t69 = t200 * t114 - t180 * t235;
t40 = t265 * t342 + t69 * t266;
t301 = m(6) * ((-t226 * t69 - t228 * t342) * t222 + t40);
t299 = m(6) * (t69 * t202 + t342 * t204);
t298 = m(6) * t40;
t297 = m(6) * (t342 * t226 - t69 * t228);
t43 = t76 * t265 + t75 * t266;
t296 = m(6) * t43;
t74 = t76 * t226;
t295 = m(6) * (-t75 * t228 + t74);
t294 = m(6) * (-t204 * t128 - t202 * t129);
t293 = m(6) * (-t128 * t228 - t129 * t226) * t222;
t292 = m(6) * (-t226 * t128 + t228 * t129);
t282 = m(3) * (-(-rSges(3,2) * t266 - t228 * rSges(3,3) + t226 * pkin(1) - t217) * t228 + (-rSges(3,2) * t265 + t226 * rSges(3,3) + t250) * t226);
t139 = t203 * rSges(4,1) - t202 * rSges(4,2) - t217 + ((rSges(4,3) + qJ(3)) * t222 + t244) * t226;
t230 = t205 * rSges(4,1) - t204 * rSges(4,2) + rSges(4,3) * t265 + t240;
t100 = t139 * t266 + t230 * t265;
t281 = m(4) * t100;
t280 = m(4) * (-t139 * t228 + t230 * t226);
t276 = m(6) * qJD(5);
t262 = Icges(6,1) * t151 - t107 - t273;
t261 = Icges(6,1) * t155 + t108 + t272;
t260 = -Icges(6,2) * t152 + t110 + t147;
t259 = Icges(6,2) * t154 - t111 + t148;
t258 = -t133 + t144;
t256 = t336 / 0.2e1;
t254 = -t336 / 0.2e1;
t253 = t321 / 0.2e1;
t252 = -t321 / 0.2e1;
t130 = (t308 + t306 + t305) * t311;
t249 = t130 * qJD(1);
t247 = t34 * t316 + (t341 + t346) * t317 + (t331 + t338) * t284;
t65 = t76 * t204;
t99 = t119 * t204;
t26 = m(5) * (t118 * t202 + t99) + m(6) * (t75 * t202 + t65);
t122 = Icges(6,5) * t151 - Icges(6,6) * t152;
t123 = Icges(6,5) * t155 + Icges(6,6) * t154;
t237 = t180 * t122 + t184 * t123;
t50 = t180 * (Icges(6,5) * t179 + Icges(6,6) * t178 + Icges(6,3) * t200) + t151 * t133 + t152 * t134;
t49 = t50 * t200;
t233 = (t34 * t180 + t35 * t184 + t49) * t284;
t232 = t262 * t180 + t261 * t184;
t231 = t260 * t180 + t259 * t184;
t145 = t178 * rSges(6,1) - t179 * rSges(6,2);
t131 = (m(4) / 0.4e1 - t248) * t311 - (m(4) + t333) * t206 / 0.2e1;
t120 = (rSges(5,3) + qJ(4)) * t204 + t238 + t319;
t88 = -t200 * t129 + t184 * t145;
t87 = t200 * t128 - t180 * t145;
t85 = t292 / 0.2e1;
t83 = t293 / 0.2e1;
t81 = t294 / 0.2e1;
t80 = t172 * t318 + t252;
t79 = t252 + t253;
t78 = -0.2e1 * t248 * t172 + t253;
t63 = -t184 * t128 + t180 * t129;
t58 = t254 + t337;
t57 = t254 + t256;
t56 = t256 - t337;
t45 = (t258 * t179 + t270 + t313) * t200;
t41 = t297 / 0.2e1;
t38 = t298 / 0.2e1;
t32 = t299 / 0.2e1;
t31 = t184 * t142 - t154 * t258 + t257 * t155;
t30 = t180 * t142 + t257 * t151 + t258 * t152;
t28 = t200 * t123 + t259 * t178 + t261 * t179;
t27 = t200 * t122 + t260 * t178 + t262 * t179;
t25 = t281 + t296 + t304;
t24 = t280 + t282 + t295 + t303;
t21 = (t134 / 0.2e1 + t143 / 0.2e1) * t178 + t320;
t19 = t301 / 0.2e1;
t13 = t41 - t292 / 0.2e1;
t12 = t85 + t41;
t11 = t85 - t297 / 0.2e1;
t9 = t38 + t19 - t293 / 0.2e1;
t8 = t83 + t38 - t301 / 0.2e1;
t7 = t83 + t19 - t298 / 0.2e1;
t6 = t32 - t294 / 0.2e1;
t5 = t81 + t32;
t4 = t81 - t299 / 0.2e1;
t2 = t49 + (t331 + t345) * t180 + t341 * t184;
t1 = t247 * t180 + t2 * t284 - t233;
t3 = [(t75 * t347 / 0.4e1 + m(5) * (-t119 + t120) * t118 / 0.4e1) * t309 + t24 * qJD(2) + t25 * qJD(3) + t26 * qJD(4) + t21 * qJD(5), t24 * qJD(1) + t131 * qJD(3) + t79 * qJD(4) + t12 * qJD(5), t25 * qJD(1) + t131 * qJD(2) + t57 * qJD(4) + t8 * qJD(5), t26 * qJD(1) + t79 * qJD(2) + t57 * qJD(3) + t5 * qJD(5), t21 * qJD(1) + t12 * qJD(2) + t8 * qJD(3) + t5 * qJD(4) + (t45 + m(6) * (t69 * t128 - t129 * t342 + t87 * t75 + t88 * t76) + t233 + (t28 / 0.2e1 + t31 / 0.2e1 - t2 / 0.2e1) * t184 + (t27 / 0.2e1 + t30 / 0.2e1 - t247) * t180) * qJD(5); -t130 * qJD(3) + t80 * qJD(4) + t11 * qJD(5) + ((-t226 * t120 + t103) * t306 + (-t226 * t77 + t74) * t305) * t310 + (-t282 / 0.4e1 - t280 / 0.4e1 - t303 / 0.4e1 - t295 / 0.4e1) * t309, 0, -t249, t80 * qJD(1), t11 * qJD(1) + (-t87 * t226 - t88 * t228) * t276; t130 * qJD(2) + t56 * qJD(4) + t7 * qJD(5) + (-t281 / 0.4e1 - t304 / 0.4e1 - t296 / 0.4e1) * t309 + (((-t139 * t226 - t228 * t230) * t222 + t100) * t308 + ((-t118 * t226 - t120 * t228) * t222 + t59) * t306 + ((-t226 * t75 - t228 * t77) * t222 + t43) * t305) * t310, t249, 0, t56 * qJD(1), t7 * qJD(1) + (-t63 * t224 + (t226 * t88 - t228 * t87) * t222) * t276; (m(5) * (-t204 * t120 + t99) + m(6) * (-t204 * t77 + t65) - t26) * qJD(1) + t78 * qJD(2) + t58 * qJD(3) + t4 * qJD(5), t78 * qJD(1), t58 * qJD(1), 0, t4 * qJD(1) + (t88 * t202 - t87 * t204 + t63 * t267) * t276; t13 * qJD(2) + t9 * qJD(3) + t6 * qJD(4) + t1 * qJD(5) + (t69 * t347 - t313 / 0.2e1 + (t316 + t284) * (t200 * t104 + t178 * t107 + t179 * t110 + t50) - t320) * qJD(1), t13 * qJD(1), t9 * qJD(1), t6 * qJD(1), t1 * qJD(1) + (m(6) * (t342 * t88 + (-t184 * t114 + t180 * t116) * t63 + t69 * t87) + t200 * (t27 * t180 + t28 * t184 + t45) / 0.2e1 + (t231 * t151 + t232 * t152 + t237 * t180 + t30 * t200) * t317 + (-t154 * t232 + t231 * t155 + t237 * t184 + t31 * t200) * t284) * qJD(5);];
Cq = t3;
