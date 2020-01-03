% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2020-01-03 11:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPPRR3_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR3_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR3_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR3_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR3_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:27:29
% EndTime: 2020-01-03 11:27:44
% DurationCPUTime: 3.75s
% Computational Cost: add. (18915->284), mult. (11629->392), div. (0->0), fcn. (10604->9), ass. (0->191)
t212 = pkin(9) + qJ(4);
t204 = sin(t212);
t206 = cos(t212);
t177 = rSges(5,1) * t204 + rSges(5,2) * t206;
t213 = qJ(1) + pkin(8);
t205 = sin(t213);
t202 = t205 ^ 2;
t207 = cos(t213);
t203 = t207 ^ 2;
t248 = t202 + t203;
t337 = t248 * t177;
t346 = -m(5) / 0.2e1;
t327 = m(6) / 0.2e1;
t208 = qJ(5) + t212;
t199 = sin(t208);
t200 = cos(t208);
t170 = rSges(6,1) * t199 + rSges(6,2) * t200;
t306 = pkin(4) * t204;
t244 = t170 + t306;
t342 = t244 * t207;
t343 = t244 * t205;
t352 = (t205 * t343 + t207 * t342) * t327;
t293 = t337 * t346 - t352;
t160 = t177 * t205;
t161 = t177 * t207;
t95 = -t160 * t205 - t161 * t207;
t302 = t95 * t346 + t352;
t25 = t302 - t293;
t353 = t25 * qJD(1);
t289 = Icges(6,4) * t199;
t169 = Icges(6,1) * t200 - t289;
t127 = -Icges(6,5) * t207 + t169 * t205;
t273 = t200 * t205;
t102 = t127 * t273;
t272 = t200 * t207;
t276 = t199 * t207;
t124 = Icges(6,5) * t272 - Icges(6,6) * t276 + Icges(6,3) * t205;
t181 = Icges(6,4) * t276;
t128 = Icges(6,1) * t272 + Icges(6,5) * t205 - t181;
t126 = Icges(6,4) * t272 - Icges(6,2) * t276 + Icges(6,6) * t205;
t285 = t126 * t199;
t228 = t128 * t200 - t285;
t351 = -t205 * t124 - t228 * t207 - t102;
t165 = Icges(6,5) * t200 - Icges(6,6) * t199;
t283 = t165 * t205;
t123 = -Icges(6,3) * t207 + t283;
t350 = t205 * t123 + t127 * t272;
t349 = t170 * t248;
t193 = Icges(6,4) * t200;
t167 = -Icges(6,2) * t199 + t193;
t336 = Icges(6,1) * t199 + t193;
t348 = t167 + t336;
t312 = -t205 / 0.2e1;
t344 = t205 / 0.2e1;
t310 = -t207 / 0.2e1;
t290 = Icges(5,4) * t204;
t173 = Icges(5,2) * t206 + t290;
t176 = Icges(5,1) * t206 - t290;
t341 = (t176 / 0.2e1 - t173 / 0.2e1) * t204;
t340 = m(6) * qJD(5);
t144 = t170 * t205;
t145 = t170 * t207;
t240 = t205 * t144 + t207 * t145;
t299 = rSges(6,1) * t200;
t171 = -rSges(6,2) * t199 + t299;
t280 = t171 * t207;
t281 = t171 * t205;
t230 = Icges(6,5) * t199 + Icges(6,6) * t200;
t138 = t205 * t230;
t139 = t230 * t207;
t257 = Icges(6,2) * t272 - t128 + t181;
t166 = Icges(6,2) * t200 + t289;
t258 = -t166 * t205 + t127;
t259 = t207 * t336 + t126;
t125 = -Icges(6,6) * t207 + t167 * t205;
t260 = -t205 * t336 - t125;
t333 = t199 * (-t205 * t257 - t207 * t258) + t200 * (t205 * t259 + t207 * t260);
t303 = (-t203 * t138 + (t207 * t139 - t333) * t205) * t310 + (t202 * t139 + (-t205 * t138 + t333) * t207) * t312;
t277 = t199 * t205;
t250 = -rSges(6,2) * t277 - t207 * rSges(6,3);
t114 = t205 * (rSges(6,1) * t273 + t250);
t236 = rSges(6,1) * t272 - rSges(6,2) * t276;
t129 = rSges(6,3) * t205 + t236;
t201 = cos(pkin(9)) * pkin(3) + pkin(2);
t305 = pkin(4) * t206;
t184 = t201 + t305;
t164 = t207 * t184;
t48 = t114 + (t184 - t201) * t202 + (-t201 * t207 + t129 + t164) * t207;
t6 = t303 + m(6) * (-t240 * t48 + t280 * t342 + t281 * t343);
t339 = t6 * qJD(5);
t338 = rSges(4,3) + qJ(3);
t195 = Icges(5,4) * t206;
t174 = -Icges(5,2) * t204 + t195;
t335 = Icges(5,1) * t204 + t195;
t332 = (t166 - t169) * t200 + t348 * t199;
t268 = t204 * t207;
t187 = Icges(5,4) * t268;
t264 = t206 * t207;
t136 = Icges(5,1) * t264 + Icges(5,5) * t205 - t187;
t253 = Icges(5,2) * t264 - t136 + t187;
t135 = -Icges(5,5) * t207 + t176 * t205;
t254 = -t173 * t205 + t135;
t134 = Icges(5,4) * t264 - Icges(5,2) * t268 + Icges(5,6) * t205;
t255 = t207 * t335 + t134;
t133 = -Icges(5,6) * t207 + t174 * t205;
t256 = -t205 * t335 - t133;
t331 = t204 * (-t205 * t253 - t207 * t254) + t206 * (t205 * t255 + t207 * t256);
t239 = t348 * t200 / 0.2e1 + (t169 / 0.2e1 - t166 / 0.2e1) * t199;
t103 = t128 * t273;
t284 = t127 * t200;
t286 = t125 * t199;
t57 = -t207 * t123 - t125 * t277 + t102;
t58 = t207 * t124 + t126 * t277 - t103;
t242 = ((t103 + (t123 - t285) * t205 - t350) * t205 + ((t123 + t228) * t207 + (t284 + t286) * t205 + t351) * t207) * t312 + (-t205 * t58 - t207 * t57) * t344 + ((-t58 + (t124 - t284) * t207 + t350) * t207 + (t57 + (t124 + t286) * t205 + t351) * t205) * t310;
t330 = 0.4e1 * qJD(1);
t329 = 2 * qJD(4);
t328 = m(5) / 0.2e1;
t210 = cos(qJ(1)) * pkin(1);
t307 = sin(qJ(1)) * pkin(1);
t326 = m(4) * (-(-t338 * t207 + t307) * t207 + (t338 * t205 + t210) * t205);
t215 = -pkin(6) - qJ(3);
t300 = rSges(5,1) * t206;
t241 = t201 + t300;
t269 = t204 * t205;
t249 = -rSges(5,2) * t269 - t207 * rSges(5,3);
t98 = t205 * t241 + t207 * t215 + t249 + t307;
t189 = rSges(5,2) * t268;
t99 = -t189 + t210 + t241 * t207 + (rSges(5,3) - t215) * t205;
t325 = m(5) * (-t160 * t98 - t161 * t99);
t61 = t205 * t99 - t207 * t98;
t324 = m(5) * t61;
t211 = -pkin(7) + t215;
t91 = t307 + t207 * t211 + (t184 + t299) * t205 + t250;
t92 = t164 + t210 + (rSges(6,3) - t211) * t205 + t236;
t237 = t91 * t280 - t281 * t92;
t319 = m(6) * (-t144 * t342 + t145 * t343 + t237);
t318 = m(6) * ((t205 * t342 - t207 * t343) * t170 + t237);
t317 = m(6) * (-t342 * t92 - t343 * t91);
t316 = m(6) * (-t144 * t91 - t145 * t92);
t315 = m(6) * (t205 * t92 - t207 * t91);
t309 = t207 / 0.2e1;
t301 = m(6) * qJD(4);
t271 = t204 * t133;
t270 = t204 * t134;
t266 = t206 * t135;
t265 = t206 * t136;
t226 = -m(6) * t349 / 0.2e1;
t235 = m(6) * t240;
t66 = -t235 / 0.2e1 + t226;
t262 = t66 * qJD(1);
t243 = t171 + t305;
t232 = Icges(5,5) * t206 - Icges(5,6) * t204;
t231 = Icges(5,5) * t204 + Icges(5,6) * t206;
t227 = t265 - t270;
t224 = -t242 + (t199 * t259 + t200 * t257 + t332 * t207 - t283) * t312 + (-t165 * t207 + t199 * t260 + t200 * t258 - t332 * t205) * t310;
t222 = -t239 + (t309 + t310) * (t126 * t200 + t128 * t199);
t178 = -rSges(5,2) * t204 + t300;
t155 = t231 * t207;
t154 = t231 * t205;
t132 = Icges(5,5) * t264 - Icges(5,6) * t268 + Icges(5,3) * t205;
t131 = -Icges(5,3) * t207 + t205 * t232;
t122 = t243 * t207;
t120 = t243 * t205;
t110 = t133 * t268;
t109 = t205 * t265;
t108 = t205 * t266;
t88 = t240 * t340;
t84 = t207 * t129 + t114;
t72 = -t248 * t306 - t240;
t67 = t235 / 0.2e1 + t226;
t65 = t205 * t132 + t227 * t207;
t64 = -t205 * t131 - t135 * t264 + t110;
t63 = t207 * t132 + t134 * t269 - t109;
t62 = -t207 * t131 - t133 * t269 + t108;
t42 = -t205 * t65 - t207 * t64;
t41 = -t205 * t63 - t207 * t62;
t32 = t239 + t316;
t29 = t318 / 0.2e1;
t28 = t319 / 0.2e1;
t26 = t293 + t302;
t23 = t315 + t324 + t326;
t19 = (t335 / 0.2e1 + t174 / 0.2e1) * t206 + t341 + t325 + t317 + t239;
t16 = (t110 - t63 + (t132 - t266) * t207) * t207 + (-t108 + t62 + (t132 + t271) * t205) * t205;
t15 = (t64 + t109 - t110 + (t131 - t270) * t205) * t205 + (-t108 - t65 + (t131 + t227) * t207 + (t266 + t271) * t205) * t207;
t8 = m(6) * (t171 * t349 - t84 * t240) + t303;
t7 = t8 * qJD(5);
t4 = t28 - t318 / 0.2e1 + t242;
t3 = t29 - t319 / 0.2e1 + t242;
t2 = t28 + t29 + t224;
t1 = (-t16 / 0.2e1 - t42 / 0.2e1) * t207 + (t41 / 0.2e1 - t15 / 0.2e1) * t205 + t242;
t5 = [qJD(3) * t23 + qJD(4) * t19 + qJD(5) * t32, 0, qJD(1) * t23 + qJD(4) * t26 + qJD(5) * t67, t19 * qJD(1) + t26 * qJD(3) + t2 * qJD(5) + ((-t120 * t92 + t122 * t91) * t327 + (-t61 * t178 + (-t160 * t207 + t161 * t205) * t177) * t328) * t329 + (t15 * t344 + (t204 * t256 + t206 * t254) * t310 + t224 + (t203 / 0.2e1 + t202 / 0.2e1) * t232 + (t204 * t255 + t206 * t253 + t41) * t312 + (t16 + t42) * t309) * qJD(4), t32 * qJD(1) + t67 * qJD(3) + t2 * qJD(4) + t224 * qJD(5) + ((-t144 * t207 + t145 * t205) * t170 + t237) * t340; 0, 0, 0, (t327 * t72 + t328 * t95) * t329 - t88, -t240 * t301 - t88; t25 * qJD(4) - t66 * qJD(5) + (-t326 / 0.4e1 - t324 / 0.4e1 - t315 / 0.4e1) * t330, 0, 0, t353 + (t120 * t207 - t122 * t205) * t301, -t262; (t222 - (t174 + t335) * t206 / 0.2e1 - t341) * qJD(1) - t25 * qJD(3) + t1 * qJD(4) + t4 * qJD(5) + (-t325 / 0.4e1 - t317 / 0.4e1) * t330, 0, -t353, t1 * qJD(1) + (m(5) * (t178 * t337 + (t207 * (rSges(5,1) * t264 + rSges(5,3) * t205 - t189) + t205 * (t205 * t300 + t249)) * t95) + (-t203 * t154 + (t207 * t155 - t331) * t205) * t310 + (t202 * t155 + (-t205 * t154 + t331) * t207) * t312 + m(6) * (t120 * t343 + t122 * t342 + t48 * t72) + t303) * qJD(4) + t339, t4 * qJD(1) + t6 * qJD(4) + t339; (t222 - t316) * qJD(1) + t66 * qJD(3) + t3 * qJD(4) + t242 * qJD(5), 0, t262, t3 * qJD(1) + ((t72 * t84 + (t120 * t205 + t122 * t207) * t170) * m(6) + t303) * qJD(4) + t7, qJD(1) * t242 + qJD(4) * t8 + t7;];
Cq = t5;
