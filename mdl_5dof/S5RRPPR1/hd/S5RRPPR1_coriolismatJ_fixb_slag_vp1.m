% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPPR1
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
% Datum: 2019-12-05 18:19
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPPR1_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR1_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR1_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR1_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR1_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:18:10
% EndTime: 2019-12-05 18:18:16
% DurationCPUTime: 1.85s
% Computational Cost: add. (14256->202), mult. (7992->260), div. (0->0), fcn. (6878->10), ass. (0->140)
t193 = qJ(1) + qJ(2);
t189 = pkin(8) + t193;
t184 = sin(t189);
t185 = cos(t189);
t195 = -pkin(7) - qJ(4);
t192 = pkin(9) + qJ(5);
t187 = sin(t192);
t234 = t185 * t187;
t212 = rSges(6,2) * t234 - t184 * rSges(6,3);
t194 = cos(pkin(9));
t188 = cos(t192);
t253 = rSges(6,1) * t188;
t213 = pkin(4) * t194 + pkin(3) + t253;
t191 = cos(t193);
t258 = pkin(2) * t191;
t105 = t184 * t195 - t213 * t185 + t212 - t258;
t260 = cos(qJ(1)) * pkin(1);
t100 = t105 - t260;
t190 = sin(t193);
t259 = pkin(2) * t190;
t306 = rSges(5,2) * sin(pkin(9)) - rSges(5,1) * t194 - pkin(3);
t307 = -rSges(5,3) - qJ(4);
t111 = -t306 * t184 + t307 * t185 + t259;
t261 = sin(qJ(1)) * pkin(1);
t108 = t111 + t261;
t103 = t185 * t108;
t112 = t307 * t184 + t306 * t185 - t258;
t109 = t112 - t260;
t235 = t185 * t111;
t238 = t184 * t187;
t219 = -rSges(6,2) * t238 - t185 * rSges(6,3);
t104 = t213 * t184 + t185 * t195 + t219 + t259;
t236 = t185 * t104;
t308 = m(6) / 0.2e1;
t309 = m(5) / 0.2e1;
t99 = t104 + t261;
t95 = t185 * t99;
t256 = (-t95 - t236 + (-t100 - t105) * t184) * t308 + (-t103 - t235 + (-t109 - t112) * t184) * t309;
t228 = t100 - t105;
t257 = (-t228 * t184 + t236 - t95) * t308 + (-t103 + t235 + (-t109 + t112) * t184) * t309;
t3 = t257 - t256;
t310 = t3 * qJD(1);
t284 = m(4) * (-t260 * (rSges(4,1) * t184 + rSges(4,2) * t185 + t259) - (-t185 * rSges(4,1) + t184 * rSges(4,2) - t258) * t261);
t280 = m(5) * (-t112 * t108 + t109 * t111);
t272 = m(6) * (t100 * t104 - t105 * t99);
t263 = m(3) * (-t260 * (rSges(3,1) * t190 + rSges(3,2) * t191) - (-t191 * rSges(3,1) + t190 * rSges(3,2)) * t261);
t163 = rSges(6,1) * t187 + rSges(6,2) * t188;
t144 = t163 * t184;
t145 = t163 * t185;
t203 = t184 * t144 + t185 * t145;
t201 = m(6) * t203;
t246 = t104 * t144;
t41 = t105 * t145 - t246;
t305 = m(6) * t41;
t164 = -rSges(6,2) * t187 + t253;
t304 = m(6) * t164;
t218 = qJD(1) + qJD(2);
t292 = t185 ^ 2;
t293 = t184 ^ 2;
t301 = t163 * (t292 + t293);
t200 = -m(6) * t301 / 0.2e1;
t68 = -t201 / 0.2e1 + t200;
t303 = t218 * t68;
t302 = -t100 * t184 - t95;
t233 = t185 * t188;
t128 = Icges(6,4) * t233 - Icges(6,2) * t234 + t184 * Icges(6,6);
t170 = Icges(6,4) * t234;
t130 = Icges(6,1) * t233 + t184 * Icges(6,5) - t170;
t300 = (t128 * t187 - t130 * t188) * t185;
t182 = Icges(6,4) * t188;
t299 = Icges(6,1) * t187 + t182;
t298 = t105 * t184 + t236;
t297 = Icges(6,2) * t187 - t182;
t249 = Icges(6,4) * t187;
t159 = Icges(6,2) * t188 + t249;
t162 = Icges(6,1) * t188 - t249;
t296 = (-t297 + t299) * t187 + (t159 - t162) * t188;
t222 = -Icges(6,2) * t233 + t130 - t170;
t224 = t185 * t299 + t128;
t295 = t222 * t187 + t224 * t188;
t169 = Icges(6,4) * t238;
t237 = t184 * t188;
t129 = -Icges(6,1) * t237 + Icges(6,5) * t185 + t169;
t223 = Icges(6,2) * t237 + t129 + t169;
t127 = Icges(6,6) * t185 + t297 * t184;
t225 = -t184 * t299 + t127;
t294 = -t223 * t187 - t225 * t188;
t210 = (-t297 / 0.2e1 + t299 / 0.2e1) * t188 + (-t159 / 0.2e1 + t162 / 0.2e1) * t187;
t291 = 0.4e1 * qJD(1);
t278 = m(5) * (-t109 * t184 - t103);
t277 = m(5) * (-t112 * t184 - t235);
t251 = t99 * t144;
t40 = t100 * t145 - t251;
t276 = m(6) * (t41 + t40);
t275 = m(6) * (t228 * t145 + t246 - t251);
t270 = m(6) * t40;
t268 = m(6) * t302;
t267 = m(6) * t298;
t266 = t184 / 0.2e1;
t265 = -t185 / 0.2e1;
t264 = t185 / 0.2e1;
t69 = t201 / 0.2e1 + t200;
t255 = t69 * qJD(4);
t242 = t127 * t187;
t158 = Icges(6,5) * t188 - Icges(6,6) * t187;
t240 = t158 * t184;
t125 = Icges(6,3) * t185 - t240;
t227 = t185 * t125 + t127 * t238;
t226 = t184 * t125 + t129 * t233;
t126 = Icges(6,5) * t233 - Icges(6,6) * t234 + Icges(6,3) * t184;
t211 = -t129 * t188 - t126;
t51 = t185 * t126 + t128 * t238 - t130 * t237;
t52 = -t127 * t234 + t226;
t53 = t126 * t184 - t300;
t11 = (t227 + t53 + t300) * t185 + (-t52 + (t211 - t242) * t185 + t51 + t226) * t184;
t50 = -t129 * t237 + t227;
t12 = (t51 + (-t126 + t242) * t185 - t226) * t185 + (t211 * t184 + t227 - t50) * t184;
t27 = t184 * t51 + t185 * t50;
t28 = t184 * t53 + t185 * t52;
t2 = (t12 / 0.2e1 + t28 / 0.2e1) * t185 + (-t27 / 0.2e1 + t11 / 0.2e1) * t184;
t217 = t68 * qJD(4) + t2 * qJD(5);
t208 = t276 / 0.2e1 + t210;
t206 = Icges(6,5) * t187 + Icges(6,6) * t188;
t138 = t184 * t206;
t199 = -t184 * t11 / 0.2e1 + (t12 + t28) * t265 + (t158 * t185 + t296 * t184 - t225 * t187 + t223 * t188) * t264 + (-t296 * t185 - t224 * t187 + t222 * t188 + t240 + t27) * t266;
t198 = -t210 + (t264 + t265) * (t188 * t128 + t187 * t130);
t139 = t206 * t185;
t65 = t68 * qJD(5);
t63 = t69 * qJD(5);
t31 = -t267 + t277;
t30 = t210 + t305;
t29 = t210 + t270;
t26 = t268 + t278;
t18 = t275 / 0.2e1;
t13 = t263 + t272 + t280 + t284;
t8 = -t275 / 0.2e1 + t208;
t7 = t18 + t208;
t6 = t18 - t276 / 0.2e1 + t198;
t5 = t256 + t257;
t1 = [t13 * qJD(2) + t26 * qJD(4) + t29 * qJD(5), t13 * qJD(1) + t5 * qJD(4) + t7 * qJD(5) + 0.2e1 * (t272 / 0.2e1 + t263 / 0.2e1 + t280 / 0.2e1 + t284 / 0.2e1) * qJD(2), 0, qJD(1) * t26 + qJD(2) * t5 + t63, t29 * qJD(1) + t7 * qJD(2) + (-t302 * t304 + t199) * qJD(5) + t255; -t3 * qJD(4) + t8 * qJD(5) + (-t263 / 0.4e1 - t284 / 0.4e1 - t280 / 0.4e1 - t272 / 0.4e1) * t291, qJD(4) * t31 + qJD(5) * t30, 0, qJD(2) * t31 - t310 + t63, t8 * qJD(1) + t30 * qJD(2) + (t298 * t304 + t199) * qJD(5) + t255; 0, 0, 0, 0, -qJD(5) * t201; t3 * qJD(2) - t65 + (-t278 / 0.4e1 - t268 / 0.4e1) * t291, t310 - t65 + 0.4e1 * (t267 / 0.4e1 - t277 / 0.4e1) * qJD(2), 0, 0, -t303; (t198 - t270) * qJD(1) + t6 * qJD(2) + t217, t6 * qJD(1) + (t198 - t305) * qJD(2) + t217, 0, t303, (m(6) * (t164 * t301 - t203 * (t185 * (rSges(6,1) * t233 - t212) - t184 * (-rSges(6,1) * t237 - t219))) + (t292 * t138 + (t295 * t184 + (-t139 - t294) * t185) * t184) * t264 + (-t293 * t139 + (t294 * t185 + (t138 - t295) * t184) * t185) * t266) * qJD(5) + t218 * t2;];
Cq = t1;
