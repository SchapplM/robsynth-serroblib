% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
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
% Datum: 2019-12-31 17:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPPRP6_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP6_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP6_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP6_coriolismatJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP6_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP6_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRP6_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:54:59
% EndTime: 2019-12-31 17:55:04
% DurationCPUTime: 3.27s
% Computational Cost: add. (5677->263), mult. (7235->367), div. (0->0), fcn. (6245->6), ass. (0->172)
t210 = sin(qJ(1));
t211 = cos(qJ(1));
t205 = pkin(7) + qJ(4);
t193 = sin(t205);
t194 = cos(t205);
t283 = rSges(6,1) + pkin(4);
t329 = rSges(6,3) + qJ(5);
t333 = t329 * t193 + t283 * t194;
t338 = t333 * t211;
t342 = t210 * t338;
t87 = t338 * t211;
t262 = t194 * t211;
t264 = t193 * t211;
t114 = Icges(6,5) * t264 - Icges(6,6) * t210 - Icges(6,3) * t262;
t270 = Icges(5,4) * t193;
t225 = Icges(5,2) * t194 + t270;
t120 = -Icges(5,6) * t210 + t211 * t225;
t341 = t114 - t120;
t179 = Icges(6,5) * t262;
t122 = Icges(6,1) * t264 - Icges(6,4) * t210 - t179;
t269 = Icges(5,4) * t194;
t226 = Icges(5,1) * t193 + t269;
t124 = -Icges(5,5) * t210 + t211 * t226;
t340 = -t122 - t124;
t164 = rSges(5,1) * t194 - rSges(5,2) * t193;
t206 = t210 ^ 2;
t207 = t211 ^ 2;
t188 = t206 + t207;
t325 = -m(5) / 0.2e1;
t331 = m(6) / 0.2e1;
t95 = t333 * t210;
t278 = (-t210 * t95 - t87) * t331 + t188 * t164 * t325;
t263 = t194 * t210;
t265 = t193 * t210;
t88 = t283 * t263 + t329 * t265;
t275 = t210 * t88;
t332 = m(5) / 0.2e1;
t142 = t164 * t210;
t144 = t164 * t211;
t78 = t210 * t142 + t144 * t211;
t279 = (t275 + t87) * t331 + t78 * t332;
t13 = t279 - t278;
t339 = t13 * qJD(1);
t221 = Icges(5,5) * t193 + Icges(5,6) * t194;
t116 = -Icges(5,3) * t210 + t211 * t221;
t118 = Icges(6,4) * t264 - Icges(6,2) * t210 - Icges(6,6) * t262;
t337 = -t116 - t118;
t276 = m(6) * qJD(5);
t90 = (0.1e1 - t188) * t194 * t193;
t336 = t90 * t276;
t243 = t283 * t193;
t335 = (Icges(6,4) + Icges(5,5)) * t194 + (-Icges(5,6) + Icges(6,6)) * t193;
t334 = (t340 * t193 + t341 * t194) * t211;
t330 = t337 * t211 + t341 * t263 + t340 * t265;
t190 = Icges(6,5) * t194;
t271 = Icges(6,1) * t193;
t121 = Icges(6,4) * t211 + (-t190 + t271) * t210;
t223 = Icges(6,4) * t193 - Icges(6,6) * t194;
t259 = t211 * (Icges(6,2) * t211 + t210 * t223) + t121 * t265;
t328 = t337 * t210 + t259 - t334;
t119 = Icges(5,6) * t211 + t210 * t225;
t180 = Icges(5,4) * t263;
t123 = Icges(5,1) * t265 + Icges(5,5) * t211 + t180;
t219 = -t119 * t194 - t123 * t193;
t326 = -t121 * t264 + t211 * t219;
t323 = t210 / 0.2e1;
t321 = t211 / 0.2e1;
t195 = t211 * qJ(2);
t281 = -pkin(6) - qJ(3);
t208 = sin(pkin(7));
t282 = pkin(3) * t208;
t217 = t211 * t282 + t195 + (-pkin(1) + t281) * t210;
t244 = -t210 * rSges(6,2) - t329 * t262;
t67 = t211 * t243 + t217 + t244;
t192 = t211 * t281;
t241 = qJ(2) + t282;
t247 = t329 * t263;
t68 = -t192 + (rSges(6,2) + pkin(1)) * t211 + (t243 + t241) * t210 - t247;
t319 = t210 * t68 + t211 * t67;
t318 = -t211 * t95 + t342;
t317 = t188 * t194;
t314 = t335 * t210;
t313 = t335 * t211;
t50 = -t211 * t88 + t342;
t280 = t50 * t331 + (-t142 * t211 + t210 * t144) * t332;
t154 = -Icges(5,2) * t193 + t269;
t268 = Icges(6,5) * t193;
t156 = Icges(6,1) * t194 + t268;
t158 = Icges(5,1) * t194 - t270;
t267 = Icges(6,3) * t193;
t312 = t193 * (t158 / 0.2e1 - t225 / 0.2e1 + t156 / 0.2e1 - Icges(6,3) * t194 / 0.2e1 + t268 / 0.2e1) - t194 * (-t226 / 0.2e1 - t154 / 0.2e1 + t190 - t271 / 0.2e1 + t267 / 0.2e1);
t310 = 0.2e1 * t188;
t309 = 0.4e1 * qJD(1);
t308 = m(3) * ((rSges(3,3) * t211 + t195) * t211 + (rSges(3,3) + qJ(2)) * t206);
t229 = rSges(4,1) * t208 + rSges(4,2) * cos(pkin(7));
t246 = rSges(4,3) + pkin(1) + qJ(3);
t98 = -t246 * t210 + t229 * t211 + t195;
t99 = t246 * t211 + (qJ(2) + t229) * t210;
t307 = m(4) * (-t210 * t98 + t211 * t99);
t306 = m(4) * (t210 * t99 + t211 * t98);
t228 = rSges(5,1) * t193 + rSges(5,2) * t194;
t212 = -t210 * rSges(5,3) + t211 * t228;
t79 = t212 + t217;
t80 = -t192 + (rSges(5,3) + pkin(1)) * t211 + (t228 + t241) * t210;
t305 = m(5) * (t142 * t80 + t144 * t79);
t304 = m(5) * (-t210 * t79 + t211 * t80);
t303 = m(5) * (t210 * t80 + t211 * t79);
t230 = -t68 * t264 + t67 * t265;
t298 = m(6) * (-t194 * t50 + t230);
t297 = m(6) * (t318 * t194 + t230);
t296 = m(6) * (t338 * t67 + t68 * t88);
t294 = m(6) * (-t210 * t67 + t211 * t68);
t293 = m(6) * t319;
t289 = m(6) * t318;
t277 = m(6) * qJD(4);
t266 = t121 * t193;
t64 = (-m(6) / 0.2e1 + t325 - m(4) / 0.2e1) * t310;
t261 = t64 * qJD(1);
t84 = m(6) * t317;
t260 = t84 * qJD(1);
t178 = Icges(6,5) * t265;
t113 = Icges(6,6) * t211 - Icges(6,3) * t263 + t178;
t257 = Icges(6,1) * t263 + t113 + t178;
t256 = t156 * t211 + t114;
t255 = -t158 * t210 + t119;
t254 = -t158 * t211 + t120;
t253 = t121 - (t190 + t267) * t210;
t252 = -Icges(6,3) * t264 + t122 - t179;
t251 = -Icges(5,2) * t265 + t123 + t180;
t250 = t154 * t211 + t124;
t249 = t329 * t194 - t243;
t31 = t319 * t194;
t245 = m(6) * t31 * qJD(1);
t41 = t211 * (Icges(5,3) * t211 + t210 * t221) + t119 * t263 + t123 * t265;
t242 = -t221 / 0.2e1 - t223 / 0.2e1;
t240 = t257 * t211;
t239 = t256 * t210;
t238 = t255 * t211;
t237 = t254 * t210;
t236 = t253 * t211;
t235 = t252 * t210;
t234 = t251 * t211;
t233 = t250 * t210;
t232 = t188 * t228;
t231 = -t113 * t194 + t118;
t94 = t249 * t210;
t96 = t249 * t211;
t227 = t94 * t210 + t211 * t96;
t39 = -t113 * t263 + t259;
t216 = -t330 * t210 / 0.2e1 + ((t231 - t266) * t211 - t326 + t330) * t323 - (t41 + t39) * t211 / 0.2e1 + ((t116 + t219) * t210 + t41 + t328 + t334) * t321;
t215 = t206 * t116 / 0.2e1 + (t210 * t231 + t328 - t39) * t323 + ((-t219 + t266 - t337) * t211 + t326 + t330) * t321;
t147 = t188 * t193;
t82 = t95 * t265;
t63 = (-m(6) / 0.4e1 - m(5) / 0.4e1 - m(4) / 0.4e1) * t310 + (m(4) + m(5) + m(6)) * t188 / 0.2e1;
t53 = -t289 / 0.2e1;
t47 = -t207 * t333 - t275;
t36 = (-t283 * t264 - t244) * t211 + (-rSges(6,2) * t211 - t283 * t265 + t247) * t210;
t24 = t264 * t338 + t317 * t36 + t82;
t21 = t297 / 0.2e1;
t19 = t298 / 0.2e1;
t18 = t53 - t280;
t17 = t289 / 0.2e1 + t280;
t16 = t53 + t280;
t15 = t294 + t304 + t307;
t12 = t278 + t279;
t10 = t293 + t303 + t306 + t308;
t9 = t296 + t305 - t312;
t4 = t21 - t298 / 0.2e1;
t3 = t21 + t19;
t2 = t19 - t297 / 0.2e1;
t1 = t210 * t216 + t211 * t215;
t5 = [t10 * qJD(2) + t15 * qJD(3) + t9 * qJD(4) - t31 * t276, qJD(1) * t10 + qJD(3) * t63 + qJD(4) * t16, qJD(1) * t15 + qJD(2) * t63 + qJD(4) * t12, t9 * qJD(1) + t16 * qJD(2) + t12 * qJD(3) + t3 * qJD(5) + (m(6) * (t67 * t94 - t68 * t96 + (-t88 + t95) * t338) + (m(5) * (-t142 * t164 + t228 * t80) + t242 * t211 - t215) * t211 + (m(5) * (t144 * t164 - t228 * t79) + t242 * t210 - t216) * t210 + (-t234 / 0.2e1 - t236 / 0.2e1 + t233 / 0.2e1 + t235 / 0.2e1) * t193 + (-t238 / 0.2e1 + t240 / 0.2e1 + t237 / 0.2e1 - t239 / 0.2e1) * t194) * qJD(4), t3 * qJD(4) - t245; t64 * qJD(3) + t17 * qJD(4) + (-t293 / 0.4e1 - t303 / 0.4e1 - t306 / 0.4e1 - t308 / 0.4e1) * t309, 0, t261, t17 * qJD(1) + 0.2e1 * (t227 * t331 + t232 * t325) * qJD(4) + t147 * t276, t147 * t277; -t64 * qJD(2) + t13 * qJD(4) - t84 * qJD(5) + (-t294 / 0.4e1 - t304 / 0.4e1 - t307 / 0.4e1) * t309, -t261, 0, t339 + (-t96 * t210 + t211 * t94) * t277, -t260; t18 * qJD(2) - t13 * qJD(3) + t1 * qJD(4) + t4 * qJD(5) + (-t296 / 0.4e1 - t305 / 0.4e1) * t309 + t312 * qJD(1), t18 * qJD(1), -t339, t1 * qJD(1) + (m(5) * (-t164 * t232 - (-t211 * t212 + (-t211 * rSges(5,3) - t210 * t228) * t210) * t78) + m(6) * (t338 * t96 + t36 * t47 + t94 * t95) + ((t314 * t210 + (t233 - t234 + t235 - t236) * t194 + ((t255 - t257) * t211 + (-t254 + t256) * t210) * t193) * t211 - t313 * t206) * t323 + ((-t313 * t211 + (-t238 + t237 + t240 - t239) * t193 + ((-t250 - t252) * t210 + (t251 + t253) * t211) * t194) * t210 + t314 * t207) * t321) * qJD(4) + t24 * t276, t4 * qJD(1) + t24 * t277 - t336; t84 * qJD(3) + t2 * qJD(4) + t245, 0, t260, t2 * qJD(1) + (t82 + (t47 + t87) * t193 + (-t227 + t36) * t194 - t24) * t277 + t336, t90 * t277;];
Cq = t5;
