% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:16
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRRPR1_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR1_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR1_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR1_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR1_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPR1_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPR1_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:15:35
% EndTime: 2019-12-05 16:15:39
% DurationCPUTime: 1.68s
% Computational Cost: add. (13474->181), mult. (7436->249), div. (0->0), fcn. (6390->8), ass. (0->131)
t267 = m(6) / 0.2e1;
t279 = m(5) / 0.2e1;
t184 = pkin(8) + qJ(2);
t182 = qJ(3) + t184;
t175 = sin(t182);
t176 = cos(t182);
t185 = cos(pkin(9));
t183 = pkin(9) + qJ(5);
t180 = cos(t183);
t237 = t180 * rSges(6,1);
t205 = pkin(4) * t185 + pkin(3) + t237;
t178 = sin(t183);
t225 = t175 * t178;
t211 = rSges(6,2) * t225 + t176 * rSges(6,3);
t243 = -pkin(7) - qJ(4);
t104 = -t205 * t175 - t176 * t243 + t211;
t245 = pkin(2) * sin(t184);
t100 = t104 - t245;
t223 = t176 * t178;
t204 = -rSges(6,2) * t223 + rSges(6,3) * t175;
t105 = -t175 * t243 + t205 * t176 + t204;
t244 = pkin(2) * cos(t184);
t101 = t105 + t244;
t46 = t100 * t176 + t101 * t175;
t55 = t104 * t176 + t105 * t175;
t276 = rSges(5,2) * sin(pkin(9)) - rSges(5,1) * t185 - pkin(3);
t277 = rSges(5,3) + qJ(4);
t115 = t276 * t175 + t277 * t176;
t109 = t115 - t245;
t116 = t277 * t175 - t276 * t176;
t110 = t116 + t244;
t62 = t109 * t176 + t110 * t175;
t69 = t115 * t176 + t116 * t175;
t241 = (t55 + t46) * t267 + (t69 + t62) * t279;
t242 = (t46 - t55) * t267 + (t62 - t69) * t279;
t4 = t242 - t241;
t281 = t4 * qJD(2);
t174 = Icges(6,4) * t180;
t152 = -Icges(6,2) * t178 + t174;
t153 = Icges(6,1) * t178 + t174;
t280 = t152 + t153;
t262 = m(5) * (-t116 * t109 + t110 * t115);
t266 = m(4) * (t244 * (-rSges(4,1) * t175 - rSges(4,2) * t176) + (rSges(4,1) * t176 - t175 * rSges(4,2)) * t245);
t210 = qJD(2) + qJD(3);
t155 = rSges(6,1) * t178 + rSges(6,2) * t180;
t136 = t155 * t175;
t137 = t155 * t176;
t195 = t175 * t136 + t176 * t137;
t189 = t195 * t267;
t272 = t176 ^ 2;
t273 = t175 ^ 2;
t274 = t155 * (t272 + t273);
t194 = m(6) * t274;
t64 = t189 + t194 / 0.2e1;
t275 = t210 * t64;
t235 = Icges(6,4) * t178;
t151 = Icges(6,2) * t180 + t235;
t154 = Icges(6,1) * t180 - t235;
t201 = t280 * t180 / 0.2e1 + (-t151 / 0.2e1 + t154 / 0.2e1) * t178;
t271 = 0.4e1 * qJD(2);
t260 = m(5) * t62;
t259 = m(5) * t69;
t40 = t100 * t136 - t101 * t137;
t41 = t104 * t136 - t105 * t137;
t258 = m(6) * (t41 + t40);
t257 = m(6) * ((-t101 + t105) * t176 + (t100 - t104) * t175) * t155;
t32 = -t105 * t100 + t101 * t104;
t254 = m(6) * t32;
t252 = m(6) * t40;
t251 = m(6) * t41;
t250 = m(6) * t46;
t249 = m(6) * t55;
t248 = -t175 / 0.2e1;
t247 = t175 / 0.2e1;
t246 = -t176 / 0.2e1;
t63 = t189 - t194 / 0.2e1;
t240 = t63 * qJD(4);
t224 = t175 * t180;
t123 = Icges(6,4) * t224 - Icges(6,2) * t225 - Icges(6,6) * t176;
t228 = t123 * t178;
t150 = Icges(6,5) * t180 - Icges(6,6) * t178;
t226 = t150 * t176;
t222 = t176 * t180;
t121 = Icges(6,5) * t224 - Icges(6,6) * t225 - Icges(6,3) * t176;
t161 = Icges(6,4) * t225;
t125 = Icges(6,1) * t224 - Icges(6,5) * t176 - t161;
t217 = -t175 * t121 - t125 * t222;
t122 = Icges(6,3) * t175 + t226;
t126 = Icges(6,5) * t175 + t154 * t176;
t216 = t175 * t122 + t126 * t222;
t215 = t153 * t175 + t123;
t124 = Icges(6,6) * t175 + t152 * t176;
t214 = -t153 * t176 - t124;
t213 = -Icges(6,2) * t224 + t125 - t161;
t212 = -t151 * t176 + t126;
t202 = t124 * t178 - t121;
t112 = t126 * t224;
t203 = t122 * t176 - t112;
t49 = -t123 * t223 - t217;
t50 = -t124 * t223 + t216;
t11 = (t202 * t176 - t216 + t50) * t176 + (t202 * t175 + t203 + t49) * t175;
t48 = -t124 * t225 - t203;
t12 = (t48 - t112 + (t122 + t228) * t176 + t217) * t176 + t216 * t175;
t22 = t175 * t48 - t176 * (-(-t125 * t180 + t228) * t175 - t121 * t176);
t23 = t175 * t50 - t176 * t49;
t2 = (t23 / 0.2e1 - t12 / 0.2e1) * t176 + (t11 / 0.2e1 + t22 / 0.2e1) * t175;
t209 = -t64 * qJD(4) + t2 * qJD(5);
t200 = t258 / 0.2e1 + t201;
t198 = Icges(6,5) * t178 + Icges(6,6) * t180;
t192 = (-t136 * t176 + t137 * t175) * t155;
t186 = (-t151 + t154) * t180 - t280 * t178;
t191 = t176 * t12 / 0.2e1 + (t11 + t22) * t248 + (t150 * t175 + t186 * t176 + t214 * t178 + t212 * t180) * t247 + (t186 * t175 - t215 * t178 + t213 * t180 - t226 + t23) * t246;
t190 = -t201 + (t247 + t248) * (t180 * t123 + t178 * t125);
t188 = t213 * t178 + t215 * t180;
t187 = -t212 * t178 + t214 * t180;
t157 = -rSges(6,2) * t178 + t237;
t131 = t176 * t198;
t130 = t198 * t175;
t60 = t64 * qJD(5);
t58 = t63 * qJD(5);
t31 = t249 + t259;
t30 = t201 + t251;
t27 = t201 + t252;
t26 = t250 + t260;
t18 = t257 / 0.2e1;
t13 = t254 + t262 + t266;
t8 = -t257 / 0.2e1 + t200;
t7 = t18 + t200;
t6 = t241 + t242;
t3 = t18 - t258 / 0.2e1 + t190;
t1 = [0, 0, 0, 0, -m(6) * t195 * qJD(5); 0, t13 * qJD(3) + t26 * qJD(4) + t27 * qJD(5), t13 * qJD(2) + t6 * qJD(4) + t7 * qJD(5) + 0.2e1 * (t262 / 0.2e1 + t32 * t267 + t266 / 0.2e1) * qJD(3), qJD(2) * t26 + qJD(3) * t6 + t58, t27 * qJD(2) + t7 * qJD(3) + ((-t46 * t157 + t192) * m(6) + t191) * qJD(5) + t240; 0, -t4 * qJD(4) + t8 * qJD(5) + (-t262 / 0.4e1 - t254 / 0.4e1 - t266 / 0.4e1) * t271, qJD(4) * t31 + qJD(5) * t30, qJD(3) * t31 - t281 + t58, t8 * qJD(2) + t30 * qJD(3) + ((-t55 * t157 + t192) * m(6) + t191) * qJD(5) + t240; 0, t4 * qJD(3) + t60 + (-t260 / 0.4e1 - t250 / 0.4e1) * t271, t281 + t60 + 0.4e1 * (-t249 / 0.4e1 - t259 / 0.4e1) * qJD(3), 0, t275; 0, (t190 - t252) * qJD(2) + t3 * qJD(3) + t209, t3 * qJD(2) + (t190 - t251) * qJD(3) + t209, -t275, (m(6) * (t157 * t274 - (t175 * (rSges(6,1) * t224 - t211) + t176 * (rSges(6,1) * t222 + t204)) * t195) + (-t273 * t131 + (t188 * t176 + (t130 + t187) * t175) * t176) * t247 + (-t272 * t130 + (t187 * t175 + (t131 + t188) * t176) * t175) * t246) * qJD(5) + t210 * t2;];
Cq = t1;
