% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-05 17:50
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRPR2_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR2_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR2_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR2_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR2_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:49:21
% EndTime: 2019-12-05 17:49:25
% DurationCPUTime: 1.73s
% Computational Cost: add. (13630->193), mult. (7600->251), div. (0->0), fcn. (6548->10), ass. (0->135)
t184 = qJ(1) + pkin(8);
t182 = qJ(3) + t184;
t175 = sin(t182);
t176 = cos(t182);
t185 = cos(pkin(9));
t290 = rSges(5,2) * sin(pkin(9)) - rSges(5,1) * t185 - pkin(3);
t291 = -rSges(5,3) - qJ(4);
t114 = -t290 * t175 + t291 * t176;
t202 = sin(qJ(1)) * pkin(1) + pkin(2) * sin(t184);
t103 = t114 + t202;
t101 = t176 * t103;
t115 = t291 * t175 + t290 * t176;
t201 = -cos(qJ(1)) * pkin(1) - pkin(2) * cos(t184);
t104 = t115 + t201;
t186 = -pkin(7) - qJ(4);
t183 = pkin(9) + qJ(5);
t178 = sin(t183);
t226 = t176 * t178;
t205 = rSges(6,2) * t226 - t175 * rSges(6,3);
t180 = cos(t183);
t245 = rSges(6,1) * t180;
t206 = pkin(4) * t185 + pkin(3) + t245;
t106 = t175 * t186 - t176 * t206 + t205;
t227 = t176 * t114;
t230 = t175 * t178;
t212 = -rSges(6,2) * t230 - t176 * rSges(6,3);
t105 = t175 * t206 + t176 * t186 + t212;
t228 = t176 * t105;
t292 = m(6) / 0.2e1;
t293 = m(5) / 0.2e1;
t97 = t105 + t202;
t93 = t176 * t97;
t98 = t106 + t201;
t248 = (-t228 - t93 + (-t106 - t98) * t175) * t292 + (-t101 - t227 + (-t104 - t115) * t175) * t293;
t241 = t106 - t98;
t249 = (t241 * t175 + t228 - t93) * t292 + (-t101 + t227 + (-t104 + t115) * t175) * t293;
t4 = t249 - t248;
t294 = t4 * qJD(1);
t266 = m(5) * (-t115 * t103 + t104 * t114);
t258 = m(6) * (t98 * t105 - t106 * t97);
t157 = rSges(6,1) * t178 + rSges(6,2) * t180;
t138 = t157 * t175;
t139 = t157 * t176;
t194 = t175 * t138 + t176 * t139;
t192 = m(6) * t194;
t270 = m(4) * (t201 * (rSges(4,1) * t175 + rSges(4,2) * t176) - (-t176 * rSges(4,1) + t175 * rSges(4,2)) * t202);
t237 = t105 * t138;
t41 = t106 * t139 - t237;
t289 = m(6) * t41;
t158 = -rSges(6,2) * t178 + t245;
t288 = m(6) * t158;
t211 = qJD(1) + qJD(3);
t276 = t176 ^ 2;
t277 = t175 ^ 2;
t285 = t157 * (t276 + t277);
t191 = -m(6) * t285 / 0.2e1;
t63 = -t192 / 0.2e1 + t191;
t287 = t211 * t63;
t286 = t175 * t98 + t93;
t225 = t176 * t180;
t124 = Icges(6,4) * t225 - Icges(6,2) * t226 + t175 * Icges(6,6);
t162 = Icges(6,4) * t226;
t126 = Icges(6,1) * t225 + t175 * Icges(6,5) - t162;
t284 = (t124 * t178 - t126 * t180) * t176;
t174 = Icges(6,4) * t180;
t283 = Icges(6,1) * t178 + t174;
t282 = t106 * t175 + t228;
t281 = Icges(6,2) * t178 - t174;
t239 = Icges(6,4) * t178;
t153 = Icges(6,2) * t180 + t239;
t156 = Icges(6,1) * t180 - t239;
t280 = (t153 - t156) * t180 + (-t281 + t283) * t178;
t215 = -Icges(6,2) * t225 + t126 - t162;
t217 = t176 * t283 + t124;
t279 = t215 * t178 + t217 * t180;
t161 = Icges(6,4) * t230;
t229 = t175 * t180;
t125 = -Icges(6,1) * t229 + Icges(6,5) * t176 + t161;
t216 = Icges(6,2) * t229 + t125 + t161;
t123 = Icges(6,6) * t176 + t281 * t175;
t218 = -t175 * t283 + t123;
t278 = -t216 * t178 - t218 * t180;
t203 = (-t281 / 0.2e1 + t283 / 0.2e1) * t180 + (-t153 / 0.2e1 + t156 / 0.2e1) * t178;
t275 = 0.4e1 * qJD(1);
t264 = m(5) * (-t104 * t175 - t101);
t263 = m(5) * (-t115 * t175 - t227);
t242 = t97 * t138;
t39 = t98 * t139 - t242;
t262 = m(6) * (t41 + t39);
t261 = m(6) * (-t241 * t139 + t237 - t242);
t256 = m(6) * t39;
t254 = m(6) * t286;
t253 = m(6) * t282;
t252 = t175 / 0.2e1;
t251 = -t176 / 0.2e1;
t250 = t176 / 0.2e1;
t64 = t192 / 0.2e1 + t191;
t247 = t64 * qJD(4);
t233 = t123 * t178;
t152 = Icges(6,5) * t180 - Icges(6,6) * t178;
t231 = t152 * t175;
t121 = Icges(6,3) * t176 - t231;
t220 = t176 * t121 + t123 * t230;
t219 = t175 * t121 + t125 * t225;
t122 = Icges(6,5) * t225 - Icges(6,6) * t226 + Icges(6,3) * t175;
t204 = -t125 * t180 - t122;
t48 = t176 * t122 + t124 * t230 - t126 * t229;
t49 = -t123 * t226 + t219;
t50 = t122 * t175 - t284;
t11 = (t220 + t50 + t284) * t176 + (-t49 + (t204 - t233) * t176 + t48 + t219) * t175;
t47 = -t125 * t229 + t220;
t12 = (t48 + (-t122 + t233) * t176 - t219) * t176 + (t175 * t204 + t220 - t47) * t175;
t23 = t175 * t48 + t176 * t47;
t24 = t175 * t50 + t176 * t49;
t2 = (t12 / 0.2e1 + t24 / 0.2e1) * t176 + (-t23 / 0.2e1 + t11 / 0.2e1) * t175;
t210 = t63 * qJD(4) + t2 * qJD(5);
t199 = t262 / 0.2e1 + t203;
t197 = Icges(6,5) * t178 + Icges(6,6) * t180;
t132 = t175 * t197;
t190 = -t175 * t11 / 0.2e1 + (t12 + t24) * t251 + (t152 * t176 + t280 * t175 - t218 * t178 + t216 * t180) * t250 + (-t280 * t176 - t217 * t178 + t215 * t180 + t23 + t231) * t252;
t189 = -t203 + (t250 + t251) * (t180 * t124 + t178 * t126);
t133 = t197 * t176;
t62 = t63 * qJD(5);
t60 = t64 * qJD(5);
t32 = -t253 + t263;
t30 = t203 + t289;
t27 = t203 + t256;
t22 = -t254 + t264;
t18 = t261 / 0.2e1;
t13 = t258 + t266 + t270;
t8 = -t261 / 0.2e1 + t199;
t7 = t18 + t199;
t5 = t248 + t249;
t3 = t18 - t262 / 0.2e1 + t189;
t1 = [t13 * qJD(3) + t22 * qJD(4) + t27 * qJD(5), 0, t13 * qJD(1) + t5 * qJD(4) + t7 * qJD(5) + 0.2e1 * (t270 / 0.2e1 + t266 / 0.2e1 + t258 / 0.2e1) * qJD(3), qJD(1) * t22 + qJD(3) * t5 + t60, t27 * qJD(1) + t7 * qJD(3) + (t286 * t288 + t190) * qJD(5) + t247; 0, 0, 0, 0, -qJD(5) * t192; -t4 * qJD(4) + t8 * qJD(5) + (-t270 / 0.4e1 - t266 / 0.4e1 - t258 / 0.4e1) * t275, 0, qJD(4) * t32 + qJD(5) * t30, qJD(3) * t32 - t294 + t60, t8 * qJD(1) + t30 * qJD(3) + (t282 * t288 + t190) * qJD(5) + t247; t4 * qJD(3) - t62 + (-t264 / 0.4e1 + t254 / 0.4e1) * t275, 0, t294 - t62 + 0.4e1 * (-t263 / 0.4e1 + t253 / 0.4e1) * qJD(3), 0, -t287; (t189 - t256) * qJD(1) + t3 * qJD(3) + t210, 0, t3 * qJD(1) + (t189 - t289) * qJD(3) + t210, t287, (m(6) * (t158 * t285 - (t176 * (rSges(6,1) * t225 - t205) - t175 * (-rSges(6,1) * t229 - t212)) * t194) + (t276 * t132 + (t279 * t175 + (-t133 - t278) * t176) * t175) * t250 + (-t277 * t133 + (t278 * t176 + (t132 - t279) * t175) * t176) * t252) * qJD(5) + t211 * t2;];
Cq = t1;
