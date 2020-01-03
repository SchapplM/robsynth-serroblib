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
% Datum: 2020-01-03 11:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:55:44
% EndTime: 2020-01-03 11:55:48
% DurationCPUTime: 1.82s
% Computational Cost: add. (14256->202), mult. (7992->261), div. (0->0), fcn. (6878->10), ass. (0->137)
t190 = qJ(1) + qJ(2);
t185 = pkin(8) + t190;
t179 = sin(t185);
t180 = cos(t185);
t187 = cos(t190);
t182 = pkin(2) * t187;
t191 = cos(pkin(9));
t287 = -rSges(5,2) * sin(pkin(9)) + pkin(3) + rSges(5,1) * t191;
t288 = -qJ(4) - rSges(5,3);
t108 = -t288 * t179 + t287 * t180 + t182;
t188 = cos(qJ(1)) * pkin(1);
t105 = t188 + t108;
t100 = t105 * t179;
t181 = pkin(4) * t191 + pkin(3);
t192 = -pkin(7) - qJ(4);
t189 = pkin(9) + qJ(5);
t183 = sin(t189);
t225 = t179 * t183;
t211 = -rSges(6,2) * t225 - t180 * rSges(6,3);
t184 = cos(t189);
t241 = rSges(6,1) * t184;
t186 = sin(t190);
t248 = pkin(2) * t186;
t101 = t248 + t180 * t192 + (t181 + t241) * t179 + t211;
t103 = t108 * t179;
t107 = t287 * t179 + t288 * t180 + t248;
t249 = sin(qJ(1)) * pkin(1);
t104 = t107 + t249;
t290 = m(6) / 0.2e1;
t291 = m(5) / 0.2e1;
t222 = t180 * t184;
t223 = t180 * t183;
t205 = rSges(6,1) * t222 - rSges(6,2) * t223;
t102 = t180 * t181 + t182 + (rSges(6,3) - t192) * t179 + t205;
t97 = t188 + t102;
t93 = t97 * t179;
t95 = t102 * t179;
t96 = t101 + t249;
t246 = (t93 + t95 + (-t101 - t96) * t180) * t290 + (t100 + t103 + (-t104 - t107) * t180) * t291;
t238 = t101 - t96;
t247 = (t238 * t180 + t93 - t95) * t290 + (t100 - t103 + (-t104 + t107) * t180) * t291;
t3 = t247 - t246;
t293 = t3 * qJD(1);
t177 = Icges(6,4) * t184;
t153 = -Icges(6,2) * t183 + t177;
t282 = Icges(6,1) * t183 + t177;
t292 = t153 + t282;
t272 = m(3) * (-t188 * (rSges(3,1) * t186 + rSges(3,2) * t187) + t249 * (t187 * rSges(3,1) - rSges(3,2) * t186));
t271 = m(4) * (-t188 * (rSges(4,1) * t179 + rSges(4,2) * t180 + t248) + t249 * (t180 * rSges(4,1) - rSges(4,2) * t179 + t182));
t254 = m(6) * (-t101 * t180 + t95);
t255 = m(6) * (-t180 * t96 + t93);
t267 = m(5) * (t104 * t108 - t107 * t105);
t259 = m(6) * (-t101 * t97 + t96 * t102);
t156 = rSges(6,1) * t183 + rSges(6,2) * t184;
t137 = t156 * t179;
t138 = t156 * t180;
t98 = -t179 * t137 - t180 * t138;
t198 = m(6) * t98;
t210 = qJD(1) + qJD(2);
t277 = t180 ^ 2;
t278 = t179 ^ 2;
t283 = t156 * (t277 + t278);
t197 = -m(6) * t283 / 0.2e1;
t68 = t198 / 0.2e1 + t197;
t286 = t210 * t68;
t236 = Icges(6,4) * t183;
t152 = Icges(6,2) * t184 + t236;
t155 = Icges(6,1) * t184 - t236;
t281 = (t152 - t155) * t184 + t292 * t183;
t162 = Icges(6,4) * t223;
t123 = Icges(6,1) * t222 + t179 * Icges(6,5) - t162;
t214 = -Icges(6,2) * t222 + t123 - t162;
t121 = Icges(6,4) * t222 - Icges(6,2) * t223 + t179 * Icges(6,6);
t216 = t180 * t282 + t121;
t280 = -t214 * t183 - t216 * t184;
t122 = -Icges(6,5) * t180 + t155 * t179;
t215 = -t152 * t179 + t122;
t120 = -Icges(6,6) * t180 + t153 * t179;
t217 = t179 * t282 + t120;
t279 = -t215 * t183 - t217 * t184;
t206 = t292 * t184 / 0.2e1 + (-t152 / 0.2e1 + t155 / 0.2e1) * t183;
t276 = 0.4e1 * qJD(1);
t265 = m(5) * (-t104 * t180 + t100);
t264 = m(5) * (-t107 * t180 + t103);
t40 = -t96 * t137 - t97 * t138;
t41 = -t101 * t137 - t102 * t138;
t263 = m(6) * (t41 + t40);
t262 = m(6) * ((t102 - t97) * t180 + t238 * t179) * t156;
t257 = m(6) * t40;
t256 = m(6) * t41;
t253 = -t179 / 0.2e1;
t252 = -t180 / 0.2e1;
t251 = t180 / 0.2e1;
t69 = -t198 / 0.2e1 + t197;
t244 = t69 * qJD(4);
t231 = t120 * t183;
t230 = t121 * t183;
t229 = t122 * t184;
t151 = Icges(6,5) * t184 - Icges(6,6) * t183;
t227 = t151 * t179;
t224 = t179 * t184;
t113 = t122 * t224;
t114 = t123 * t224;
t115 = t120 * t223;
t118 = -Icges(6,3) * t180 + t227;
t199 = t123 * t184 - t230;
t52 = -t118 * t179 - t122 * t222 + t115;
t119 = Icges(6,5) * t222 - Icges(6,6) * t223 + Icges(6,3) * t179;
t53 = t119 * t179 + t199 * t180;
t11 = (t52 + t114 - t115 + (t118 - t230) * t179) * t179 + (-t113 - t53 + (t118 + t199) * t180 + (t229 + t231) * t179) * t180;
t50 = -t118 * t180 - t120 * t225 + t113;
t51 = t119 * t180 + t121 * t225 - t114;
t12 = (t115 - t51 + (t119 - t229) * t180) * t180 + (-t113 + t50 + (t119 + t231) * t179) * t179;
t27 = -t179 * t51 - t180 * t50;
t28 = -t179 * t53 - t180 * t52;
t2 = (-t12 / 0.2e1 - t28 / 0.2e1) * t180 + (t27 / 0.2e1 - t11 / 0.2e1) * t179;
t209 = t68 * qJD(4) + t2 * qJD(5);
t203 = t263 / 0.2e1 + t206;
t201 = Icges(6,5) * t183 + Icges(6,6) * t184;
t196 = t179 * t11 / 0.2e1 + (-t151 * t180 - t281 * t179 - t217 * t183 + t215 * t184) * t252 + (t12 + t28) * t251 + (t281 * t180 + t216 * t183 - t214 * t184 - t227 + t27) * t253;
t195 = -t206 + (t251 + t252) * (t184 * t121 + t183 * t123);
t157 = -rSges(6,2) * t183 + t241;
t132 = t201 * t180;
t131 = t179 * t201;
t65 = t68 * qJD(5);
t63 = t69 * qJD(5);
t31 = t254 + t264;
t30 = t206 + t256;
t29 = t206 + t257;
t26 = t255 + t265;
t18 = t262 / 0.2e1;
t13 = t259 + t267 + t271 + t272;
t8 = -t262 / 0.2e1 + t203;
t7 = t18 + t203;
t6 = t18 - t263 / 0.2e1 + t195;
t5 = t246 + t247;
t1 = [t13 * qJD(2) + t26 * qJD(4) + t29 * qJD(5), t13 * qJD(1) + t5 * qJD(4) + t7 * qJD(5) + 0.2e1 * (t272 / 0.2e1 + t271 / 0.2e1 + t267 / 0.2e1 + t259 / 0.2e1) * qJD(2), 0, qJD(1) * t26 + qJD(2) * t5 + t63, t29 * qJD(1) + t7 * qJD(2) + (-t157 * t255 + t196) * qJD(5) + t244; -t3 * qJD(4) + t8 * qJD(5) + (-t272 / 0.4e1 - t271 / 0.4e1 - t267 / 0.4e1 - t259 / 0.4e1) * t276, qJD(4) * t31 + qJD(5) * t30, 0, qJD(2) * t31 - t293 + t63, t8 * qJD(1) + t30 * qJD(2) + (-t157 * t254 + t196) * qJD(5) + t244; 0, 0, 0, 0, qJD(5) * t198; (-t265 / 0.4e1 - t255 / 0.4e1) * t276 + t3 * qJD(2) - t65, t293 - t65 + 0.4e1 * (-t254 / 0.4e1 - t264 / 0.4e1) * qJD(2), 0, 0, -t286; (t195 - t257) * qJD(1) + t6 * qJD(2) + t209, t6 * qJD(1) + (t195 - t256) * qJD(2) + t209, 0, t286, (m(6) * (t157 * t283 + (t180 * (rSges(6,3) * t179 + t205) + t179 * (rSges(6,1) * t224 + t211)) * t98) + (-t277 * t131 + (t280 * t179 + (t132 - t279) * t180) * t179) * t252 + (t278 * t132 + (t279 * t180 + (-t131 - t280) * t179) * t180) * t253) * qJD(5) + t210 * t2;];
Cq = t1;
