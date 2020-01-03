% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
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
% Datum: 2019-12-31 17:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPPRP2_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP2_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP2_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP2_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP2_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP2_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRP2_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:49:06
% EndTime: 2019-12-31 17:49:11
% DurationCPUTime: 3.07s
% Computational Cost: add. (8302->230), mult. (6434->317), div. (0->0), fcn. (5671->7), ass. (0->158)
t185 = qJ(1) + pkin(7);
t181 = sin(t185);
t184 = pkin(8) + qJ(4);
t182 = cos(t184);
t180 = sin(t184);
t266 = rSges(6,1) + pkin(4);
t219 = t266 * t180;
t250 = rSges(6,3) + qJ(5);
t308 = -t250 * t182 + t219;
t311 = t308 * t181;
t312 = t181 * t311;
t183 = cos(t185);
t171 = Icges(6,5) * t180;
t248 = Icges(6,1) * t182;
t200 = t171 + t248;
t111 = Icges(6,4) * t181 + t200 * t183;
t247 = Icges(5,4) * t180;
t147 = Icges(5,1) * t182 - t247;
t113 = Icges(5,5) * t181 + t147 * t183;
t310 = t111 + t113;
t283 = m(6) / 0.2e1;
t150 = rSges(5,1) * t180 + rSges(5,2) * t182;
t178 = t181 ^ 2;
t179 = t183 ^ 2;
t222 = t178 + t179;
t294 = t222 * t150;
t94 = t308 * t183;
t261 = (-t183 * t94 - t312) * t283 - m(5) * t294 / 0.2e1;
t236 = t182 * t183;
t76 = -t183 * t219 + t250 * t236;
t252 = t183 * t76;
t130 = t150 * t181;
t132 = t150 * t183;
t65 = t130 * t181 + t132 * t183;
t303 = m(5) * t65;
t262 = (-t252 + t312) * t283 + t303 / 0.2e1;
t11 = t262 - t261;
t309 = t11 * qJD(1);
t290 = t250 * t180 + t266 * t182;
t157 = Icges(6,5) * t236;
t238 = t180 * t183;
t103 = Icges(6,6) * t181 + Icges(6,3) * t238 + t157;
t140 = Icges(5,5) * t182 - Icges(5,6) * t180;
t105 = Icges(5,3) * t181 + t140 * t183;
t141 = Icges(6,4) * t182 + Icges(6,6) * t180;
t107 = Icges(6,2) * t181 + t141 * t183;
t307 = t103 * t238 + t310 * t236 + (t105 + t107) * t181;
t92 = t181 * t290;
t305 = (-Icges(5,6) + Icges(6,6)) * t182 + (-Icges(6,4) - Icges(5,5)) * t180;
t142 = Icges(5,2) * t182 + t247;
t244 = Icges(6,3) * t182;
t197 = t244 - t171;
t304 = (-t142 - t197) * t183 + t310;
t237 = t181 * t182;
t240 = t180 * t181;
t108 = Icges(5,4) * t237 - Icges(5,2) * t240 - Icges(5,6) * t183;
t174 = Icges(5,4) * t182;
t245 = Icges(5,2) * t180;
t109 = Icges(5,6) * t181 + (t174 - t245) * t183;
t84 = t113 * t237;
t212 = t105 * t183 - t84;
t104 = Icges(5,5) * t237 - Icges(5,6) * t240 - Icges(5,3) * t183;
t158 = Icges(5,4) * t240;
t112 = Icges(5,1) * t237 - Icges(5,5) * t183 - t158;
t258 = -t181 * t104 - t112 * t236;
t302 = -t108 * t238 - t109 * t240 - t212 - t258;
t301 = -t109 * t238 + t307;
t242 = (-Icges(6,2) * t183 + t141 * t181) * t183;
t300 = t242 + t307;
t298 = -t181 / 0.2e1;
t268 = t181 / 0.2e1;
t297 = -t183 / 0.2e1;
t256 = m(6) * qJD(4);
t246 = Icges(6,5) * t182;
t139 = Icges(6,3) * t180 + t246;
t102 = -Icges(6,6) * t183 + t139 * t181;
t110 = -Icges(6,4) * t183 + t200 * t181;
t295 = (t102 * t180 + t110 * t182) * t181;
t95 = t290 * t183;
t292 = t305 * t181;
t291 = t305 * t183;
t289 = t304 * t181;
t249 = Icges(5,1) * t180;
t201 = -t174 - t249;
t208 = (t201 * t183 - t109) * t181;
t210 = (-Icges(6,1) * t238 + t103 + t157) * t181;
t288 = t208 + t210;
t144 = Icges(6,1) * t180 - t246;
t286 = -(t147 / 0.2e1 - t142 / 0.2e1 + t171 + t248 / 0.2e1 - t244 / 0.2e1) * t180 - (t174 + t249 / 0.2e1 - t245 / 0.2e1 + t144 / 0.2e1 - t139 / 0.2e1) * t182;
t177 = cos(pkin(8)) * pkin(3) + pkin(2);
t285 = t177 + t290;
t284 = 0.4e1 * qJD(1);
t251 = rSges(4,3) + qJ(3);
t264 = sin(qJ(1)) * pkin(1);
t265 = cos(qJ(1)) * pkin(1);
t282 = m(4) * (t183 * (t251 * t183 - t264) + (t251 * t181 + t265) * t181);
t263 = -pkin(6) - qJ(3);
t191 = -t183 * t263 - t264;
t254 = rSges(5,1) * t182;
t214 = t177 + t254;
t223 = rSges(5,2) * t240 + t183 * rSges(5,3);
t67 = -t214 * t181 + t191 + t223;
t213 = -rSges(5,2) * t238 + t181 * rSges(5,3);
t215 = -t181 * t263 + t265;
t68 = t214 * t183 + t213 + t215;
t281 = m(5) * (t130 * t67 - t132 * t68);
t280 = m(5) * (t68 * t181 + t183 * t67);
t176 = t183 * rSges(6,2);
t56 = -t181 * t285 + t176 + t191;
t253 = t181 * rSges(6,2);
t57 = t183 * t285 + t215 + t253;
t260 = t56 * t236 + t57 * t237;
t276 = m(6) * ((t181 * t76 + t183 * t311) * t180 + t260);
t275 = m(6) * (-t238 * t311 + t94 * t240 + t260);
t274 = m(6) * (t311 * t56 + t57 * t76);
t272 = m(6) * (t57 * t181 + t183 * t56);
t133 = t222 * t180;
t69 = m(6) * t133;
t259 = -t94 * t236 - t237 * t311;
t255 = m(6) * qJD(5);
t241 = t108 * t180;
t239 = t180 * t182;
t233 = t69 * qJD(1);
t232 = -t144 * t181 + t102;
t231 = -t201 * t181 + t108;
t230 = -t197 * t181 + t110;
t228 = -Icges(5,2) * t237 + t112 - t158;
t226 = t222 * t239;
t101 = (t133 / 0.2e1 - t180 / 0.2e1) * m(6);
t221 = t101 * qJD(2);
t29 = t57 * t238 - t56 * t240;
t220 = m(6) * t29 * qJD(1);
t218 = t140 / 0.2e1 + t141 / 0.2e1;
t211 = t232 * t183;
t209 = t231 * t183;
t207 = t230 * t183;
t205 = t228 * t183;
t203 = t109 * t180 - t104;
t202 = -t103 * t240 + t107 * t183 - t111 * t237;
t36 = -t242 + t295;
t193 = (t36 - t295 + t300) * t298 + t301 * t268 + (-t84 + (t105 + t241) * t183 + t258 + t302) * t297;
t192 = t202 * t298 + (-(-t112 * t182 + t241) * t181 - t104 * t183 + t36) * t297 + (t203 * t183 - t300 + t301) * t183 / 0.2e1 + (t102 * t238 + t110 * t236 + t203 * t181 + t202 + t212 + t302) * t268;
t153 = -rSges(5,2) * t180 + t254;
t100 = t69 / 0.2e1 + t180 * t283;
t81 = t226 - t239;
t44 = -t178 * t308 + t252;
t33 = (t253 + t95) * t183 + (-t176 + t92) * t181;
t19 = t133 * t33 + t259;
t17 = t275 / 0.2e1;
t15 = t276 / 0.2e1;
t14 = t272 + t280 + t282;
t12 = t261 + t262;
t9 = t274 + t281 - t286;
t4 = t17 - t276 / 0.2e1;
t3 = t17 + t15;
t2 = t15 - t275 / 0.2e1;
t1 = t192 * t181 + t193 * t183;
t5 = [t14 * qJD(3) + t9 * qJD(4) + t29 * t255, 0, qJD(1) * t14 + qJD(4) * t12, t9 * qJD(1) + t12 * qJD(3) + (-t56 * t95 - t57 * t92 + (-t76 - t94) * t311) * t256 + t3 * qJD(5) + ((m(5) * (-t130 * t150 - t153 * t67) + t218 * t183 - t193) * t183 + (m(5) * (t132 * t150 - t153 * t68) + t218 * t181 - t192) * t181 + (t208 / 0.2e1 + t210 / 0.2e1 + t209 / 0.2e1 - t211 / 0.2e1) * t180 + (-t205 / 0.2e1 - t207 / 0.2e1 + t304 * t268) * t182) * qJD(4), t3 * qJD(4) + t220; 0, 0, 0, 0.2e1 * (-t303 / 0.2e1 + t44 * t283) * qJD(4) + t100 * qJD(5), t100 * qJD(4); t11 * qJD(4) - t69 * qJD(5) + (-t282 / 0.4e1 - t280 / 0.4e1 - t272 / 0.4e1) * t284, 0, 0, t309 + (-t181 * t95 + t183 * t92) * t256, -t233; -t11 * qJD(3) + t1 * qJD(4) + t4 * qJD(5) + (-t281 / 0.4e1 - t274 / 0.4e1) * t284 + t286 * qJD(1), qJD(5) * t101, -t309, t1 * qJD(1) + (m(5) * (t153 * t294 - (t181 * (rSges(5,1) * t237 - t223) + t183 * (rSges(5,1) * t236 + t213)) * t65) + m(6) * (t311 * t92 + t33 * t44 + t94 * t95) + ((-t292 * t181 + (t205 + t207 - t289) * t180 + ((t231 - t232) * t183 + t288) * t182) * t183 + t291 * t178) * t268 + ((-t291 * t183 + (t209 - t211 + t288) * t182 + ((t228 + t230) * t183 - t289) * t180) * t181 + t292 * t179) * t297) * qJD(4) + t19 * t255, t4 * qJD(1) + t221 + t19 * t256 + (-t133 * t182 + t226 - t81) * t255; t69 * qJD(3) + t2 * qJD(4) - t220, -t101 * qJD(4), t233, t2 * qJD(1) - t221 + (-t182 * t44 + (-t181 * t92 - t183 * t95 + t33) * t180 - t19 + t259) * t256 + t81 * t255, t81 * t256;];
Cq = t5;
