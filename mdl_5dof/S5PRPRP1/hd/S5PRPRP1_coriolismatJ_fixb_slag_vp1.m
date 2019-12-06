% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRPRP1_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP1_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP1_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP1_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP1_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRP1_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRP1_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:28:19
% EndTime: 2019-12-05 15:28:27
% DurationCPUTime: 3.06s
% Computational Cost: add. (8202->225), mult. (6326->312), div. (0->0), fcn. (5569->5), ass. (0->155)
t183 = pkin(7) + qJ(2);
t179 = sin(t183);
t182 = pkin(8) + qJ(4);
t180 = cos(t182);
t178 = sin(t182);
t259 = rSges(6,1) + pkin(4);
t214 = t259 * t178;
t245 = rSges(6,3) + qJ(5);
t301 = -t245 * t180 + t214;
t304 = t301 * t179;
t305 = t179 * t304;
t181 = cos(t183);
t169 = Icges(6,5) * t178;
t243 = Icges(6,1) * t180;
t195 = t169 + t243;
t111 = Icges(6,4) * t179 + t195 * t181;
t242 = Icges(5,4) * t178;
t145 = Icges(5,1) * t180 - t242;
t113 = Icges(5,5) * t179 + t145 * t181;
t303 = t111 + t113;
t276 = m(6) / 0.2e1;
t148 = rSges(5,1) * t178 + rSges(5,2) * t180;
t176 = t179 ^ 2;
t177 = t181 ^ 2;
t217 = t176 + t177;
t287 = t217 * t148;
t94 = t301 * t181;
t256 = (-t181 * t94 - t305) * t276 - m(5) * t287 / 0.2e1;
t231 = t180 * t181;
t75 = -t181 * t214 + t245 * t231;
t247 = t181 * t75;
t130 = t148 * t179;
t132 = t148 * t181;
t64 = t130 * t179 + t132 * t181;
t296 = m(5) * t64;
t257 = (-t247 + t305) * t276 + t296 / 0.2e1;
t11 = t257 - t256;
t302 = t11 * qJD(2);
t283 = t245 * t178 + t259 * t180;
t157 = Icges(6,5) * t231;
t233 = t178 * t181;
t103 = Icges(6,6) * t179 + Icges(6,3) * t233 + t157;
t138 = Icges(5,5) * t180 - Icges(5,6) * t178;
t105 = Icges(5,3) * t179 + t138 * t181;
t139 = Icges(6,4) * t180 + Icges(6,6) * t178;
t107 = Icges(6,2) * t179 + t139 * t181;
t300 = t103 * t233 + t303 * t231 + (t105 + t107) * t179;
t92 = t179 * t283;
t298 = (-Icges(5,6) + Icges(6,6)) * t180 + (-Icges(6,4) - Icges(5,5)) * t178;
t140 = Icges(5,2) * t180 + t242;
t239 = Icges(6,3) * t180;
t192 = t239 - t169;
t297 = (-t140 - t192) * t181 + t303;
t232 = t180 * t179;
t235 = t178 * t179;
t108 = Icges(5,4) * t232 - Icges(5,2) * t235 - Icges(5,6) * t181;
t172 = Icges(5,4) * t180;
t240 = Icges(5,2) * t178;
t109 = Icges(5,6) * t179 + (t172 - t240) * t181;
t82 = t113 * t232;
t207 = t105 * t181 - t82;
t104 = Icges(5,5) * t232 - Icges(5,6) * t235 - Icges(5,3) * t181;
t158 = Icges(5,4) * t235;
t112 = Icges(5,1) * t232 - Icges(5,5) * t181 - t158;
t253 = -t179 * t104 - t112 * t231;
t295 = -t108 * t233 - t109 * t235 - t207 - t253;
t294 = -t109 * t233 + t300;
t237 = (-Icges(6,2) * t181 + t139 * t179) * t181;
t293 = t237 + t300;
t291 = -t179 / 0.2e1;
t261 = t179 / 0.2e1;
t290 = -t181 / 0.2e1;
t251 = m(6) * qJD(4);
t241 = Icges(6,5) * t180;
t137 = Icges(6,3) * t178 + t241;
t102 = -Icges(6,6) * t181 + t137 * t179;
t110 = -Icges(6,4) * t181 + t195 * t179;
t288 = (t102 * t178 + t110 * t180) * t179;
t95 = t283 * t181;
t285 = t298 * t179;
t284 = t298 * t181;
t282 = t297 * t179;
t244 = Icges(5,1) * t178;
t196 = -t172 - t244;
t203 = (t196 * t181 - t109) * t179;
t205 = (-Icges(6,1) * t233 + t103 + t157) * t179;
t281 = t203 + t205;
t142 = Icges(6,1) * t178 - t241;
t279 = -t178 * (t145 / 0.2e1 - t140 / 0.2e1 + t169 + t243 / 0.2e1 - t239 / 0.2e1) - t180 * (t172 + t244 / 0.2e1 - t240 / 0.2e1 + t142 / 0.2e1 - t137 / 0.2e1);
t175 = cos(pkin(8)) * pkin(3) + pkin(2);
t278 = t175 + t283;
t277 = 0.4e1 * qJD(2);
t275 = m(4) * t217 * (rSges(4,3) + qJ(3));
t249 = rSges(5,1) * t180;
t209 = t175 + t249;
t258 = -pkin(6) - qJ(3);
t211 = t181 * t258;
t218 = rSges(5,2) * t235 + t181 * rSges(5,3);
t69 = -t209 * t179 - t211 + t218;
t163 = t179 * t258;
t208 = -rSges(5,2) * t233 + rSges(5,3) * t179;
t70 = t209 * t181 - t163 + t208;
t274 = m(5) * (t130 * t69 - t132 * t70);
t273 = m(5) * (t70 * t179 + t181 * t69);
t174 = t181 * rSges(6,2);
t57 = -t179 * t278 + t174 - t211;
t248 = rSges(6,2) * t179;
t58 = t181 * t278 - t163 + t248;
t255 = t57 * t231 + t58 * t232;
t269 = m(6) * ((t179 * t75 + t181 * t304) * t178 + t255);
t268 = m(6) * (-t233 * t304 + t235 * t94 + t255);
t267 = m(6) * (t304 * t57 + t58 * t75);
t265 = m(6) * (t58 * t179 + t181 * t57);
t133 = t217 * t178;
t67 = m(6) * t133;
t254 = -t94 * t231 - t232 * t304;
t250 = m(6) * qJD(5);
t236 = t108 * t178;
t234 = t178 * t180;
t228 = t67 * qJD(2);
t227 = -t142 * t179 + t102;
t226 = -t196 * t179 + t108;
t225 = -t192 * t179 + t110;
t223 = -Icges(5,2) * t232 + t112 - t158;
t221 = t217 * t234;
t101 = (t133 / 0.2e1 - t178 / 0.2e1) * m(6);
t216 = t101 * qJD(1);
t29 = t58 * t233 - t235 * t57;
t215 = m(6) * t29 * qJD(2);
t213 = t138 / 0.2e1 + t139 / 0.2e1;
t206 = t227 * t181;
t204 = t226 * t181;
t202 = t225 * t181;
t200 = t223 * t181;
t198 = t109 * t178 - t104;
t197 = -t103 * t235 + t107 * t181 - t111 * t232;
t35 = -t237 + t288;
t188 = (t35 - t288 + t293) * t291 + t294 * t261 + (-t82 + (t105 + t236) * t181 + t253 + t295) * t290;
t187 = t197 * t291 + (-t179 * (-t112 * t180 + t236) - t104 * t181 + t35) * t290 + (t181 * t198 - t293 + t294) * t181 / 0.2e1 + (t102 * t233 + t110 * t231 + t179 * t198 + t197 + t207 + t295) * t261;
t152 = -rSges(5,2) * t178 + t249;
t100 = t67 / 0.2e1 + t178 * t276;
t79 = t221 - t234;
t43 = -t176 * t301 + t247;
t33 = (t248 + t95) * t181 + (-t174 + t92) * t179;
t19 = t133 * t33 + t254;
t17 = t268 / 0.2e1;
t15 = t269 / 0.2e1;
t14 = t265 + t273 + t275;
t12 = t256 + t257;
t9 = t267 + t274 - t279;
t4 = t17 - t269 / 0.2e1;
t3 = t17 + t15;
t2 = t15 - t268 / 0.2e1;
t1 = t179 * t187 + t181 * t188;
t5 = [0, 0, 0, 0.2e1 * (-t296 / 0.2e1 + t43 * t276) * qJD(4) + t100 * qJD(5), t100 * qJD(4); 0, t14 * qJD(3) + t9 * qJD(4) + t250 * t29, qJD(2) * t14 + qJD(4) * t12, t9 * qJD(2) + t12 * qJD(3) + (-t57 * t95 - t58 * t92 + (-t75 - t94) * t304) * t251 + t3 * qJD(5) + ((m(5) * (-t130 * t148 - t152 * t69) + t213 * t181 - t188) * t181 + (m(5) * (t132 * t148 - t152 * t70) + t213 * t179 - t187) * t179 + (t203 / 0.2e1 + t205 / 0.2e1 + t204 / 0.2e1 - t206 / 0.2e1) * t178 + (-t200 / 0.2e1 - t202 / 0.2e1 + t297 * t261) * t180) * qJD(4), t3 * qJD(4) + t215; 0, t11 * qJD(4) - t67 * qJD(5) + (-t275 / 0.4e1 - t273 / 0.4e1 - t265 / 0.4e1) * t277, 0, t302 + (-t179 * t95 + t181 * t92) * t251, -t228; qJD(5) * t101, -t11 * qJD(3) + t1 * qJD(4) + t4 * qJD(5) + (-t274 / 0.4e1 - t267 / 0.4e1) * t277 + t279 * qJD(2), -t302, t1 * qJD(2) + (m(5) * (t152 * t287 - (t179 * (rSges(5,1) * t232 - t218) + t181 * (rSges(5,1) * t231 + t208)) * t64) + m(6) * (t304 * t92 + t33 * t43 + t94 * t95) + ((-t285 * t179 + (t200 + t202 - t282) * t178 + ((t226 - t227) * t181 + t281) * t180) * t181 + t284 * t176) * t261 + ((-t284 * t181 + (t204 - t206 + t281) * t180 + ((t223 + t225) * t181 - t282) * t178) * t179 + t285 * t177) * t290) * qJD(4) + t19 * t250, t216 + t4 * qJD(2) + t19 * t251 + (-t133 * t180 + t221 - t79) * t250; -t101 * qJD(4), t67 * qJD(3) + t2 * qJD(4) - t215, t228, -t216 + t2 * qJD(2) + (-t180 * t43 + (-t179 * t92 - t181 * t95 + t33) * t178 - t19 + t254) * t251 + t79 * t250, t79 * t251;];
Cq = t5;
