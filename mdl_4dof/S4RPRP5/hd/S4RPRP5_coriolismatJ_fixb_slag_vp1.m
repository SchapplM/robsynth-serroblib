% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RPRP5_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP5_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP5_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP5_coriolismatJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP5_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRP5_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRP5_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:44:46
% EndTime: 2019-12-31 16:44:51
% DurationCPUTime: 2.73s
% Computational Cost: add. (4868->217), mult. (5949->298), div. (0->0), fcn. (5305->5), ass. (0->150)
t174 = sin(qJ(1));
t170 = pkin(6) + qJ(3);
t161 = cos(t170);
t160 = sin(t170);
t247 = rSges(5,1) + pkin(3);
t205 = t247 * t160;
t231 = rSges(5,3) + qJ(4);
t288 = -t231 * t161 + t205;
t291 = t288 * t174;
t292 = t174 * t291;
t175 = cos(qJ(1));
t156 = Icges(5,5) * t160;
t228 = Icges(5,1) * t161;
t184 = t156 + t228;
t104 = Icges(5,4) * t174 + t184 * t175;
t227 = Icges(4,4) * t160;
t136 = Icges(4,1) * t161 - t227;
t106 = Icges(4,5) * t174 + t136 * t175;
t290 = t104 + t106;
t139 = rSges(4,1) * t160 + rSges(4,2) * t161;
t171 = t174 ^ 2;
t172 = t175 ^ 2;
t207 = t171 + t172;
t274 = t207 * t139;
t283 = m(5) / 0.2e1;
t77 = t288 * t175;
t244 = (-t175 * t77 - t292) * t283 - m(4) * t274 / 0.2e1;
t218 = t161 * t175;
t72 = -t175 * t205 + t231 * t218;
t235 = t175 * t72;
t123 = t139 * t174;
t125 = t139 * t175;
t60 = t174 * t123 + t125 * t175;
t245 = (-t235 + t292) * t283 + m(4) * t60 / 0.2e1;
t11 = t245 - t244;
t289 = t11 * qJD(1);
t270 = t231 * t160 + t247 * t161;
t130 = Icges(5,4) * t161 + Icges(5,6) * t160;
t100 = Icges(5,2) * t174 + t130 * t175;
t220 = t160 * t175;
t150 = Icges(5,5) * t218;
t96 = Icges(5,6) * t174 + Icges(5,3) * t220 + t150;
t129 = Icges(4,5) * t161 - Icges(4,6) * t160;
t98 = Icges(4,3) * t174 + t129 * t175;
t287 = t96 * t220 + t290 * t218 + (t100 + t98) * t174;
t75 = t270 * t174;
t285 = (-Icges(4,6) + Icges(5,6)) * t161 + (-Icges(5,4) - Icges(4,5)) * t160;
t131 = Icges(4,2) * t161 + t227;
t224 = Icges(5,3) * t161;
t181 = t224 - t156;
t284 = (-t131 - t181) * t175 + t290;
t219 = t161 * t174;
t221 = t160 * t174;
t101 = Icges(4,4) * t219 - Icges(4,2) * t221 - Icges(4,6) * t175;
t157 = Icges(4,4) * t161;
t225 = Icges(4,2) * t160;
t102 = Icges(4,6) * t174 + (t157 - t225) * t175;
t81 = t106 * t219;
t200 = t175 * t98 - t81;
t151 = Icges(4,4) * t221;
t105 = Icges(4,1) * t219 - Icges(4,5) * t175 - t151;
t97 = Icges(4,5) * t219 - Icges(4,6) * t221 - Icges(4,3) * t175;
t241 = -t105 * t218 - t174 * t97;
t282 = -t101 * t220 - t102 * t221 - t200 - t241;
t281 = -t102 * t220 + t287;
t233 = t175 * (-Icges(5,2) * t175 + t174 * t130);
t280 = t233 + t287;
t278 = -t174 / 0.2e1;
t249 = t174 / 0.2e1;
t277 = -t175 / 0.2e1;
t239 = m(5) * qJD(3);
t103 = -Icges(5,4) * t175 + t184 * t174;
t226 = Icges(5,5) * t161;
t128 = Icges(5,3) * t160 + t226;
t95 = -Icges(5,6) * t175 + t128 * t174;
t275 = (t103 * t161 + t160 * t95) * t174;
t78 = t270 * t175;
t272 = t285 * t174;
t271 = t285 * t175;
t269 = t284 * t174;
t229 = Icges(4,1) * t160;
t185 = -t157 - t229;
t193 = (t185 * t175 - t102) * t174;
t198 = (-Icges(5,1) * t220 + t150 + t96) * t174;
t268 = t193 + t198;
t133 = Icges(5,1) * t160 - t226;
t266 = -(t136 / 0.2e1 - t131 / 0.2e1 + t156 + t228 / 0.2e1 - t224 / 0.2e1) * t160 - (t157 + t229 / 0.2e1 - t225 / 0.2e1 + t133 / 0.2e1 - t128 / 0.2e1) * t161;
t158 = cos(pkin(6)) * pkin(2) + pkin(1);
t265 = t158 + t270;
t264 = 0.4e1 * qJD(1);
t263 = m(3) * t207 * (rSges(3,3) + qJ(2));
t237 = rSges(4,1) * t161;
t197 = t158 + t237;
t246 = -pkin(5) - qJ(2);
t202 = t175 * t246;
t208 = rSges(4,2) * t221 + t175 * rSges(4,3);
t68 = -t197 * t174 - t202 + t208;
t159 = t174 * t246;
t196 = -rSges(4,2) * t220 + t174 * rSges(4,3);
t69 = t197 * t175 - t159 + t196;
t262 = m(4) * (t123 * t68 - t125 * t69);
t261 = m(4) * (t69 * t174 + t175 * t68);
t169 = t175 * rSges(5,2);
t55 = -t265 * t174 + t169 - t202;
t236 = t174 * rSges(5,2);
t56 = t265 * t175 - t159 + t236;
t243 = t55 * t218 + t56 * t219;
t257 = m(5) * ((t174 * t72 + t175 * t291) * t160 + t243);
t256 = m(5) * (-t220 * t291 + t77 * t221 + t243);
t254 = m(5) * (t291 * t55 + t56 * t72);
t253 = m(5) * (t56 * t174 + t175 * t55);
t242 = -t77 * t218 - t219 * t291;
t238 = m(5) * qJD(4);
t230 = -t133 * t174 + t95;
t223 = t101 * t160;
t222 = t160 * t161;
t126 = t207 * t160;
t64 = m(5) * t126;
t217 = t64 * qJD(1);
t216 = -t185 * t174 + t101;
t215 = -t181 * t174 + t103;
t213 = -Icges(4,2) * t219 + t105 - t151;
t209 = t207 * t222;
t27 = t56 * t220 - t55 * t221;
t206 = m(5) * t27 * qJD(1);
t204 = t129 / 0.2e1 + t130 / 0.2e1;
t199 = t230 * t175;
t195 = t102 * t160 - t97;
t194 = t216 * t175;
t192 = t215 * t175;
t190 = t213 * t175;
t188 = t100 * t175 - t104 * t219 - t96 * t221;
t32 = -t233 + t275;
t179 = (t32 - t275 + t280) * t278 + t281 * t249 + (-t81 + (t98 + t223) * t175 + t241 + t282) * t277;
t178 = t188 * t278 + (-(-t105 * t161 + t223) * t174 - t175 * t97 + t32) * t277 + (t195 * t175 - t280 + t281) * t175 / 0.2e1 + (t103 * t218 + t195 * t174 + t95 * t220 + t188 + t200 + t282) * t249;
t142 = -rSges(4,2) * t160 + t237;
t73 = t209 - t222;
t40 = -t171 * t288 + t235;
t30 = (t236 + t78) * t175 + (-t169 + t75) * t174;
t19 = t126 * t30 + t242;
t17 = t256 / 0.2e1;
t15 = t257 / 0.2e1;
t14 = t253 + t261 + t263;
t12 = t244 + t245;
t9 = t254 + t262 - t266;
t4 = t17 - t257 / 0.2e1;
t3 = t17 + t15;
t2 = t15 - t256 / 0.2e1;
t1 = t178 * t174 + t179 * t175;
t5 = [t14 * qJD(2) + t9 * qJD(3) + t27 * t238, qJD(1) * t14 + qJD(3) * t12, t9 * qJD(1) + t12 * qJD(2) + (-t55 * t78 - t56 * t75 + (-t72 - t77) * t291) * t239 + t3 * qJD(4) + ((m(4) * (-t123 * t139 - t142 * t68) + t204 * t175 - t179) * t175 + (m(4) * (t125 * t139 - t142 * t69) + t204 * t174 - t178) * t174 + (t193 / 0.2e1 + t198 / 0.2e1 + t194 / 0.2e1 - t199 / 0.2e1) * t160 + (-t190 / 0.2e1 - t192 / 0.2e1 + t284 * t249) * t161) * qJD(3), t3 * qJD(3) + t206; t11 * qJD(3) - t64 * qJD(4) + (-t263 / 0.4e1 - t261 / 0.4e1 - t253 / 0.4e1) * t264, 0, t289 + (-t78 * t174 + t175 * t75) * t239, -t217; -t11 * qJD(2) + t1 * qJD(3) + t4 * qJD(4) + (-t262 / 0.4e1 - t254 / 0.4e1) * t264 + t266 * qJD(1), -t289, t1 * qJD(1) + (m(4) * (t142 * t274 - (t174 * (rSges(4,1) * t219 - t208) + t175 * (rSges(4,1) * t218 + t196)) * t60) + m(5) * (t291 * t75 + t30 * t40 + t77 * t78) + ((-t272 * t174 + (t190 + t192 - t269) * t160 + ((t216 - t230) * t175 + t268) * t161) * t175 + t271 * t171) * t249 + ((-t271 * t175 + (t194 - t199 + t268) * t161 + ((t213 + t215) * t175 - t269) * t160) * t174 + t272 * t172) * t277) * qJD(3) + t19 * t238, t4 * qJD(1) + t19 * t239 + (-t126 * t161 + t209 - t73) * t238; t64 * qJD(2) + t2 * qJD(3) - t206, t217, t2 * qJD(1) + (-t161 * t40 + (-t174 * t75 - t175 * t78 + t30) * t160 - t19 + t242) * t239 + t73 * t238, t73 * t239;];
Cq = t5;
