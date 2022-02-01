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
% m [6x1]
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
% Datum: 2022-01-23 09:19
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-23 09:18:49
% EndTime: 2022-01-23 09:18:52
% DurationCPUTime: 1.66s
% Computational Cost: add. (13630->183), mult. (7600->251), div. (0->0), fcn. (6548->10), ass. (0->131)
t271 = m(6) / 0.2e1;
t186 = qJ(1) + pkin(8);
t184 = qJ(3) + t186;
t177 = sin(t184);
t178 = cos(t184);
t187 = cos(pkin(9));
t185 = pkin(9) + qJ(5);
t182 = cos(t185);
t243 = t182 * rSges(6,1);
t211 = pkin(4) * t187 + pkin(3) + t243;
t180 = sin(t185);
t231 = t177 * t180;
t217 = rSges(6,2) * t231 + t178 * rSges(6,3);
t249 = -pkin(7) - qJ(4);
t106 = -t211 * t177 - t178 * t249 + t217;
t229 = t178 * t180;
t210 = -rSges(6,2) * t229 + t177 * rSges(6,3);
t107 = -t177 * t249 + t211 * t178 + t210;
t278 = -t106 * t178 - t107 * t177;
t284 = m(5) / 0.2e1;
t206 = -sin(qJ(1)) * pkin(1) - pkin(2) * sin(t186);
t96 = t106 + t206;
t205 = cos(qJ(1)) * pkin(1) + pkin(2) * cos(t186);
t97 = t107 + t205;
t46 = t97 * t177 + t96 * t178;
t281 = rSges(5,2) * sin(pkin(9)) - rSges(5,1) * t187 - pkin(3);
t282 = rSges(5,3) + qJ(4);
t115 = t281 * t177 + t282 * t178;
t104 = t115 + t206;
t116 = t282 * t177 - t281 * t178;
t105 = t116 + t205;
t55 = t104 * t178 + t105 * t177;
t69 = t115 * t178 + t116 * t177;
t247 = (-t278 + t46) * t271 + (t69 + t55) * t284;
t248 = (t46 + t278) * t271 + (t55 - t69) * t284;
t4 = t248 - t247;
t286 = t4 * qJD(1);
t176 = Icges(6,4) * t182;
t154 = -Icges(6,2) * t180 + t176;
t155 = Icges(6,1) * t180 + t176;
t285 = t154 + t155;
t266 = m(5) * (-t116 * t104 + t105 * t115);
t270 = m(4) * (t205 * (-rSges(4,1) * t177 - rSges(4,2) * t178) - (t178 * rSges(4,1) - t177 * rSges(4,2)) * t206);
t216 = qJD(1) + qJD(3);
t157 = rSges(6,1) * t180 + rSges(6,2) * t182;
t138 = t157 * t177;
t139 = t157 * t178;
t199 = t177 * t138 + t178 * t139;
t193 = t199 * t271;
t276 = t178 ^ 2;
t277 = t177 ^ 2;
t279 = t157 * (t276 + t277);
t198 = m(6) * t279;
t64 = t193 + t198 / 0.2e1;
t280 = t216 * t64;
t241 = Icges(6,4) * t180;
t153 = Icges(6,2) * t182 + t241;
t156 = Icges(6,1) * t182 - t241;
t207 = t285 * t182 / 0.2e1 + (-t153 / 0.2e1 + t156 / 0.2e1) * t180;
t275 = 0.4e1 * qJD(1);
t264 = m(5) * t55;
t263 = m(5) * t69;
t39 = t96 * t138 - t97 * t139;
t41 = t106 * t138 - t107 * t139;
t262 = m(6) * (t41 + t39);
t261 = m(6) * ((t107 - t97) * t178 + (-t106 + t96) * t177) * t157;
t31 = t97 * t106 - t107 * t96;
t258 = m(6) * t31;
t256 = m(6) * t39;
t255 = m(6) * t41;
t254 = m(6) * t46;
t253 = m(6) * t278;
t252 = -t177 / 0.2e1;
t251 = t177 / 0.2e1;
t250 = -t178 / 0.2e1;
t63 = t193 - t198 / 0.2e1;
t246 = t63 * qJD(4);
t230 = t177 * t182;
t123 = Icges(6,4) * t230 - Icges(6,2) * t231 - Icges(6,6) * t178;
t234 = t123 * t180;
t152 = Icges(6,5) * t182 - Icges(6,6) * t180;
t232 = t152 * t178;
t228 = t178 * t182;
t121 = Icges(6,5) * t230 - Icges(6,6) * t231 - Icges(6,3) * t178;
t161 = Icges(6,4) * t231;
t125 = Icges(6,1) * t230 - Icges(6,5) * t178 - t161;
t223 = -t177 * t121 - t125 * t228;
t122 = Icges(6,3) * t177 + t232;
t126 = Icges(6,5) * t177 + t156 * t178;
t222 = t177 * t122 + t126 * t228;
t221 = t155 * t177 + t123;
t124 = Icges(6,6) * t177 + t154 * t178;
t220 = -t155 * t178 - t124;
t219 = -Icges(6,2) * t230 + t125 - t161;
t218 = -t153 * t178 + t126;
t208 = t124 * t180 - t121;
t112 = t126 * t230;
t209 = t122 * t178 - t112;
t49 = -t123 * t229 - t223;
t50 = -t124 * t229 + t222;
t11 = (t208 * t178 - t222 + t50) * t178 + (t208 * t177 + t209 + t49) * t177;
t48 = -t124 * t231 - t209;
t12 = (t48 - t112 + (t122 + t234) * t178 + t223) * t178 + t222 * t177;
t23 = t177 * t48 - t178 * (-(-t125 * t182 + t234) * t177 - t121 * t178);
t24 = t177 * t50 - t178 * t49;
t2 = (t24 / 0.2e1 - t12 / 0.2e1) * t178 + (t11 / 0.2e1 + t23 / 0.2e1) * t177;
t215 = -t64 * qJD(4) + t2 * qJD(5);
t204 = t262 / 0.2e1 + t207;
t202 = Icges(6,5) * t180 + Icges(6,6) * t182;
t196 = (-t138 * t178 + t139 * t177) * t157;
t190 = (-t153 + t156) * t182 - t285 * t180;
t195 = t178 * t12 / 0.2e1 + (t11 + t23) * t252 + (t152 * t177 + t190 * t178 + t220 * t180 + t218 * t182) * t251 + (t190 * t177 - t221 * t180 + t219 * t182 - t232 + t24) * t250;
t194 = -t207 + (t251 + t252) * (t182 * t123 + t180 * t125);
t192 = t219 * t180 + t221 * t182;
t191 = -t218 * t180 + t220 * t182;
t158 = -rSges(6,2) * t180 + t243;
t133 = t178 * t202;
t132 = t202 * t177;
t61 = t64 * qJD(5);
t59 = t63 * qJD(5);
t32 = -t253 + t263;
t30 = t207 + t255;
t27 = t207 + t256;
t22 = t254 + t264;
t18 = t261 / 0.2e1;
t13 = t258 + t266 + t270;
t8 = -t261 / 0.2e1 + t204;
t7 = t18 + t204;
t5 = t247 + t248;
t3 = t18 - t262 / 0.2e1 + t194;
t1 = [qJD(3) * t13 + qJD(4) * t22 + qJD(5) * t27, 0, t13 * qJD(1) + t5 * qJD(4) + t7 * qJD(5) + 0.2e1 * (t31 * t271 + t266 / 0.2e1 + t270 / 0.2e1) * qJD(3), qJD(1) * t22 + qJD(3) * t5 + t59, t27 * qJD(1) + t7 * qJD(3) + ((-t46 * t158 + t196) * m(6) + t195) * qJD(5) + t246; 0, 0, 0, 0, -m(6) * t199 * qJD(5); -t4 * qJD(4) + t8 * qJD(5) + (-t258 / 0.4e1 - t266 / 0.4e1 - t270 / 0.4e1) * t275, 0, qJD(4) * t32 + qJD(5) * t30, qJD(3) * t32 - t286 + t59, t8 * qJD(1) + t30 * qJD(3) + ((t278 * t158 + t196) * m(6) + t195) * qJD(5) + t246; t4 * qJD(3) + t61 + (-t254 / 0.4e1 - t264 / 0.4e1) * t275, 0, t286 + t61 + 0.4e1 * (t253 / 0.4e1 - t263 / 0.4e1) * qJD(3), 0, t280; (t194 - t256) * qJD(1) + t3 * qJD(3) + t215, 0, t3 * qJD(1) + (t194 - t255) * qJD(3) + t215, -t280, (m(6) * (t158 * t279 - (t177 * (rSges(6,1) * t230 - t217) + t178 * (rSges(6,1) * t228 + t210)) * t199) + (-t277 * t133 + (t192 * t178 + (t132 + t191) * t177) * t178) * t251 + (-t276 * t132 + (t191 * t177 + (t133 + t192) * t178) * t177) * t250) * qJD(5) + t216 * t2;];
Cq = t1;
