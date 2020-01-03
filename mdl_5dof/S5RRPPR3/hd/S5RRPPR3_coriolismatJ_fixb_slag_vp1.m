% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPPR3_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR3_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR3_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR3_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR3_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR3_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR3_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:26:28
% EndTime: 2019-12-31 19:26:31
% DurationCPUTime: 1.60s
% Computational Cost: add. (12020->183), mult. (7489->255), div. (0->0), fcn. (6440->8), ass. (0->126)
t286 = m(6) / 0.2e1;
t287 = m(5) / 0.2e1;
t188 = qJ(1) + qJ(2);
t185 = pkin(8) + t188;
t183 = sin(t185);
t184 = cos(t185);
t243 = sin(qJ(1)) * pkin(1);
t189 = sin(qJ(5));
t191 = cos(qJ(5));
t204 = rSges(6,1) * t189 + rSges(6,2) * t191;
t193 = -t183 * rSges(6,3) + t184 * t204;
t186 = sin(t188);
t241 = pkin(2) * t186;
t210 = t184 * qJ(4) - t241;
t268 = -pkin(3) - pkin(7);
t99 = t268 * t183 + t193 + t210;
t97 = t99 - t243;
t187 = cos(t188);
t240 = pkin(2) * t187;
t100 = t240 + (rSges(6,3) - t268) * t184 + (qJ(4) + t204) * t183;
t242 = cos(qJ(1)) * pkin(1);
t98 = t100 + t242;
t44 = t98 * t183 + t97 * t184;
t47 = t100 * t183 + t99 * t184;
t284 = pkin(3) - rSges(5,2);
t120 = t184 * rSges(5,3) - t284 * t183 + t210;
t114 = t120 - t243;
t121 = t240 + (rSges(5,3) + qJ(4)) * t183 + t284 * t184;
t115 = t121 + t242;
t61 = t114 * t184 + t115 * t183;
t65 = t120 * t184 + t121 * t183;
t238 = (t47 + t44) * t286 + (t65 + t61) * t287;
t239 = (t44 - t47) * t286 + (t61 - t65) * t287;
t3 = t239 - t238;
t289 = t3 * qJD(1);
t232 = Icges(6,4) * t191;
t164 = -Icges(6,2) * t189 + t232;
t203 = Icges(6,1) * t189 + t232;
t288 = t164 + t203;
t267 = m(3) * (t242 * (-rSges(3,1) * t186 - rSges(3,2) * t187) + (t187 * rSges(3,1) - t186 * rSges(3,2)) * t243);
t266 = m(4) * (t242 * (-rSges(4,1) * t183 - rSges(4,2) * t184 - t241) + (t184 * rSges(4,1) - t183 * rSges(4,2) + t240) * t243);
t262 = m(5) * (-t121 * t114 + t115 * t120);
t254 = m(6) * (-t100 * t97 + t98 * t99);
t200 = Icges(6,5) * t189 + Icges(6,6) * t191;
t283 = t183 * t200;
t282 = t184 * t200;
t169 = rSges(6,1) * t191 - rSges(6,2) * t189;
t148 = t169 * t184;
t147 = t169 * t183;
t231 = t100 * t147;
t39 = t99 * t148 + t231;
t281 = m(6) * t39;
t280 = m(6) * t204;
t233 = Icges(6,4) * t189;
t202 = Icges(6,2) * t191 + t233;
t131 = -Icges(6,6) * t183 + t184 * t202;
t133 = -Icges(6,5) * t183 + t184 * t203;
t279 = (t131 * t191 + t133 * t189) * t184;
t166 = Icges(6,1) * t191 - t233;
t278 = t288 * t189 - (-t202 + t166) * t191;
t215 = t164 * t184 + t133;
t217 = -t166 * t184 + t131;
t277 = t217 * t189 - t215 * t191;
t224 = t183 * t191;
t161 = Icges(6,4) * t224;
t225 = t183 * t189;
t132 = Icges(6,1) * t225 + Icges(6,5) * t184 + t161;
t216 = -Icges(6,2) * t225 + t132 + t161;
t130 = Icges(6,6) * t184 + t183 * t202;
t218 = -t166 * t183 + t130;
t276 = t218 * t189 - t216 * t191;
t207 = -t288 * t191 / 0.2e1 + (t202 / 0.2e1 - t166 / 0.2e1) * t189;
t180 = t183 ^ 2;
t181 = t184 ^ 2;
t275 = 0.4e1 * qJD(1);
t260 = m(5) * t61;
t259 = m(5) * t65;
t234 = t98 * t147;
t38 = t97 * t148 + t234;
t258 = m(6) * (t39 + t38);
t257 = m(6) * (t234 - t231 + (t97 - t99) * t148);
t252 = m(6) * t38;
t250 = m(6) * t44;
t249 = m(6) * t47;
t248 = t183 / 0.2e1;
t247 = -t184 / 0.2e1;
t246 = t184 / 0.2e1;
t101 = -t147 * t184 + t148 * t183;
t244 = m(6) * t101;
t211 = t244 / 0.2e1;
t237 = qJD(4) * t211;
t128 = Icges(6,3) * t184 + t283;
t123 = t183 * t128;
t129 = -Icges(6,3) * t183 + t282;
t199 = -t130 * t191 - t132 * t189;
t54 = t184 * t128 + t130 * t224 + t132 * t225;
t55 = -t184 * t129 - t131 * t224 - t133 * t225;
t56 = t184 * t199 + t123;
t57 = -t129 * t183 + t279;
t11 = (-t56 + t123 + t55) * t183 + (t57 - t279 + (t129 + t199) * t183 + t54) * t184;
t12 = t129 * t180 + (-t123 + t55 + (t129 - t199) * t184) * t184;
t27 = t183 * t55 + t184 * t54;
t28 = t183 * t57 + t184 * t56;
t2 = (t12 / 0.2e1 + t28 / 0.2e1) * t184 + (-t27 / 0.2e1 + t11 / 0.2e1) * t183;
t212 = -qJD(4) * t244 / 0.2e1 + t2 * qJD(5);
t209 = qJD(1) / 0.4e1 + qJD(2) / 0.4e1;
t208 = (t180 + t181) * t204;
t205 = t258 / 0.2e1 + t207;
t201 = Icges(6,5) * t191 - Icges(6,6) * t189;
t141 = t183 * t201;
t195 = -t183 * t11 / 0.2e1 + (t12 + t28) * t247 + (-t278 * t183 - t216 * t189 - t218 * t191 - t282) * t246 + (t278 * t184 + t215 * t189 + t217 * t191 + t27 - t283) * t248;
t194 = -t207 + (t246 + t247) * (t189 * t131 - t191 * t133);
t142 = t201 * t184;
t102 = -t147 * t183 - t148 * t184;
t91 = qJD(5) * t211;
t33 = t207 + t281;
t32 = t207 + t252;
t31 = t249 + t259;
t26 = t250 + t260;
t18 = t257 / 0.2e1;
t13 = t254 + t262 + t266 + t267;
t8 = -t257 / 0.2e1 + t205;
t7 = t18 + t205;
t6 = t18 - t258 / 0.2e1 + t194;
t4 = t238 + t239;
t1 = [t13 * qJD(2) + t26 * qJD(4) + t32 * qJD(5), t13 * qJD(1) + t4 * qJD(4) + t7 * qJD(5) + 0.2e1 * (t254 / 0.2e1 + t262 / 0.2e1 + t266 / 0.2e1 + t267 / 0.2e1) * qJD(2), 0, t26 * qJD(1) + t4 * qJD(2) + t91, t32 * qJD(1) + t7 * qJD(2) + (-(t183 * t97 - t184 * t98) * t280 + t195) * qJD(5) + t237; -t3 * qJD(4) + t8 * qJD(5) + (-t254 / 0.4e1 - t262 / 0.4e1 - t267 / 0.4e1 - t266 / 0.4e1) * t275, qJD(4) * t31 + qJD(5) * t33, 0, t31 * qJD(2) - t289 + t91, t8 * qJD(1) + t33 * qJD(2) + (-(-t100 * t184 + t183 * t99) * t280 + t195) * qJD(5) + t237; 0, 0, 0, 0, m(6) * t102 * qJD(5); t3 * qJD(2) + t91 + (-t250 / 0.4e1 - t260 / 0.4e1) * t275, t289 + t91 + 0.4e1 * (-t259 / 0.4e1 - t249 / 0.4e1) * qJD(2), 0, 0, 0.2e1 * (-qJD(5) * t208 / 0.2e1 + t209 * t101) * m(6); (t194 - t252) * qJD(1) + t6 * qJD(2) + t212, t6 * qJD(1) + (t194 - t281) * qJD(2) + t212, 0, -0.2e1 * t209 * t244, (m(6) * (t102 * (-t184 * t193 + (-t184 * rSges(6,3) - t183 * t204) * t183) - t169 * t208) + (t181 * t141 + (t277 * t183 + (-t142 - t276) * t184) * t183) * t246 + (-t180 * t142 + (t276 * t184 + (t141 - t277) * t183) * t184) * t248) * qJD(5) + (qJD(1) + qJD(2)) * t2;];
Cq = t1;
