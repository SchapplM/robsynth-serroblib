% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
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
% Datum: 2019-12-31 17:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPPRP3_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP3_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP3_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP3_coriolismatJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP3_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP3_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRP3_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:50:44
% EndTime: 2019-12-31 17:50:49
% DurationCPUTime: 2.48s
% Computational Cost: add. (4980->176), mult. (5197->257), div. (0->0), fcn. (4570->6), ass. (0->116)
t261 = Icges(5,5) + Icges(6,5);
t260 = Icges(5,6) + Icges(6,6);
t163 = cos(qJ(4));
t161 = sin(qJ(4));
t200 = Icges(6,4) * t161;
t175 = Icges(6,2) * t163 + t200;
t202 = Icges(5,4) * t161;
t176 = Icges(5,2) * t163 + t202;
t262 = t175 + t176;
t160 = qJ(1) + pkin(7);
t158 = sin(t160);
t197 = t158 * t163;
t259 = (-Icges(5,4) - Icges(6,4)) * t197;
t159 = cos(t160);
t255 = t262 * t158 + t260 * t159;
t198 = t158 * t161;
t253 = -t259 + (Icges(5,1) + Icges(6,1)) * t198 + t261 * t159;
t258 = Icges(6,3) + Icges(5,3);
t212 = rSges(6,1) + pkin(4);
t251 = t163 * rSges(6,2) + t212 * t161;
t104 = t251 * t159;
t257 = t158 * rSges(6,3) - t104;
t171 = Icges(6,5) * t161 + Icges(6,6) * t163;
t173 = Icges(5,5) * t161 + Icges(5,6) * t163;
t256 = t171 + t173;
t254 = t260 * t158 - t262 * t159;
t199 = Icges(6,4) * t163;
t177 = Icges(6,1) * t161 + t199;
t201 = Icges(5,4) * t163;
t178 = Icges(5,1) * t161 + t201;
t252 = (-t177 - t178) * t159 + t261 * t158;
t245 = -t258 * t158 + t256 * t159;
t243 = t253 * t161 + t255 * t163;
t250 = (t252 * t161 + t254 * t163) * t159;
t229 = m(5) / 0.2e1;
t228 = m(6) / 0.2e1;
t188 = -sin(qJ(1)) * pkin(1) + t159 * qJ(3);
t208 = -qJ(5) - pkin(6);
t55 = (-pkin(2) + t208) * t158 + t188 - t257;
t145 = t159 * t208;
t211 = cos(qJ(1)) * pkin(1);
t56 = t211 - t145 + (rSges(6,3) + pkin(2)) * t159 + (qJ(3) + t251) * t158;
t249 = m(6) * (-t158 * t55 + t159 * t56);
t248 = t253 * t198 + t255 * t197 + (t256 * t158 + t258 * t159) * t159;
t247 = -t245 * t159 + t254 * t197 + t252 * t198;
t246 = -t245 * t158 - t250;
t244 = -t260 * t161 + t261 * t163;
t242 = t243 * t159;
t156 = t158 ^ 2;
t241 = t156 / 0.2e1;
t214 = t158 / 0.2e1;
t239 = -t159 / 0.2e1;
t213 = t159 / 0.2e1;
t237 = qJD(1) * t249;
t102 = t251 * t158;
t132 = Icges(6,1) * t163 - t200;
t134 = Icges(5,1) * t163 - t202;
t234 = t132 + t134;
t192 = qJD(4) * t158;
t139 = rSges(5,1) * t163 - rSges(5,2) * t161;
t119 = t139 * t158;
t120 = t139 * t159;
t144 = pkin(4) * t197;
t185 = rSges(6,1) * t197 - rSges(6,2) * t198;
t86 = t144 + t185;
t204 = rSges(6,2) * t161;
t87 = (t212 * t163 - t204) * t159;
t207 = (t158 * t87 - t159 * t86) * t228 + (-t119 * t159 + t120 * t158) * t229;
t128 = -Icges(6,2) * t161 + t199;
t130 = -Icges(5,2) * t161 + t201;
t233 = (t134 / 0.2e1 - t176 / 0.2e1 + t132 / 0.2e1 - t175 / 0.2e1) * t161 - (-t178 / 0.2e1 - t130 / 0.2e1 - t177 / 0.2e1 - t128 / 0.2e1) * t163;
t157 = t159 ^ 2;
t231 = 0.4e1 * qJD(1);
t230 = 0.2e1 * qJD(4);
t227 = pkin(2) + pkin(6);
t226 = m(4) * (t159 * (rSges(4,3) * t159 + t188) + (t211 + (rSges(4,3) + qJ(3)) * t158) * t158);
t184 = rSges(5,1) * t161 + rSges(5,2) * t163;
t165 = -t158 * rSges(5,3) + t184 * t159;
t63 = -t227 * t158 + t165 + t188;
t64 = t211 + (rSges(5,3) + t227) * t159 + (qJ(3) + t184) * t158;
t225 = m(5) * (t119 * t64 + t120 * t63);
t224 = m(5) * (t64 * t158 + t159 * t63);
t222 = m(6) * (t55 * t87 + t56 * t86);
t221 = m(6) * (t158 * t56 + t159 * t55);
t219 = m(6) * (t158 * t86 + t159 * t87);
t138 = rSges(6,1) * t163 - t204;
t103 = t138 * t158 + t144;
t209 = pkin(4) * t163;
t105 = (t138 + t209) * t159;
t218 = m(6) * (t103 * t159 - t105 * t158);
t216 = m(6) * (-t103 * t158 - t159 * t105);
t205 = m(6) * qJD(4);
t74 = 0.2e1 * (t241 + t157 / 0.2e1) * m(6);
t194 = t74 * qJD(1);
t193 = t156 + t157;
t191 = t244 * t214;
t190 = t244 * t239;
t189 = -t173 / 0.2e1 - t171 / 0.2e1;
t186 = t193 * t184;
t170 = -t247 * t158 / 0.2e1 + t248 * t239 + (t242 + t247) * t214 + ((-t243 + t245) * t158 + t246 + t248 + t250) * t213;
t169 = t245 * t241 + t246 * t214 + ((t243 + t245) * t159 - t242 + t247) * t213;
t62 = -t119 * t158 - t120 * t159;
t45 = t216 / 0.2e1;
t43 = t218 / 0.2e1;
t39 = -t138 * t157 - t158 * t185 - t193 * t209;
t38 = t219 / 0.2e1;
t18 = t45 - t219 / 0.2e1;
t17 = t45 + t38;
t16 = t38 - t216 / 0.2e1;
t11 = -t218 / 0.2e1 + t207;
t10 = t43 + t207;
t9 = t43 - t207;
t8 = t221 + t224 + t226;
t6 = t222 + t225 - t233;
t1 = t170 * t158 + t169 * t159;
t2 = [t8 * qJD(3) + t6 * qJD(4) + qJD(5) * t249, 0, qJD(1) * t8 + qJD(4) * t10, t6 * qJD(1) + t10 * qJD(3) + (-t102 * t55 + t103 * t87 + t104 * t56 - t105 * t86) * t205 + t17 * qJD(5) + (m(5) * (t120 * t139 - t184 * t63) + t189 * t158 - t170) * t192 + ((m(5) * (-t119 * t139 + t184 * t64) + t189 * t159 - t169) * t159 + (((t128 + t130) * t159 - t252) * t214 + ((Icges(5,2) + Icges(6,2)) * t198 - t253 + t259) * t213) * t161 + ((-t234 * t159 - t254) * t214 + (t234 * t158 - t255) * t213) * t163) * qJD(4), t17 * qJD(4) + t237; 0, 0, 0, (t39 * t228 + t62 * t229) * t230, 0; t11 * qJD(4) - t74 * qJD(5) + (-t226 / 0.4e1 - t224 / 0.4e1 - t221 / 0.4e1) * t231, 0, 0, t11 * qJD(1) + ((-t102 * t158 - t104 * t159) * t228 - t186 * t229) * t230, -t194; t9 * qJD(3) + t1 * qJD(4) + t18 * qJD(5) + (-t225 / 0.4e1 - t222 / 0.4e1) * t231 + t233 * qJD(1), 0, t9 * qJD(1), t1 * qJD(1) + (t190 * t156 + (t191 * t158 + t190 * t159) * t159) * t192 + (m(5) * (-t139 * t186 + (-t159 * t165 + (-t159 * rSges(5,3) - t184 * t158) * t158) * t62) + m(6) * (-t102 * t103 - t104 * t105 + ((t145 - t102) * t158 + ((-rSges(6,3) - t208) * t158 + t257) * t159) * t39) + t191 * t157 * t159) * qJD(4), t18 * qJD(1); t74 * qJD(3) + t16 * qJD(4) - t237, 0, t194, t16 * qJD(1) + (-t102 * t159 + t104 * t158) * t205, 0;];
Cq = t2;
