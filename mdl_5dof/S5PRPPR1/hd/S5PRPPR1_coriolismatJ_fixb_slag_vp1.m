% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
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
% Datum: 2019-12-05 15:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRPPR1_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR1_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR1_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR1_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR1_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPPR1_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPPR1_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:21:47
% EndTime: 2019-12-05 15:21:55
% DurationCPUTime: 2.59s
% Computational Cost: add. (8948->194), mult. (8085->312), div. (0->0), fcn. (8558->8), ass. (0->114)
t145 = pkin(9) + qJ(5);
t141 = sin(t145);
t143 = cos(t145);
t148 = sin(pkin(8));
t150 = cos(pkin(8));
t180 = Icges(6,4) * t141;
t115 = -Icges(6,5) * t150 + (Icges(6,1) * t143 - t180) * t148;
t128 = (-Icges(6,2) * t143 - t180) * t148;
t168 = t115 + t128;
t179 = Icges(6,4) * t143;
t114 = -Icges(6,6) * t150 + (-Icges(6,2) * t141 + t179) * t148;
t129 = (-Icges(6,1) * t141 - t179) * t148;
t169 = -t114 + t129;
t127 = (-Icges(6,5) * t141 - Icges(6,6) * t143) * t148;
t173 = t150 * t127;
t228 = (-t173 + (-t168 * t141 + t169 * t143) * t148) * t150;
t146 = pkin(7) + qJ(2);
t144 = cos(t146);
t142 = sin(t146);
t177 = t142 * t150;
t117 = t141 * t177 + t144 * t143;
t118 = -t144 * t141 + t143 * t177;
t102 = -rSges(6,1) * t117 - rSges(6,2) * t118;
t174 = t144 * t150;
t119 = -t141 * t174 + t142 * t143;
t120 = t142 * t141 + t143 * t174;
t103 = rSges(6,1) * t119 - rSges(6,2) * t120;
t139 = t144 * qJ(3);
t147 = sin(pkin(9));
t176 = t144 * t147;
t215 = -t118 * rSges(6,1) + t117 * rSges(6,2);
t149 = cos(pkin(9));
t219 = (pkin(4) * t149 + pkin(3)) * t150 + pkin(2) + (pkin(6) + qJ(4) + rSges(6,3)) * t148;
t220 = pkin(4) * t176 - t219 * t142 + t139 + t215;
t157 = t120 * rSges(6,1) + t119 * rSges(6,2);
t57 = (pkin(4) * t147 + qJ(3)) * t142 + t157 + t219 * t144;
t16 = (-(t115 / 0.2e1 + t128 / 0.2e1) * t141 + (t129 / 0.2e1 - t114 / 0.2e1) * t143) * t148 + m(6) * (-t102 * t220 + t103 * t57) - t173 / 0.2e1;
t227 = t16 * qJD(2);
t175 = t144 * t148;
t116 = -t150 * rSges(6,3) + (rSges(6,1) * t143 - rSges(6,2) * t141) * t148;
t178 = t142 * t148;
t93 = rSges(6,3) * t178 - t215;
t225 = t116 * t178 + t150 * t93;
t109 = Icges(6,4) * t118;
t87 = -Icges(6,2) * t117 + Icges(6,6) * t178 + t109;
t108 = Icges(6,4) * t117;
t91 = -Icges(6,1) * t118 - Icges(6,5) * t178 + t108;
t224 = -t119 * t87 + t120 * t91;
t84 = Icges(6,5) * t118 - Icges(6,6) * t117 + Icges(6,3) * t178;
t221 = t84 * t175;
t218 = t148 / 0.2e1;
t192 = -t150 / 0.2e1;
t26 = t221 - t224;
t216 = t26 - t221;
t48 = (t102 * t144 - t103 * t142) * t148;
t213 = pkin(2) + (rSges(5,3) + qJ(4)) * t148;
t131 = (-t142 ^ 2 - t144 ^ 2) * t148;
t211 = 0.2e1 * t131;
t209 = 0.4e1 * qJD(2);
t208 = m(5) / 0.2e1;
t206 = m(6) / 0.2e1;
t204 = m(4) * ((rSges(4,2) * t178 + t144 * rSges(4,3) + t139) * t144 + (-rSges(4,2) * t175 + (rSges(4,3) + qJ(3)) * t142) * t142);
t151 = -(t149 * t177 - t176) * rSges(5,1) + (t144 * t149 + t147 * t177) * rSges(5,2) + t139 + (-t150 * pkin(3) - t213) * t142;
t78 = (rSges(5,1) * t147 + rSges(5,2) * t149 + qJ(3)) * t142 + ((rSges(5,1) * t149 - rSges(5,2) * t147 + pkin(3)) * t150 + t213) * t144;
t64 = t78 * t175;
t203 = m(5) * (-t151 * t178 + t64);
t202 = m(5) * (t78 * t142 + t151 * t144);
t51 = t57 * t175;
t199 = m(6) * (-t178 * t220 + t51);
t198 = m(6) * (t57 * t142 + t220 * t144);
t95 = rSges(6,3) * t175 + t157;
t63 = t116 * t175 + t150 * t95;
t197 = m(6) * (-t63 * t175 - t178 * t225);
t196 = m(6) * (-t63 * t142 + t225 * t144);
t195 = m(6) * t48;
t194 = m(6) * (-t102 * t142 - t103 * t144);
t189 = -Icges(6,2) * t118 - t108 - t91;
t110 = Icges(6,4) * t119;
t92 = Icges(6,1) * t120 + Icges(6,5) * t175 + t110;
t188 = -Icges(6,2) * t120 + t110 + t92;
t187 = m(6) * qJD(5);
t183 = -Icges(6,1) * t117 - t109 - t87;
t181 = Icges(6,4) * t120;
t89 = Icges(6,2) * t119 + Icges(6,6) * t175 + t181;
t182 = Icges(6,1) * t119 - t181 - t89;
t75 = (t208 + t206) * t211;
t171 = t75 * qJD(2);
t27 = t119 * t89 + t120 * t92 + (Icges(6,5) * t120 + Icges(6,6) * t119 + Icges(6,3) * t175) * t175;
t166 = -(t142 * t26 + t144 * t27) * t148 / 0.2e1 + ((t84 * t178 + t27) * t144 + t216 * t142) * t218;
t165 = (t216 + t224) * t218 * t144 + (t150 / 0.2e1 + t192) * ((-Icges(6,3) * t150 + (Icges(6,5) * t143 - Icges(6,6) * t141) * t148) * t178 - t114 * t117 + t115 * t118);
t130 = (-rSges(6,1) * t141 - rSges(6,2) * t143) * t148;
t97 = Icges(6,5) * t119 - Icges(6,6) * t120;
t96 = -Icges(6,5) * t117 - Icges(6,6) * t118;
t80 = -t103 * t150 - t130 * t175;
t79 = t102 * t150 + t130 * t178;
t74 = (m(5) / 0.4e1 + m(6) / 0.4e1) * t211 - (m(5) + m(6)) * t131 / 0.2e1;
t50 = t194 / 0.2e1;
t47 = -t195 / 0.2e1;
t34 = t196 / 0.2e1;
t32 = t168 * t119 + t169 * t120 + t127 * t175;
t31 = -t168 * t117 + t169 * t118 + t127 * t178;
t30 = t197 / 0.2e1;
t23 = -t150 * t97 + (-t188 * t141 + t182 * t143) * t148;
t22 = -t150 * t96 + (-t189 * t141 + t183 * t143) * t148;
t17 = t199 + t203;
t15 = t198 + t202 + t204;
t9 = t34 - t194 / 0.2e1;
t8 = t50 + t34;
t7 = t50 - t196 / 0.2e1;
t6 = t30 + t195 / 0.2e1;
t5 = t47 + t30;
t4 = t47 - t197 / 0.2e1;
t1 = (t166 * t142 + t165 * t144) * t148;
t2 = [0, 0, 0, 0, t48 * t187; 0, t15 * qJD(3) + t17 * qJD(4) + t16 * qJD(5), qJD(2) * t15 + qJD(4) * t74 + qJD(5) * t8, qJD(2) * t17 + qJD(3) * t74 + qJD(5) * t5, t227 + t8 * qJD(3) + t5 * qJD(4) + (-t228 + (-t102 * t225 - t63 * t103 + t220 * t79 + t80 * t57) * m(6) + ((t23 / 0.2e1 + t32 / 0.2e1 - t165) * t144 + (t22 / 0.2e1 + t31 / 0.2e1 - t166) * t142) * t148) * qJD(5); 0, t75 * qJD(4) + t7 * qJD(5) + (-t198 / 0.4e1 - t202 / 0.4e1 - t204 / 0.4e1) * t209, 0, t171, t7 * qJD(2) + (t142 * t79 - t144 * t80) * t187; 0, -t75 * qJD(3) + t4 * qJD(5) + (-t199 / 0.4e1 - t203 / 0.4e1) * t209 + 0.2e1 * (t51 * t206 + t64 * t208 + (-t57 * t206 - t78 * t208) * t175) * qJD(2), -t171, 0, t4 * qJD(2) + (-t48 * t150 + (t142 * t80 + t144 * t79) * t148) * t187; 0, t9 * qJD(3) + t6 * qJD(4) + t1 * qJD(5) - t227, t9 * qJD(2), t6 * qJD(2), t1 * qJD(2) + (m(6) * (t225 * t79 + (-t142 * t95 + t144 * t93) * t148 * t48 - t63 * t80) + ((t188 * t119 + t182 * t120 + t97 * t175) * t175 + (t189 * t119 + t183 * t120 + t96 * t175) * t178 - t32 * t150) * t175 / 0.2e1 + ((-t188 * t117 + t182 * t118 + t97 * t178) * t175 + (-t189 * t117 + t183 * t118 + t96 * t178) * t178 - t31 * t150) * t178 / 0.2e1 + (-t228 + (t142 * t22 + t144 * t23) * t148) * t192) * qJD(5);];
Cq = t2;
