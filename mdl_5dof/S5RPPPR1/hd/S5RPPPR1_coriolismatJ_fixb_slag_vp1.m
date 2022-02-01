% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% Datum: 2022-01-20 09:13
% Revision: 008671b0a00594318b890887636eaaff83cd5e2f (2021-12-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPPPR1_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR1_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR1_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR1_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPPR1_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:12:14
% EndTime: 2022-01-20 09:12:21
% DurationCPUTime: 2.70s
% Computational Cost: add. (9030->198), mult. (8175->314), div. (0->0), fcn. (8642->10), ass. (0->115)
t147 = pkin(9) + qJ(5);
t143 = sin(t147);
t145 = cos(t147);
t150 = sin(pkin(8));
t152 = cos(pkin(8));
t185 = Icges(6,4) * t143;
t115 = -Icges(6,5) * t152 + (Icges(6,1) * t145 - t185) * t150;
t128 = (-Icges(6,2) * t145 - t185) * t150;
t173 = t115 + t128;
t184 = Icges(6,4) * t145;
t114 = -Icges(6,6) * t152 + (-Icges(6,2) * t143 + t184) * t150;
t129 = (-Icges(6,1) * t143 - t184) * t150;
t174 = -t114 + t129;
t127 = (-Icges(6,5) * t143 - Icges(6,6) * t145) * t150;
t178 = t152 * t127;
t234 = (-t178 + (-t173 * t143 + t174 * t145) * t150) * t152;
t148 = qJ(1) + pkin(7);
t146 = cos(t148);
t144 = sin(t148);
t182 = t144 * t152;
t117 = t143 * t182 + t146 * t145;
t118 = -t146 * t143 + t145 * t182;
t102 = -t117 * rSges(6,1) - t118 * rSges(6,2);
t179 = t146 * t152;
t119 = -t143 * t179 + t144 * t145;
t120 = t144 * t143 + t145 * t179;
t103 = t119 * rSges(6,1) - t120 * rSges(6,2);
t163 = -sin(qJ(1)) * pkin(1) + t146 * qJ(3);
t149 = sin(pkin(9));
t181 = t146 * t149;
t221 = -t118 * rSges(6,1) + t117 * rSges(6,2);
t151 = cos(pkin(9));
t225 = (t151 * pkin(4) + pkin(3)) * t152 + pkin(2) + (pkin(6) + qJ(4) + rSges(6,3)) * t150;
t226 = pkin(4) * t181 - t225 * t144 + t163 + t221;
t161 = t120 * rSges(6,1) + t119 * rSges(6,2);
t197 = cos(qJ(1)) * pkin(1);
t57 = t197 + (pkin(4) * t149 + qJ(3)) * t144 + t161 + t225 * t146;
t16 = (-(t115 / 0.2e1 + t128 / 0.2e1) * t143 + (t129 / 0.2e1 - t114 / 0.2e1) * t145) * t150 - t178 / 0.2e1 + m(6) * (-t102 * t226 + t57 * t103);
t233 = t16 * qJD(1);
t180 = t146 * t150;
t116 = -t152 * rSges(6,3) + (rSges(6,1) * t145 - rSges(6,2) * t143) * t150;
t183 = t144 * t150;
t93 = rSges(6,3) * t183 - t221;
t231 = t116 * t183 + t152 * t93;
t109 = Icges(6,4) * t118;
t87 = -Icges(6,2) * t117 + Icges(6,6) * t183 + t109;
t108 = Icges(6,4) * t117;
t91 = -Icges(6,1) * t118 - Icges(6,5) * t183 + t108;
t230 = -t119 * t87 + t120 * t91;
t84 = Icges(6,5) * t118 - Icges(6,6) * t117 + Icges(6,3) * t183;
t227 = t84 * t180;
t224 = t150 / 0.2e1;
t198 = -t152 / 0.2e1;
t26 = t227 - t230;
t222 = t26 - t227;
t48 = (t102 * t146 - t103 * t144) * t150;
t219 = pkin(2) + (rSges(5,3) + qJ(4)) * t150;
t133 = (-t144 ^ 2 - t146 ^ 2) * t150;
t217 = 0.2e1 * t133;
t215 = 0.4e1 * qJD(1);
t214 = m(5) / 0.2e1;
t212 = m(6) / 0.2e1;
t210 = m(4) * ((rSges(4,2) * t183 + t146 * rSges(4,3) + t163) * t146 + (t197 - rSges(4,2) * t180 + (rSges(4,3) + qJ(3)) * t144) * t144);
t155 = -(t151 * t182 - t181) * rSges(5,1) + (t146 * t151 + t149 * t182) * rSges(5,2) + (-t152 * pkin(3) - t219) * t144 + t163;
t70 = t197 + (t149 * rSges(5,1) + t151 * rSges(5,2) + qJ(3)) * t144 + ((rSges(5,1) * t151 - rSges(5,2) * t149 + pkin(3)) * t152 + t219) * t146;
t64 = t70 * t180;
t209 = m(5) * (-t155 * t183 + t64);
t208 = m(5) * (t70 * t144 + t155 * t146);
t50 = t57 * t180;
t205 = m(6) * (-t183 * t226 + t50);
t204 = m(6) * (t57 * t144 + t226 * t146);
t95 = rSges(6,3) * t180 + t161;
t63 = t116 * t180 + t152 * t95;
t203 = m(6) * (-t63 * t180 - t183 * t231);
t202 = m(6) * (-t63 * t144 + t231 * t146);
t201 = m(6) * t48;
t200 = m(6) * (-t144 * t102 - t146 * t103);
t194 = -Icges(6,2) * t118 - t108 - t91;
t110 = Icges(6,4) * t119;
t92 = Icges(6,1) * t120 + Icges(6,5) * t180 + t110;
t193 = -Icges(6,2) * t120 + t110 + t92;
t192 = m(6) * qJD(5);
t188 = -Icges(6,1) * t117 - t109 - t87;
t186 = Icges(6,4) * t120;
t89 = Icges(6,2) * t119 + Icges(6,6) * t180 + t186;
t187 = Icges(6,1) * t119 - t186 - t89;
t78 = (t212 + t214) * t217;
t176 = t78 * qJD(1);
t27 = t119 * t89 + t120 * t92 + (Icges(6,5) * t120 + Icges(6,6) * t119 + Icges(6,3) * t180) * t180;
t171 = -(t144 * t26 + t146 * t27) * t150 / 0.2e1 + ((t84 * t183 + t27) * t146 + t222 * t144) * t224;
t170 = (t222 + t230) * t224 * t146 + (t152 / 0.2e1 + t198) * ((-Icges(6,3) * t152 + (Icges(6,5) * t145 - Icges(6,6) * t143) * t150) * t183 - t117 * t114 + t118 * t115);
t130 = (-rSges(6,1) * t143 - rSges(6,2) * t145) * t150;
t97 = Icges(6,5) * t119 - Icges(6,6) * t120;
t96 = -Icges(6,5) * t117 - Icges(6,6) * t118;
t80 = -t152 * t103 - t130 * t180;
t79 = t152 * t102 + t130 * t183;
t77 = (m(6) / 0.4e1 + m(5) / 0.4e1) * t217 - (m(6) + m(5)) * t133 / 0.2e1;
t51 = t200 / 0.2e1;
t47 = -t201 / 0.2e1;
t34 = t202 / 0.2e1;
t32 = t173 * t119 + t174 * t120 + t127 * t180;
t31 = -t173 * t117 + t174 * t118 + t127 * t183;
t30 = t203 / 0.2e1;
t23 = -t152 * t97 + (-t193 * t143 + t187 * t145) * t150;
t22 = -t152 * t96 + (-t194 * t143 + t188 * t145) * t150;
t17 = t205 + t209;
t15 = t204 + t208 + t210;
t9 = t34 - t200 / 0.2e1;
t8 = t51 + t34;
t7 = t51 - t202 / 0.2e1;
t6 = t30 + t201 / 0.2e1;
t5 = t47 + t30;
t4 = t47 - t203 / 0.2e1;
t1 = (t171 * t144 + t170 * t146) * t150;
t2 = [t15 * qJD(3) + t17 * qJD(4) + t16 * qJD(5), 0, t15 * qJD(1) + t77 * qJD(4) + t8 * qJD(5), t17 * qJD(1) + t77 * qJD(3) + t5 * qJD(5), t233 + t8 * qJD(3) + t5 * qJD(4) + (-t234 + (-t102 * t231 - t63 * t103 + t226 * t79 + t80 * t57) * m(6) + ((t23 / 0.2e1 + t32 / 0.2e1 - t170) * t146 + (t22 / 0.2e1 + t31 / 0.2e1 - t171) * t144) * t150) * qJD(5); 0, 0, 0, 0, t48 * t192; t78 * qJD(4) + t7 * qJD(5) + (-t204 / 0.4e1 - t208 / 0.4e1 - t210 / 0.4e1) * t215, 0, 0, t176, t7 * qJD(1) + (t79 * t144 - t80 * t146) * t192; -t78 * qJD(3) + t4 * qJD(5) + (-t205 / 0.4e1 - t209 / 0.4e1) * t215 + 0.2e1 * (t50 * t212 + t64 * t214 + (-t57 * t212 - t70 * t214) * t180) * qJD(1), 0, -t176, 0, t4 * qJD(1) + (-t48 * t152 + (t144 * t80 + t146 * t79) * t150) * t192; t9 * qJD(3) + t6 * qJD(4) + t1 * qJD(5) - t233, 0, t9 * qJD(1), t6 * qJD(1), t1 * qJD(1) + (m(6) * (t231 * t79 + (-t144 * t95 + t146 * t93) * t150 * t48 - t63 * t80) + ((t193 * t119 + t187 * t120 + t97 * t180) * t180 + (t194 * t119 + t188 * t120 + t96 * t180) * t183 - t32 * t152) * t180 / 0.2e1 + ((-t193 * t117 + t187 * t118 + t97 * t183) * t180 + (-t194 * t117 + t188 * t118 + t96 * t183) * t183 - t31 * t152) * t183 / 0.2e1 + (-t234 + (t144 * t22 + t146 * t23) * t150) * t198) * qJD(5);];
Cq = t2;
