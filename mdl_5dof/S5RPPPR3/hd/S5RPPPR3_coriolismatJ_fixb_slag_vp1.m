% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3]';
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
% Datum: 2019-12-31 17:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPPPR3_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR3_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR3_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR3_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR3_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR3_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPPR3_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:43:52
% EndTime: 2019-12-31 17:43:54
% DurationCPUTime: 1.23s
% Computational Cost: add. (4614->158), mult. (5704->236), div. (0->0), fcn. (6067->8), ass. (0->99)
t110 = qJ(1) + pkin(7);
t109 = cos(t110);
t111 = sin(pkin(8));
t139 = t109 * t111;
t108 = sin(t110);
t112 = cos(pkin(8));
t160 = sin(qJ(5));
t161 = cos(qJ(5));
t117 = t111 * t160 + t112 * t161;
t85 = t117 * t108;
t80 = Icges(6,4) * t85;
t97 = t111 * t161 - t112 * t160;
t84 = t97 * t108;
t54 = -Icges(6,2) * t84 - Icges(6,6) * t109 - t80;
t79 = Icges(6,4) * t84;
t56 = Icges(6,1) * t85 + Icges(6,5) * t109 + t79;
t189 = -t84 * t54 + t56 * t85;
t86 = t97 * t109;
t87 = t117 * t109;
t156 = -t54 * t86 + t56 * t87;
t145 = Icges(6,4) * t87;
t55 = Icges(6,2) * t86 - Icges(6,6) * t108 + t145;
t81 = Icges(6,4) * t86;
t58 = Icges(6,1) * t87 - Icges(6,5) * t108 + t81;
t155 = t55 * t86 + t58 * t87;
t52 = Icges(6,5) * t87 + Icges(6,6) * t86 - Icges(6,3) * t108;
t125 = t108 * t52 - t155;
t50 = Icges(6,5) * t85 + Icges(6,6) * t84 + Icges(6,3) * t109;
t14 = t109 * t50 + t189;
t187 = t125 - t14;
t140 = t108 * t111;
t129 = -sin(qJ(1)) * pkin(1) + t109 * qJ(3);
t181 = qJ(4) * t111 + pkin(2) + (pkin(3) + pkin(4)) * t112;
t183 = -t85 * rSges(6,1) - t84 * rSges(6,2);
t186 = (-rSges(6,3) - pkin(6)) * t109 - t181 * t108 + t129 + t183;
t120 = t87 * rSges(6,1) + t86 * rSges(6,2) - rSges(6,3) * t108;
t159 = cos(qJ(1)) * pkin(1);
t34 = t159 + t181 * t109 + t120 + (qJ(3) - pkin(6)) * t108;
t66 = rSges(6,1) * t84 - rSges(6,2) * t85;
t67 = rSges(6,1) * t86 - rSges(6,2) * t87;
t144 = Icges(6,4) * t97;
t74 = -Icges(6,2) * t117 + t144;
t75 = -Icges(6,1) * t117 - t144;
t185 = (t75 / 0.2e1 - t74 / 0.2e1) * t97 + m(6) * (-t186 * t66 + t34 * t67);
t184 = t55 * t84 + t58 * t85;
t95 = Icges(6,4) * t117;
t73 = -Icges(6,2) * t97 - t95;
t76 = Icges(6,1) * t97 - t95;
t150 = t73 + t76;
t136 = t108 ^ 2 + t109 ^ 2;
t94 = t136 * t111;
t151 = -Icges(6,2) * t87 + t58 + t81;
t152 = -Icges(6,2) * t85 + t56 + t79;
t180 = t108 * t151 - t109 * t152;
t179 = pkin(2) + (rSges(5,1) + pkin(3)) * t112 + (rSges(5,3) + qJ(4)) * t111;
t178 = -0.2e1 * t94;
t177 = 4 * qJD(1);
t176 = m(6) / 0.2e1;
t170 = m(4) * ((rSges(4,2) * t140 + t109 * rSges(4,3) + t129) * t109 + (t159 - rSges(4,2) * t139 + (rSges(4,3) + qJ(3)) * t108) * t108);
t47 = t109 * rSges(5,2) - t108 * t179 + t129;
t48 = t159 + (rSges(5,2) + qJ(3)) * t108 + t179 * t109;
t169 = m(5) * (t139 * t48 - t140 * t47);
t168 = m(5) * (t108 * t48 + t109 * t47);
t167 = m(6) * (t34 * t139 - t140 * t186);
t166 = (t108 * t34 + t186 * t109) * m(6);
t165 = -t108 / 0.2e1;
t163 = t109 / 0.2e1;
t154 = Icges(6,1) * t84 + t54 - t80;
t153 = Icges(6,1) * t86 - t145 - t55;
t149 = -t74 + t75;
t146 = m(6) * qJD(5);
t119 = (t108 * t67 - t109 * t66) * t111;
t10 = -t119 * m(6) / 0.2e1;
t141 = t10 * qJD(4);
t31 = -t108 * t66 - t109 * t67;
t118 = t31 * t176;
t78 = rSges(6,1) * t97 - rSges(6,2) * t117;
t127 = t136 * t78;
t126 = m(6) * t127;
t19 = t118 - t126 / 0.2e1;
t138 = t19 * qJD(1);
t37 = (m(5) / 0.2e1 + t176) * t178;
t137 = t37 * qJD(1);
t123 = t108 * (Icges(6,5) * t86 - Icges(6,6) * t87) - t109 * (Icges(6,5) * t84 - Icges(6,6) * t85);
t116 = -t108 * t153 + t109 * t154;
t77 = -rSges(6,1) * t117 - rSges(6,2) * t97;
t71 = -Icges(6,5) * t117 - Icges(6,6) * t97;
t36 = (m(5) / 0.4e1 + m(6) / 0.4e1) * t178 + (m(5) + m(6)) * t94 / 0.2e1;
t18 = t118 + t126 / 0.2e1;
t11 = t119 * t176;
t8 = t167 + t169;
t7 = t108 * t125 + (-t108 * t50 + t156) * t109;
t6 = -(t109 * t52 + t184) * t108 + t14 * t109;
t5 = t166 + t168 + t170;
t4 = -(t76 / 0.2e1 + t73 / 0.2e1) * t117 + t185;
t3 = t156 * t109 + (t187 + t189) * t108;
t2 = (t155 + t187) * t109 + t184 * t108;
t1 = (-t7 / 0.2e1 + t3 / 0.2e1) * t109 + (-t2 / 0.2e1 - t6 / 0.2e1) * t108;
t9 = [t5 * qJD(3) + t8 * qJD(4) + t4 * qJD(5), 0, qJD(1) * t5 + qJD(4) * t36 + qJD(5) * t18, qJD(1) * t8 + qJD(3) * t36 - qJD(5) * t10, t4 * qJD(1) + t18 * qJD(3) - t141 + (-t109 * t3 / 0.2e1 + (-t108 * t71 - t117 * t151 + t149 * t87 + t150 * t86 + t153 * t97) * t165 + (t2 + t6) * t108 / 0.2e1 + ((t186 * t77 - t66 * t78) * t109 + (t34 * t77 + t67 * t78) * t108) * m(6) + (t109 * t71 - t117 * t152 + t149 * t85 + t150 * t84 + t154 * t97 + t7) * t163) * qJD(5); 0, 0, 0, 0, t31 * t146; t37 * qJD(4) + t19 * qJD(5) + (-t170 / 0.4e1 - t168 / 0.4e1 - t166 / 0.4e1) * t177, 0, 0, t137, t138; -t37 * qJD(3) + t11 * qJD(5) + (-t169 / 0.4e1 - t167 / 0.4e1) * t177, 0, -t137, 0, t11 * qJD(1) + (-t31 * t112 + t77 * t94) * t146; -t19 * qJD(3) + t1 * qJD(5) + t141 + (t150 * t117 / 0.2e1 - t185) * qJD(1), 0, -t138, t10 * qJD(1), t1 * qJD(1) + (m(6) * (t127 * t77 + (-t108 * (t109 * rSges(6,3) - t183) - t109 * t120) * t31) + (t123 * t108 + t116 * t87 - t180 * t86) * t165 + (-t123 * t109 + t116 * t85 - t180 * t84) * t163) * qJD(5);];
Cq = t9;
