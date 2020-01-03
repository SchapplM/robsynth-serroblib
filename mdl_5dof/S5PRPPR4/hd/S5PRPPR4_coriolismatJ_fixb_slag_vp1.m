% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
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
% Datum: 2019-12-31 17:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRPPR4_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR4_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR4_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR4_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR4_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPPR4_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPPR4_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:36:47
% EndTime: 2019-12-31 17:36:49
% DurationCPUTime: 1.19s
% Computational Cost: add. (4536->154), mult. (5618->234), div. (0->0), fcn. (5987->6), ass. (0->98)
t108 = pkin(7) + qJ(2);
t107 = cos(t108);
t109 = sin(pkin(8));
t134 = t107 * t109;
t106 = sin(t108);
t110 = cos(pkin(8));
t154 = sin(qJ(5));
t155 = cos(qJ(5));
t113 = t109 * t154 + t110 * t155;
t85 = t113 * t106;
t80 = Icges(6,4) * t85;
t97 = t109 * t155 - t110 * t154;
t84 = t97 * t106;
t54 = -Icges(6,2) * t84 - Icges(6,6) * t107 - t80;
t79 = Icges(6,4) * t84;
t56 = Icges(6,1) * t85 + Icges(6,5) * t107 + t79;
t183 = -t84 * t54 + t85 * t56;
t86 = t97 * t107;
t87 = t113 * t107;
t151 = -t86 * t54 + t87 * t56;
t140 = Icges(6,4) * t87;
t55 = Icges(6,2) * t86 - Icges(6,6) * t106 + t140;
t81 = Icges(6,4) * t86;
t58 = Icges(6,1) * t87 - Icges(6,5) * t106 + t81;
t150 = t86 * t55 + t87 * t58;
t52 = Icges(6,5) * t87 + Icges(6,6) * t86 - Icges(6,3) * t106;
t121 = t106 * t52 - t150;
t50 = Icges(6,5) * t85 + Icges(6,6) * t84 + Icges(6,3) * t107;
t13 = t107 * t50 + t183;
t181 = t121 - t13;
t135 = t106 * t109;
t100 = t107 * qJ(3);
t175 = qJ(4) * t109 + pkin(2) + (pkin(3) + pkin(4)) * t110;
t177 = -t85 * rSges(6,1) - t84 * rSges(6,2);
t180 = t100 + (-rSges(6,3) - pkin(6)) * t107 - t175 * t106 + t177;
t116 = t87 * rSges(6,1) + t86 * rSges(6,2) - t106 * rSges(6,3);
t34 = t175 * t107 + t116 + (qJ(3) - pkin(6)) * t106;
t66 = rSges(6,1) * t84 - rSges(6,2) * t85;
t67 = rSges(6,1) * t86 - rSges(6,2) * t87;
t139 = Icges(6,4) * t97;
t74 = -Icges(6,2) * t113 + t139;
t75 = -Icges(6,1) * t113 - t139;
t179 = (t75 / 0.2e1 - t74 / 0.2e1) * t97 + m(6) * (-t180 * t66 + t34 * t67);
t178 = t84 * t55 + t85 * t58;
t93 = Icges(6,4) * t113;
t73 = -Icges(6,2) * t97 - t93;
t76 = Icges(6,1) * t97 - t93;
t145 = t73 + t76;
t131 = t106 ^ 2 + t107 ^ 2;
t92 = t131 * t109;
t146 = -Icges(6,2) * t87 + t58 + t81;
t147 = -Icges(6,2) * t85 + t56 + t79;
t174 = t146 * t106 - t147 * t107;
t173 = pkin(2) + (rSges(5,1) + pkin(3)) * t110 + (rSges(5,3) + qJ(4)) * t109;
t172 = -0.2e1 * t92;
t171 = 4 * qJD(2);
t170 = m(6) / 0.2e1;
t164 = m(4) * ((rSges(4,2) * t135 + t107 * rSges(4,3) + t100) * t107 + (-rSges(4,2) * t134 + (rSges(4,3) + qJ(3)) * t106) * t106);
t48 = t107 * rSges(5,2) - t173 * t106 + t100;
t49 = (rSges(5,2) + qJ(3)) * t106 + t173 * t107;
t163 = m(5) * (t49 * t134 - t48 * t135);
t162 = m(5) * (t49 * t106 + t48 * t107);
t161 = m(6) * (t34 * t134 - t135 * t180);
t160 = (t34 * t106 + t180 * t107) * m(6);
t159 = -t106 / 0.2e1;
t157 = t107 / 0.2e1;
t149 = Icges(6,1) * t84 + t54 - t80;
t148 = Icges(6,1) * t86 - t140 - t55;
t144 = -t74 + t75;
t141 = m(6) * qJD(5);
t115 = (t106 * t67 - t107 * t66) * t109;
t10 = -t115 * m(6) / 0.2e1;
t136 = t10 * qJD(4);
t29 = -t106 * t66 - t107 * t67;
t114 = t29 * t170;
t78 = t97 * rSges(6,1) - rSges(6,2) * t113;
t123 = t131 * t78;
t122 = m(6) * t123;
t19 = t114 - t122 / 0.2e1;
t133 = t19 * qJD(2);
t37 = (m(5) / 0.2e1 + t170) * t172;
t132 = t37 * qJD(2);
t119 = t106 * (Icges(6,5) * t86 - Icges(6,6) * t87) - t107 * (Icges(6,5) * t84 - Icges(6,6) * t85);
t112 = -t148 * t106 + t149 * t107;
t77 = -rSges(6,1) * t113 - t97 * rSges(6,2);
t71 = -Icges(6,5) * t113 - Icges(6,6) * t97;
t36 = (m(5) / 0.4e1 + m(6) / 0.4e1) * t172 + (m(5) + m(6)) * t92 / 0.2e1;
t18 = t114 + t122 / 0.2e1;
t11 = t115 * t170;
t8 = t161 + t163;
t7 = t160 + t162 + t164;
t6 = t106 * t121 + (-t106 * t50 + t151) * t107;
t5 = -(t107 * t52 + t178) * t106 + t13 * t107;
t4 = -(t76 / 0.2e1 + t73 / 0.2e1) * t113 + t179;
t3 = t151 * t107 + (t181 + t183) * t106;
t2 = (t150 + t181) * t107 + t178 * t106;
t1 = (-t6 / 0.2e1 + t3 / 0.2e1) * t107 + (-t2 / 0.2e1 - t5 / 0.2e1) * t106;
t9 = [0, 0, 0, 0, t29 * t141; 0, t7 * qJD(3) + t8 * qJD(4) + t4 * qJD(5), qJD(2) * t7 + qJD(4) * t36 + qJD(5) * t18, qJD(2) * t8 + qJD(3) * t36 - qJD(5) * t10, t4 * qJD(2) + t18 * qJD(3) - t136 + (-t107 * t3 / 0.2e1 + (-t106 * t71 - t113 * t146 + t144 * t87 + t145 * t86 + t148 * t97) * t159 + (t2 + t5) * t106 / 0.2e1 + ((t180 * t77 - t66 * t78) * t107 + (t34 * t77 + t67 * t78) * t106) * m(6) + (t107 * t71 - t113 * t147 + t144 * t85 + t145 * t84 + t149 * t97 + t6) * t157) * qJD(5); 0, t37 * qJD(4) + t19 * qJD(5) + (-t164 / 0.4e1 - t162 / 0.4e1 - t160 / 0.4e1) * t171, 0, t132, t133; 0, -t37 * qJD(3) + t11 * qJD(5) + (-t163 / 0.4e1 - t161 / 0.4e1) * t171, -t132, 0, t11 * qJD(2) + (-t29 * t110 + t77 * t92) * t141; 0, -t19 * qJD(3) + t1 * qJD(5) + t136 + (t145 * t113 / 0.2e1 - t179) * qJD(2), -t133, t10 * qJD(2), t1 * qJD(2) + (m(6) * (t77 * t123 + (-t106 * (t107 * rSges(6,3) - t177) - t107 * t116) * t29) + (t119 * t106 + t112 * t87 - t174 * t86) * t159 + (-t119 * t107 + t112 * t85 - t174 * t84) * t157) * qJD(5);];
Cq = t9;
