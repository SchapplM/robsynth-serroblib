% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
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
% Datum: 2019-12-31 16:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RPPR6_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR6_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR6_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR6_coriolismatJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR6_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPR6_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPR6_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:40:38
% EndTime: 2019-12-31 16:40:39
% DurationCPUTime: 1.08s
% Computational Cost: add. (2112->153), mult. (5417->233), div. (0->0), fcn. (5793->6), ass. (0->96)
t105 = sin(pkin(6));
t108 = cos(qJ(1));
t131 = t105 * t108;
t107 = sin(qJ(1));
t106 = cos(pkin(6));
t149 = sin(qJ(4));
t150 = cos(qJ(4));
t111 = t105 * t149 + t106 * t150;
t84 = t111 * t107;
t79 = Icges(5,4) * t84;
t93 = t105 * t150 - t106 * t149;
t83 = t93 * t107;
t51 = -Icges(5,2) * t83 - Icges(5,6) * t108 - t79;
t78 = Icges(5,4) * t83;
t53 = Icges(5,1) * t84 + Icges(5,5) * t108 + t78;
t178 = -t83 * t51 + t84 * t53;
t85 = t93 * t108;
t86 = t111 * t108;
t146 = -t51 * t85 + t86 * t53;
t136 = Icges(5,4) * t86;
t52 = Icges(5,2) * t85 - Icges(5,6) * t107 + t136;
t80 = Icges(5,4) * t85;
t55 = Icges(5,1) * t86 - Icges(5,5) * t107 + t80;
t145 = t85 * t52 + t86 * t55;
t49 = Icges(5,5) * t86 + Icges(5,6) * t85 - Icges(5,3) * t107;
t119 = t107 * t49 - t145;
t47 = Icges(5,5) * t84 + Icges(5,6) * t83 + Icges(5,3) * t108;
t12 = t47 * t108 + t178;
t176 = t119 - t12;
t132 = t105 * t107;
t100 = t108 * qJ(2);
t169 = qJ(3) * t105 + pkin(1) + (pkin(2) + pkin(3)) * t106;
t172 = -t84 * rSges(5,1) - t83 * rSges(5,2);
t175 = t100 + (-rSges(5,3) - pkin(5)) * t108 - t169 * t107 + t172;
t114 = t86 * rSges(5,1) + t85 * rSges(5,2) - t107 * rSges(5,3);
t33 = t169 * t108 + t114 + (qJ(2) - pkin(5)) * t107;
t65 = rSges(5,1) * t83 - rSges(5,2) * t84;
t66 = rSges(5,1) * t85 - rSges(5,2) * t86;
t135 = Icges(5,4) * t93;
t71 = -Icges(5,2) * t111 + t135;
t72 = -Icges(5,1) * t111 - t135;
t174 = (t72 / 0.2e1 - t71 / 0.2e1) * t93 + m(5) * (-t175 * t65 + t33 * t66);
t173 = t83 * t52 + t84 * t55;
t91 = Icges(5,4) * t111;
t70 = -Icges(5,2) * t93 - t91;
t73 = Icges(5,1) * t93 - t91;
t140 = t70 + t73;
t128 = t107 ^ 2 + t108 ^ 2;
t94 = t128 * t105;
t141 = -Icges(5,2) * t86 + t55 + t80;
t142 = -Icges(5,2) * t84 + t53 + t78;
t168 = t107 * t141 - t108 * t142;
t167 = pkin(1) + (rSges(4,1) + pkin(2)) * t106 + (rSges(4,3) + qJ(3)) * t105;
t166 = -0.2e1 * t94;
t165 = 4 * qJD(1);
t164 = m(5) / 0.2e1;
t159 = m(3) * (t108 * (rSges(3,2) * t132 + rSges(3,3) * t108 + t100) + (-rSges(3,2) * t131 + (rSges(3,3) + qJ(2)) * t107) * t107);
t57 = rSges(4,2) * t108 - t167 * t107 + t100;
t58 = (rSges(4,2) + qJ(2)) * t107 + t167 * t108;
t158 = m(4) * (t58 * t131 - t132 * t57);
t157 = m(4) * (t58 * t107 + t108 * t57);
t156 = m(5) * (t33 * t131 - t132 * t175);
t155 = (t33 * t107 + t175 * t108) * m(5);
t154 = -t107 / 0.2e1;
t152 = t108 / 0.2e1;
t144 = Icges(5,1) * t83 + t51 - t79;
t143 = Icges(5,1) * t85 - t136 - t52;
t139 = -t71 + t72;
t113 = (t107 * t66 - t108 * t65) * t105;
t10 = -t113 * m(5) / 0.2e1;
t133 = t10 * qJD(3);
t27 = -t107 * t65 - t108 * t66;
t112 = t27 * t164;
t75 = rSges(5,1) * t93 - rSges(5,2) * t111;
t121 = t128 * t75;
t120 = m(5) * t121;
t19 = t112 - t120 / 0.2e1;
t130 = t19 * qJD(1);
t36 = (m(4) / 0.2e1 + t164) * t166;
t129 = t36 * qJD(1);
t117 = t107 * (Icges(5,5) * t85 - Icges(5,6) * t86) - t108 * (Icges(5,5) * t83 - Icges(5,6) * t84);
t110 = -t107 * t143 + t108 * t144;
t74 = -rSges(5,1) * t111 - rSges(5,2) * t93;
t68 = -Icges(5,5) * t111 - Icges(5,6) * t93;
t35 = (m(4) / 0.4e1 + m(5) / 0.4e1) * t166 + (m(4) + m(5)) * t94 / 0.2e1;
t18 = t112 + t120 / 0.2e1;
t11 = t113 * t164;
t8 = t156 + t158;
t7 = t155 + t157 + t159;
t6 = t107 * t119 + t108 * (-t107 * t47 + t146);
t5 = -(t108 * t49 + t173) * t107 + t108 * t12;
t4 = -(t73 / 0.2e1 + t70 / 0.2e1) * t111 + t174;
t3 = t146 * t108 + (t176 + t178) * t107;
t2 = (t145 + t176) * t108 + t173 * t107;
t1 = (-t6 / 0.2e1 + t3 / 0.2e1) * t108 + (-t2 / 0.2e1 - t5 / 0.2e1) * t107;
t9 = [t7 * qJD(2) + t8 * qJD(3) + t4 * qJD(4), qJD(1) * t7 + qJD(3) * t35 + qJD(4) * t18, qJD(1) * t8 + qJD(2) * t35 - qJD(4) * t10, t4 * qJD(1) + t18 * qJD(2) - t133 + (-t108 * t3 / 0.2e1 + ((t175 * t74 - t65 * t75) * t108 + (t33 * t74 + t66 * t75) * t107) * m(5) + (t108 * t68 - t111 * t142 + t139 * t84 + t140 * t83 + t144 * t93 + t6) * t152 + (-t107 * t68 - t111 * t141 + t139 * t86 + t140 * t85 + t143 * t93) * t154 + (t2 + t5) * t107 / 0.2e1) * qJD(4); t36 * qJD(3) + t19 * qJD(4) + (-t159 / 0.4e1 - t157 / 0.4e1 - t155 / 0.4e1) * t165, 0, t129, t130; -t36 * qJD(2) + t11 * qJD(4) + (-t158 / 0.4e1 - t156 / 0.4e1) * t165, -t129, 0, t11 * qJD(1) + m(5) * (-t106 * t27 + t74 * t94) * qJD(4); -t19 * qJD(2) + t1 * qJD(4) + t133 + (t140 * t111 / 0.2e1 - t174) * qJD(1), -t130, t10 * qJD(1), t1 * qJD(1) + (m(5) * (t121 * t74 + (-t107 * (rSges(5,3) * t108 - t172) - t108 * t114) * t27) + (t117 * t107 + t110 * t86 - t168 * t85) * t154 + (-t117 * t108 + t110 * t84 - t168 * t83) * t152) * qJD(4);];
Cq = t9;
