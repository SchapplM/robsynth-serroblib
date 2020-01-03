% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-31 17:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PPRPR4_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR4_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR4_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR4_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR4_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRPR4_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRPR4_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:32:17
% EndTime: 2019-12-31 17:32:19
% DurationCPUTime: 1.17s
% Computational Cost: add. (3264->133), mult. (5265->223), div. (0->0), fcn. (6208->8), ass. (0->101)
t94 = pkin(8) + qJ(5);
t92 = sin(t94);
t93 = cos(t94);
t104 = Icges(6,5) * t93 - Icges(6,6) * t92;
t116 = sin(pkin(7));
t117 = cos(pkin(7));
t137 = sin(qJ(3));
t138 = cos(qJ(3));
t85 = -t116 * t137 - t117 * t138;
t86 = -t116 * t138 + t117 * t137;
t49 = -Icges(6,3) * t85 - t104 * t86;
t161 = t86 * t49;
t121 = Icges(6,4) * t93;
t106 = -Icges(6,2) * t92 + t121;
t52 = -Icges(6,6) * t85 - t106 * t86;
t134 = t52 * t92;
t122 = Icges(6,4) * t92;
t108 = Icges(6,1) * t93 - t122;
t55 = -Icges(6,5) * t85 - t108 * t86;
t159 = -t55 * t93 + t134;
t51 = -Icges(6,6) * t86 + t106 * t85;
t54 = -Icges(6,5) * t86 + t108 * t85;
t156 = -t51 * t92 + t54 * t93;
t96 = -pkin(6) - qJ(4);
t128 = -rSges(6,3) + t96;
t136 = rSges(6,1) * t93;
t130 = t86 * t92;
t72 = rSges(6,2) * t130;
t95 = cos(pkin(8));
t91 = pkin(4) * t95 + pkin(3);
t34 = t72 + (-t91 - t136) * t86 + t128 * t85;
t135 = rSges(6,2) * t92;
t84 = t135 - t136;
t35 = t128 * t86 + (-t84 + t91) * t85;
t109 = rSges(6,1) * t92 + rSges(6,2) * t93;
t66 = t109 * t86;
t67 = t109 * t85;
t157 = (t108 / 0.2e1 - Icges(6,2) * t93 / 0.2e1 - t122 / 0.2e1) * t92 + m(6) * (t34 * t66 - t35 * t67);
t28 = t85 * t34;
t14 = -t35 * t86 - t28;
t155 = -rSges(5,3) - qJ(4);
t154 = sin(pkin(8)) * rSges(5,2) - rSges(5,1) * t95 - pkin(3);
t48 = -Icges(6,3) * t86 + t104 * t85;
t152 = t85 ^ 2;
t151 = t86 ^ 2;
t150 = 2 * qJD(3);
t149 = m(5) / 0.2e1;
t148 = m(6) / 0.2e1;
t131 = t85 * t93;
t75 = t86 * rSges(6,3);
t124 = rSges(6,1) * t131 - t75;
t132 = t85 * t92;
t57 = t84 * t85 + t75;
t21 = (-rSges(6,2) * t132 + t124 + t57) * t86;
t146 = -t21 / 0.2e1;
t143 = t86 / 0.2e1;
t37 = t154 * t86 + t155 * t85;
t97 = -t154 * t85 + t155 * t86;
t141 = m(5) * (-t85 * t37 - t86 * t97);
t140 = m(6) * t14;
t139 = m(6) * t109;
t133 = t85 * t86;
t129 = t86 * t93;
t33 = t86 * t96 + (t91 - t135) * t85 + t124;
t127 = t33 - t35;
t126 = t55 * t129 + t85 * t49;
t125 = -t55 * t131 + t161;
t123 = m(6) * qJD(5);
t26 = t86 * t66 + t67 * t85;
t101 = t26 * t148;
t110 = (t151 + t152) * t109;
t102 = m(6) * t110;
t17 = t101 + t102 / 0.2e1;
t115 = t17 * qJD(3);
t98 = -(-t116 * t86 - t117 * t85) * t139 / 0.2e1;
t99 = m(6) * (t66 * t116 + t67 * t117);
t18 = -t99 / 0.2e1 + t98;
t114 = t18 * qJD(2);
t113 = t21 * t148;
t103 = Icges(6,5) * t92 + Icges(6,6) * t93;
t11 = t129 * t54 - t130 * t51 + t48 * t85;
t100 = -t85 * t116 + t86 * t117;
t81 = Icges(6,1) * t92 + t121;
t61 = t103 * t85;
t60 = t103 * t86;
t39 = 0.2e1 * (t149 + t148) * t100;
t22 = t86 * (-rSges(6,1) * t129 - t85 * rSges(6,3) + t72) + t85 * t57;
t20 = qJD(3) * t113;
t19 = t99 / 0.2e1 + t98;
t16 = t101 - t102 / 0.2e1;
t12 = t52 * t132 + t125;
t10 = t52 * t130 - t126;
t8 = (t81 / 0.2e1 + t106 / 0.2e1) * t93 + t157;
t7 = t22 * t21;
t6 = t140 + t141;
t5 = -t12 * t85 + (t156 * t85 - t48 * t86) * t86;
t4 = -t10 * t85 + t11 * t86;
t3 = (t11 - t12 + t125) * t86 + t126 * t85;
t2 = (-t10 + (t134 + t48) * t86 - t126) * t86 + (-t11 + (t48 + t159) * t85 + t161) * t85;
t1 = m(6) * t7 + (t3 / 0.2e1 - t4 / 0.2e1) * t86 + (-t5 / 0.2e1 - t2 / 0.2e1) * t85;
t9 = [0, 0, qJD(5) * t113, 0, t26 * t123 + t20; 0, 0, t39 * qJD(4) + t19 * qJD(5) + (m(4) * ((rSges(4,1) * t85 + rSges(4,2) * t86) * t116 + (-rSges(4,1) * t86 + rSges(4,2) * t85) * t117) / 0.2e1 + (t116 * t97 + t37 * t117) * t149 + (t33 * t116 + t34 * t117) * t148) * t150, t39 * qJD(3), t100 * t84 * t123 + t19 * qJD(3); t123 * t146, -t18 * qJD(5), m(6) * t127 * t34 * qJD(3) + t6 * qJD(4) + t8 * qJD(5), t6 * qJD(3) + t16 * qJD(5), qJD(1) * t146 * m(6) + t8 * qJD(3) + t16 * qJD(4) - t114 + (t4 * t143 + (t14 * t84 - (-t66 * t85 + t67 * t86) * t109 - t7) * m(6) - (-t152 / 0.2e1 - t151 / 0.2e1) * t104 - (t156 + t3) * t86 / 0.2e1 + (t2 + t5 + t159) * t85 / 0.2e1) * qJD(5); 0, 0, t17 * qJD(5) + 0.4e1 * (-t140 / 0.4e1 - t141 / 0.4e1) * qJD(3) + (t33 * t86 + t14 + t28) * t148 * t150, 0, t115; t20, t18 * qJD(3), qJD(1) * t113 - t17 * qJD(4) + t1 * qJD(5) + t114 + ((-t106 - t81) * t93 / 0.2e1 + t127 * t139 * t85 - t157) * qJD(3), -t115, t1 * qJD(3) + (m(6) * (-t84 * t110 + t22 * t26) + (-t60 * t133 + t151 * t61) * t143 - t85 * (-t61 * t133 + t152 * t60) / 0.2e1) * qJD(5);];
Cq = t9;
