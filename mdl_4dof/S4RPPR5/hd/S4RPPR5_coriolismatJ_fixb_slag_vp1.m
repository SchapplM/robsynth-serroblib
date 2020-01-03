% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
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
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RPPR5_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR5_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR5_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR5_coriolismatJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR5_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPR5_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPR5_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:39:43
% EndTime: 2019-12-31 16:39:45
% DurationCPUTime: 1.06s
% Computational Cost: add. (2238->114), mult. (4795->187), div. (0->0), fcn. (5618->6), ass. (0->74)
t92 = cos(qJ(4));
t114 = Icges(5,4) * t92;
t91 = sin(qJ(4));
t101 = -Icges(5,2) * t91 + t114;
t109 = sin(pkin(6));
t110 = cos(pkin(6));
t128 = sin(qJ(1));
t129 = cos(qJ(1));
t73 = -t128 * t109 - t129 * t110;
t74 = t129 * t109 - t128 * t110;
t46 = Icges(5,6) * t73 + t101 * t74;
t121 = t46 * t91;
t99 = Icges(5,5) * t92 - Icges(5,6) * t91;
t45 = Icges(5,3) * t74 - t99 * t73;
t148 = t74 * t45;
t48 = Icges(5,6) * t74 - t101 * t73;
t120 = t91 * t48;
t43 = Icges(5,3) * t73 + t99 * t74;
t147 = -t120 + t43;
t115 = Icges(5,4) * t91;
t100 = Icges(5,2) * t92 + t115;
t103 = Icges(5,1) * t92 - t115;
t122 = t91 * rSges(5,2);
t107 = t74 * rSges(5,3) + t73 * t122;
t127 = rSges(5,1) * t92;
t97 = t129 * pkin(1) + t128 * qJ(2);
t95 = t129 * pkin(2) + t97;
t27 = t74 * pkin(5) + (-pkin(3) - t127) * t73 + t95 + t107;
t104 = t91 * rSges(5,1) + rSges(5,2) * t92;
t65 = t104 * t74;
t66 = t104 * t73;
t105 = -t122 + t127;
t106 = -t128 * pkin(1) + t129 * qJ(2);
t96 = -t128 * pkin(2) + t106;
t93 = (rSges(5,3) + pkin(5)) * t73 + (pkin(3) + t105) * t74 + t96;
t146 = (t103 / 0.2e1 - t100 / 0.2e1) * t91 + m(5) * (t27 * t66 - t65 * t93);
t49 = Icges(5,5) * t73 + t103 * t74;
t144 = t73 ^ 2;
t143 = t74 ^ 2;
t139 = -t73 / 0.2e1;
t137 = t74 / 0.2e1;
t135 = m(3) * ((t129 * rSges(3,3) + t106) * t129 + (t128 * rSges(3,3) + t97) * t128);
t134 = m(4) * ((-t73 * rSges(4,1) - t74 * rSges(4,2) + t95) * t128 + (t74 * rSges(4,1) - t73 * rSges(4,2) + t96) * t129);
t133 = m(5) * (t27 * t128 + t93 * t129);
t132 = m(5) * (-t128 * t65 - t129 * t66);
t130 = m(5) * (t128 * t74 + t129 * t73) * t104;
t51 = Icges(5,5) * t74 - t103 * t73;
t126 = t51 * t92;
t124 = t73 * t92;
t123 = t74 * t73;
t118 = -t74 * t126 - t73 * t45;
t117 = -t51 * t124 + t148;
t116 = m(5) * qJD(4);
t102 = Icges(5,1) * t91 + t114;
t98 = Icges(5,5) * t91 + Icges(5,6) * t92;
t12 = -t73 * t121 + t124 * t49 - t74 * t43;
t60 = t98 * t73;
t59 = t98 * t74;
t40 = t130 / 0.2e1;
t29 = t132 / 0.2e1;
t22 = t65 * t74 + t66 * t73;
t16 = t29 + t40;
t15 = t29 - t130 / 0.2e1;
t14 = t40 - t132 / 0.2e1;
t13 = t73 * t120 + t117;
t11 = t74 * t120 + t118;
t8 = (t102 / 0.2e1 + t101 / 0.2e1) * t92 + t146;
t6 = t133 + t134 + t135;
t5 = -t12 * t73 + t13 * t74;
t4 = -(-(-t49 * t92 + t121) * t74 + t73 * t43) * t73 + t11 * t74;
t3 = (t11 - t12 - t118) * t73 + t117 * t74;
t2 = (t12 + (t126 + t147) * t74) * t74 + (t147 * t73 - t117 + t13 + t148) * t73;
t1 = (t2 / 0.2e1 + t4 / 0.2e1) * t74 + (t5 / 0.2e1 - t3 / 0.2e1) * t73;
t7 = [t6 * qJD(2) + t8 * qJD(4), qJD(1) * t6 + qJD(4) * t16, 0, t8 * qJD(1) + t16 * qJD(2) + (((-t100 * t73 - t51) * t92 + (-t102 * t73 + t48) * t91) * t137 + t73 * t3 / 0.2e1 + ((t104 * t66 + t105 * t27) * t74 + (-t104 * t65 + t105 * t93) * t73) * m(5) - (t143 / 0.2e1 + t144 / 0.2e1) * t99 + ((-t100 * t74 + t49) * t92 + (-t102 * t74 - t46) * t91 + t5) * t139 - (t2 + t4) * t74 / 0.2e1) * qJD(4); t15 * qJD(4) + 0.4e1 * (-t135 / 0.4e1 - t134 / 0.4e1 - t133 / 0.4e1) * qJD(1), 0, 0, t15 * qJD(1) - (-t128 * t73 + t129 * t74) * t105 * t116; 0, 0, 0, -t22 * t116; t14 * qJD(2) + t1 * qJD(4) + ((-t101 - t102) * t92 / 0.2e1 - t146) * qJD(1), t14 * qJD(1), 0, t1 * qJD(1) + (m(5) * ((t74 * (-t73 * rSges(5,3) - t105 * t74) + t73 * (-rSges(5,1) * t124 + t107)) * t22 + (t143 + t144) * t105 * t104) + (-t59 * t123 + t143 * t60) * t137 + (-t60 * t123 + t144 * t59) * t139) * qJD(4);];
Cq = t7;
