% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RPPR7
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
% Datum: 2019-12-31 16:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RPPR7_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR7_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR7_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR7_coriolismatJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR7_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPR7_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPR7_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:41:36
% EndTime: 2019-12-31 16:41:37
% DurationCPUTime: 0.77s
% Computational Cost: add. (1910->127), mult. (2637->202), div. (0->0), fcn. (2230->6), ass. (0->89)
t92 = sin(qJ(1));
t88 = t92 ^ 2;
t141 = pkin(1) + qJ(3);
t87 = pkin(6) + qJ(4);
t80 = cos(t87);
t112 = Icges(5,4) * t80;
t79 = sin(t87);
t114 = Icges(5,1) * t79;
t100 = t112 + t114;
t138 = -pkin(5) - t141;
t93 = cos(qJ(1));
t81 = t93 * qJ(2);
t85 = t92 * rSges(5,3);
t104 = t79 * rSges(5,1) + t80 * rSges(5,2);
t90 = sin(pkin(6));
t94 = t90 * pkin(3) + t104;
t30 = t138 * t92 + t94 * t93 + t81 - t85;
t31 = (qJ(2) + t94) * t92 + (rSges(5,3) - t138) * t93;
t69 = rSges(5,1) * t80 - rSges(5,2) * t79;
t59 = t69 * t92;
t60 = t69 * t93;
t111 = Icges(5,2) * t79;
t65 = -t111 + t112;
t140 = (t100 / 0.2e1 + t65 / 0.2e1) * t80 - m(5) * (t30 * t60 + t31 * t59);
t29 = -t59 * t93 + t92 * t60;
t121 = m(5) * t29;
t113 = Icges(5,4) * t79;
t99 = Icges(5,2) * t80 + t113;
t48 = -Icges(5,6) * t92 + t99 * t93;
t50 = -Icges(5,5) * t92 + t100 * t93;
t137 = (t80 * t48 + t79 * t50) * t93;
t136 = -t30 * t92 + t93 * t31;
t89 = t93 ^ 2;
t77 = t88 + t89;
t134 = 0.2e1 * t77;
t133 = 4 * qJD(1);
t131 = -t77 / 0.2e1;
t129 = t92 / 0.2e1;
t127 = t93 / 0.2e1;
t126 = m(3) * ((rSges(3,3) * t93 + t81) * t93 + (rSges(3,3) + qJ(2)) * t88);
t105 = rSges(4,1) * t90 + rSges(4,2) * cos(pkin(6));
t108 = rSges(4,3) + t141;
t36 = t105 * t93 - t108 * t92 + t81;
t37 = t108 * t93 + (qJ(2) + t105) * t92;
t125 = m(4) * (-t36 * t92 + t93 * t37);
t124 = m(4) * (t36 * t93 + t92 * t37);
t123 = m(5) * t136;
t122 = m(5) * (t30 * t93 + t92 * t31);
t119 = t79 * t92;
t118 = t80 * t92;
t117 = t92 * t93;
t116 = qJD(4) / 0.2e1;
t101 = t92 * t59 + t93 * t60;
t95 = m(5) * t101;
t96 = m(5) * t69 * t131;
t16 = -t95 / 0.2e1 + t96;
t110 = t16 * qJD(1);
t27 = (m(4) / 0.2e1 + m(5) / 0.2e1) * t134;
t109 = t27 * qJD(1);
t97 = Icges(5,5) * t79 + Icges(5,6) * t80;
t45 = Icges(5,3) * t93 + t97 * t92;
t47 = Icges(5,6) * t93 + t99 * t92;
t76 = t92 * t112;
t49 = Icges(5,5) * t93 + t92 * t114 + t76;
t12 = t47 * t118 + t49 * t119 + t93 * t45;
t46 = -Icges(5,3) * t92 + t97 * t93;
t13 = -t48 * t118 - t50 * t119 - t93 * t46;
t107 = -t121 / 0.2e1;
t106 = t77 * t104;
t103 = -t80 * t47 - t79 * t49;
t67 = Icges(5,1) * t80 - t113;
t98 = Icges(5,5) * t80 - Icges(5,6) * t79;
t54 = t98 * t93;
t53 = t92 * t98;
t42 = t92 * t45;
t28 = (m(4) / 0.4e1 + m(5) / 0.4e1) * t134 + (m(4) + m(5)) * t131;
t23 = t116 * t121;
t17 = t95 / 0.2e1 + t96;
t15 = -t92 * t46 + t137;
t14 = t103 * t93 + t42;
t8 = t123 + t125;
t7 = t14 * t93 + t15 * t92;
t6 = t12 * t93 + t13 * t92;
t5 = (-t67 / 0.2e1 + t99 / 0.2e1) * t79 - t140;
t4 = t122 + t124 + t126;
t3 = t88 * t46 + (t13 - t42 + (-t103 + t46) * t93) * t93;
t2 = (-t14 + t42 + t13) * t92 + (t15 - t137 + (t103 + t46) * t92 + t12) * t93;
t1 = (t3 / 0.2e1 + t7 / 0.2e1) * t93 + (-t6 / 0.2e1 + t2 / 0.2e1) * t92;
t9 = [qJD(2) * t4 + qJD(3) * t8 + qJD(4) * t5, qJD(1) * t4 + qJD(3) * t28 + t23, qJD(1) * t8 + qJD(2) * t28 + qJD(4) * t17, t5 * qJD(1) + t17 * qJD(3) + qJD(2) * t121 / 0.2e1 + (((t67 * t92 - t47) * t80 + (t92 * t111 - t49 - t76) * t79) * t127 - t92 * t2 / 0.2e1 + (t104 * t136 + t29 * t69) * m(5) - (t89 / 0.2e1 + t88 / 0.2e1) * t97 + ((-t67 * t93 + t48) * t80 + (t65 * t93 + t50) * t79 + t6) * t129 - (t3 + t7) * t93 / 0.2e1) * qJD(4); -t27 * qJD(3) + t23 + (-t126 / 0.4e1 - t124 / 0.4e1 - t122 / 0.4e1) * t133, 0, -t109, 0.2e1 * (t29 * qJD(1) / 0.4e1 - t106 * t116) * m(5); t27 * qJD(2) - t16 * qJD(4) + (-t125 / 0.4e1 - t123 / 0.4e1) * t133, t109, 0, -t110; ((-t99 + t67) * t79 / 0.2e1 + t140) * qJD(1) + qJD(2) * t107 + t16 * qJD(3) + t1 * qJD(4), qJD(1) * t107, t110, t1 * qJD(1) + (m(5) * (-(-t93 * (t104 * t93 - t85) + (-t93 * rSges(5,3) - t104 * t92) * t92) * t101 - t69 * t106) + (-t54 * t117 + t89 * t53) * t127 + (t53 * t117 - t88 * t54) * t129) * qJD(4);];
Cq = t9;
