% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
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
% Datum: 2019-12-31 17:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPPPR4_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR4_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR4_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR4_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR4_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR4_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPPR4_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:45:07
% EndTime: 2019-12-31 17:45:09
% DurationCPUTime: 0.83s
% Computational Cost: add. (3500->133), mult. (2834->207), div. (0->0), fcn. (2386->8), ass. (0->91)
t150 = pkin(2) + qJ(4);
t95 = pkin(8) + qJ(5);
t93 = cos(t95);
t120 = Icges(6,4) * t93;
t91 = sin(t95);
t122 = Icges(6,1) * t91;
t107 = t120 + t122;
t111 = t91 * rSges(6,1) + t93 * rSges(6,2);
t97 = sin(pkin(8));
t101 = t97 * pkin(4) + t111;
t96 = qJ(1) + pkin(7);
t94 = cos(t96);
t114 = -sin(qJ(1)) * pkin(1) + t94 * qJ(3);
t147 = -pkin(6) - t150;
t92 = sin(t96);
t87 = t92 * rSges(6,3);
t27 = t101 * t94 + t147 * t92 + t114 - t87;
t129 = cos(qJ(1)) * pkin(1);
t28 = t129 + (qJ(3) + t101) * t92 + (rSges(6,3) - t147) * t94;
t77 = rSges(6,1) * t93 - rSges(6,2) * t91;
t63 = t77 * t92;
t64 = t77 * t94;
t119 = Icges(6,2) * t91;
t73 = -t119 + t120;
t149 = (t107 / 0.2e1 + t73 / 0.2e1) * t93 - m(6) * (t27 * t64 + t28 * t63);
t121 = Icges(6,4) * t91;
t106 = Icges(6,2) * t93 + t121;
t53 = -Icges(6,6) * t92 + t106 * t94;
t55 = -Icges(6,5) * t92 + t107 * t94;
t146 = (t53 * t93 + t55 * t91) * t94;
t145 = -t27 * t92 + t94 * t28;
t89 = t92 ^ 2;
t90 = t94 ^ 2;
t78 = t89 + t90;
t143 = 0.2e1 * t78;
t142 = 4 * qJD(1);
t33 = -t63 * t94 + t64 * t92;
t140 = t33 / 0.2e1;
t139 = -t78 / 0.2e1;
t137 = t92 / 0.2e1;
t135 = t94 / 0.2e1;
t134 = m(4) * ((rSges(4,3) * t94 + t114) * t94 + (t129 + (rSges(4,3) + qJ(3)) * t92) * t92);
t112 = rSges(5,1) * t97 + rSges(5,2) * cos(pkin(8));
t116 = rSges(5,3) + t150;
t37 = t112 * t94 - t116 * t92 + t114;
t38 = t129 + t116 * t94 + (qJ(3) + t112) * t92;
t133 = m(5) * (-t37 * t92 + t94 * t38);
t132 = m(5) * (t37 * t94 + t92 * t38);
t131 = m(6) * t145;
t130 = m(6) * (t27 * t94 + t92 * t28);
t127 = t91 * t92;
t126 = t92 * t93;
t125 = t92 * t94;
t123 = m(6) * qJD(5);
t108 = t92 * t63 + t94 * t64;
t102 = m(6) * t108;
t103 = m(6) * t77 * t139;
t18 = -t102 / 0.2e1 + t103;
t118 = t18 * qJD(1);
t30 = (m(5) / 0.2e1 + m(6) / 0.2e1) * t143;
t117 = t30 * qJD(1);
t104 = Icges(6,5) * t91 + Icges(6,6) * t93;
t50 = Icges(6,3) * t94 + t104 * t92;
t52 = Icges(6,6) * t94 + t106 * t92;
t79 = t92 * t120;
t54 = Icges(6,5) * t94 + t122 * t92 + t79;
t12 = t52 * t126 + t54 * t127 + t94 * t50;
t51 = -Icges(6,3) * t92 + t104 * t94;
t13 = -t53 * t126 - t55 * t127 - t94 * t51;
t115 = -m(6) * t33 / 0.2e1;
t113 = t78 * t111;
t110 = -t52 * t93 - t54 * t91;
t75 = Icges(6,1) * t93 - t121;
t105 = Icges(6,5) * t93 - Icges(6,6) * t91;
t58 = t105 * t94;
t57 = t92 * t105;
t45 = t92 * t50;
t31 = (m(5) / 0.4e1 + m(6) / 0.4e1) * t143 + (m(5) + m(6)) * t139;
t29 = t123 * t140;
t19 = t102 / 0.2e1 + t103;
t15 = -t51 * t92 + t146;
t14 = t110 * t94 + t45;
t8 = t14 * t94 + t15 * t92;
t7 = t12 * t94 + t13 * t92;
t6 = t131 + t133;
t5 = (-t75 / 0.2e1 + t106 / 0.2e1) * t91 - t149;
t4 = t130 + t132 + t134;
t3 = t51 * t89 + (t13 - t45 + (-t110 + t51) * t94) * t94;
t2 = (-t14 + t45 + t13) * t92 + (t15 - t146 + (t110 + t51) * t92 + t12) * t94;
t1 = (t3 / 0.2e1 + t8 / 0.2e1) * t94 + (-t7 / 0.2e1 + t2 / 0.2e1) * t92;
t9 = [qJD(3) * t4 + qJD(4) * t6 + qJD(5) * t5, 0, qJD(1) * t4 + qJD(4) * t31 + t29, qJD(1) * t6 + qJD(3) * t31 + qJD(5) * t19, qJD(3) * t140 * m(6) + t5 * qJD(1) + t19 * qJD(4) + (((t75 * t92 - t52) * t93 + (t119 * t92 - t54 - t79) * t91) * t135 - t92 * t2 / 0.2e1 + (t111 * t145 + t33 * t77) * m(6) - (t90 / 0.2e1 + t89 / 0.2e1) * t104 + ((-t75 * t94 + t53) * t93 + (t73 * t94 + t55) * t91 + t7) * t137 - (t3 + t8) * t94 / 0.2e1) * qJD(5); 0, 0, 0, 0, -t123 * t108; -t30 * qJD(4) + t29 + (-t134 / 0.4e1 - t132 / 0.4e1 - t130 / 0.4e1) * t142, 0, 0, -t117, 0.2e1 * (t33 * qJD(1) / 0.4e1 - qJD(5) * t113 / 0.2e1) * m(6); t30 * qJD(3) - t18 * qJD(5) + (-t133 / 0.4e1 - t131 / 0.4e1) * t142, 0, t117, 0, -t118; ((-t106 + t75) * t91 / 0.2e1 + t149) * qJD(1) + t18 * qJD(4) + t1 * qJD(5) + qJD(3) * t115, 0, qJD(1) * t115, t118, t1 * qJD(1) + (m(6) * (-t113 * t77 - (-t94 * (t111 * t94 - t87) + (-t94 * rSges(6,3) - t111 * t92) * t92) * t108) + (-t125 * t58 + t57 * t90) * t135 + (t125 * t57 - t58 * t89) * t137) * qJD(5);];
Cq = t9;
