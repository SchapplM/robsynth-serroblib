% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:46
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Cq = S4RPPP1_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPP1_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_coriolismatJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPP1_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPP1_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPP1_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:45:25
% EndTime: 2018-11-14 13:45:26
% DurationCPUTime: 0.58s
% Computational Cost: add. (3013->100), mult. (3926->146), div. (0->0), fcn. (3592->9), ass. (0->76)
t107 = sin(pkin(6));
t92 = sin(qJ(1));
t93 = cos(qJ(1));
t100 = pkin(4) + pkin(6);
t101 = pkin(4) - pkin(6);
t94 = cos(t101) / 0.2e1 + cos(t100) / 0.2e1;
t75 = t92 * t107 - t93 * t94;
t108 = cos(pkin(6));
t81 = sin(t100) / 0.2e1 - sin(t101) / 0.2e1;
t76 = t92 * t108 + t93 * t81;
t77 = t93 * t107 + t92 * t94;
t78 = t93 * t108 - t92 * t81;
t146 = m(5) * (-t75 * t78 + t76 * t77);
t37 = -t146 / 0.2e1;
t38 = t146 / 0.2e1;
t91 = sin(pkin(4));
t114 = t91 * t93;
t115 = t91 * t92;
t141 = m(5) * (t78 * t114 + t76 * t115);
t50 = -t141 / 0.2e1;
t53 = t141 / 0.2e1;
t61 = t77 * t114 + t75 * t115;
t145 = (m(4) + m(5)) * t61;
t144 = t61 * (-m(4) / 0.2e1 - m(5) / 0.2e1);
t142 = rSges(5,1) + pkin(3);
t109 = -rSges(5,3) - qJ(4);
t139 = -pkin(2) + t109;
t137 = 2 * qJD(1);
t136 = 4 * qJD(1);
t135 = m(4) / 0.2e1;
t134 = m(5) / 0.2e1;
t119 = t93 * pkin(1);
t99 = t77 * qJ(3) + t119;
t29 = t77 * rSges(5,2) + (qJ(2) + t142) * t115 + t99 - t139 * t78;
t22 = t29 * t77;
t71 = t75 * qJ(3);
t98 = -t92 * pkin(1) + qJ(2) * t114;
t96 = t142 * t114 + t98;
t28 = -t75 * rSges(5,2) + t139 * t76 - t71 + t96;
t9 = -t28 * t75 + t22;
t133 = m(5) * t9;
t47 = t78 * rSges(3,1) - t77 * rSges(3,2) + t119 + (rSges(3,3) + qJ(2)) * t115;
t95 = -t76 * rSges(3,1) + t75 * rSges(3,2) + rSges(3,3) * t114 + t98;
t26 = t114 * t95 + t47 * t115;
t129 = m(3) * t26;
t118 = rSges(4,2) - pkin(2);
t35 = t77 * rSges(4,3) + (rSges(4,1) + qJ(2)) * t115 - t118 * t78 + t99;
t31 = t35 * t77;
t97 = rSges(4,1) * t114 + t98;
t34 = -t75 * rSges(4,3) + t118 * t76 - t71 + t97;
t12 = -t34 * t75 + t31;
t128 = m(4) * t12;
t13 = t34 * t114 + t35 * t115;
t127 = m(4) * t13;
t11 = t28 * t114 + t29 * t115;
t126 = m(5) * t11;
t23 = t29 * t78;
t10 = -t28 * t76 + t23;
t113 = -t145 / 0.2e1;
t111 = t145 / 0.2e1;
t110 = m(5) * qJD(1);
t74 = t76 * pkin(2);
t36 = t76 * rSges(4,2) - t74 + (-rSges(4,3) - qJ(3)) * t75 + t97;
t27 = -t74 + t109 * t76 + (-rSges(5,2) - qJ(3)) * t75 + t96;
t18 = 0.2e1 * t53;
t17 = t50 + t53;
t16 = 0.2e1 * t50;
t8 = t111 - t144;
t7 = t111 + t113;
t6 = t113 + t144;
t5 = 0.2e1 * t38;
t4 = t37 + t38;
t3 = 0.2e1 * t37;
t2 = t128 + t133;
t1 = t126 + t127 + t129;
t14 = [(m(4) * (-t34 + t36) * t35 / 0.4e1 + m(5) * (t27 - t28) * t29 / 0.4e1) * t136 + t1 * qJD(2) + t2 * qJD(3) + m(5) * t10 * qJD(4), t1 * qJD(1) + t7 * qJD(3) + t17 * qJD(4), t2 * qJD(1) + t7 * qJD(2) + t4 * qJD(4), t17 * qJD(2) + t4 * qJD(3) + t10 * t110; t6 * qJD(3) + t16 * qJD(4) + (-t129 / 0.4e1 - t127 / 0.4e1 - t126 / 0.4e1) * t136 + (m(3) * ((-t47 * t92 - t93 * t95) * t91 + t26) / 0.2e1 + ((-t35 * t92 - t36 * t93) * t91 + t13) * t135 + ((-t27 * t93 - t29 * t92) * t91 + t11) * t134) * t137, 0, t6 * qJD(1), t16 * qJD(1); t8 * qJD(2) + t3 * qJD(4) + (-t128 / 0.4e1 - t133 / 0.4e1) * t136 + ((t75 * t36 + t12 - t31) * t135 + (t75 * t27 - t22 + t9) * t134) * t137, t8 * qJD(1), 0, t3 * qJD(1); (t76 * t27 - t23) * t110 + t18 * qJD(2) + t5 * qJD(3), t18 * qJD(1), t5 * qJD(1), 0;];
Cq  = t14;
