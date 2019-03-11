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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

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
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPP1_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPP1_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPP1_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:26:20
% EndTime: 2019-03-08 18:26:21
% DurationCPUTime: 0.51s
% Computational Cost: add. (1285->96), mult. (3350->144), div. (0->0), fcn. (3592->6), ass. (0->75)
t105 = sin(pkin(6));
t91 = sin(qJ(1));
t92 = cos(qJ(1));
t106 = cos(pkin(6));
t107 = cos(pkin(4));
t95 = t107 * t106;
t75 = t105 * t91 - t92 * t95;
t94 = t107 * t105;
t76 = t106 * t91 + t92 * t94;
t77 = t105 * t92 + t91 * t95;
t78 = t106 * t92 - t91 * t94;
t145 = m(5) * (-t75 * t78 + t76 * t77);
t37 = -t145 / 0.2e1;
t38 = t145 / 0.2e1;
t90 = sin(pkin(4));
t113 = t90 * t92;
t114 = t90 * t91;
t140 = m(5) * (t78 * t113 + t76 * t114);
t50 = -t140 / 0.2e1;
t53 = t140 / 0.2e1;
t61 = t77 * t113 + t75 * t114;
t144 = (m(4) + m(5)) * t61;
t143 = t61 * (-m(4) / 0.2e1 - m(5) / 0.2e1);
t141 = rSges(5,1) + pkin(3);
t108 = -rSges(5,3) - qJ(4);
t138 = -pkin(2) + t108;
t136 = 2 * qJD(1);
t135 = 4 * qJD(1);
t134 = m(4) / 0.2e1;
t133 = m(5) / 0.2e1;
t118 = t92 * pkin(1);
t99 = t77 * qJ(3) + t118;
t29 = t77 * rSges(5,2) + (qJ(2) + t141) * t114 + t99 - t138 * t78;
t21 = t29 * t77;
t71 = t75 * qJ(3);
t98 = -t91 * pkin(1) + qJ(2) * t113;
t96 = t141 * t113 + t98;
t28 = -t75 * rSges(5,2) + t138 * t76 - t71 + t96;
t9 = -t28 * t75 + t21;
t132 = m(5) * t9;
t47 = t78 * rSges(3,1) - t77 * rSges(3,2) + t118 + (rSges(3,3) + qJ(2)) * t114;
t93 = -t76 * rSges(3,1) + t75 * rSges(3,2) + rSges(3,3) * t113 + t98;
t24 = t113 * t93 + t47 * t114;
t128 = m(3) * t24;
t117 = rSges(4,2) - pkin(2);
t35 = t77 * rSges(4,3) + (rSges(4,1) + qJ(2)) * t114 - t117 * t78 + t99;
t31 = t35 * t77;
t97 = rSges(4,1) * t113 + t98;
t34 = -t75 * rSges(4,3) + t117 * t76 - t71 + t97;
t12 = -t34 * t75 + t31;
t127 = m(4) * t12;
t13 = t34 * t113 + t35 * t114;
t126 = m(4) * t13;
t11 = t28 * t113 + t29 * t114;
t125 = m(5) * t11;
t22 = t29 * t78;
t10 = -t28 * t76 + t22;
t112 = -t144 / 0.2e1;
t110 = t144 / 0.2e1;
t109 = m(5) * qJD(1);
t74 = t76 * pkin(2);
t36 = t76 * rSges(4,2) - t74 + (-rSges(4,3) - qJ(3)) * t75 + t97;
t27 = -t74 + t108 * t76 + (-rSges(5,2) - qJ(3)) * t75 + t96;
t16 = 0.2e1 * t53;
t15 = t50 + t53;
t14 = 0.2e1 * t50;
t8 = 0.2e1 * t38;
t7 = t37 + t38;
t6 = 0.2e1 * t37;
t5 = t110 - t143;
t4 = t110 + t112;
t3 = t112 + t143;
t2 = t127 + t132;
t1 = t125 + t126 + t128;
t17 = [(m(4) * (-t34 + t36) * t35 / 0.4e1 + m(5) * (t27 - t28) * t29 / 0.4e1) * t135 + t1 * qJD(2) + t2 * qJD(3) + m(5) * t10 * qJD(4), t1 * qJD(1) + t4 * qJD(3) + t15 * qJD(4), t2 * qJD(1) + t4 * qJD(2) + t7 * qJD(4), t15 * qJD(2) + t7 * qJD(3) + t10 * t109; t3 * qJD(3) + t14 * qJD(4) + (-t128 / 0.4e1 - t126 / 0.4e1 - t125 / 0.4e1) * t135 + (m(3) * ((-t47 * t91 - t92 * t93) * t90 + t24) / 0.2e1 + ((-t35 * t91 - t36 * t92) * t90 + t13) * t134 + ((-t27 * t92 - t29 * t91) * t90 + t11) * t133) * t136, 0, t3 * qJD(1), t14 * qJD(1); t5 * qJD(2) + t6 * qJD(4) + (-t127 / 0.4e1 - t132 / 0.4e1) * t135 + ((t75 * t36 + t12 - t31) * t134 + (t75 * t27 - t21 + t9) * t133) * t136, t5 * qJD(1), 0, t6 * qJD(1); (t76 * t27 - t22) * t109 + t16 * qJD(2) + t8 * qJD(3), t16 * qJD(1), t8 * qJD(1), 0;];
Cq  = t17;
