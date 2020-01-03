% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2,theta3]';
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
% Datum: 2019-12-31 16:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RPPR3_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR3_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR3_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPPR3_coriolismatJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR3_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPR3_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPR3_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:37:50
% EndTime: 2019-12-31 16:37:51
% DurationCPUTime: 0.66s
% Computational Cost: add. (2848->106), mult. (2224->163), div. (0->0), fcn. (1972->7), ass. (0->78)
t78 = pkin(7) + qJ(4);
t74 = sin(t78);
t79 = qJ(1) + pkin(6);
t75 = sin(t79);
t115 = t74 * t75;
t77 = cos(t79);
t107 = rSges(5,2) * t115 + t77 * rSges(5,3);
t110 = pkin(5) + qJ(3);
t119 = sin(qJ(1)) * pkin(1);
t76 = cos(t78);
t112 = t76 * rSges(5,1);
t94 = cos(pkin(7)) * pkin(3) + pkin(2) + t112;
t24 = t77 * t110 - t94 * t75 + t107 - t119;
t120 = cos(qJ(1)) * pkin(1);
t114 = t74 * t77;
t93 = -rSges(5,2) * t114 + t75 * rSges(5,3);
t25 = t75 * t110 + t94 * t77 + t120 + t93;
t57 = rSges(5,1) * t74 + rSges(5,2) * t76;
t48 = t57 * t75;
t49 = t57 * t77;
t100 = Icges(5,2) * t76;
t103 = Icges(5,4) * t74;
t53 = t100 + t103;
t104 = Icges(5,1) * t76;
t56 = -t103 + t104;
t134 = -(t56 / 0.2e1 - t53 / 0.2e1) * t74 - m(5) * (t24 * t48 - t25 * t49);
t128 = t77 ^ 2;
t129 = t75 ^ 2;
t133 = t57 * (t128 + t129);
t131 = t24 * t77 + t25 * t75;
t71 = Icges(5,4) * t76;
t54 = -Icges(5,2) * t74 + t71;
t55 = Icges(5,1) * t74 + t71;
t126 = m(5) * t131;
t124 = t75 / 0.2e1;
t122 = -t77 / 0.2e1;
t98 = rSges(4,3) + qJ(3);
t121 = m(4) * ((t98 * t77 - t119) * t77 + (t98 * t75 + t120) * t75);
t37 = -Icges(5,6) * t77 + t54 * t75;
t116 = t74 * t37;
t113 = t75 * t77;
t111 = t76 * t77;
t102 = Icges(5,5) * t76;
t99 = Icges(5,6) * t75;
t35 = -Icges(5,3) * t77 + t75 * t102 - t74 * t99;
t61 = t75 * t103;
t39 = -Icges(5,5) * t77 + t75 * t104 - t61;
t109 = -t39 * t111 - t75 * t35;
t86 = -Icges(5,6) * t74 + t102;
t36 = Icges(5,3) * t75 + t86 * t77;
t40 = Icges(5,5) * t75 + t56 * t77;
t108 = t40 * t111 + t75 * t36;
t106 = m(5) * qJD(4);
t88 = t75 * t48 + t77 * t49;
t83 = m(5) * t88 / 0.2e1;
t84 = m(5) * t133;
t15 = t83 + t84 / 0.2e1;
t97 = t15 * qJD(1);
t30 = t75 * t76 * t40;
t92 = t36 * t77 - t30;
t38 = t54 * t77 + t99;
t91 = t74 * t38 - t35;
t85 = -Icges(5,5) * t74 - Icges(5,6) * t76;
t58 = -rSges(5,2) * t74 + t112;
t43 = t85 * t77;
t42 = t85 * t75;
t14 = t83 - t84 / 0.2e1;
t13 = -t38 * t114 + t108;
t12 = -t37 * t114 - t109;
t11 = -t38 * t115 - t92;
t7 = t121 + t126;
t6 = (t55 / 0.2e1 + t54 / 0.2e1) * t76 - t134;
t5 = -t12 * t77 + t13 * t75;
t4 = -(-(-t76 * t39 + t116) * t75 - t35 * t77) * t77 + t11 * t75;
t3 = (t11 - t30 + (t36 + t116) * t77 + t109) * t77 + t108 * t75;
t2 = (t91 * t77 - t108 + t13) * t77 + (t91 * t75 + t12 + t92) * t75;
t1 = (t5 / 0.2e1 - t3 / 0.2e1) * t77 + (t2 / 0.2e1 + t4 / 0.2e1) * t75;
t8 = [qJD(3) * t7 + qJD(4) * t6, 0, qJD(1) * t7 + qJD(4) * t14, t6 * qJD(1) + t14 * qJD(3) + (-t131 * t58 + (-t48 * t77 + t49 * t75) * t57) * t106 + (((-t53 * t77 + t40) * t76 + (-t55 * t77 - t38) * t74) * t124 + t77 * t3 / 0.2e1 + (t129 / 0.2e1 + t128 / 0.2e1) * t86 - (t2 + t4) * t75 / 0.2e1 + ((-t75 * t100 + t39 - t61) * t76 + (-t55 * t75 - t37) * t74 + t5) * t122) * qJD(4); 0, 0, 0, -t88 * t106; t15 * qJD(4) + 0.4e1 * (-t121 / 0.4e1 - t126 / 0.4e1) * qJD(1), 0, 0, t97; (-(t55 + t54) * t76 / 0.2e1 + t134) * qJD(1) - t15 * qJD(3) + t1 * qJD(4), 0, -t97, t1 * qJD(1) + (m(5) * (t58 * t133 - (t75 * (t75 * t112 - t107) + t77 * (rSges(5,1) * t111 + t93)) * t88) + (-t42 * t113 + t129 * t43) * t124 + (-t43 * t113 + t128 * t42) * t122) * qJD(4);];
Cq = t8;
