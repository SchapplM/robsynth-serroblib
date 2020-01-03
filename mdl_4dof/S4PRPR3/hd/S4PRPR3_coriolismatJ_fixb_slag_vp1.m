% Calculate matrix of centrifugal and coriolis load on the joints for
% S4PRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-31 16:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4PRPR3_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR3_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR3_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR3_coriolismatJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR3_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPR3_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPR3_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:20:48
% EndTime: 2019-12-31 16:20:50
% DurationCPUTime: 0.74s
% Computational Cost: add. (2804->101), mult. (2172->158), div. (0->0), fcn. (1926->5), ass. (0->76)
t77 = pkin(6) + qJ(2);
t75 = cos(t77);
t122 = t75 ^ 2;
t73 = sin(t77);
t123 = t73 ^ 2;
t129 = t122 + t123;
t76 = pkin(7) + qJ(4);
t72 = sin(t76);
t112 = t72 * t73;
t103 = rSges(5,2) * t112 + t75 * rSges(5,3);
t106 = pkin(5) + qJ(3);
t74 = cos(t76);
t109 = t74 * rSges(5,1);
t90 = cos(pkin(7)) * pkin(3) + pkin(2) + t109;
t24 = t75 * t106 - t90 * t73 + t103;
t111 = t72 * t75;
t89 = -rSges(5,2) * t111 + rSges(5,3) * t73;
t25 = t73 * t106 + t90 * t75 + t89;
t55 = rSges(5,1) * t72 + rSges(5,2) * t74;
t48 = t55 * t73;
t49 = t55 * t75;
t96 = Icges(5,2) * t74;
t99 = Icges(5,4) * t72;
t51 = t96 + t99;
t100 = Icges(5,1) * t74;
t54 = -t99 + t100;
t128 = -(t54 / 0.2e1 - t51 / 0.2e1) * t72 - m(5) * (t24 * t48 - t25 * t49);
t127 = t55 * t129;
t69 = Icges(5,4) * t74;
t52 = -Icges(5,2) * t72 + t69;
t125 = t24 * t75 + t25 * t73;
t53 = Icges(5,1) * t72 + t69;
t119 = t73 / 0.2e1;
t117 = -t75 / 0.2e1;
t116 = m(4) * t129 * (rSges(4,3) + qJ(3));
t115 = m(5) * t125;
t37 = -Icges(5,6) * t75 + t52 * t73;
t113 = t72 * t37;
t110 = t73 * t75;
t108 = t74 * t75;
t95 = Icges(5,6) * t73;
t98 = Icges(5,5) * t74;
t35 = -Icges(5,3) * t75 - t72 * t95 + t73 * t98;
t61 = t73 * t99;
t39 = -Icges(5,5) * t75 + t73 * t100 - t61;
t105 = -t39 * t108 - t73 * t35;
t82 = -Icges(5,6) * t72 + t98;
t36 = Icges(5,3) * t73 + t82 * t75;
t40 = Icges(5,5) * t73 + t54 * t75;
t104 = t40 * t108 + t73 * t36;
t102 = m(5) * qJD(4);
t84 = t73 * t48 + t75 * t49;
t79 = m(5) * t84 / 0.2e1;
t80 = m(5) * t127;
t15 = t79 + t80 / 0.2e1;
t93 = t15 * qJD(2);
t28 = t73 * t74 * t40;
t88 = t75 * t36 - t28;
t38 = t52 * t75 + t95;
t87 = t72 * t38 - t35;
t81 = -Icges(5,5) * t72 - Icges(5,6) * t74;
t57 = -rSges(5,2) * t72 + t109;
t43 = t81 * t75;
t42 = t81 * t73;
t14 = t79 - t80 / 0.2e1;
t12 = -t38 * t111 + t104;
t11 = -t37 * t111 - t105;
t10 = -t38 * t112 - t88;
t7 = t115 + t116;
t6 = (t53 / 0.2e1 + t52 / 0.2e1) * t74 - t128;
t5 = -t11 * t75 + t12 * t73;
t4 = t10 * t73 - t75 * (-(-t74 * t39 + t113) * t73 - t75 * t35);
t3 = (t10 - t28 + (t36 + t113) * t75 + t105) * t75 + t104 * t73;
t2 = (t87 * t75 - t104 + t12) * t75 + (t87 * t73 + t11 + t88) * t73;
t1 = (t5 / 0.2e1 - t3 / 0.2e1) * t75 + (t2 / 0.2e1 + t4 / 0.2e1) * t73;
t8 = [0, 0, 0, -t84 * t102; 0, qJD(3) * t7 + qJD(4) * t6, qJD(2) * t7 + qJD(4) * t14, t6 * qJD(2) + t14 * qJD(3) + (-t125 * t57 + (-t48 * t75 + t49 * t73) * t55) * t102 + (((-t51 * t75 + t40) * t74 + (-t53 * t75 - t38) * t72) * t119 + t75 * t3 / 0.2e1 + (t123 / 0.2e1 + t122 / 0.2e1) * t82 - (t2 + t4) * t73 / 0.2e1 + ((-t73 * t96 + t39 - t61) * t74 + (-t53 * t73 - t37) * t72 + t5) * t117) * qJD(4); 0, t15 * qJD(4) + 0.4e1 * (-t116 / 0.4e1 - t115 / 0.4e1) * qJD(2), 0, t93; 0, (-(t53 + t52) * t74 / 0.2e1 + t128) * qJD(2) - t15 * qJD(3) + t1 * qJD(4), -t93, t1 * qJD(2) + (m(5) * (t57 * t127 - (t73 * (t73 * t109 - t103) + t75 * (rSges(5,1) * t108 + t89)) * t84) + (-t42 * t110 + t123 * t43) * t119 + (-t43 * t110 + t122 * t42) * t117) * qJD(4);];
Cq = t8;
