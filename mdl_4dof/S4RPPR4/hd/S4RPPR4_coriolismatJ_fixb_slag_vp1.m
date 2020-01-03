% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RPPR4
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
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RPPR4_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR4_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR4_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR4_coriolismatJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR4_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPR4_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPR4_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:38:49
% EndTime: 2019-12-31 16:38:50
% DurationCPUTime: 0.64s
% Computational Cost: add. (2050->97), mult. (2129->159), div. (0->0), fcn. (1906->6), ass. (0->65)
t72 = qJ(1) + pkin(6);
t70 = sin(t72);
t71 = cos(t72);
t75 = cos(qJ(4));
t73 = sin(qJ(4));
t90 = Icges(5,4) * t73;
t80 = Icges(5,2) * t75 + t90;
t38 = -Icges(5,6) * t70 + t80 * t71;
t89 = Icges(5,4) * t75;
t91 = Icges(5,1) * t73;
t81 = t89 + t91;
t40 = -Icges(5,5) * t70 + t81 * t71;
t113 = t75 * t38 + t73 * t40;
t37 = Icges(5,6) * t71 + t80 * t70;
t61 = t70 * t89;
t39 = Icges(5,5) * t71 + t70 * t91 + t61;
t112 = t75 * t37 + t73 * t39;
t100 = -pkin(2) - pkin(5);
t84 = rSges(5,1) * t73 + rSges(5,2) * t75;
t77 = -t70 * rSges(5,3) + t84 * t71;
t86 = -sin(qJ(1)) * pkin(1) + t71 * qJ(3);
t23 = t100 * t70 + t77 + t86;
t98 = cos(qJ(1)) * pkin(1);
t24 = t98 + (rSges(5,3) - t100) * t71 + (qJ(3) + t84) * t70;
t59 = rSges(5,1) * t75 - rSges(5,2) * t73;
t48 = t59 * t70;
t49 = t59 * t71;
t88 = Icges(5,2) * t73;
t54 = -t88 + t89;
t111 = (t81 / 0.2e1 + t54 / 0.2e1) * t75 - m(5) * (t23 * t49 + t24 * t48);
t109 = t113 * t71;
t68 = t70 ^ 2;
t69 = t71 ^ 2;
t106 = m(5) * (t23 * t71 + t24 * t70);
t21 = -t48 * t71 + t49 * t70;
t105 = t21 / 0.2e1;
t104 = t70 / 0.2e1;
t102 = t71 / 0.2e1;
t99 = m(4) * ((rSges(4,3) * t71 + t86) * t71 + (t98 + (rSges(4,3) + qJ(3)) * t70) * t70);
t97 = t70 * t71;
t92 = m(5) * qJD(4);
t78 = Icges(5,5) * t73 + Icges(5,6) * t75;
t35 = Icges(5,3) * t71 + t78 * t70;
t10 = t112 * t70 + t71 * t35;
t36 = -Icges(5,3) * t70 + t78 * t71;
t11 = -t113 * t70 - t71 * t36;
t87 = -m(5) * t21 / 0.2e1;
t85 = (t68 + t69) * t84;
t56 = Icges(5,1) * t75 - t90;
t79 = Icges(5,5) * t75 - Icges(5,6) * t73;
t43 = t79 * t71;
t42 = t70 * t79;
t30 = t70 * t35;
t22 = -t48 * t70 - t49 * t71;
t18 = t92 * t105;
t13 = -t70 * t36 + t109;
t12 = -t112 * t71 + t30;
t7 = t99 + t106;
t6 = (-t56 / 0.2e1 + t80 / 0.2e1) * t73 - t111;
t5 = t12 * t71 + t13 * t70;
t4 = t10 * t71 + t11 * t70;
t3 = t68 * t36 + (t11 - t30 + (t36 + t112) * t71) * t71;
t2 = (-t12 + t30 + t11) * t70 + (t13 - t109 + (t36 - t112) * t70 + t10) * t71;
t1 = (t3 / 0.2e1 + t5 / 0.2e1) * t71 + (-t4 / 0.2e1 + t2 / 0.2e1) * t70;
t8 = [qJD(3) * t7 + qJD(4) * t6, 0, qJD(1) * t7 + t18, qJD(3) * t105 * m(5) + t6 * qJD(1) + (((t56 * t70 - t37) * t75 + (t70 * t88 - t39 - t61) * t73) * t102 - t70 * t2 / 0.2e1 + (t21 * t59 - (t23 * t70 - t24 * t71) * t84) * m(5) - (t69 / 0.2e1 + t68 / 0.2e1) * t78 + ((-t56 * t71 + t38) * t75 + (t54 * t71 + t40) * t73 + t4) * t104 - (t3 + t5) * t71 / 0.2e1) * qJD(4); 0, 0, 0, t22 * t92; t18 + 0.4e1 * (-t99 / 0.4e1 - t106 / 0.4e1) * qJD(1), 0, 0, 0.2e1 * (t21 * qJD(1) / 0.4e1 - qJD(4) * t85 / 0.2e1) * m(5); ((-t80 + t56) * t73 / 0.2e1 + t111) * qJD(1) + t1 * qJD(4) + qJD(3) * t87, 0, qJD(1) * t87, t1 * qJD(1) + (m(5) * ((-t71 * t77 + (-t71 * rSges(5,3) - t84 * t70) * t70) * t22 - t59 * t85) + (t69 * t42 - t43 * t97) * t102 + (t42 * t97 - t68 * t43) * t104) * qJD(4);];
Cq = t8;
