% Calculate matrix of centrifugal and coriolis load on the joints for
% S4PPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4PPRR3_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR3_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR3_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR3_coriolismatJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR3_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRR3_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPRR3_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:17:22
% EndTime: 2019-12-31 16:17:24
% DurationCPUTime: 0.97s
% Computational Cost: add. (1918->107), mult. (4257->184), div. (0->0), fcn. (5094->6), ass. (0->74)
t104 = sin(qJ(3));
t105 = cos(qJ(3));
t85 = sin(pkin(6));
t86 = cos(pkin(6));
t56 = -t104 * t85 - t105 * t86;
t57 = t104 * t86 - t105 * t85;
t72 = cos(qJ(4));
t71 = sin(qJ(4));
t91 = Icges(5,4) * t71;
t80 = Icges(5,1) * t72 - t91;
t38 = -Icges(5,5) * t57 + t56 * t80;
t90 = Icges(5,4) * t72;
t78 = -Icges(5,2) * t71 + t90;
t97 = (-Icges(5,6) * t57 + t56 * t78) * t71;
t121 = t38 * t72 - t97;
t76 = Icges(5,5) * t72 - Icges(5,6) * t71;
t33 = -Icges(5,3) * t56 - t57 * t76;
t122 = t57 * t33;
t39 = -Icges(5,5) * t56 - t57 * t80;
t98 = t71 * (-Icges(5,6) * t56 - t57 * t78);
t120 = -t39 * t72 + t98;
t103 = rSges(5,1) * t72;
t99 = t71 * rSges(5,2);
t82 = -t56 * rSges(5,3) + t57 * t99;
t22 = -t56 * pkin(5) + (-pkin(3) - t103) * t57 + t82;
t67 = t99 - t103;
t23 = (-rSges(5,3) - pkin(5)) * t57 + (pkin(3) - t67) * t56;
t81 = t71 * rSges(5,1) + rSges(5,2) * t72;
t48 = t81 * t57;
t49 = t81 * t56;
t117 = (t80 / 0.2e1 - Icges(5,2) * t72 / 0.2e1 - t91 / 0.2e1) * t71 + m(5) * (t22 * t48 - t23 * t49);
t32 = -Icges(5,3) * t57 + t56 * t76;
t115 = t56 ^ 2;
t114 = t57 ^ 2;
t113 = m(5) / 0.2e1;
t54 = t57 * rSges(5,3);
t41 = t56 * t67 + t54;
t101 = t56 * t72;
t93 = rSges(5,1) * t101 - t54;
t16 = (-t56 * t99 + t41 + t93) * t57;
t111 = -t16 / 0.2e1;
t108 = t57 / 0.2e1;
t106 = m(5) * t81;
t102 = t56 * t57;
t100 = t57 * t72;
t21 = -t57 * pkin(5) + (pkin(3) - t99) * t56 + t93;
t96 = t21 - t23;
t95 = t39 * t100 + t56 * t33;
t94 = -t39 * t101 + t122;
t92 = m(5) * qJD(4);
t73 = -(-t56 * t86 - t57 * t85) * t106 / 0.2e1;
t74 = m(5) * (t48 * t85 + t49 * t86);
t13 = -t74 / 0.2e1 + t73;
t84 = t13 * qJD(2);
t83 = t16 * t113;
t75 = Icges(5,5) * t71 + Icges(5,6) * t72;
t10 = t100 * t38 + t56 * t32 - t57 * t97;
t64 = Icges(5,1) * t71 + t90;
t43 = t75 * t56;
t42 = t75 * t57;
t19 = t48 * t57 + t49 * t56;
t17 = t57 * (-rSges(5,1) * t100 + t82) + t56 * t41;
t15 = qJD(3) * t83;
t14 = t74 / 0.2e1 + t73;
t11 = t56 * t98 + t94;
t9 = t57 * t98 - t95;
t7 = (t64 / 0.2e1 + t78 / 0.2e1) * t72 + t117;
t6 = t17 * t16;
t5 = -t11 * t56 + (t121 * t56 - t57 * t32) * t57;
t4 = t10 * t57 - t56 * t9;
t3 = (t10 - t11 + t94) * t57 + t95 * t56;
t2 = (-t9 + (t32 + t98) * t57 - t95) * t57 + (-t10 + (t32 + t120) * t56 + t122) * t56;
t1 = m(5) * t6 + (t3 / 0.2e1 - t4 / 0.2e1) * t57 + (-t5 / 0.2e1 - t2 / 0.2e1) * t56;
t8 = [0, 0, qJD(4) * t83, t19 * t92 + t15; 0, 0, t14 * qJD(4) + 0.2e1 * (m(4) * ((rSges(4,1) * t56 + rSges(4,2) * t57) * t85 + (-rSges(4,1) * t57 + rSges(4,2) * t56) * t86) / 0.2e1 + (t21 * t85 + t22 * t86) * t113) * qJD(3), t14 * qJD(3) + (-t56 * t85 + t57 * t86) * t67 * t92; t92 * t111, -t13 * qJD(4), m(5) * qJD(3) * t22 * t96 + t7 * qJD(4), qJD(1) * t111 * m(5) + t7 * qJD(3) - t84 + (t4 * t108 + ((-t23 * t67 - t49 * t81) * t57 + (-t22 * t67 + t48 * t81) * t56 - t6) * m(5) - (-t114 / 0.2e1 - t115 / 0.2e1) * t76 - (t3 + t121) * t57 / 0.2e1 + (t2 + t5 + t120) * t56 / 0.2e1) * qJD(4); t15, t13 * qJD(3), qJD(1) * t83 + t1 * qJD(4) + t84 + ((-t78 - t64) * t72 / 0.2e1 + t96 * t106 * t56 - t117) * qJD(3), t1 * qJD(3) + (m(5) * (t17 * t19 - (t114 + t115) * t67 * t81) + (-t102 * t42 + t114 * t43) * t108 - t56 * (-t102 * t43 + t115 * t42) / 0.2e1) * qJD(4);];
Cq = t8;
