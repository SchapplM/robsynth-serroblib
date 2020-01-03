% Calculate matrix of centrifugal and coriolis load on the joints for
% S4PRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
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
% Datum: 2019-12-31 16:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4PRPR4_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR4_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR4_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR4_coriolismatJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR4_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPR4_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPR4_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:21:56
% EndTime: 2019-12-31 16:21:57
% DurationCPUTime: 0.64s
% Computational Cost: add. (2006->94), mult. (2077->156), div. (0->0), fcn. (1860->4), ass. (0->64)
t70 = pkin(6) + qJ(2);
t68 = sin(t70);
t66 = t68 ^ 2;
t69 = cos(t70);
t72 = cos(qJ(4));
t71 = sin(qJ(4));
t85 = Icges(5,4) * t71;
t76 = Icges(5,2) * t72 + t85;
t38 = -Icges(5,6) * t68 + t76 * t69;
t84 = Icges(5,4) * t72;
t86 = Icges(5,1) * t71;
t77 = t84 + t86;
t40 = -Icges(5,5) * t68 + t77 * t69;
t107 = t72 * t38 + t71 * t40;
t37 = Icges(5,6) * t69 + t76 * t68;
t59 = t68 * t84;
t39 = Icges(5,5) * t69 + t68 * t86 + t59;
t106 = t72 * t37 + t71 * t39;
t60 = t69 * qJ(3);
t80 = t71 * rSges(5,1) + rSges(5,2) * t72;
t73 = -t68 * rSges(5,3) + t80 * t69;
t94 = -pkin(2) - pkin(5);
t23 = t94 * t68 + t60 + t73;
t24 = (rSges(5,3) - t94) * t69 + (qJ(3) + t80) * t68;
t58 = rSges(5,1) * t72 - t71 * rSges(5,2);
t48 = t58 * t68;
t49 = t58 * t69;
t83 = Icges(5,2) * t71;
t54 = -t83 + t84;
t105 = (t77 / 0.2e1 + t54 / 0.2e1) * t72 - m(5) * (t23 * t49 + t24 * t48);
t103 = t107 * t69;
t67 = t69 ^ 2;
t100 = m(5) * (t23 * t69 + t24 * t68);
t21 = -t48 * t69 + t49 * t68;
t99 = t21 / 0.2e1;
t98 = t68 / 0.2e1;
t96 = t69 / 0.2e1;
t93 = m(4) * ((rSges(4,3) * t69 + t60) * t69 + (rSges(4,3) + qJ(3)) * t66);
t92 = t68 * t69;
t87 = m(5) * qJD(4);
t74 = Icges(5,5) * t71 + Icges(5,6) * t72;
t35 = Icges(5,3) * t69 + t74 * t68;
t10 = t106 * t68 + t69 * t35;
t36 = -Icges(5,3) * t68 + t74 * t69;
t11 = -t107 * t68 - t69 * t36;
t82 = -m(5) * t21 / 0.2e1;
t81 = (t66 + t67) * t80;
t56 = Icges(5,1) * t72 - t85;
t75 = Icges(5,5) * t72 - Icges(5,6) * t71;
t43 = t75 * t69;
t42 = t68 * t75;
t30 = t68 * t35;
t22 = -t48 * t68 - t49 * t69;
t18 = t87 * t99;
t13 = -t68 * t36 + t103;
t12 = -t106 * t69 + t30;
t7 = t93 + t100;
t6 = (-t56 / 0.2e1 + t76 / 0.2e1) * t71 - t105;
t5 = t12 * t69 + t13 * t68;
t4 = t10 * t69 + t11 * t68;
t3 = t66 * t36 + (t11 - t30 + (t36 + t106) * t69) * t69;
t2 = (-t12 + t30 + t11) * t68 + (t13 - t103 + (t36 - t106) * t68 + t10) * t69;
t1 = (t3 / 0.2e1 + t5 / 0.2e1) * t69 + (-t4 / 0.2e1 + t2 / 0.2e1) * t68;
t8 = [0, 0, 0, t22 * t87; 0, qJD(3) * t7 + qJD(4) * t6, qJD(2) * t7 + t18, qJD(3) * t99 * m(5) + t6 * qJD(2) + (((t56 * t68 - t37) * t72 + (t68 * t83 - t39 - t59) * t71) * t96 - t68 * t2 / 0.2e1 + (t21 * t58 - (t23 * t68 - t24 * t69) * t80) * m(5) - (t67 / 0.2e1 + t66 / 0.2e1) * t74 + ((-t56 * t69 + t38) * t72 + (t54 * t69 + t40) * t71 + t4) * t98 - (t3 + t5) * t69 / 0.2e1) * qJD(4); 0, t18 + 0.4e1 * (-t93 / 0.4e1 - t100 / 0.4e1) * qJD(2), 0, 0.2e1 * (t21 * qJD(2) / 0.4e1 - qJD(4) * t81 / 0.2e1) * m(5); 0, ((-t76 + t56) * t71 / 0.2e1 + t105) * qJD(2) + t1 * qJD(4) + qJD(3) * t82, qJD(2) * t82, t1 * qJD(2) + (m(5) * ((-t69 * t73 + (-t69 * rSges(5,3) - t80 * t68) * t68) * t22 - t58 * t81) + (t67 * t42 - t43 * t92) * t96 + (t42 * t92 - t66 * t43) * t98) * qJD(4);];
Cq = t8;
