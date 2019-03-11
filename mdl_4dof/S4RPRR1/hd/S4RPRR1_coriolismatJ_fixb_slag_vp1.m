% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2019-03-08 18:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RPRR1_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR1_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR1_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR1_coriolismatJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR1_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR1_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR1_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:31:39
% EndTime: 2019-03-08 18:31:40
% DurationCPUTime: 0.24s
% Computational Cost: add. (2985->35), mult. (1421->47), div. (0->0), fcn. (948->8), ass. (0->34)
t60 = qJ(1) + pkin(7);
t59 = qJ(3) + t60;
t66 = qJ(4) + t59;
t54 = cos(t66);
t63 = sin(t66);
t44 = -t63 * rSges(5,1) - t54 * rSges(5,2);
t55 = sin(t59);
t42 = -pkin(3) * t55 + t44;
t45 = t54 * rSges(5,1) - t63 * rSges(5,2);
t56 = cos(t59);
t43 = pkin(3) * t56 + t45;
t70 = t45 * t42 - t43 * t44;
t93 = m(5) * qJD(3) * t70;
t65 = -sin(qJ(1)) * pkin(1) - pkin(2) * sin(t60);
t36 = t42 + t65;
t64 = cos(qJ(1)) * pkin(1) + pkin(2) * cos(t60);
t37 = t43 + t64;
t86 = t45 * t36 - t37 * t44;
t92 = t86 * m(5) * qJD(1);
t91 = m(4) * (t64 * (-t55 * rSges(4,1) - t56 * rSges(4,2)) - (t56 * rSges(4,1) - t55 * rSges(4,2)) * t65);
t90 = m(5) * (-t43 * t36 + t37 * t42);
t67 = m(5) * qJD(4);
t89 = t86 * t67;
t88 = t70 * t67;
t4 = t90 + t91;
t87 = t4 * qJD(1);
t83 = m(5) * (-t70 - t86);
t82 = m(5) * (t70 - t86);
t6 = t82 / 0.2e1;
t5 = t83 / 0.2e1;
t3 = t6 - t83 / 0.2e1;
t2 = t5 - t82 / 0.2e1;
t1 = t5 + t6;
t7 = [t4 * qJD(3) - t89, 0, t87 + t1 * qJD(4) + 0.2e1 * (t90 / 0.2e1 + t91 / 0.2e1) * qJD(3), t1 * qJD(3) - t89 - t92; 0, 0, 0, 0; t2 * qJD(4) - t87, 0, -t88, t2 * qJD(1) - t88 - t93; t3 * qJD(3) + t92, 0, t3 * qJD(1) + t93, 0;];
Cq  = t7;
