% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
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
% Datum: 2018-11-14 13:52
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Cq = S4RRPP2_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP2_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP2_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP2_coriolismatJ_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP2_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPP2_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPP2_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:52:25
% EndTime: 2018-11-14 13:52:25
% DurationCPUTime: 0.31s
% Computational Cost: add. (2614->49), mult. (1970->70), div. (0->0), fcn. (1484->4), ass. (0->39)
t108 = m(5) / 0.2e1;
t109 = m(4) / 0.2e1;
t73 = qJ(1) + qJ(2);
t72 = cos(t73);
t67 = t72 * qJ(3);
t71 = sin(t73);
t77 = rSges(5,1) + pkin(2) + pkin(3);
t48 = t72 * rSges(5,2) - t77 * t71 + t67;
t87 = sin(qJ(1)) * pkin(1);
t44 = t48 - t87;
t49 = (rSges(5,2) + qJ(3)) * t71 + t77 * t72;
t86 = cos(qJ(1)) * pkin(1);
t45 = t49 + t86;
t21 = t44 * t72 + t45 * t71;
t24 = t48 * t72 + t49 * t71;
t85 = rSges(4,1) + pkin(2);
t59 = t72 * rSges(4,3) - t85 * t71 + t67;
t55 = t59 - t87;
t60 = t85 * t72 + (rSges(4,3) + qJ(3)) * t71;
t56 = t60 + t86;
t25 = t55 * t72 + t56 * t71;
t33 = t59 * t72 + t60 * t71;
t83 = (t33 + t25) * t109 + (t24 + t21) * t108;
t84 = (t25 - t33) * t109 + (t21 - t24) * t108;
t1 = t84 - t83;
t110 = t1 * qJD(1);
t101 = m(3) * ((-t71 * rSges(3,1) - t72 * rSges(3,2)) * t86 + (t72 * rSges(3,1) - t71 * rSges(3,2)) * t87);
t97 = m(4) * (-t60 * t55 + t56 * t59);
t91 = m(5) * (-t49 * t44 + t45 * t48);
t107 = 0.4e1 * qJD(1);
t95 = m(4) * t25;
t94 = m(4) * t33;
t89 = m(5) * t21;
t88 = m(5) * t24;
t12 = t88 + t94;
t7 = t89 + t95;
t4 = t91 + t97 + t101;
t2 = t83 + t84;
t3 = [t4 * qJD(2) + t7 * qJD(3), t4 * qJD(1) + t2 * qJD(3) + 0.2e1 * (t91 / 0.2e1 + t97 / 0.2e1 + t101 / 0.2e1) * qJD(2), t7 * qJD(1) + t2 * qJD(2), 0; -t1 * qJD(3) + (-t101 / 0.4e1 - t97 / 0.4e1 - t91 / 0.4e1) * t107, t12 * qJD(3), t12 * qJD(2) - t110, 0; t1 * qJD(2) + (-t95 / 0.4e1 - t89 / 0.4e1) * t107, t110 + 0.4e1 * (-t94 / 0.4e1 - t88 / 0.4e1) * qJD(2), 0, 0; 0, 0, 0, 0;];
Cq  = t3;
