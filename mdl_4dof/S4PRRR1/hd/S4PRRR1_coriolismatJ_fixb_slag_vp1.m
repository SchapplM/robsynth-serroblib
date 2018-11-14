% Calculate matrix of centrifugal and coriolis load on the joints for
% S4PRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2018-11-14 13:45
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Cq = S4PRRR1_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR1_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR1_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR1_coriolismatJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR1_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRR1_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRR1_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:44:18
% EndTime: 2018-11-14 13:44:19
% DurationCPUTime: 0.27s
% Computational Cost: add. (2907->33), mult. (1335->45), div. (0->0), fcn. (868->6), ass. (0->34)
t58 = pkin(7) + qJ(2);
t57 = qJ(3) + t58;
t60 = qJ(4) + t57;
t52 = cos(t60);
t59 = sin(t60);
t44 = -t59 * rSges(5,1) - t52 * rSges(5,2);
t53 = sin(t57);
t40 = -pkin(3) * t53 + t44;
t45 = t52 * rSges(5,1) - t59 * rSges(5,2);
t54 = cos(t57);
t41 = pkin(3) * t54 + t45;
t64 = t45 * t40 - t41 * t44;
t89 = m(5) * qJD(3) * t64;
t71 = pkin(2) * sin(t58);
t38 = t40 - t71;
t70 = pkin(2) * cos(t58);
t39 = t41 + t70;
t82 = t45 * t38 - t39 * t44;
t88 = t82 * m(5) * qJD(2);
t87 = m(4) * (t70 * (-t53 * rSges(4,1) - t54 * rSges(4,2)) + (t54 * rSges(4,1) - t53 * rSges(4,2)) * t71);
t86 = m(5) * (-t41 * t38 + t39 * t40);
t61 = m(5) * qJD(4);
t85 = t82 * t61;
t84 = t64 * t61;
t4 = t86 + t87;
t83 = t4 * qJD(2);
t79 = m(5) * (-t64 - t82);
t78 = m(5) * (t64 - t82);
t6 = t78 / 0.2e1;
t5 = t79 / 0.2e1;
t3 = t6 - t79 / 0.2e1;
t2 = t5 - t78 / 0.2e1;
t1 = t5 + t6;
t7 = [0, 0, 0, 0; 0, t4 * qJD(3) - t85, t83 + t1 * qJD(4) + 0.2e1 * (t86 / 0.2e1 + t87 / 0.2e1) * qJD(3), t1 * qJD(3) - t85 - t88; 0, t2 * qJD(4) - t83, -t84, t2 * qJD(2) - t84 - t89; 0, t3 * qJD(3) + t88, t3 * qJD(2) + t89, 0;];
Cq  = t7;
