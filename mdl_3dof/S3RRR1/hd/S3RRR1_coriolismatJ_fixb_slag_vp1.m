% Calculate matrix of centrifugal and coriolis load on the joints for
% S3RRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2,d3]';
% m_mdh [4x1]
%   mass of all robot links (including the base)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [4x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Cq [3x3]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:16
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Cq = S3RRR1_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(5,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRR1_coriolismatJ_fixb_slag_vp1: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RRR1_coriolismatJ_fixb_slag_vp1: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S3RRR1_coriolismatJ_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [4 1]), ...
  'S3RRR1_coriolismatJ_fixb_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'S3RRR1_coriolismatJ_fixb_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'S3RRR1_coriolismatJ_fixb_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:15:54
% EndTime: 2018-11-14 10:15:54
% DurationCPUTime: 0.24s
% Computational Cost: add. (2039->32), mult. (1335->45), div. (0->0), fcn. (868->6), ass. (0->33)
t55 = qJ(1) + qJ(2);
t59 = qJ(3) + t55;
t52 = cos(t59);
t58 = sin(t59);
t44 = -t58 * rSges(4,1) - t52 * rSges(4,2);
t53 = sin(t55);
t40 = -pkin(2) * t53 + t44;
t45 = t52 * rSges(4,1) - t58 * rSges(4,2);
t54 = cos(t55);
t41 = pkin(2) * t54 + t45;
t63 = t45 * t40 - t41 * t44;
t88 = m(4) * qJD(2) * t63;
t70 = sin(qJ(1)) * pkin(1);
t38 = t40 - t70;
t69 = cos(qJ(1)) * pkin(1);
t39 = t41 + t69;
t81 = t45 * t38 - t39 * t44;
t87 = t81 * m(4) * qJD(1);
t86 = m(3) * (t69 * (-t53 * rSges(3,1) - t54 * rSges(3,2)) + (t54 * rSges(3,1) - t53 * rSges(3,2)) * t70);
t85 = m(4) * (-t41 * t38 + t39 * t40);
t60 = m(4) * qJD(3);
t84 = t81 * t60;
t83 = t63 * t60;
t4 = t85 + t86;
t82 = t4 * qJD(1);
t78 = m(4) * (-t63 - t81);
t77 = m(4) * (t63 - t81);
t6 = t77 / 0.2e1;
t5 = t78 / 0.2e1;
t3 = t6 - t78 / 0.2e1;
t2 = t5 - t77 / 0.2e1;
t1 = t5 + t6;
t7 = [t4 * qJD(2) - t84, t82 + t1 * qJD(3) + 0.2e1 * (t85 / 0.2e1 + t86 / 0.2e1) * qJD(2), t1 * qJD(2) - t84 - t87; t2 * qJD(3) - t82, -t83, t2 * qJD(1) - t83 - t88; t3 * qJD(2) + t87, t3 * qJD(1) + t88, 0;];
Cq  = t7;
