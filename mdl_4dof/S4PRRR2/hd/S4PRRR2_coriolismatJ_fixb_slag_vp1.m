% Calculate matrix of centrifugal and coriolis load on the joints for
% S4PRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
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
% Datum: 2019-07-18 13:27
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4PRRR2_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(2,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR2_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR2_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S4PRRR2_coriolismatJ_fixb_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR2_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRR2_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRR2_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:27:22
% EndTime: 2019-07-18 13:27:22
% DurationCPUTime: 0.24s
% Computational Cost: add. (2039->31), mult. (1335->45), div. (0->0), fcn. (868->6), ass. (0->33)
t57 = qJ(2) + qJ(3);
t56 = qJ(4) + t57;
t52 = sin(t56);
t53 = cos(t56);
t44 = t52 * rSges(5,1) + t53 * rSges(5,2);
t54 = sin(t57);
t40 = pkin(2) * t54 + t44;
t45 = -t53 * rSges(5,1) + t52 * rSges(5,2);
t55 = cos(t57);
t41 = -pkin(2) * t55 + t45;
t16 = -t45 * t40 + t41 * t44;
t87 = t16 * m(5) * qJD(3);
t69 = cos(qJ(2)) * pkin(1);
t70 = sin(qJ(2)) * pkin(1);
t86 = m(4) * (-t69 * (t54 * rSges(4,1) + t55 * rSges(4,2)) - (-t55 * rSges(4,1) + t54 * rSges(4,2)) * t70);
t38 = t40 + t70;
t39 = t41 - t69;
t85 = m(5) * (-t41 * t38 + t39 * t40);
t15 = -t45 * t38 + t39 * t44;
t60 = m(5) * qJD(4);
t84 = t15 * t60;
t83 = t16 * t60;
t4 = t85 + t86;
t82 = t4 * qJD(2);
t81 = m(5) * t15 * qJD(2);
t78 = m(5) * (t16 + t15);
t77 = m(5) * (-t16 + t15);
t6 = t77 / 0.2e1;
t5 = t78 / 0.2e1;
t3 = t6 - t78 / 0.2e1;
t2 = t5 - t77 / 0.2e1;
t1 = t5 + t6;
t7 = [0, 0, 0, 0; 0, t4 * qJD(3) + t84, t82 + t1 * qJD(4) + 0.2e1 * (t85 / 0.2e1 + t86 / 0.2e1) * qJD(3), t1 * qJD(3) + t81 + t84; 0, t2 * qJD(4) - t82, t83, t2 * qJD(2) + t83 + t87; 0, t3 * qJD(3) - t81, t3 * qJD(2) - t87, 0;];
Cq  = t7;
