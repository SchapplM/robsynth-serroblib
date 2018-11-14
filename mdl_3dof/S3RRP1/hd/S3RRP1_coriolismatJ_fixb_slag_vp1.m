% Calculate matrix of centrifugal and coriolis load on the joints for
% S3RRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2]';
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
% Datum: 2018-11-14 10:15
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Cq = S3RRP1_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRP1_coriolismatJ_fixb_slag_vp1: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RRP1_coriolismatJ_fixb_slag_vp1: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RRP1_coriolismatJ_fixb_slag_vp1: pkin has to be [4x1] (double)');
assert( isreal(m) && all(size(m) == [4 1]), ...
  'S3RRP1_coriolismatJ_fixb_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'S3RRP1_coriolismatJ_fixb_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'S3RRP1_coriolismatJ_fixb_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:15:07
% EndTime: 2018-11-14 10:15:08
% DurationCPUTime: 0.17s
% Computational Cost: add. (1274->28), mult. (1030->43), div. (0->0), fcn. (740->4), ass. (0->28)
t47 = qJ(1) + qJ(2);
t45 = sin(t47);
t46 = cos(t47);
t57 = cos(qJ(1)) * pkin(1);
t58 = sin(qJ(1)) * pkin(1);
t73 = m(3) * (t57 * (-t45 * rSges(3,1) - t46 * rSges(3,2)) + (t46 * rSges(3,1) - t45 * rSges(3,2)) * t58);
t56 = rSges(4,1) + pkin(2);
t69 = rSges(4,3) + qJ(3);
t34 = -t56 * t45 + t69 * t46;
t30 = t34 - t58;
t35 = t69 * t45 + t56 * t46;
t31 = t35 + t57;
t72 = m(4) * (-t35 * t30 + t31 * t34);
t17 = t34 * t46 + t35 * t45;
t71 = m(4) * t17 * qJD(2);
t4 = t72 + t73;
t70 = t4 * qJD(1);
t13 = t30 * t46 + t31 * t45;
t68 = m(4) * t13 * qJD(1);
t65 = m(4) * (t13 - t17);
t62 = m(4) * (t17 + t13);
t50 = m(4) * qJD(3);
t8 = t62 / 0.2e1;
t7 = t65 / 0.2e1;
t3 = t8 - t65 / 0.2e1;
t2 = t8 + t7;
t1 = t7 - t62 / 0.2e1;
t5 = [t4 * qJD(2) + t13 * t50, t70 + t2 * qJD(3) + 0.2e1 * (t72 / 0.2e1 + t73 / 0.2e1) * qJD(2), t2 * qJD(2) + t68; t3 * qJD(3) - t70, t17 * t50, t3 * qJD(1) + t71; t1 * qJD(2) - t68, t1 * qJD(1) - t71, 0;];
Cq  = t5;
