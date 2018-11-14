% Calculate matrix of centrifugal and coriolis load on the joints for
% S3RPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d3]';
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
% Datum: 2018-11-14 10:14
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Cq = S3RPR1_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPR1_coriolismatJ_fixb_slag_vp1: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RPR1_coriolismatJ_fixb_slag_vp1: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RPR1_coriolismatJ_fixb_slag_vp1: pkin has to be [4x1] (double)');
assert( isreal(m) && all(size(m) == [4 1]), ...
  'S3RPR1_coriolismatJ_fixb_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'S3RPR1_coriolismatJ_fixb_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'S3RPR1_coriolismatJ_fixb_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:14:28
% EndTime: 2018-11-14 10:14:28
% DurationCPUTime: 0.17s
% Computational Cost: add. (530->24), mult. (1070->41), div. (0->0), fcn. (1066->4), ass. (0->26)
t37 = sin(qJ(1));
t38 = cos(qJ(1));
t50 = sin(qJ(3));
t51 = cos(qJ(3));
t29 = -t37 * t50 - t38 * t51;
t30 = -t37 * t51 + t38 * t50;
t22 = t29 * rSges(4,1) + t30 * rSges(4,2);
t40 = -t30 * rSges(4,1) + t29 * rSges(4,2);
t60 = t22 * t37 + t38 * t40;
t64 = m(4) * t60;
t10 = -t64 / 0.2e1;
t11 = t64 / 0.2e1;
t56 = pkin(1) + pkin(2);
t19 = t37 * qJ(2) + t56 * t38 - t22;
t36 = t38 * qJ(2);
t39 = -t56 * t37 + t36 - t40;
t4 = t19 * t40 - t22 * t39;
t41 = m(4) * qJD(3);
t63 = t4 * t41;
t5 = m(3) * ((t38 * rSges(3,3) + t36) * t38 + (rSges(3,3) + qJ(2)) * t37 ^ 2) + m(4) * (t19 * t37 + t38 * t39);
t62 = t5 * qJD(1);
t61 = m(4) * t4 * qJD(1);
t3 = 0.2e1 * t11;
t2 = t10 + t11;
t1 = 0.2e1 * t10;
t6 = [t5 * qJD(2) + t63, t2 * qJD(3) + t62, t2 * qJD(2) + t61 - t63; t1 * qJD(3) - t62, 0, t1 * qJD(1) + t60 * t41; t3 * qJD(2) - t61, t3 * qJD(1), 0;];
Cq  = t6;
