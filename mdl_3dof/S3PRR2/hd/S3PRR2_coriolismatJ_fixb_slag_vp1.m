% Calculate matrix of centrifugal and coriolis load on the joints for
% S3PRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2,d3]';
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
% Datum: 2018-11-14 10:13
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Cq = S3PRR2_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRR2_coriolismatJ_fixb_slag_vp1: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PRR2_coriolismatJ_fixb_slag_vp1: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3PRR2_coriolismatJ_fixb_slag_vp1: pkin has to be [4x1] (double)');
assert( isreal(m) && all(size(m) == [4 1]), ...
  'S3PRR2_coriolismatJ_fixb_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'S3PRR2_coriolismatJ_fixb_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'S3PRR2_coriolismatJ_fixb_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:12:45
% EndTime: 2018-11-14 10:12:45
% DurationCPUTime: 0.08s
% Computational Cost: add. (203->12), mult. (238->23), div. (0->0), fcn. (131->4), ass. (0->13)
t17 = qJ(2) + qJ(3);
t16 = cos(t17);
t23 = sin(t17);
t11 = -t23 * rSges(4,1) - t16 * rSges(4,2);
t12 = t16 * rSges(4,1) - t23 * rSges(4,2);
t19 = cos(qJ(2));
t18 = sin(qJ(2));
t9 = -t18 * pkin(2) + t11;
t1 = t11 * (t19 * pkin(2) + t12) - t9 * t12;
t21 = m(4) * qJD(2);
t20 = m(4) * qJD(3);
t8 = t11 * t20;
t2 = [0, t8 + 0.2e1 * (m(3) * (-t18 * rSges(3,1) - t19 * rSges(3,2)) / 0.2e1 + m(4) * t9 / 0.2e1) * qJD(2), t11 * t21 + t8; 0, t1 * t20 (t20 + t21) * t1; 0, -t1 * t21, 0;];
Cq  = t2;
