% Calculate matrix of centrifugal and coriolis load on the joints for
% S3RPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1]';
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

function Cq = S3RPP1_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPP1_coriolismatJ_fixb_slag_vp1: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RPP1_coriolismatJ_fixb_slag_vp1: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3RPP1_coriolismatJ_fixb_slag_vp1: pkin has to be [3x1] (double)');
assert( isreal(m) && all(size(m) == [4 1]), ...
  'S3RPP1_coriolismatJ_fixb_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'S3RPP1_coriolismatJ_fixb_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'S3RPP1_coriolismatJ_fixb_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:13:33
% EndTime: 2018-11-14 10:13:33
% DurationCPUTime: 0.09s
% Computational Cost: add. (167->17), mult. (331->29), div. (0->0), fcn. (204->2), ass. (0->15)
t20 = sin(qJ(1));
t33 = t20 ^ 2;
t21 = cos(qJ(1));
t23 = rSges(4,3) + pkin(1) + qJ(3);
t10 = (rSges(4,2) + qJ(2)) * t20 + t23 * t21;
t18 = t21 * qJ(2);
t9 = t21 * rSges(4,2) - t23 * t20 + t18;
t32 = m(4) * (t21 * t10 - t20 * t9);
t31 = qJD(1) * t32;
t27 = m(3) * ((t21 * rSges(3,3) + t18) * t21 + (rSges(3,3) + qJ(2)) * t33);
t26 = m(4) * (t20 * t10 + t21 * t9);
t7 = m(4) * (t21 ^ 2 + t33);
t24 = t7 * qJD(1);
t1 = t26 + t27;
t2 = [t1 * qJD(2) + qJD(3) * t32, t1 * qJD(1), t31; -t7 * qJD(3) + 0.4e1 * (-t27 / 0.4e1 - t26 / 0.4e1) * qJD(1), 0, -t24; t7 * qJD(2) - t31, t24, 0;];
Cq  = t2;
