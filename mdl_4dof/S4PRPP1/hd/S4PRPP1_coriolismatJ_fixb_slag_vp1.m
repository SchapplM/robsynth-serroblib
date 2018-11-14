% Calculate matrix of centrifugal and coriolis load on the joints for
% S4PRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,theta1]';
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
% Datum: 2018-11-14 13:41
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Cq = S4PRPP1_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP1_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPP1_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP1_coriolismatJ_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPP1_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPP1_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPP1_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:40:59
% EndTime: 2018-11-14 13:40:59
% DurationCPUTime: 0.09s
% Computational Cost: add. (371->18), mult. (331->29), div. (0->0), fcn. (204->2), ass. (0->16)
t22 = pkin(5) + qJ(2);
t20 = sin(t22);
t34 = t20 ^ 2;
t21 = cos(t22);
t18 = t21 * qJ(3);
t24 = rSges(5,3) + pkin(2) + qJ(4);
t7 = t21 * rSges(5,2) - t24 * t20 + t18;
t8 = (rSges(5,2) + qJ(3)) * t20 + t24 * t21;
t33 = m(5) * (-t20 * t7 + t21 * t8);
t32 = qJD(2) * t33;
t28 = m(4) * ((t21 * rSges(4,3) + t18) * t21 + (rSges(4,3) + qJ(3)) * t34);
t27 = m(5) * (t20 * t8 + t21 * t7);
t9 = m(5) * (t21 ^ 2 + t34);
t25 = t9 * qJD(2);
t1 = t27 + t28;
t2 = [0, 0, 0, 0; 0, t1 * qJD(3) + qJD(4) * t33, t1 * qJD(2), t32; 0, -t9 * qJD(4) + 0.4e1 * (-t28 / 0.4e1 - t27 / 0.4e1) * qJD(2), 0, -t25; 0, t9 * qJD(3) - t32, t25, 0;];
Cq  = t2;
