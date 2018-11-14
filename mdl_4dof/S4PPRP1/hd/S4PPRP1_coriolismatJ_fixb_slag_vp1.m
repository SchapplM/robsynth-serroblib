% Calculate matrix of centrifugal and coriolis load on the joints for
% S4PPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,theta1]';
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
% Datum: 2018-11-14 13:39
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Cq = S4PPRP1_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP1_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRP1_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPRP1_coriolismatJ_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRP1_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRP1_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPRP1_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:39:14
% EndTime: 2018-11-14 13:39:14
% DurationCPUTime: 0.09s
% Computational Cost: add. (184->15), mult. (386->32), div. (0->0), fcn. (432->4), ass. (0->16)
t28 = m(5) * qJD(3);
t18 = sin(pkin(5));
t19 = cos(pkin(5));
t23 = sin(qJ(3));
t24 = cos(qJ(3));
t12 = -t18 * t23 - t19 * t24;
t13 = -t18 * t24 + t19 * t23;
t25 = rSges(5,1) + pkin(3);
t26 = -rSges(5,3) - qJ(4);
t17 = t25 * t12 + t26 * t13;
t5 = t26 * t12 - t25 * t13;
t1 = -t12 * t5 - t17 * t13;
t27 = t1 * t28;
t21 = m(5) * qJD(4);
t8 = -t12 * t18 + t13 * t19;
t2 = [0, 0, 0, 0; 0, 0, 0.2e1 * (m(4) * ((-t13 * rSges(4,1) + t12 * rSges(4,2)) * t19 + (t12 * rSges(4,1) + t13 * rSges(4,2)) * t18) / 0.2e1 + m(5) * (t17 * t18 + t5 * t19) / 0.2e1) * qJD(3) + t8 * t21, t8 * t28; 0, 0, t1 * t21, t27; 0, 0, -t27, 0;];
Cq  = t2;
