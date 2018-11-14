% Calculate matrix of centrifugal and coriolis load on the joints for
% S4PPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta2]';
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
% Datum: 2018-11-14 14:08
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Cq = S4PPRR3_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR3_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR3_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR3_coriolismatJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR3_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRR3_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPRR3_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:08:16
% EndTime: 2018-11-14 14:08:16
% DurationCPUTime: 0.08s
% Computational Cost: add. (334->13), mult. (238->23), div. (0->0), fcn. (131->4), ass. (0->14)
t20 = pkin(6) + qJ(3);
t19 = qJ(4) + t20;
t16 = cos(t19);
t24 = sin(t19);
t11 = -t24 * rSges(5,1) - t16 * rSges(5,2);
t12 = t16 * rSges(5,1) - t24 * rSges(5,2);
t18 = cos(t20);
t17 = sin(t20);
t8 = -pkin(3) * t17 + t11;
t1 = t11 * (pkin(3) * t18 + t12) - t8 * t12;
t22 = m(5) * qJD(3);
t21 = m(5) * qJD(4);
t10 = t11 * t21;
t2 = [0, 0, t10 + 0.2e1 * (m(4) * (-t17 * rSges(4,1) - t18 * rSges(4,2)) / 0.2e1 + m(5) * t8 / 0.2e1) * qJD(3), t11 * t22 + t10; 0, 0, 0, 0; 0, 0, t1 * t21 (t21 + t22) * t1; 0, 0, -t1 * t22, 0;];
Cq  = t2;
