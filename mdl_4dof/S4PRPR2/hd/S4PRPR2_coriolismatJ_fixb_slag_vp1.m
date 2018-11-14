% Calculate matrix of centrifugal and coriolis load on the joints for
% S4PRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta3]';
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
% Datum: 2018-11-14 14:02
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Cq = S4PRPR2_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR2_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR2_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR2_coriolismatJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR2_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPR2_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPR2_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:02:11
% EndTime: 2018-11-14 14:02:11
% DurationCPUTime: 0.08s
% Computational Cost: add. (364->18), mult. (282->29), div. (0->0), fcn. (163->6), ass. (0->17)
t22 = qJ(2) + pkin(6);
t21 = qJ(4) + t22;
t18 = cos(t21);
t29 = sin(t21);
t13 = -t29 * rSges(5,1) - t18 * rSges(5,2);
t14 = t18 * rSges(5,1) - t29 * rSges(5,2);
t20 = cos(t22);
t24 = cos(qJ(2));
t19 = sin(t22);
t23 = sin(qJ(2));
t28 = t23 * pkin(2);
t8 = -pkin(3) * t19 + t13 - t28;
t1 = (t24 * pkin(2) + pkin(3) * t20 + t14) * t13 - t14 * t8;
t26 = m(5) * qJD(2);
t25 = m(5) * qJD(4);
t10 = t13 * t25;
t2 = [0, t10 + 0.2e1 * (m(3) * (-t23 * rSges(3,1) - t24 * rSges(3,2)) / 0.2e1 + m(4) * (-t19 * rSges(4,1) - t20 * rSges(4,2) - t28) / 0.2e1 + m(5) * t8 / 0.2e1) * qJD(2), 0, t13 * t26 + t10; 0, t1 * t25, 0 (t25 + t26) * t1; 0, 0, 0, 0; 0, -t1 * t26, 0, 0;];
Cq  = t2;
