% Calculate matrix of centrifugal and coriolis load on the joints for
% S4PPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
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
% Datum: 2018-11-14 13:40
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Cq = S4PPRR1_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR1_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR1_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR1_coriolismatJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR1_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRR1_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPRR1_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:40:13
% EndTime: 2018-11-14 13:40:13
% DurationCPUTime: 0.17s
% Computational Cost: add. (500->26), mult. (572->49), div. (0->0), fcn. (590->6), ass. (0->26)
t30 = sin(pkin(6));
t31 = cos(pkin(6));
t37 = qJ(3) + qJ(4);
t35 = sin(t37);
t36 = cos(t37);
t21 = -t30 * t35 - t31 * t36;
t22 = -t30 * t36 + t31 * t35;
t43 = t21 * rSges(5,1) + t22 * rSges(5,2);
t34 = t22 * rSges(5,1) - t21 * rSges(5,2);
t32 = sin(qJ(3));
t40 = t31 * t32;
t33 = cos(qJ(3));
t42 = t33 * pkin(3);
t10 = -pkin(3) * t40 + t42 * t30 - t34;
t41 = t30 * t32;
t11 = -pkin(3) * t41 - t42 * t31 + t43;
t1 = t10 * t43 + t11 * t34;
t39 = m(5) * qJD(3);
t38 = m(5) * qJD(4);
t24 = -t30 * t33 + t40;
t23 = -t31 * t33 - t41;
t13 = t24 * pkin(3) + t34;
t12 = t23 * pkin(3) + t43;
t7 = t30 * t43 - t31 * t34;
t6 = t7 * t38;
t2 = [0, 0, 0, 0; 0, 0, t6 + 0.2e1 * (m(4) * ((t23 * rSges(4,1) + t24 * rSges(4,2)) * t30 + (-t24 * rSges(4,1) + t23 * rSges(4,2)) * t31) / 0.2e1 + m(5) * (t12 * t30 - t13 * t31) / 0.2e1) * qJD(3), t7 * t39 + t6; 0, 0 (t10 * t12 + t11 * t13) * t39 + t1 * t38 (t38 + t39) * t1; 0, 0 (-t12 * t34 + t13 * t43) * t39, 0;];
Cq  = t2;
