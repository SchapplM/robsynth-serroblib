% Calculate matrix of centrifugal and coriolis load on the joints for
% S4PRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3]';
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
% Datum: 2018-11-14 14:03
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Cq = S4PRRP2_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP2_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP2_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRRP2_coriolismatJ_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP2_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRP2_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRP2_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:03:19
% EndTime: 2018-11-14 14:03:19
% DurationCPUTime: 0.12s
% Computational Cost: add. (508->22), mult. (500->36), div. (0->0), fcn. (304->4), ass. (0->23)
t32 = qJ(2) + qJ(3);
t30 = sin(t32);
t31 = cos(t32);
t37 = rSges(5,1) + pkin(3);
t20 = -t31 * rSges(5,2) - t37 * t30;
t33 = sin(qJ(2));
t39 = t33 * pkin(2);
t18 = t20 - t39;
t21 = -t30 * rSges(5,2) + t37 * t31;
t34 = cos(qJ(2));
t38 = t34 * pkin(2);
t4 = -t21 * t18 + (t21 + t38) * t20;
t24 = -t30 * rSges(4,1) - t31 * rSges(4,2);
t22 = t24 - t39;
t25 = t31 * rSges(4,1) - t30 * rSges(4,2);
t6 = -t25 * t22 + (t25 + t38) * t24;
t1 = m(4) * t6 + m(5) * t4;
t46 = t1 * qJD(2);
t45 = m(4) / 0.2e1;
t44 = m(5) / 0.2e1;
t11 = m(4) * t24 + m(5) * t20;
t8 = t11 * qJD(3);
t2 = [0, t8 + 0.2e1 * (m(3) * (-t33 * rSges(3,1) - t34 * rSges(3,2)) / 0.2e1 + t22 * t45 + t18 * t44) * qJD(2), t11 * qJD(2) + t8, 0; 0, t1 * qJD(3), t46 + 0.2e1 * (t4 * t44 + t6 * t45) * qJD(3), 0; 0, -t46, 0, 0; 0, 0, 0, 0;];
Cq  = t2;
