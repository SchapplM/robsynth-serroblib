% Calculate matrix of centrifugal and coriolis load on the joints for
% S4PRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
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
% Datum: 2018-11-14 13:43
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Cq = S4PRRP1_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP1_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP1_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP1_coriolismatJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP1_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRP1_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRP1_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:43:27
% EndTime: 2018-11-14 13:43:27
% DurationCPUTime: 0.19s
% Computational Cost: add. (2014->29), mult. (1030->43), div. (0->0), fcn. (740->4), ass. (0->29)
t50 = pkin(6) + qJ(2);
t49 = qJ(3) + t50;
t45 = sin(t49);
t46 = cos(t49);
t58 = pkin(2) * cos(t50);
t59 = pkin(2) * sin(t50);
t74 = m(4) * (t58 * (-t45 * rSges(4,1) - t46 * rSges(4,2)) + (t46 * rSges(4,1) - t45 * rSges(4,2)) * t59);
t57 = rSges(5,1) + pkin(3);
t70 = rSges(5,3) + qJ(4);
t34 = -t57 * t45 + t70 * t46;
t30 = t34 - t59;
t35 = t70 * t45 + t57 * t46;
t31 = t35 + t58;
t73 = m(5) * (-t35 * t30 + t31 * t34);
t17 = t34 * t46 + t35 * t45;
t72 = m(5) * t17 * qJD(3);
t4 = t73 + t74;
t71 = t4 * qJD(2);
t13 = t30 * t46 + t31 * t45;
t69 = m(5) * t13 * qJD(2);
t66 = m(5) * (t13 - t17);
t63 = m(5) * (t17 + t13);
t51 = m(5) * qJD(4);
t8 = t63 / 0.2e1;
t7 = t66 / 0.2e1;
t3 = t8 - t66 / 0.2e1;
t2 = t8 + t7;
t1 = t7 - t63 / 0.2e1;
t5 = [0, 0, 0, 0; 0, t4 * qJD(3) + t13 * t51, t71 + t2 * qJD(4) + 0.2e1 * (t73 / 0.2e1 + t74 / 0.2e1) * qJD(3), t2 * qJD(3) + t69; 0, t3 * qJD(4) - t71, t17 * t51, t3 * qJD(2) + t72; 0, t1 * qJD(3) - t69, t1 * qJD(2) - t72, 0;];
Cq  = t5;
