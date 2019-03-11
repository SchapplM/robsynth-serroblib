% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RPRP1_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP1_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP1_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP1_coriolismatJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP1_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRP1_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRP1_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:29:42
% EndTime: 2019-03-08 18:29:42
% DurationCPUTime: 0.18s
% Computational Cost: add. (2088->31), mult. (1112->45), div. (0->0), fcn. (816->6), ass. (0->29)
t52 = qJ(1) + pkin(6);
t51 = qJ(3) + t52;
t47 = sin(t51);
t48 = cos(t51);
t55 = cos(qJ(1)) * pkin(1) + pkin(2) * cos(t52);
t56 = -sin(qJ(1)) * pkin(1) - pkin(2) * sin(t52);
t78 = m(4) * (t55 * (-t47 * rSges(4,1) - t48 * rSges(4,2)) - (t48 * rSges(4,1) - t47 * rSges(4,2)) * t56);
t63 = rSges(5,1) + pkin(3);
t74 = rSges(5,3) + qJ(4);
t34 = -t63 * t47 + t74 * t48;
t27 = t34 + t56;
t35 = t74 * t47 + t63 * t48;
t28 = t35 + t55;
t77 = m(5) * (-t35 * t27 + t28 * t34);
t17 = t34 * t48 + t35 * t47;
t76 = t17 * m(5) * qJD(3);
t4 = t77 + t78;
t75 = t4 * qJD(1);
t13 = t27 * t48 + t28 * t47;
t73 = t13 * m(5) * qJD(1);
t70 = m(5) * (t13 - t17);
t67 = m(5) * (t17 + t13);
t57 = m(5) * qJD(4);
t8 = t67 / 0.2e1;
t7 = t70 / 0.2e1;
t3 = t8 - t70 / 0.2e1;
t2 = t8 + t7;
t1 = t7 - t67 / 0.2e1;
t5 = [t4 * qJD(3) + t13 * t57, 0, t75 + t2 * qJD(4) + 0.2e1 * (t77 / 0.2e1 + t78 / 0.2e1) * qJD(3), t2 * qJD(3) + t73; 0, 0, 0, 0; t3 * qJD(4) - t75, 0, t17 * t57, t3 * qJD(1) + t76; t1 * qJD(3) - t73, 0, t1 * qJD(1) - t76, 0;];
Cq  = t5;
