% Calculate vector of inverse dynamics joint torques for
% S4PPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
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
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:40
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau = S4PPRR1_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR1_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR1_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR1_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR1_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR1_invdynJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR1_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRR1_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPRR1_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:40:08
% EndTime: 2018-11-14 13:40:09
% DurationCPUTime: 0.39s
% Computational Cost: add. (541->80), mult. (723->105), div. (0->0), fcn. (610->6), ass. (0->55)
t64 = qJ(3) + qJ(4);
t46 = sin(t64);
t49 = sin(pkin(6));
t50 = cos(pkin(6));
t58 = cos(t64);
t29 = -t49 * t46 - t50 * t58;
t30 = t50 * t46 - t49 * t58;
t11 = t29 * rSges(5,1) + t30 * rSges(5,2);
t48 = qJD(3) + qJD(4);
t77 = t48 * t11;
t52 = cos(qJ(3));
t74 = t52 * pkin(3);
t76 = t74 * t50;
t59 = t74 * t49;
t51 = sin(qJ(3));
t67 = t50 * t51;
t61 = pkin(3) * t67;
t23 = t59 - t61;
t45 = qJD(2) * t49;
t12 = -t30 * rSges(5,1) + t29 * rSges(5,2);
t71 = t48 * t12;
t5 = qJD(3) * t23 + t45 + t71;
t21 = t30 * t48;
t22 = t29 * t48;
t55 = t21 * rSges(5,1) - t22 * rSges(5,2);
t70 = t49 * t51;
t24 = pkin(3) * t70 + t76;
t65 = qJD(2) * t50;
t6 = -qJD(3) * t24 - t65 + t77;
t8 = t22 * rSges(5,1) + t21 * rSges(5,2);
t75 = t5 * t8 + t6 * t55;
t73 = qJD(3) ^ 2;
t69 = t49 * t52;
t47 = -qJDD(3) - qJDD(4);
t66 = Icges(5,3) * t47;
t63 = qJDD(2) * t50;
t62 = -m(3) - m(4) - m(5);
t57 = pkin(3) * t69 - t61;
t33 = -t50 * t52 - t70;
t54 = -t67 + t69;
t16 = t33 * rSges(4,1) - rSges(4,2) * t54;
t53 = t33 * pkin(3);
t32 = qJD(3) * t33;
t44 = qJDD(2) * t49;
t31 = t54 * qJD(3);
t17 = rSges(4,1) * t54 + t33 * rSges(4,2);
t15 = t32 * rSges(4,1) - t31 * rSges(4,2);
t14 = t31 * rSges(4,1) + t32 * rSges(4,2);
t10 = qJD(3) * t16 - t65;
t9 = qJD(3) * t17 + t45;
t4 = -qJD(3) * t14 + qJDD(3) * t16 - t63;
t3 = qJD(3) * t15 + qJDD(3) * t17 + t44;
t2 = -t73 * t54 * pkin(3) - qJDD(3) * t24 - t11 * t47 + t48 * t55 - t63;
t1 = qJDD(3) * t23 - t47 * t12 + t48 * t8 + t73 * t53 + t44;
t7 = [(-g(3) + qJDD(1)) * (m(2) - t62); t62 * (g(1) * t49 - g(2) * t50) + m(4) * (t3 * t49 - t4 * t50) + m(5) * (t1 * t49 - t2 * t50) + m(3) * (t49 ^ 2 + t50 ^ 2) * qJDD(2); Icges(4,3) * qJDD(3) - t66 + (-g(1) * (t57 + t12) - g(2) * (t53 + t11) + t1 * (t59 + t12) + t2 * (t11 - t76) + ((-t1 * t50 - t2 * t49) * t51 + (t5 * t33 - t6 * t54) * qJD(3)) * pkin(3) - t5 * (pkin(3) * t32 + t77) - t6 * (-qJD(3) * t57 - t71) + t75) * m(5) + (-t10 * t14 + t9 * t15 + (qJD(3) * t10 - g(1) + t3) * t17 + (-qJD(3) * t9 - g(2) + t4) * t16) * m(4); -t66 + ((t48 * t6 - g(1) + t1) * t12 + t75 + (-t48 * t5 - g(2) + t2) * t11) * m(5);];
tau  = t7;
