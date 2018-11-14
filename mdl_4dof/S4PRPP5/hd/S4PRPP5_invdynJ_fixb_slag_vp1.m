% Calculate vector of inverse dynamics joint torques for
% S4PRPP5
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2]';
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
% Datum: 2018-11-14 14:10
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau = S4PRPP5_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(4,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP5_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPP5_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPP5_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPP5_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S4PRPP5_invdynJ_fixb_slag_vp1: pkin has to be [4x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPP5_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPP5_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPP5_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:10:12
% EndTime: 2018-11-14 14:10:13
% DurationCPUTime: 0.48s
% Computational Cost: add. (334->102), mult. (657->106), div. (0->0), fcn. (340->2), ass. (0->49)
t52 = sin(qJ(2));
t53 = cos(qJ(2));
t32 = t53 * rSges(4,1) + t52 * rSges(4,3);
t30 = t53 * pkin(2) + t52 * qJ(3);
t42 = qJD(3) * t53;
t77 = qJD(2) * t30 - t42;
t56 = -qJD(2) * t32 - t77;
t79 = -rSges(5,1) - pkin(3);
t60 = -pkin(2) + t79;
t76 = t60 * t52;
t69 = -rSges(4,1) - pkin(2);
t75 = t69 * t52;
t31 = t53 * rSges(5,1) + t52 * rSges(5,2);
t70 = t53 * pkin(3);
t58 = t31 + t70;
t44 = t53 * qJ(3);
t26 = t52 * pkin(2) - t44;
t41 = qJD(3) * t52;
t62 = qJD(2) * t53;
t64 = qJ(3) * t62 + t41;
t74 = qJD(2) * t26 - t41 + t64;
t54 = qJD(2) ^ 2;
t73 = t52 / 0.2e1;
t72 = -t53 / 0.2e1;
t71 = -m(4) - m(5);
t47 = t53 * rSges(4,3);
t28 = t52 * rSges(4,1) - t47;
t66 = -t26 - t28;
t61 = qJD(2) * qJD(3);
t65 = qJDD(3) * t52 + t53 * t61;
t63 = qJD(2) * t52;
t48 = t53 * rSges(5,2);
t59 = t79 * t52 + t48;
t57 = -t26 + t59;
t33 = t53 * rSges(3,1) - t52 * rSges(3,2);
t29 = t52 * rSges(3,1) + t53 * rSges(3,2);
t55 = -qJDD(3) * t53 + qJDD(1) + qJD(2) * (-pkin(2) * t63 + t64) + qJDD(2) * t30 + t52 * t61;
t39 = rSges(5,2) * t62;
t38 = rSges(4,3) * t62;
t21 = qJD(2) * t29;
t20 = qJD(2) * t33 + qJD(1);
t13 = -qJD(2) * t21 + qJDD(2) * t33 + qJDD(1);
t9 = qJD(2) * t57 + t41;
t7 = t58 * qJD(2) + qJD(1) + t77;
t4 = t56 * qJD(2) + t66 * qJDD(2) + t65;
t3 = qJDD(2) * t32 + qJD(2) * (-rSges(4,1) * t63 + t38) + t55;
t2 = -t54 * t70 + t57 * qJDD(2) + (-t31 * qJD(2) - t77) * qJD(2) + t65;
t1 = qJDD(2) * t31 + qJD(2) * (-rSges(5,1) * t63 + t39) + (qJDD(2) * t53 - t54 * t52) * pkin(3) + t55;
t5 = [m(2) * qJDD(1) + m(3) * t13 + m(4) * t3 + m(5) * t1 + (-m(2) - m(3) + t71) * g(1); (Icges(3,3) + Icges(4,2) + Icges(5,3)) * qJDD(2) + ((t77 + t42) * t9 + (t39 + t74) * t7 + (-g(1) + t1) * (t30 + t58) + (-g(2) + t2) * (t44 + t48 + t76) + ((-t59 + t76) * t7 + (t58 + t60 * t53 + (-rSges(5,2) - qJ(3)) * t52) * t9) * qJD(2)) * m(5) + ((-g(1) + t3) * (t30 + t32) + (-g(2) + t4) * (t44 + t47 + t75) + (t38 + (t28 + t75) * qJD(2) + t74) * (qJD(1) - t56) + (-t56 + t42 + (t69 * t53 + (-rSges(4,3) - qJ(3)) * t52) * qJD(2)) * (qJD(2) * t66 + t41)) * m(4) + (-t20 * t21 + (t13 - g(1)) * t33 + (qJD(2) * t20 + qJDD(2) * t29 + t33 * t54 + g(2)) * t29) * m(3); t71 * (-g(1) * t53 + g(2) * t52) + 0.2e1 * (t1 * t72 + t2 * t73) * m(5) + 0.2e1 * (t3 * t72 + t4 * t73) * m(4); (-g(3) + qJDD(4)) * m(5);];
tau  = t5;
