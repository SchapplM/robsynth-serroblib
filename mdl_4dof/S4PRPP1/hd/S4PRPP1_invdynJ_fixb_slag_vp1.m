% Calculate vector of inverse dynamics joint torques for
% S4PRPP1
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
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:41
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau = S4PRPP1_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP1_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPP1_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPP1_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPP1_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP1_invdynJ_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPP1_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPP1_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPP1_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:40:59
% EndTime: 2018-11-14 13:41:00
% DurationCPUTime: 0.37s
% Computational Cost: add. (767->101), mult. (729->114), div. (0->0), fcn. (398->2), ass. (0->54)
t52 = pkin(5) + qJ(2);
t50 = sin(t52);
t51 = cos(t52);
t28 = t51 * pkin(2) + t50 * qJ(3);
t42 = qJD(3) * t51;
t19 = t28 * qJD(2) - t42;
t25 = t50 * rSges(5,2) + t51 * rSges(5,3);
t44 = t51 * qJ(3);
t24 = t50 * pkin(2) - t44;
t29 = t51 * rSges(5,2) - t50 * rSges(5,3);
t57 = -t50 * qJ(4) + t29;
t54 = -t24 + t57;
t53 = qJD(2) ^ 2;
t55 = -t53 * qJ(4) + qJDD(4);
t62 = qJD(4) * t50;
t61 = qJD(2) * qJD(3);
t70 = qJDD(3) * t50 + t51 * t61;
t1 = t55 * t51 + t54 * qJDD(2) + (-t25 * qJD(2) - t19 - 0.2e1 * t62) * qJD(2) + t70;
t82 = -g(1) + t1;
t63 = qJD(2) * t51;
t38 = rSges(5,2) * t63;
t64 = qJD(2) * t50;
t35 = qJ(3) * t63;
t41 = qJD(3) * t50;
t69 = t35 + t41;
t59 = qJD(2) * (-pkin(2) * t64 + t69) + qJDD(2) * t28 + t50 * t61;
t2 = qJD(2) * t38 + qJDD(2) * t25 + (-t53 * rSges(5,3) + t55) * t50 + (qJDD(2) * qJ(4) + 0.2e1 * qJD(4) * qJD(2) - qJDD(3)) * t51 + t59;
t81 = -g(2) + t2;
t80 = t51 * rSges(4,2) - t50 * rSges(4,3);
t79 = t51 * qJ(4) + t25 + t28;
t78 = t50 / 0.2e1;
t77 = -t51 / 0.2e1;
t76 = -m(4) - m(5);
t75 = rSges(4,2) - pkin(2);
t73 = t51 * rSges(4,3);
t72 = -pkin(2) - qJ(4);
t26 = t50 * rSges(4,2) + t73;
t71 = -t24 + t26;
t17 = t28 - t80;
t68 = rSges(4,2) * t64 + rSges(4,3) * t63;
t40 = qJD(4) * t51;
t67 = t40 + t41;
t66 = -qJD(2) * t24 + t41;
t60 = -rSges(5,3) + t72;
t56 = t42 - t62;
t31 = t51 * rSges(3,1) - t50 * rSges(3,2);
t27 = t50 * rSges(3,1) + t51 * rSges(3,2);
t13 = t17 * qJD(2) - t42;
t12 = t71 * qJD(2) + t41;
t9 = t79 * qJD(2) - t56;
t8 = t54 * qJD(2) + t67;
t4 = qJD(2) * t68 - qJDD(2) * t80 - qJDD(3) * t51 + t59;
t3 = t71 * qJDD(2) + (t80 * qJD(2) - t19) * qJD(2) + t70;
t5 = [(-g(3) + qJDD(1)) * (m(2) + m(3) - t76); -m(3) * (-g(1) * t27 + g(2) * t31) + (-(t57 * qJD(2) + t40 + t66 - t8) * t9 + t8 * t56 + t9 * (t35 + t38 + t67) + (t8 * t60 * t51 + (t8 * (-rSges(5,2) - qJ(3)) + t9 * t60) * t50) * qJD(2) + t81 * t79 + t82 * (t72 * t50 + t29 + t44)) * m(5) + (-(qJD(2) * t26 - t12 + t66) * t13 + t12 * t42 + t13 * (t68 + t69) + (t12 * t75 * t51 + (t12 * (-rSges(4,3) - qJ(3)) - t13 * pkin(2)) * t50) * qJD(2) + (-g(2) + t4) * t17 + (-g(1) + t3) * (t75 * t50 + t44 + t73)) * m(4) + (Icges(3,3) + Icges(4,1) + Icges(5,1) + m(3) * (t27 ^ 2 + t31 ^ 2)) * qJDD(2); t76 * (g(1) * t50 - g(2) * t51) + 0.2e1 * (t1 * t78 + t2 * t77) * m(5) + 0.2e1 * (t3 * t78 + t4 * t77) * m(4); (t81 * t50 + t82 * t51) * m(5);];
tau  = t5;
