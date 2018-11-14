% Calculate vector of inverse dynamics joint torques for
% S4PRRP3
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
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:12
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau = S4PRRP3_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP3_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP3_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP3_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP3_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRRP3_invdynJ_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP3_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRP3_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRP3_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:12:18
% EndTime: 2018-11-14 14:12:19
% DurationCPUTime: 0.39s
% Computational Cost: add. (539->88), mult. (573->92), div. (0->0), fcn. (258->4), ass. (0->51)
t43 = sin(qJ(2));
t44 = cos(qJ(2));
t45 = qJD(2) ^ 2;
t67 = (-qJDD(2) * t43 - t45 * t44) * pkin(2) - g(2);
t55 = pkin(2) * qJD(2);
t53 = t43 * t55;
t42 = qJ(2) + qJ(3);
t37 = sin(t42);
t38 = cos(t42);
t21 = t37 * rSges(4,1) + t38 * rSges(4,2);
t41 = qJD(2) + qJD(3);
t57 = t41 * t21;
t9 = -t53 - t57;
t59 = t38 * rSges(5,2);
t20 = t37 * rSges(5,1) + t59;
t60 = t37 * t41;
t26 = rSges(5,2) * t60;
t40 = -qJDD(2) - qJDD(3);
t58 = t38 * t41;
t61 = -rSges(5,1) - pkin(3);
t66 = (pkin(3) * t37 + t20) * t40 + (t58 * t61 + t26) * t41 + t67;
t13 = rSges(4,1) * t58 - rSges(4,2) * t60;
t65 = -t41 * t13 + t40 * t21 + t67;
t52 = t61 * t37;
t14 = t52 - t59;
t22 = t38 * rSges(5,1) - t37 * rSges(5,2);
t15 = pkin(3) * t38 + t22;
t39 = t44 * pkin(2);
t62 = t43 * pkin(2);
t48 = qJDD(2) * t39 - t45 * t62 + qJDD(1);
t1 = t14 * t41 ^ 2 - t15 * t40 + t48;
t64 = t1 - g(1);
t23 = t38 * rSges(4,1) - t37 * rSges(4,2);
t3 = -t40 * t23 - t41 * t57 + t48;
t63 = t3 - g(1);
t56 = -pkin(3) * t58 - t41 * t22;
t36 = t44 * t55;
t54 = -t36 + t56;
t51 = -pkin(3) * t60 - t41 * t20;
t50 = -t41 * t23 - t36;
t31 = t44 * rSges(3,1) - t43 * rSges(3,2);
t30 = t43 * rSges(3,1) + t44 * rSges(3,2);
t49 = (-Icges(4,3) - Icges(5,3)) * t40;
t6 = t51 - t53;
t5 = qJD(1) - t54;
t46 = (t5 * t52 + (-t5 * rSges(5,2) + t6 * t61) * t38) * t41;
t24 = qJD(2) * t30;
t19 = qJD(2) * t31 + qJD(1);
t8 = qJD(1) - t50;
t7 = -qJD(2) * t24 + qJDD(2) * t31 + qJDD(1);
t2 = [m(2) * qJDD(1) + m(3) * t7 + m(4) * t3 + m(5) * t1 + (-m(2) - m(3) - m(4) - m(5)) * g(1); Icges(3,3) * qJDD(2) + t49 + (-t5 * t53 + t46 + t64 * (t15 + t39) + t66 * (t14 - t62) + (-t5 - t54 + t26 - t36) * t6) * m(5) + (t63 * (t23 + t39) + t65 * (-t21 - t62) + (-t50 - t13 - t36) * t9) * m(4) + (-t19 * t24 + (t7 - g(1)) * t31 + (qJD(2) * t19 + qJDD(2) * t30 + t31 * t45 + g(2)) * t30) * m(3); t49 + (-t5 * t51 + t46 + (t26 - t56) * t6 + t64 * t15 + t66 * t14) * m(5) + (-t8 * t57 - t9 * t13 + (t41 * t9 + t63) * t23 + (t41 * t8 - t65) * t21) * m(4); (g(3) + qJDD(4)) * m(5);];
tau  = t2;
