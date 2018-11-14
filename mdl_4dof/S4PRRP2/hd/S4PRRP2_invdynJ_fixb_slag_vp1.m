% Calculate vector of inverse dynamics joint torques for
% S4PRRP2
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
% Datum: 2018-11-14 14:03
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau = S4PRRP2_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP2_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP2_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP2_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP2_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRRP2_invdynJ_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP2_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRP2_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRP2_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:03:19
% EndTime: 2018-11-14 14:03:19
% DurationCPUTime: 0.41s
% Computational Cost: add. (534->89), mult. (573->101), div. (0->0), fcn. (258->4), ass. (0->52)
t44 = sin(qJ(2));
t45 = cos(qJ(2));
t46 = qJD(2) ^ 2;
t69 = (-qJDD(2) * t44 - t46 * t45) * pkin(2) - g(1);
t54 = pkin(2) * qJD(2);
t53 = t44 * t54;
t43 = qJ(2) + qJ(3);
t37 = sin(t43);
t38 = cos(t43);
t21 = t37 * rSges(4,1) + t38 * rSges(4,2);
t42 = qJD(2) + qJD(3);
t56 = t42 * t21;
t9 = -t53 - t56;
t57 = t38 * t42;
t59 = t37 * t42;
t13 = rSges(4,1) * t57 - rSges(4,2) * t59;
t41 = qJDD(2) + qJDD(3);
t68 = -t42 * t13 - t41 * t21 + t69;
t58 = t38 * rSges(5,2);
t20 = t37 * rSges(5,1) + t58;
t22 = t38 * rSges(5,1) - t37 * rSges(5,2);
t40 = t42 ^ 2;
t39 = t45 * pkin(2);
t61 = t44 * pkin(2);
t49 = qJDD(2) * t39 - t46 * t61 + qJDD(1);
t1 = t41 * t22 - t20 * t40 + (-t40 * t37 + t41 * t38) * pkin(3) + t49;
t67 = t1 - g(2);
t26 = rSges(5,2) * t59;
t66 = -t41 * t20 - t42 * (rSges(5,1) * t57 - t26) + (-t41 * t37 - t40 * t38) * pkin(3) + t69;
t23 = t38 * rSges(4,1) - t37 * rSges(4,2);
t3 = t41 * t23 - t42 * t56 + t49;
t65 = t3 - g(2);
t63 = pkin(3) * t38 + t22;
t64 = t42 * t63;
t62 = pkin(3) * t37;
t60 = -rSges(5,1) - pkin(3);
t55 = (Icges(4,3) + Icges(5,3)) * t41;
t52 = t60 * t37;
t36 = t45 * t54;
t50 = -t42 * t23 - t36;
t29 = t45 * rSges(3,1) - t44 * rSges(3,2);
t28 = t44 * rSges(3,1) + t45 * rSges(3,2);
t14 = t52 - t58;
t5 = qJD(1) + t36 + t64;
t6 = -t53 + (-t20 - t62) * t42;
t47 = (t5 * t52 + (-t5 * rSges(5,2) + t6 * t60) * t38) * t42;
t24 = qJD(2) * t28;
t19 = qJD(2) * t29 + qJD(1);
t18 = t42 * t20;
t8 = qJD(1) - t50;
t7 = -qJD(2) * t24 + qJDD(2) * t29 + qJDD(1);
t2 = [m(2) * qJDD(1) + m(3) * t7 + m(4) * t3 + m(5) * t1 + (-m(2) - m(3) - m(4) - m(5)) * g(2); Icges(3,3) * qJDD(2) + t55 + (t47 + t67 * (t63 + t39) + t66 * (t14 - t61) + (t26 + t64) * t6 + (pkin(3) * t59 + t18) * t5) * m(5) + (t65 * (t23 + t39) + t68 * (-t21 - t61) + (-t50 - t13 - t36) * t9) * m(4) + (-t19 * t24 + (t7 - g(2)) * t29 + (qJD(2) * t19 + qJDD(2) * t28 + t29 * t46 + g(1)) * t28) * m(3); t55 + (t6 * t26 + t47 + t5 * t18 - (-t5 * t62 - t6 * t63) * t42 + t67 * t63 + t66 * t14) * m(5) + (-t8 * t56 - t9 * t13 + (t42 * t9 + t65) * t23 + (t42 * t8 - t68) * t21) * m(4); (-g(3) + qJDD(4)) * m(5);];
tau  = t2;
