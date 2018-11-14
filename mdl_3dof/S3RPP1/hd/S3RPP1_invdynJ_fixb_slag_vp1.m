% Calculate vector of inverse dynamics joint torques for
% S3RPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% qJDD [3x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1]';
% m_mdh [4x1]
%   mass of all robot links (including the base)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [4x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tau [3x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:14
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau = S3RPP1_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPP1_invdynJ_fixb_slag_vp1: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RPP1_invdynJ_fixb_slag_vp1: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3RPP1_invdynJ_fixb_slag_vp1: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RPP1_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3RPP1_invdynJ_fixb_slag_vp1: pkin has to be [3x1] (double)');
assert( isreal(m) && all(size(m) == [4 1]), ...
  'S3RPP1_invdynJ_fixb_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'S3RPP1_invdynJ_fixb_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'S3RPP1_invdynJ_fixb_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:13:33
% EndTime: 2018-11-14 10:13:34
% DurationCPUTime: 0.38s
% Computational Cost: add. (362->97), mult. (721->113), div. (0->0), fcn. (398->2), ass. (0->52)
t50 = sin(qJ(1));
t51 = cos(qJ(1));
t28 = t51 * pkin(1) + t50 * qJ(2);
t42 = qJD(2) * t51;
t19 = qJD(1) * t28 - t42;
t25 = t50 * rSges(4,2) + t51 * rSges(4,3);
t44 = t51 * qJ(2);
t24 = t50 * pkin(1) - t44;
t29 = t51 * rSges(4,2) - t50 * rSges(4,3);
t56 = -t50 * qJ(3) + t29;
t53 = -t24 + t56;
t52 = qJD(1) ^ 2;
t54 = -t52 * qJ(3) + qJDD(3);
t60 = qJD(3) * t50;
t59 = qJD(1) * qJD(2);
t68 = qJDD(2) * t50 + t51 * t59;
t1 = t54 * t51 + t53 * qJDD(1) + (-qJD(1) * t25 - t19 - 0.2e1 * t60) * qJD(1) + t68;
t79 = -g(1) + t1;
t61 = qJD(1) * t51;
t38 = rSges(4,2) * t61;
t62 = qJD(1) * t50;
t35 = qJ(2) * t61;
t41 = qJD(2) * t50;
t67 = t35 + t41;
t57 = qJD(1) * (-pkin(1) * t62 + t67) + qJDD(1) * t28 + t50 * t59;
t2 = qJD(1) * t38 + qJDD(1) * t25 + (-t52 * rSges(4,3) + t54) * t50 + (qJDD(1) * qJ(3) + 0.2e1 * qJD(3) * qJD(1) - qJDD(2)) * t51 + t57;
t78 = -g(2) + t2;
t77 = t51 * rSges(3,2) - t50 * rSges(3,3);
t76 = t51 * qJ(3) + t25 + t28;
t75 = t50 / 0.2e1;
t74 = -t51 / 0.2e1;
t73 = rSges(3,2) - pkin(1);
t71 = t51 * rSges(3,3);
t70 = -pkin(1) - qJ(3);
t26 = t50 * rSges(3,2) + t71;
t69 = -t24 + t26;
t18 = t28 - t77;
t66 = rSges(3,2) * t62 + rSges(3,3) * t61;
t40 = qJD(3) * t51;
t65 = t40 + t41;
t64 = -qJD(1) * t24 + t41;
t58 = -rSges(4,3) + t70;
t55 = t42 - t60;
t31 = t51 * rSges(2,1) - t50 * rSges(2,2);
t27 = t50 * rSges(2,1) + t51 * rSges(2,2);
t13 = t18 * qJD(1) - t42;
t12 = t69 * qJD(1) + t41;
t9 = qJD(1) * t76 - t55;
t8 = t53 * qJD(1) + t65;
t4 = qJD(1) * t66 - qJDD(1) * t77 - qJDD(2) * t51 + t57;
t3 = t69 * qJDD(1) + (qJD(1) * t77 - t19) * qJD(1) + t68;
t5 = [-m(2) * (-g(1) * t27 + g(2) * t31) + (-(t56 * qJD(1) + t40 + t64 - t8) * t9 + t8 * t55 + t9 * (t35 + t38 + t65) + (t8 * t58 * t51 + (t8 * (-rSges(4,2) - qJ(2)) + t9 * t58) * t50) * qJD(1) + t78 * t76 + t79 * (t70 * t50 + t29 + t44)) * m(4) + (-(qJD(1) * t26 - t12 + t64) * t13 + t12 * t42 + t13 * (t66 + t67) + (t12 * t73 * t51 + (t12 * (-rSges(3,3) - qJ(2)) - t13 * pkin(1)) * t50) * qJD(1) + (-g(2) + t4) * t18 + (-g(1) + t3) * (t73 * t50 + t44 + t71)) * m(3) + (Icges(2,3) + Icges(3,1) + Icges(4,1) + m(2) * (t27 ^ 2 + t31 ^ 2)) * qJDD(1); (-m(3) - m(4)) * (g(1) * t50 - g(2) * t51) + 0.2e1 * (t1 * t75 + t2 * t74) * m(4) + 0.2e1 * (t3 * t75 + t4 * t74) * m(3); (t78 * t50 + t79 * t51) * m(4);];
tau  = t5;
