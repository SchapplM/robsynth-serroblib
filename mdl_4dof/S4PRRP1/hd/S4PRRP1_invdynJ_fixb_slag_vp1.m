% Calculate vector of inverse dynamics joint torques for
% S4PRRP1
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
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PRRP1_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP1_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP1_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP1_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP1_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP1_invdynJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP1_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRP1_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRP1_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:22:49
% EndTime: 2019-03-08 18:22:49
% DurationCPUTime: 0.47s
% Computational Cost: add. (1105->81), mult. (695->88), div. (0->0), fcn. (358->4), ass. (0->49)
t53 = pkin(6) + qJ(2);
t51 = qJ(3) + t53;
t45 = sin(t51);
t34 = qJD(4) * t45;
t54 = qJD(2) + qJD(3);
t46 = cos(t51);
t38 = t46 * qJ(4);
t84 = rSges(5,1) + pkin(3);
t85 = -t46 * rSges(5,3) + t45 * t84 - t38;
t76 = t85 * t54;
t88 = t34 - t76;
t87 = t84 * t46;
t50 = cos(t53);
t44 = pkin(2) * t50;
t55 = qJD(2) ^ 2;
t49 = sin(t53);
t74 = pkin(2) * t49;
t83 = -qJDD(2) * t44 + t55 * t74 + g(2);
t82 = (-qJDD(2) * t49 - t55 * t50) * pkin(2) - g(1);
t67 = pkin(2) * qJD(2);
t64 = t49 * t67;
t25 = t45 * rSges(4,1) + t46 * rSges(4,2);
t70 = t54 * t25;
t13 = -t64 - t70;
t81 = rSges(5,3) + qJ(4);
t71 = t46 * t54;
t72 = t45 * t54;
t16 = rSges(4,1) * t71 - rSges(4,2) * t72;
t52 = qJDD(2) + qJDD(3);
t80 = -t54 * t16 - t52 * t25 + t82;
t28 = t46 * rSges(4,1) - t45 * rSges(4,2);
t79 = t52 * t28 - t54 * t70 - t83;
t12 = t81 * t45 + t87;
t58 = -t54 * t84 + qJD(4);
t75 = rSges(5,3) * t71 + t54 * t38;
t65 = t34 + t75;
t78 = qJDD(4) * t46 - t12 * t52 - (t45 * t58 + t65) * t54 + t83;
t35 = qJD(4) * t46;
t77 = qJDD(4) * t45 - t85 * t52 + (t46 * t58 - t72 * t81 + t35) * t54 + t82;
t68 = (Icges(5,2) + Icges(4,3)) * t52;
t63 = t50 * t67;
t59 = t35 - t63;
t30 = t50 * rSges(3,1) - t49 * rSges(3,2);
t29 = t49 * rSges(3,1) + t50 * rSges(3,2);
t7 = -t64 + t88;
t8 = t12 * t54 - t59;
t56 = (-t7 * t87 + (-t7 * t81 - t8 * t84) * t45) * t54;
t14 = t54 * t28 + t63;
t1 = [(-g(3) + qJDD(1)) * (m(2) + m(3) + m(4) + m(5)); Icges(3,3) * qJDD(2) + t68 + (t79 * (t28 + t44) + t80 * (-t25 - t74) + (-t16 - t63 + t14) * t13) * m(4) + (g(1) * t29 - g(2) * t30 + (t29 ^ 2 + t30 ^ 2) * qJDD(2)) * m(3) + (t7 * t59 + t56 + t77 * (-t85 - t74) - t78 * (t44 + t12) + (t7 + t75 + t76) * t8) * m(5); t68 + (t56 + (t65 - t88) * t8 - t77 * t85 + (t54 * t7 - t78) * t12) * m(5) + (-t13 * t16 - t14 * t70 + (t13 * t54 + t79) * t28 + (t14 * t54 - t80) * t25) * m(4); (t77 * t45 + t78 * t46) * m(5);];
tau  = t1;
