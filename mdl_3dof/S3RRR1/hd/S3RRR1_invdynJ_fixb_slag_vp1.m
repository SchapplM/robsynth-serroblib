% Calculate vector of inverse dynamics joint torques for
% S3RRR1
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2,d3]';
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S3RRR1_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(5,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRR1_invdynJ_fixb_slag_vp1: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RRR1_invdynJ_fixb_slag_vp1: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3RRR1_invdynJ_fixb_slag_vp1: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RRR1_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S3RRR1_invdynJ_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'S3RRR1_invdynJ_fixb_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'S3RRR1_invdynJ_fixb_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'S3RRR1_invdynJ_fixb_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:07:57
% EndTime: 2019-03-08 18:07:58
% DurationCPUTime: 0.37s
% Computational Cost: add. (735->68), mult. (612->83), div. (0->0), fcn. (288->6), ass. (0->46)
t49 = cos(qJ(1));
t57 = pkin(1) * qJD(1);
t55 = t49 * t57;
t47 = qJ(1) + qJ(2);
t41 = cos(t47);
t46 = qJD(1) + qJD(2);
t60 = t41 * t46;
t42 = qJ(3) + t47;
t34 = sin(t42);
t35 = cos(t42);
t21 = t35 * rSges(4,1) - t34 * rSges(4,2);
t39 = qJD(3) + t46;
t72 = t39 * t21;
t6 = pkin(2) * t60 + t55 + t72;
t48 = sin(qJ(1));
t50 = qJD(1) ^ 2;
t74 = (-qJDD(1) * t48 - t49 * t50) * pkin(1) - g(1);
t43 = t49 * pkin(1);
t65 = t48 * pkin(1);
t73 = qJDD(1) * t43 - t50 * t65 - g(2);
t56 = t48 * t57;
t40 = sin(t47);
t22 = t40 * rSges(3,1) + t41 * rSges(3,2);
t59 = t46 * t22;
t11 = -t56 - t59;
t20 = t34 * rSges(4,1) + t35 * rSges(4,2);
t45 = qJDD(1) + qJDD(2);
t38 = qJDD(3) + t45;
t44 = t46 ^ 2;
t71 = -t39 * t72 - t38 * t20 + (-t40 * t45 - t41 * t44) * pkin(2) + t74;
t61 = t40 * t46;
t17 = rSges(3,1) * t60 - rSges(3,2) * t61;
t70 = -t46 * t17 - t45 * t22 + t74;
t15 = t39 * t20;
t69 = t38 * t21 - t39 * t15 + (-t40 * t44 + t41 * t45) * pkin(2) + t73;
t23 = t41 * rSges(3,1) - t40 * rSges(3,2);
t68 = t45 * t23 - t46 * t59 + t73;
t31 = Icges(4,3) * t38;
t58 = Icges(3,3) * t45 + t31;
t14 = pkin(2) * t41 + t21;
t29 = t49 * rSges(2,1) - t48 * rSges(2,2);
t28 = t48 * rSges(2,1) + t49 * rSges(2,2);
t13 = -pkin(2) * t40 - t20;
t12 = t46 * t23 + t55;
t5 = -pkin(2) * t61 - t15 - t56;
t1 = [Icges(2,3) * qJDD(1) + t58 + (t68 * (t23 + t43) + t70 * (-t22 - t65) + (-t17 - t55 + t12) * t11) * m(3) + (g(1) * t28 - g(2) * t29 + (t28 ^ 2 + t29 ^ 2) * qJDD(1)) * m(2) + (t69 * (t14 + t43) + t71 * (t13 - t65)) * m(4); t58 + (t71 * t13 + t69 * t14) * m(4) + (-t11 * t17 - t12 * t59 + (t11 * t46 + t68) * t23 + (t12 * t46 - t70) * t22) * m(3); t31 + ((t39 * t5 + t69) * t21 + (t39 * t6 - t71) * t20 - t5 * t72 - t6 * t15) * m(4);];
tau  = t1;
