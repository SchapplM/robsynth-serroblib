% Calculate vector of inverse dynamics joint torques for
% S4PRPR2
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
%   pkin=[a2,a3,a4,d2,d4,theta3]';
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
% Datum: 2019-03-08 18:22
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PRPR2_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR2_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR2_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR2_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR2_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR2_invdynJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR2_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPR2_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPR2_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:21:48
% EndTime: 2019-03-08 18:21:49
% DurationCPUTime: 0.42s
% Computational Cost: add. (450->67), mult. (475->72), div. (0->0), fcn. (210->6), ass. (0->37)
t38 = qJ(2) + pkin(6);
t33 = cos(t38);
t40 = cos(qJ(2));
t35 = t40 * pkin(2);
t64 = pkin(3) * t33 + t35;
t41 = qJD(2) ^ 2;
t39 = sin(qJ(2));
t63 = (-qJDD(2) * t39 - t41 * t40) * pkin(2) - g(1);
t34 = qJ(4) + t38;
t27 = sin(t34);
t28 = cos(t34);
t16 = t28 * rSges(5,1) - t27 * rSges(5,2);
t37 = qJD(2) + qJD(4);
t11 = t37 * t16;
t15 = t27 * rSges(5,1) + t28 * rSges(5,2);
t32 = sin(t38);
t36 = qJDD(2) + qJDD(4);
t62 = -t37 * t11 - t36 * t15 + (-qJDD(2) * t32 - t33 * t41) * pkin(3) + t63;
t55 = t39 * pkin(2);
t44 = qJDD(2) * t35 - t41 * t55 + qJDD(1);
t51 = t37 * t15;
t1 = -t37 * t51 + t36 * t16 + (qJDD(2) * t33 - t32 * t41) * pkin(3) + t44;
t61 = t1 - g(2);
t46 = -pkin(3) * t32 - t55;
t5 = t46 * qJD(2) - t51;
t60 = m(4) + m(5);
t18 = t33 * rSges(4,1) - t32 * rSges(4,2);
t24 = t40 * rSges(3,1) - t39 * rSges(3,2);
t23 = t39 * rSges(3,1) + t40 * rSges(3,2);
t17 = t32 * rSges(4,1) + t33 * rSges(4,2);
t30 = Icges(5,3) * t36;
t20 = qJD(2) * t23;
t19 = qJD(2) * t24 + qJD(1);
t6 = -qJD(2) * t20 + qJDD(2) * t24 + qJDD(1);
t4 = t64 * qJD(2) + qJD(1) + t11;
t3 = qJDD(2) * t18 - t17 * t41 + t44;
t2 = [m(2) * qJDD(1) + m(3) * t6 + m(4) * t3 + m(5) * t1 + (-m(2) - m(3) - t60) * g(2); t30 + (Icges(3,3) + Icges(4,3)) * qJDD(2) + (t61 * (t16 + t64) + t62 * (-t15 + t46)) * m(5) + (-t19 * t20 + (t6 - g(2)) * t24 + (qJD(2) * t19 + qJDD(2) * t23 + t24 * t41 + g(1)) * t23) * m(3) + ((t3 - g(2)) * (t18 + t35) + (-qJDD(2) * t17 - t41 * t18 + t63) * (-t17 - t55)) * m(4); t60 * (-g(3) + qJDD(3)); t30 + (-t4 * t51 - t5 * t11 + (t37 * t5 + t61) * t16 + (t37 * t4 - t62) * t15) * m(5);];
tau  = t2;
