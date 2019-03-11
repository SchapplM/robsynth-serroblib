% Calculate vector of inverse dynamics joint torques for
% S4PPRR2
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
%   pkin=[a2,a3,a4,d3,d4,theta2]';
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
% Datum: 2019-03-08 18:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PPRR2_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR2_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR2_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR2_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR2_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR2_invdynJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR2_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRR2_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPRR2_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:16:49
% EndTime: 2019-03-08 18:16:50
% DurationCPUTime: 0.20s
% Computational Cost: add. (374->48), mult. (303->54), div. (0->0), fcn. (130->4), ass. (0->27)
t26 = pkin(6) + qJ(3);
t24 = qJ(4) + t26;
t19 = sin(t24);
t20 = cos(t24);
t14 = t20 * rSges(5,1) - t19 * rSges(5,2);
t27 = qJD(3) + qJD(4);
t7 = t27 * t14;
t22 = sin(t26);
t33 = pkin(3) * qJD(3);
t13 = t19 * rSges(5,1) + t20 * rSges(5,2);
t34 = t27 * t13;
t5 = -t22 * t33 - t34;
t23 = cos(t26);
t25 = qJDD(3) + qJDD(4);
t28 = qJD(3) ^ 2;
t1 = t25 * t14 - t27 * t34 + qJDD(1) + (qJDD(3) * t23 - t22 * t28) * pkin(3);
t37 = t1 - g(2);
t36 = -t25 * t13 - t27 * t7 + (-qJDD(3) * t22 - t23 * t28) * pkin(3) - g(1);
t32 = m(3) + m(4) + m(5);
t16 = t23 * rSges(4,1) - t22 * rSges(4,2);
t15 = t22 * rSges(4,1) + t23 * rSges(4,2);
t21 = Icges(5,3) * t25;
t11 = qJD(3) * t15;
t10 = qJD(3) * t16 + qJD(1);
t4 = t23 * t33 + qJD(1) + t7;
t3 = -qJD(3) * t11 + qJDD(3) * t16 + qJDD(1);
t2 = [(m(2) + m(3)) * qJDD(1) + m(4) * t3 + m(5) * t1 + (-m(2) - t32) * g(2); (-g(3) + qJDD(2)) * t32; Icges(4,3) * qJDD(3) + t21 + (t37 * (pkin(3) * t23 + t14) + t36 * (-pkin(3) * t22 - t13)) * m(5) + (-t10 * t11 + (-g(2) + t3) * t16 + (qJD(3) * t10 + qJDD(3) * t15 + t16 * t28 + g(1)) * t15) * m(4); t21 + (-t4 * t34 - t5 * t7 + (t27 * t5 + t37) * t14 + (t27 * t4 - t36) * t13) * m(5);];
tau  = t2;
