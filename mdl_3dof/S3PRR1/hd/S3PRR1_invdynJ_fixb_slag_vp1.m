% Calculate vector of inverse dynamics joint torques for
% S3PRR1
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2,d3]';
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
% Datum: 2019-03-08 18:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S3PRR1_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRR1_invdynJ_fixb_slag_vp1: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PRR1_invdynJ_fixb_slag_vp1: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3PRR1_invdynJ_fixb_slag_vp1: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3PRR1_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3PRR1_invdynJ_fixb_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'S3PRR1_invdynJ_fixb_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'S3PRR1_invdynJ_fixb_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'S3PRR1_invdynJ_fixb_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:03:50
% EndTime: 2019-03-08 18:03:50
% DurationCPUTime: 0.23s
% Computational Cost: add. (237->45), mult. (295->53), div. (0->0), fcn. (130->4), ass. (0->25)
t24 = qJ(2) + qJ(3);
t20 = sin(t24);
t21 = cos(t24);
t12 = t21 * rSges(4,1) - t20 * rSges(4,2);
t23 = qJD(2) + qJD(3);
t7 = t23 * t12;
t25 = sin(qJ(2));
t31 = pkin(2) * qJD(2);
t11 = t20 * rSges(4,1) + t21 * rSges(4,2);
t32 = t23 * t11;
t5 = -t25 * t31 - t32;
t22 = qJDD(2) + qJDD(3);
t26 = cos(qJ(2));
t27 = qJD(2) ^ 2;
t1 = t22 * t12 - t23 * t32 + qJDD(1) + (qJDD(2) * t26 - t25 * t27) * pkin(2);
t35 = t1 - g(2);
t34 = -t22 * t11 - t23 * t7 + (-qJDD(2) * t25 - t26 * t27) * pkin(2) - g(1);
t17 = t26 * rSges(3,1) - t25 * rSges(3,2);
t16 = t25 * rSges(3,1) + t26 * rSges(3,2);
t19 = Icges(4,3) * t22;
t13 = qJD(2) * t16;
t10 = qJD(2) * t17 + qJD(1);
t4 = t26 * t31 + qJD(1) + t7;
t3 = -qJD(2) * t13 + qJDD(2) * t17 + qJDD(1);
t2 = [m(2) * qJDD(1) + m(3) * t3 + m(4) * t1 + (-m(2) - m(3) - m(4)) * g(2); Icges(3,3) * qJDD(2) + t19 + (t35 * (t26 * pkin(2) + t12) + t34 * (-t25 * pkin(2) - t11)) * m(4) + (-t10 * t13 + (t3 - g(2)) * t17 + (qJD(2) * t10 + qJDD(2) * t16 + t17 * t27 + g(1)) * t16) * m(3); t19 + (-t4 * t32 - t5 * t7 + (t23 * t5 + t35) * t12 + (t23 * t4 - t34) * t11) * m(4);];
tau  = t2;
