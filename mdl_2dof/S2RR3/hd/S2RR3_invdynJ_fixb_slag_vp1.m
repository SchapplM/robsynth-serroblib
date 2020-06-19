% Calculate vector of inverse dynamics joint torques for
% S2RR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% qJDD [2x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,d1,d2]';
% m [3x1]
%   mass of all robot links (including the base)
% rSges [3x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [3x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tau [2x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 09:14
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S2RR3_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(2,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR3_invdynJ_fixb_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR3_invdynJ_fixb_slag_vp1: qJD has to be [2x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [2 1]), ...
  'S2RR3_invdynJ_fixb_slag_vp1: qJDD has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S2RR3_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S2RR3_invdynJ_fixb_slag_vp1: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'S2RR3_invdynJ_fixb_slag_vp1: m has to be [3x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [3,3]), ...
  'S2RR3_invdynJ_fixb_slag_vp1: rSges has to be [3x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [3 6]), ...
  'S2RR3_invdynJ_fixb_slag_vp1: Icges has to be [3x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-19 09:14:25
% EndTime: 2020-06-19 09:14:25
% DurationCPUTime: 0.36s
% Computational Cost: add. (205->32), mult. (258->44), div. (0->0), fcn. (120->4), ass. (0->21)
t24 = cos(qJ(1));
t28 = pkin(1) * qJD(1);
t22 = qJ(1) + qJ(2);
t18 = sin(t22);
t19 = cos(t22);
t10 = t19 * rSges(3,1) - t18 * rSges(3,2);
t21 = qJD(1) + qJD(2);
t6 = t21 * t10;
t33 = -t24 * t28 - t6;
t23 = sin(qJ(1));
t9 = t18 * rSges(3,1) + t19 * rSges(3,2);
t30 = t21 * t9;
t3 = -t23 * t28 - t30;
t20 = qJDD(1) + qJDD(2);
t25 = qJD(1) ^ 2;
t32 = -t20 * t9 - t21 * t6 + (-qJDD(1) * t23 - t24 * t25) * pkin(1) - g(1);
t31 = t20 * t10 - t21 * t30 + (qJDD(1) * t24 - t23 * t25) * pkin(1) - g(2);
t15 = t24 * rSges(2,1) - t23 * rSges(2,2);
t14 = t23 * rSges(2,1) + t24 * rSges(2,2);
t17 = Icges(3,3) * t20;
t1 = [Icges(2,3) * qJDD(1) + t17 + (t31 * (t24 * pkin(1) + t10) + t32 * (-t23 * pkin(1) - t9)) * m(3) + ((t14 ^ 2 + t15 ^ 2) * qJDD(1) + g(1) * t14 - g(2) * t15) * m(2); t17 + (-t3 * t6 + t33 * t30 + (-t21 * t33 - t32) * t9 + (t21 * t3 + t31) * t10) * m(3);];
tau = t1;
