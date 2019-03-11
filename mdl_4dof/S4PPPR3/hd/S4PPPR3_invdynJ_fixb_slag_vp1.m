% Calculate vector of inverse dynamics joint torques for
% S4PPPR3
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
%   pkin=[a2,a3,a4,d4,theta3]';
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
% Datum: 2019-03-08 18:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PPPR3_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR3_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPPR3_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPPR3_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPPR3_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR3_invdynJ_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPPR3_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPPR3_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPPR3_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:11:11
% EndTime: 2019-03-08 18:11:11
% DurationCPUTime: 0.08s
% Computational Cost: add. (95->27), mult. (112->25), div. (0->0), fcn. (40->2), ass. (0->14)
t15 = m(3) + m(4);
t11 = pkin(5) + qJ(4);
t10 = cos(t11);
t9 = sin(t11);
t7 = rSges(5,1) * t9 + rSges(5,2) * t10;
t14 = qJD(4) * t7;
t8 = rSges(5,1) * t10 - rSges(5,2) * t9;
t6 = qJD(4) * t8;
t12 = -m(5) - t15;
t4 = qJD(1) - t14;
t3 = qJD(2) + t6;
t2 = -qJD(4) * t6 - qJDD(4) * t7 + qJDD(1);
t1 = -qJD(4) * t14 + qJDD(4) * t8 + qJDD(2);
t5 = [m(5) * t2 + (m(2) + t15) * qJDD(1) + (-m(2) + t12) * g(2); m(5) * t1 + t12 * g(1) + t15 * qJDD(2); (m(4) + m(5)) * (g(3) + qJDD(3)); Icges(5,3) * qJDD(4) + (-t3 * t14 - t4 * t6 + (qJD(4) * t4 - g(1) + t1) * t8 + (qJD(4) * t3 + g(2) - t2) * t7) * m(5);];
tau  = t5;
