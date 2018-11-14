% Calculate vector of inverse dynamics joint torques for
% S4PPPR5
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
%   pkin=[a2,a3,a4,d4,theta2]';
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
% Datum: 2018-11-14 14:05
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau = S4PPPR5_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR5_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPPR5_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPPR5_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPPR5_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR5_invdynJ_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PPPR5_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPPR5_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPPR5_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:05:01
% EndTime: 2018-11-14 14:05:01
% DurationCPUTime: 0.13s
% Computational Cost: add. (114->35), mult. (248->46), div. (0->0), fcn. (195->4), ass. (0->20)
t16 = sin(pkin(5));
t17 = cos(pkin(5));
t18 = sin(qJ(4));
t19 = cos(qJ(4));
t12 = -t16 * t18 - t17 * t19;
t13 = t16 * t19 - t17 * t18;
t7 = t12 * rSges(5,1) - t13 * rSges(5,2);
t22 = qJD(4) * t7;
t21 = -m(4) - m(5);
t20 = m(3) - t21;
t15 = -qJDD(3) * t17 + qJDD(1);
t11 = t13 * qJD(4);
t10 = t12 * qJD(4);
t8 = t13 * rSges(5,1) + t12 * rSges(5,2);
t6 = t11 * rSges(5,1) + t10 * rSges(5,2);
t5 = t10 * rSges(5,1) - t11 * rSges(5,2);
t3 = -qJD(3) * t17 + qJD(1) + t22;
t2 = qJD(4) * t5 + qJDD(3) * t16 + qJDD(4) * t8;
t1 = -qJD(4) * t6 + qJDD(4) * t7 + t15;
t4 = [(m(2) + m(3)) * qJDD(1) + m(4) * t15 + m(5) * t1 + (-m(2) - t20) * g(1); (g(3) + qJDD(2)) * t20; t21 * (-g(1) * t17 + g(2) * t16) + m(4) * (qJDD(3) * t16 ^ 2 - t15 * t17) + m(5) * (-t1 * t17 + t2 * t16); Icges(5,3) * qJDD(4) + (-t3 * t6 + (t5 - t22) * (qJD(3) * t16 + qJD(4) * t8) + (qJD(4) * t3 - g(2) + t2) * t8 + (-g(1) + t1) * t7) * m(5);];
tau  = t4;
