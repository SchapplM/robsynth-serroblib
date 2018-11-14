% Calculate vector of inverse dynamics joint torques for
% S4PPRP4
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3]';
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
% Datum: 2018-11-14 14:06
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau = S4PPRP4_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(4,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP4_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRP4_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRP4_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRP4_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S4PPRP4_invdynJ_fixb_slag_vp1: pkin has to be [4x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRP4_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRP4_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPRP4_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:06:03
% EndTime: 2018-11-14 14:06:03
% DurationCPUTime: 0.21s
% Computational Cost: add. (114->51), mult. (248->55), div. (0->0), fcn. (100->2), ass. (0->24)
t21 = sin(qJ(3));
t22 = cos(qJ(3));
t29 = t22 * rSges(5,2);
t13 = -t21 * rSges(5,1) - t29;
t24 = t21 * rSges(4,1) + t22 * rSges(4,2);
t11 = t24 * qJD(3);
t32 = -rSges(5,1) - pkin(3);
t34 = t32 * t22;
t23 = qJD(3) ^ 2;
t30 = t22 * rSges(5,1);
t16 = t22 * rSges(4,1) - t21 * rSges(4,2);
t28 = qJD(3) * t16;
t27 = qJD(3) * t21;
t26 = m(3) + m(4) + m(5);
t19 = t21 * rSges(5,2);
t17 = rSges(5,2) * t27;
t15 = -t19 + t30;
t10 = qJD(1) - t11;
t9 = -qJD(2) - t28;
t4 = -qJD(3) * t28 - qJDD(3) * t24 + qJDD(1);
t3 = -qJD(3) * t11 + qJDD(3) * t16 + qJDD(2);
t2 = qJDD(1) + qJDD(3) * t13 + qJD(3) * (-qJD(3) * t30 + t17) + (-qJDD(3) * t21 - t23 * t22) * pkin(3);
t1 = qJDD(3) * t15 + qJDD(2) + t13 * t23 + (qJDD(3) * t22 - t23 * t21) * pkin(3);
t5 = [(m(2) + m(3)) * qJDD(1) + m(4) * t4 + m(5) * t2 + (-m(2) - t26) * g(1); m(3) * qJDD(2) + m(4) * t3 + m(5) * t1 + t26 * g(2); (Icges(4,3) + Icges(5,3)) * qJDD(3) + (-t10 * t28 + t9 * t11 + (t10 * qJD(3) + g(2) + t3) * t16 - (t9 * qJD(3) - g(1) + t4) * t24) * m(4) + ((-g(2) - t1) * (t19 + t34) + (-g(1) + t2) * (t32 * t21 - t29) + (t17 + (t22 * pkin(3) + t15 + t34) * qJD(3)) * (-pkin(3) * t27 + qJD(3) * t13 + qJD(1))) * m(5); (g(3) + qJDD(4)) * m(5);];
tau  = t5;
