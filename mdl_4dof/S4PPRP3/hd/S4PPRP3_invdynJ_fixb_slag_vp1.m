% Calculate vector of inverse dynamics joint torques for
% S4PPRP3
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
% Datum: 2018-11-14 13:58
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau = S4PPRP3_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(4,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP3_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRP3_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRP3_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRP3_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S4PPRP3_invdynJ_fixb_slag_vp1: pkin has to be [4x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRP3_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRP3_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPRP3_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:58:14
% EndTime: 2018-11-14 13:58:15
% DurationCPUTime: 0.20s
% Computational Cost: add. (114->48), mult. (248->51), div. (0->0), fcn. (100->2), ass. (0->20)
t18 = sin(qJ(3));
t19 = cos(qJ(3));
t15 = t19 * rSges(5,1) - t18 * rSges(5,2);
t8 = t19 * pkin(3) + t15;
t20 = qJD(3) ^ 2;
t25 = t19 * rSges(5,2);
t14 = t18 * rSges(4,1) + t19 * rSges(4,2);
t24 = qJD(3) * t14;
t16 = t19 * rSges(4,1) - t18 * rSges(4,2);
t12 = qJD(3) * t16;
t23 = -m(3) - m(4) - m(5);
t13 = t18 * rSges(5,1) + t25;
t7 = -t25 + (-rSges(5,1) - pkin(3)) * t18;
t10 = qJD(1) - t24;
t9 = qJD(2) + t12;
t4 = -qJD(3) * t12 - qJDD(3) * t14 + qJDD(1);
t3 = -qJD(3) * t24 + qJDD(3) * t16 + qJDD(2);
t2 = -qJDD(3) * t13 + qJDD(1) - t15 * t20 + (-qJDD(3) * t18 - t20 * t19) * pkin(3);
t1 = qJDD(3) * t15 + qJDD(2) - t13 * t20 + (qJDD(3) * t19 - t20 * t18) * pkin(3);
t5 = [(m(2) + m(3)) * qJDD(1) + m(4) * t4 + m(5) * t2 + (-m(2) + t23) * g(2); m(3) * qJDD(2) + m(4) * t3 + m(5) * t1 + t23 * g(1); (Icges(4,3) + Icges(5,3)) * qJDD(3) + ((-g(1) + t1) * t8 + (-g(2) + t2) * t7 + (t18 * pkin(3) + t13 + t7) * (t8 * qJD(3) + qJD(2)) * qJD(3)) * m(5) + (-t10 * t12 - t9 * t24 + (qJD(3) * t10 - g(1) + t3) * t16 + (qJD(3) * t9 + g(2) - t4) * t14) * m(4); (g(3) + qJDD(4)) * m(5);];
tau  = t5;
