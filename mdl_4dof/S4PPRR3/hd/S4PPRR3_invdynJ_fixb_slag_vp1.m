% Calculate vector of inverse dynamics joint torques for
% S4PPRR3
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:08
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau = S4PPRR3_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR3_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR3_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR3_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR3_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR3_invdynJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR3_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRR3_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPRR3_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:08:15
% EndTime: 2018-11-14 14:08:16
% DurationCPUTime: 0.21s
% Computational Cost: add. (374->48), mult. (303->54), div. (0->0), fcn. (130->4), ass. (0->27)
t25 = pkin(6) + qJ(3);
t23 = qJ(4) + t25;
t19 = sin(t23);
t20 = cos(t23);
t14 = t20 * rSges(5,1) - t19 * rSges(5,2);
t26 = qJD(3) + qJD(4);
t7 = t26 * t14;
t21 = sin(t25);
t32 = pkin(3) * qJD(3);
t13 = t19 * rSges(5,1) + t20 * rSges(5,2);
t34 = t26 * t13;
t5 = -t21 * t32 - t34;
t22 = cos(t25);
t24 = -qJDD(3) - qJDD(4);
t27 = qJD(3) ^ 2;
t1 = -t24 * t14 - t26 * t34 + qJDD(1) + (qJDD(3) * t22 - t21 * t27) * pkin(3);
t37 = -g(1) + t1;
t36 = -g(2) + t24 * t13 - t26 * t7 + (-qJDD(3) * t21 - t22 * t27) * pkin(3);
t33 = Icges(5,3) * t24;
t31 = m(3) + m(4) + m(5);
t16 = t22 * rSges(4,1) - t21 * rSges(4,2);
t15 = t21 * rSges(4,1) + t22 * rSges(4,2);
t11 = qJD(3) * t15;
t10 = qJD(3) * t16 + qJD(1);
t4 = t22 * t32 + qJD(1) + t7;
t3 = -qJD(3) * t11 + qJDD(3) * t16 + qJDD(1);
t2 = [(m(2) + m(3)) * qJDD(1) + m(4) * t3 + m(5) * t1 + (-m(2) - t31) * g(1); (g(3) + qJDD(2)) * t31; Icges(4,3) * qJDD(3) - t33 + (t37 * (pkin(3) * t22 + t14) + t36 * (-pkin(3) * t21 - t13)) * m(5) + (-t10 * t11 + (-g(1) + t3) * t16 + (qJD(3) * t10 + qJDD(3) * t15 + t16 * t27 + g(2)) * t15) * m(4); -t33 + (-t4 * t34 - t5 * t7 + (t26 * t5 + t37) * t14 + (t26 * t4 - t36) * t13) * m(5);];
tau  = t2;
