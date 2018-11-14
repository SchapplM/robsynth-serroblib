% Calculate vector of inverse dynamics joint torques for
% S4PRPR3
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:11
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau = S4PRPR3_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR3_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR3_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR3_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR3_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR3_invdynJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR3_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPR3_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPR3_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:11:15
% EndTime: 2018-11-14 14:11:15
% DurationCPUTime: 0.45s
% Computational Cost: add. (450->67), mult. (475->72), div. (0->0), fcn. (210->6), ass. (0->37)
t37 = qJ(2) + pkin(6);
t32 = cos(t37);
t39 = cos(qJ(2));
t34 = t39 * pkin(2);
t64 = pkin(3) * t32 + t34;
t40 = qJD(2) ^ 2;
t38 = sin(qJ(2));
t63 = (-qJDD(2) * t38 - t40 * t39) * pkin(2) - g(2);
t33 = qJ(4) + t37;
t27 = sin(t33);
t28 = cos(t33);
t16 = t28 * rSges(5,1) - t27 * rSges(5,2);
t36 = qJD(2) + qJD(4);
t11 = t36 * t16;
t15 = t27 * rSges(5,1) + t28 * rSges(5,2);
t31 = sin(t37);
t35 = -qJDD(2) - qJDD(4);
t62 = -t36 * t11 + t35 * t15 + (-qJDD(2) * t31 - t32 * t40) * pkin(3) + t63;
t55 = t38 * pkin(2);
t43 = qJDD(2) * t34 - t40 * t55 + qJDD(1);
t51 = t36 * t15;
t1 = -t36 * t51 - t35 * t16 + (qJDD(2) * t32 - t31 * t40) * pkin(3) + t43;
t61 = t1 - g(1);
t45 = -pkin(3) * t31 - t55;
t5 = t45 * qJD(2) - t51;
t60 = m(4) + m(5);
t18 = t32 * rSges(4,1) - t31 * rSges(4,2);
t49 = Icges(5,3) * t35;
t24 = t39 * rSges(3,1) - t38 * rSges(3,2);
t23 = t38 * rSges(3,1) + t39 * rSges(3,2);
t17 = t31 * rSges(4,1) + t32 * rSges(4,2);
t20 = qJD(2) * t23;
t19 = qJD(2) * t24 + qJD(1);
t6 = -qJD(2) * t20 + qJDD(2) * t24 + qJDD(1);
t4 = t64 * qJD(2) + qJD(1) + t11;
t3 = qJDD(2) * t18 - t17 * t40 + t43;
t2 = [m(2) * qJDD(1) + m(3) * t6 + m(4) * t3 + m(5) * t1 + (-m(2) - m(3) - t60) * g(1); -t49 + (Icges(3,3) + Icges(4,3)) * qJDD(2) + (t61 * (t16 + t64) + t62 * (-t15 + t45)) * m(5) + (-t19 * t20 + (t6 - g(1)) * t24 + (qJD(2) * t19 + qJDD(2) * t23 + t24 * t40 + g(2)) * t23) * m(3) + ((t3 - g(1)) * (t18 + t34) + (-qJDD(2) * t17 - t40 * t18 + t63) * (-t17 - t55)) * m(4); t60 * (g(3) + qJDD(3)); -t49 + (-t4 * t51 - t5 * t11 + (t36 * t5 + t61) * t16 + (t36 * t4 - t62) * t15) * m(5);];
tau  = t2;
