% Calculate vector of inverse dynamics joint torques for
% S4PRPP2
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
%   pkin=[a2,a3,a4,d2,theta3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:00
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau = S4PRPP2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP2_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPP2_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPP2_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPP2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP2_invdynJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPP2_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPP2_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPP2_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:00:20
% EndTime: 2018-11-14 14:00:21
% DurationCPUTime: 0.32s
% Computational Cost: add. (222->86), mult. (448->105), div. (0->0), fcn. (260->6), ass. (0->41)
t44 = m(4) + m(5);
t29 = sin(pkin(5));
t30 = cos(pkin(5));
t31 = sin(qJ(2));
t32 = cos(qJ(2));
t16 = t29 * t31 - t30 * t32;
t13 = t16 * qJD(1);
t46 = qJD(4) + t13;
t45 = t44 * pkin(2) + mrSges(3,1);
t37 = qJD(1) * qJD(2);
t19 = t32 * qJDD(1) - t31 * t37;
t15 = qJDD(2) * pkin(2) + t19;
t20 = t31 * qJDD(1) + t32 * t37;
t4 = t29 * t15 + t30 * t20;
t43 = t31 * mrSges(3,2);
t42 = t32 * mrSges(3,2);
t41 = mrSges(4,1) + mrSges(5,1);
t40 = -mrSges(4,2) + mrSges(5,3);
t39 = qJD(1) * t31;
t36 = m(5) * pkin(3) + t41;
t35 = -m(5) * qJ(4) - t40;
t21 = qJD(2) * pkin(2) + t32 * qJD(1);
t8 = t29 * t21 + t30 * t39;
t34 = t31 * mrSges(3,1) + t42;
t3 = t30 * t15 - t29 * t20;
t17 = t29 * t32 + t30 * t31;
t7 = t30 * t21 - t29 * t39;
t33 = qJD(2) ^ 2;
t28 = qJ(2) + pkin(5);
t26 = cos(t28);
t25 = sin(t28);
t24 = -t30 * pkin(2) - pkin(3);
t22 = t29 * pkin(2) + qJ(4);
t12 = t17 * qJD(2);
t11 = t16 * qJD(2);
t10 = t17 * qJD(1);
t6 = qJD(2) * qJ(4) + t8;
t5 = -qJD(2) * pkin(3) + qJD(4) - t7;
t2 = -qJDD(2) * pkin(3) + qJDD(4) - t3;
t1 = qJDD(2) * qJ(4) + qJD(2) * qJD(4) + t4;
t9 = [m(2) * qJDD(1) - t34 * t33 + m(3) * (t19 * t32 + t20 * t31) + m(4) * (-t8 * t11 - t7 * t12 - t3 * t16 + t4 * t17) + m(5) * (t1 * t17 - t6 * t11 + t5 * t12 + t2 * t16) + (-t40 * t11 - t41 * t12) * qJD(2) + (-m(2) - m(3) - t44) * g(2) + (t32 * mrSges(3,1) - t41 * t16 + t40 * t17 - t43) * qJDD(2); t19 * mrSges(3,1) + t3 * mrSges(4,1) - t2 * mrSges(5,1) - t20 * mrSges(3,2) - t4 * mrSges(4,2) + t1 * mrSges(5,3) + (t35 * t25 - t36 * t26 - t45 * t32 + t43) * g(2) + (t36 * t25 + t35 * t26 + t45 * t31 + t42) * g(1) + (-t24 * mrSges(5,1) + t22 * mrSges(5,3) + Ifges(5,2) + Ifges(3,3) + Ifges(4,3) + (mrSges(4,1) * t30 - mrSges(4,2) * t29) * pkin(2)) * qJDD(2) + (-t13 * mrSges(4,2) + t46 * mrSges(5,3) + t34 * qJD(1) + t41 * t10) * qJD(2) + (t1 * t22 - t5 * t10 + t2 * t24 + t46 * t6) * m(5) + (t7 * t10 + t8 * t13 + (t29 * t4 + t3 * t30) * pkin(2)) * m(4); t44 * (-g(3) + qJDD(3)); -qJDD(2) * mrSges(5,1) - t33 * mrSges(5,3) + (-g(1) * t25 + g(2) * t26 - t6 * qJD(2) + t2) * m(5);];
tau  = t9;
