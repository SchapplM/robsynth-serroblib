% Calculate vector of inverse dynamics joint torques for
% S4RPPR1
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
%   pkin=[a2,a3,a4,d1,d4,theta2]';
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
% Datum: 2018-11-14 13:47
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau = S4RPPR1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR1_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR1_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR1_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR1_invdynJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR1_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR1_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR1_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:46:34
% EndTime: 2018-11-14 13:46:35
% DurationCPUTime: 0.37s
% Computational Cost: add. (359->99), mult. (548->116), div. (0->0), fcn. (239->8), ass. (0->40)
t48 = m(3) * pkin(1);
t47 = -m(5) - m(4);
t28 = -qJDD(1) + qJDD(4);
t46 = Ifges(5,3) * t28;
t45 = mrSges(3,1) + mrSges(4,1);
t44 = -mrSges(4,3) + mrSges(3,2);
t29 = -qJD(1) + qJD(4);
t43 = qJD(4) * t29;
t42 = m(3) - t47;
t41 = qJD(1) * qJD(3);
t32 = cos(pkin(6));
t24 = -t32 * pkin(1) - pkin(2);
t31 = sin(pkin(6));
t22 = t31 * pkin(1) + qJ(3);
t19 = -pkin(3) + t24;
t30 = qJ(1) + pkin(6);
t25 = sin(t30);
t26 = cos(t30);
t39 = -g(1) * t25 + g(2) * t26;
t12 = qJD(1) * t19 + qJD(3);
t15 = t22 * qJD(1);
t33 = sin(qJ(4));
t35 = cos(qJ(4));
t3 = t35 * t12 - t33 * t15;
t4 = t33 * t12 + t35 * t15;
t38 = -t3 * t33 + t35 * t4;
t7 = t35 * t19 - t33 * t22;
t8 = t33 * t19 + t35 * t22;
t36 = cos(qJ(1));
t34 = sin(qJ(1));
t14 = qJDD(1) * t24 + qJDD(3);
t13 = qJDD(1) * t22 + t41;
t11 = qJDD(1) * t19 + qJDD(3);
t10 = -t25 * t35 + t26 * t33;
t9 = -t25 * t33 - t26 * t35;
t6 = -t33 * qJD(3) - qJD(4) * t8;
t5 = t35 * qJD(3) + qJD(4) * t7;
t2 = -qJD(4) * t4 + t35 * t11 - t33 * t13;
t1 = qJD(4) * t3 + t33 * t11 + t35 * t13;
t16 = [-t14 * mrSges(4,1) - t46 + (t13 + t41) * mrSges(4,3) + m(5) * (t1 * t8 + t2 * t7 + t3 * t6 + t4 * t5) + m(4) * (t15 * qJD(3) + t13 * t22 + t14 * t24) + (-t8 * t28 - t5 * t29 + t1) * mrSges(5,2) + (t7 * t28 + t6 * t29 - t2) * mrSges(5,1) + (t9 * mrSges(5,1) + t10 * mrSges(5,2) + t34 * mrSges(2,2) + (-mrSges(2,1) - t48) * t36 + (-m(5) * pkin(3) - t45) * t26 + t44 * t25 + t47 * (t36 * pkin(1) + t26 * pkin(2) + t25 * qJ(3))) * g(2) + (-t10 * mrSges(5,1) + t36 * mrSges(2,2) + t9 * mrSges(5,2) + (pkin(1) * t42 + mrSges(2,1)) * t34 + (-m(5) * (-pkin(2) - pkin(3)) + m(4) * pkin(2) + t45) * t25 + (qJ(3) * t47 + t44) * t26) * g(1) + (-t24 * mrSges(4,1) + t22 * mrSges(4,3) + Ifges(4,2) + Ifges(2,3) + Ifges(3,3) + (0.2e1 * t32 * mrSges(3,1) - 0.2e1 * t31 * mrSges(3,2) + (t31 ^ 2 + t32 ^ 2) * t48) * pkin(1)) * qJDD(1); (-g(3) + qJDD(2)) * t42; -qJDD(1) * mrSges(4,1) + (-t33 * t28 - t35 * t43) * mrSges(5,2) + (t35 * t28 - t33 * t43) * mrSges(5,1) + (qJD(4) * t38 + t1 * t33 + t2 * t35 + t39) * m(5) + (t14 + t39) * m(4) + (-m(4) * t15 - m(5) * t38 - mrSges(4,3) * qJD(1) + (t33 * mrSges(5,1) + t35 * mrSges(5,2)) * t29) * qJD(1); t46 + (-g(1) * t9 - g(2) * t10 + t3 * t29 - t1) * mrSges(5,2) + (g(1) * t10 - g(2) * t9 + t4 * t29 + t2) * mrSges(5,1);];
tau  = t16;
