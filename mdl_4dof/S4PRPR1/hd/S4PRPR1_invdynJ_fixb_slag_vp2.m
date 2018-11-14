% Calculate vector of inverse dynamics joint torques for
% S4PRPR1
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
%   pkin=[a2,a3,a4,d2,d4,theta1]';
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
% Datum: 2018-11-14 13:43
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau = S4PRPR1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR1_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR1_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR1_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR1_invdynJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR1_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR1_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR1_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:42:12
% EndTime: 2018-11-14 13:42:13
% DurationCPUTime: 0.31s
% Computational Cost: add. (263->80), mult. (368->95), div. (0->0), fcn. (142->4), ass. (0->34)
t36 = qJD(2) - qJD(4);
t28 = -pkin(2) - pkin(3);
t23 = -qJDD(2) + qJDD(4);
t41 = Ifges(5,3) * t23;
t40 = -mrSges(3,1) - mrSges(4,1);
t24 = pkin(6) + qJ(2);
t21 = sin(t24);
t22 = cos(t24);
t39 = t22 * pkin(2) + t21 * qJ(3);
t38 = qJ(3) * qJD(2);
t37 = qJD(2) * qJD(3);
t16 = qJ(3) * qJDD(2) + t37;
t34 = t16 + t37;
t33 = -g(1) * t21 + g(2) * t22;
t26 = sin(qJ(4));
t27 = cos(qJ(4));
t10 = t27 * qJ(3) + t26 * t28;
t9 = -t26 * qJ(3) + t27 * t28;
t11 = t28 * qJDD(2) + qJDD(3);
t15 = t28 * qJD(2) + qJD(3);
t5 = t27 * t15 - t26 * t38;
t1 = qJD(4) * t5 + t26 * t11 + t27 * t16;
t7 = -t21 * t26 - t22 * t27;
t8 = -t21 * t27 + t22 * t26;
t31 = -g(1) * t7 - g(2) * t8 - t1;
t6 = t26 * t15 + t27 * t38;
t2 = -qJD(4) * t6 + t27 * t11 - t26 * t16;
t30 = g(1) * t8 - g(2) * t7 + t2;
t29 = qJD(2) ^ 2;
t20 = -pkin(2) * qJDD(2) + qJDD(3);
t18 = t22 * qJ(3);
t4 = -t26 * qJD(3) - t10 * qJD(4);
t3 = t27 * qJD(3) + t9 * qJD(4);
t12 = [(-g(3) + qJDD(1)) * (m(2) + m(3) + m(4) + m(5)); -t41 - t20 * mrSges(4,1) + m(4) * (-t20 * pkin(2) + t34 * qJ(3)) + (-m(4) * t39 + t21 * mrSges(3,2) + t40 * t22) * g(2) + (-m(4) * t18 + t22 * mrSges(3,2) + (m(4) * pkin(2) - t40) * t21) * g(1) + (-g(1) * t22 - g(2) * t21 + t34) * mrSges(4,3) + (-g(1) * (t28 * t21 + t18) - g(2) * (t22 * pkin(3) + t39) + t1 * t10 + t2 * t9 + t6 * t3 + t5 * t4) * m(5) + (-t10 * t23 + t3 * t36 - t31) * mrSges(5,2) + (t9 * t23 - t36 * t4 - t30) * mrSges(5,1) + (pkin(2) * mrSges(4,1) + qJ(3) * mrSges(4,3) + Ifges(4,2) + Ifges(3,3)) * qJDD(2); -qJDD(2) * mrSges(4,1) - t29 * mrSges(4,3) + (t27 * mrSges(5,1) - t26 * mrSges(5,2)) * t23 + (t1 * t26 + t2 * t27 + t33 - t36 * (-t26 * t5 + t27 * t6)) * m(5) + (-t29 * qJ(3) + t20 + t33) * m(4) - (mrSges(5,1) * t26 + mrSges(5,2) * t27) * t36 ^ 2; t41 + (-t36 * t5 + t31) * mrSges(5,2) + (-t36 * t6 + t30) * mrSges(5,1);];
tau  = t12;
