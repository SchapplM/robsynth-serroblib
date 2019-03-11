% Calculate vector of inverse dynamics joint torques for
% S4PRRP1
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
%   pkin=[a2,a3,a4,d2,d3,theta1]';
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PRRP1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP1_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP1_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP1_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP1_invdynJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP1_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP1_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRP1_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:22:49
% EndTime: 2019-03-08 18:22:49
% DurationCPUTime: 0.30s
% Computational Cost: add. (229->72), mult. (282->85), div. (0->0), fcn. (90->6), ass. (0->38)
t44 = mrSges(4,1) + mrSges(5,1);
t52 = m(5) * pkin(3) + t44;
t50 = mrSges(4,2) - mrSges(5,3);
t32 = cos(qJ(3));
t31 = sin(qJ(3));
t41 = qJD(3) * t31;
t5 = (qJD(2) * t41 - qJDD(2) * t32) * pkin(2);
t49 = t44 * t31;
t47 = m(4) + m(5);
t48 = t47 * pkin(2) + mrSges(3,1);
t28 = qJDD(2) + qJDD(3);
t46 = t28 * pkin(3);
t45 = t31 * pkin(2);
t43 = pkin(2) * qJD(2);
t40 = t32 * t43;
t6 = qJD(3) * t40 + qJDD(2) * t45;
t42 = pkin(2) * qJD(3);
t29 = pkin(6) + qJ(2);
t27 = qJ(3) + t29;
t18 = sin(t27);
t3 = qJDD(4) - t46 + t5;
t38 = g(1) * t18 - t3;
t19 = cos(t27);
t37 = t50 * t19;
t30 = qJD(2) + qJD(3);
t36 = t28 * qJ(4) + t30 * qJD(4);
t2 = t36 + t6;
t35 = -t5 * mrSges(4,1) - t6 * mrSges(4,2) + t2 * mrSges(5,3) + (Ifges(5,2) + Ifges(4,3)) * t28;
t33 = -t52 * t19 + (-m(5) * qJ(4) + t50) * t18;
t26 = cos(t29);
t25 = sin(t29);
t21 = -t32 * pkin(2) - pkin(3);
t20 = qJ(4) + t45;
t15 = t32 * t42 + qJD(4);
t10 = t19 * qJ(4);
t8 = t30 * qJ(4) + t31 * t43;
t7 = -t30 * pkin(3) + qJD(4) - t40;
t1 = [(-g(3) + qJDD(1)) * (m(2) + m(3) + t47); -t3 * mrSges(5,1) + Ifges(3,3) * qJDD(2) + m(4) * (t31 * t6 - t32 * t5) * pkin(2) + m(5) * (t7 * pkin(2) * t41 + t8 * t15 + t2 * t20 + t3 * t21) + (t25 * mrSges(3,2) - t48 * t26 + t33) * g(2) + (-m(5) * t10 + t26 * mrSges(3,2) + t52 * t18 + t48 * t25 + t37) * g(1) + (t15 * mrSges(5,3) + (-t32 * mrSges(4,2) - t49) * t42) * t30 + (-t21 * mrSges(5,1) + t20 * mrSges(5,3) + (mrSges(4,1) * t32 - mrSges(4,2) * t31) * pkin(2)) * t28 + t35; m(5) * (-t3 * pkin(3) + t2 * qJ(4) + t8 * qJD(4)) + t33 * g(2) + (t18 * mrSges(4,1) - m(5) * (-t18 * pkin(3) + t10) + t37) * g(1) + t36 * mrSges(5,3) + (t38 + t46) * mrSges(5,1) + (-m(5) * (t31 * t7 + t32 * t8) + (t50 * t32 + t49) * t30) * t43 + t35; -t30 ^ 2 * mrSges(5,3) - t28 * mrSges(5,1) + (g(2) * t19 - t8 * t30 - t38) * m(5);];
tau  = t1;
