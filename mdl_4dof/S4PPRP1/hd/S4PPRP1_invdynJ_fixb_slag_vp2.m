% Calculate vector of inverse dynamics joint torques for
% S4PPRP1
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
%   pkin=[a2,a3,a4,d3,theta1]';
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
% Datum: 2019-03-08 18:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PPRP1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP1_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRP1_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRP1_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRP1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPRP1_invdynJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRP1_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRP1_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRP1_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:12:14
% EndTime: 2019-03-08 18:12:14
% DurationCPUTime: 0.18s
% Computational Cost: add. (108->57), mult. (215->66), div. (0->0), fcn. (110->4), ass. (0->22)
t15 = sin(qJ(3));
t16 = cos(qJ(3));
t26 = mrSges(4,2) - mrSges(5,3);
t27 = mrSges(4,1) + mrSges(5,1);
t28 = t27 * t15 + t26 * t16;
t24 = qJD(2) * qJD(3);
t6 = t15 * qJDD(2) + t16 * t24;
t25 = sin(pkin(5));
t23 = qJD(3) * qJD(4);
t21 = m(5) * qJ(4) - mrSges(4,2);
t20 = m(5) * pkin(3) + t27;
t8 = qJD(3) * qJ(4) + t15 * qJD(2);
t19 = t15 * (-qJD(3) * pkin(3) - t16 * qJD(2) + qJD(4)) + t16 * t8;
t14 = cos(pkin(5));
t18 = -g(1) * t25 + g(2) * t14;
t5 = t16 * qJDD(2) - t15 * t24;
t17 = qJD(3) ^ 2;
t4 = t14 * t15 - t25 * t16;
t3 = -t14 * t16 - t25 * t15;
t2 = -qJDD(3) * pkin(3) + qJDD(4) - t5;
t1 = qJDD(3) * qJ(4) + t23 + t6;
t7 = [(-g(3) + qJDD(1)) * (m(2) + m(3) + m(4) + m(5)); (t19 * qJD(3) + t1 * t15 - t2 * t16 + t18) * m(5) + (t6 * t15 + t5 * t16 + t18) * m(4) + (qJDD(2) + t18) * m(3) - t28 * t17 + (-t26 * t15 + t27 * t16) * qJDD(3); -t6 * mrSges(4,2) + t5 * mrSges(4,1) - t2 * mrSges(5,1) + m(5) * (-t2 * pkin(3) + t1 * qJ(4) + t8 * qJD(4)) + (-t20 * t3 + t21 * t4) * g(2) + (t20 * t4 + t21 * t3) * g(1) + (g(1) * t3 + g(2) * t4 + t1 + t23) * mrSges(5,3) + (pkin(3) * mrSges(5,1) + qJ(4) * mrSges(5,3) + Ifges(5,2) + Ifges(4,3)) * qJDD(3) + (-m(5) * t19 + t28 * qJD(3)) * qJD(2); -qJDD(3) * mrSges(5,1) - t17 * mrSges(5,3) + (-g(1) * t4 + g(2) * t3 - t8 * qJD(3) + t2) * m(5);];
tau  = t7;
