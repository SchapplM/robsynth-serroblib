% Calculate vector of inverse dynamics joint torques for
% S4PRPP1
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
%   pkin=[a2,a3,a4,d2,theta1]';
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
% Datum: 2019-03-08 18:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PRPP1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP1_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPP1_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPP1_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPP1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP1_invdynJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPP1_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPP1_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPP1_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:17:54
% EndTime: 2019-03-08 18:17:55
% DurationCPUTime: 0.24s
% Computational Cost: add. (121->58), mult. (151->50), div. (0->0), fcn. (28->2), ass. (0->21)
t26 = -m(4) - m(5);
t25 = mrSges(4,2) - mrSges(5,3);
t18 = qJD(2) * qJD(3);
t4 = -qJ(3) * qJDD(2) - t18;
t24 = -mrSges(3,1) + t25;
t12 = pkin(5) + qJ(2);
t10 = sin(t12);
t11 = cos(t12);
t23 = -g(1) * t11 - g(2) * t10;
t22 = t23 - t4;
t13 = -pkin(2) - qJ(4);
t20 = pkin(2) * qJDD(2);
t17 = qJD(2) * qJD(4);
t15 = -g(1) * t10 + g(2) * t11;
t14 = qJD(2) ^ 2;
t9 = qJDD(3) - t20;
t8 = qJ(3) * qJD(2) + qJD(4);
t3 = t13 * qJD(2) + qJD(3);
t2 = qJDD(4) - t4;
t1 = t13 * qJDD(2) + qJDD(3) - t17;
t5 = [(-g(3) + qJDD(1)) * (m(2) + m(3) - t26); m(4) * (-t9 * pkin(2) + (-t4 + t18) * qJ(3)) + m(5) * (t2 * qJ(3) + t8 * qJD(3) - t3 * qJD(4) + t1 * t13) + (t17 - t1) * mrSges(5,3) + (-t20 + t9) * mrSges(4,2) + (t10 * mrSges(3,2) + t26 * (t11 * pkin(2) + t10 * qJ(3)) + (-m(5) * qJ(4) + t24) * t11) * g(2) + ((m(4) * pkin(2) - m(5) * t13 - t24) * t10 + (t26 * qJ(3) + mrSges(3,2)) * t11) * g(1) + (t22 - t4) * mrSges(4,3) + (t2 + t22) * mrSges(5,2) + (-t13 * mrSges(5,3) + Ifges(4,1) + Ifges(5,1) + Ifges(3,3)) * qJDD(2); (-mrSges(4,3) - mrSges(5,2)) * t14 + t25 * qJDD(2) + (-t8 * qJD(2) + t1 + t15) * m(5) + (-t14 * qJ(3) + t15 + t9) * m(4); qJDD(2) * mrSges(5,2) - t14 * mrSges(5,3) + (t3 * qJD(2) + t2 + t23) * m(5);];
tau  = t5;
