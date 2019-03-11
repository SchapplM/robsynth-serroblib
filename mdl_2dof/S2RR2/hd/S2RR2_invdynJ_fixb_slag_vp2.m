% Calculate vector of inverse dynamics joint torques for
% S2RR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% qJDD [2x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d2]';
% m_mdh [3x1]
%   mass of all robot links (including the base)
% mrSges [3x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [3x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [2x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S2RR2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(2,1),zeros(3,1),zeros(1,1),zeros(3,1),zeros(3,3),zeros(3,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR2_invdynJ_fixb_slag_vp2: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR2_invdynJ_fixb_slag_vp2: qJD has to be [2x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [2 1]), ...
  'S2RR2_invdynJ_fixb_slag_vp2: qJDD has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S2RR2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR2_invdynJ_fixb_slag_vp2: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'S2RR2_invdynJ_fixb_slag_vp2: m has to be [3x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [3,3]), ...
  'S2RR2_invdynJ_fixb_slag_vp2: mrSges has to be [3x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [3 6]), ...
  'S2RR2_invdynJ_fixb_slag_vp2: Ifges has to be [3x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:00:52
% EndTime: 2019-03-08 18:00:53
% DurationCPUTime: 0.22s
% Computational Cost: add. (90->56), mult. (239->92), div. (0->0), fcn. (108->4), ass. (0->22)
t10 = sin(qJ(2));
t12 = cos(qJ(2));
t21 = qJD(1) * t12;
t22 = qJD(1) * t10;
t24 = Ifges(3,2) * t12;
t26 = Ifges(3,4) * t10;
t27 = t12 / 0.2e1;
t9 = Ifges(3,4) * t21;
t28 = -(t10 * (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t21) + t12 * (qJD(2) * mrSges(3,1) - mrSges(3,3) * t22)) * pkin(1) - t10 * (Ifges(3,6) * qJD(2) + (t24 + t26) * qJD(1)) / 0.2e1 + (Ifges(3,1) * t22 + Ifges(3,5) * qJD(2) + t9) * t27;
t25 = Ifges(3,5) * t12;
t23 = Ifges(3,6) * t10;
t20 = qJD(1) * qJD(2);
t19 = m(3) * pkin(1) - mrSges(2,2) + mrSges(3,3);
t11 = sin(qJ(1));
t13 = cos(qJ(1));
t18 = g(1) * t11 + g(3) * t13;
t14 = t10 * (Ifges(3,1) * t12 - t26);
t6 = t10 * qJDD(1) + t12 * t20;
t5 = t12 * qJDD(1) - t10 * t20;
t4 = pkin(1) * t6;
t3 = pkin(1) * t5;
t1 = [Ifges(2,3) * qJDD(1) + (t13 * mrSges(2,1) + t19 * t11) * g(3) + (t11 * mrSges(2,1) - t19 * t13) * g(1) + (t3 * mrSges(3,3) + Ifges(3,4) * t6 + Ifges(3,2) * t5 + Ifges(3,6) * qJDD(2) + (m(3) * t3 - qJDD(2) * mrSges(3,2) + t5 * mrSges(3,3)) * pkin(1) + t18 * mrSges(3,1)) * t12 + (t4 * mrSges(3,3) + Ifges(3,1) * t6 + Ifges(3,4) * t5 + Ifges(3,5) * qJDD(2) + (m(3) * t4 - qJDD(2) * mrSges(3,1) + t6 * mrSges(3,3)) * pkin(1) - t18 * mrSges(3,2)) * t10 + ((t25 / 0.2e1 - t23 / 0.2e1) * qJD(2) + (t14 / 0.2e1 + (Ifges(3,4) * t12 - Ifges(3,2) * t10) * t27) * qJD(1) + t28) * qJD(2); Ifges(3,5) * t6 + Ifges(3,6) * t5 + Ifges(3,3) * qJDD(2) - t3 * mrSges(3,2) - t4 * mrSges(3,1) - g(2) * (t12 * mrSges(3,1) - t10 * mrSges(3,2)) + (-t12 * t9 / 0.2e1 - qJD(2) * (-t23 + t25) / 0.2e1 + (-t14 / 0.2e1 + t10 * t24 / 0.2e1) * qJD(1) - t28) * qJD(1) + (g(1) * t13 - g(3) * t11) * (mrSges(3,1) * t10 + mrSges(3,2) * t12);];
tau  = t1;
