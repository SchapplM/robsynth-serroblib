% Calculate vector of inverse dynamics joint torques for
% S2RR1
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
% Datum: 2019-03-08 18:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S2RR1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(2,1),zeros(3,1),zeros(1,1),zeros(3,1),zeros(3,3),zeros(3,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR1_invdynJ_fixb_slag_vp2: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR1_invdynJ_fixb_slag_vp2: qJD has to be [2x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [2 1]), ...
  'S2RR1_invdynJ_fixb_slag_vp2: qJDD has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S2RR1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR1_invdynJ_fixb_slag_vp2: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'S2RR1_invdynJ_fixb_slag_vp2: m has to be [3x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [3,3]), ...
  'S2RR1_invdynJ_fixb_slag_vp2: mrSges has to be [3x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [3 6]), ...
  'S2RR1_invdynJ_fixb_slag_vp2: Ifges has to be [3x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 17:59:55
% EndTime: 2019-03-08 17:59:56
% DurationCPUTime: 0.25s
% Computational Cost: add. (90->50), mult. (239->78), div. (0->0), fcn. (108->4), ass. (0->20)
t11 = cos(qJ(2));
t29 = t11 / 0.2e1;
t9 = sin(qJ(2));
t28 = t11 * mrSges(3,1) - t9 * mrSges(3,2);
t23 = Ifges(3,4) * t11;
t25 = Ifges(3,4) * t9;
t26 = t9 / 0.2e1;
t27 = -t28 * qJD(2) * pkin(1) + ((-Ifges(3,1) * t11 + t25) * t26 + (Ifges(3,2) * t9 - t23) * t29) * qJD(1) + (Ifges(3,5) * qJD(2) + (-Ifges(3,1) * t9 - t23) * qJD(1)) * t29 - (Ifges(3,6) * qJD(2) + (-Ifges(3,2) * t11 - t25) * qJD(1)) * t26;
t24 = Ifges(3,6) * t9;
t22 = Ifges(3,5) * t11;
t20 = qJD(1) * qJD(2);
t19 = -m(3) * pkin(1) - mrSges(2,2) - mrSges(3,3);
t10 = sin(qJ(1));
t12 = cos(qJ(1));
t17 = g(1) * t12 - g(3) * t10;
t5 = -t11 * qJDD(1) + t9 * t20;
t14 = t9 * qJDD(1) + t11 * t20;
t4 = t14 * pkin(1);
t3 = pkin(1) * t5;
t1 = [Ifges(2,3) * qJDD(1) + (-t10 * mrSges(2,1) + t19 * t12) * g(3) + (t12 * mrSges(2,1) + t19 * t10) * g(1) + (t4 * mrSges(3,3) + Ifges(3,1) * t14 - Ifges(3,4) * t5 - Ifges(3,5) * qJDD(2) + (m(3) * t4 + qJDD(2) * mrSges(3,1) + mrSges(3,3) * t14) * pkin(1) - t17 * mrSges(3,2)) * t9 + (-t3 * mrSges(3,3) + Ifges(3,4) * t14 - Ifges(3,2) * t5 - Ifges(3,6) * qJDD(2) + (-m(3) * t3 + qJDD(2) * mrSges(3,2) - t5 * mrSges(3,3)) * pkin(1) + t17 * mrSges(3,1)) * t11 + ((-t22 / 0.2e1 + t24 / 0.2e1) * qJD(2) - t27) * qJD(2); -Ifges(3,5) * t14 + Ifges(3,6) * t5 + Ifges(3,3) * qJDD(2) - t3 * mrSges(3,2) + t4 * mrSges(3,1) + g(2) * t28 + (-qJD(2) * (-t22 + t24) / 0.2e1 + t27) * qJD(1) + (-g(1) * t10 - g(3) * t12) * (mrSges(3,1) * t9 + mrSges(3,2) * t11);];
tau  = t1;
