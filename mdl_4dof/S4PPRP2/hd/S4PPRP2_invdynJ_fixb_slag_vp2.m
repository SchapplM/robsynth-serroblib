% Calculate vector of inverse dynamics joint torques for
% S4PPRP2
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
%   pkin=[a2,a3,a4,d3,theta2]';
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
% Datum: 2019-03-08 18:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PPRP2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP2_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRP2_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRP2_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRP2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPRP2_invdynJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRP2_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRP2_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRP2_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:13:15
% EndTime: 2019-03-08 18:13:15
% DurationCPUTime: 0.22s
% Computational Cost: add. (166->60), mult. (376->74), div. (0->0), fcn. (240->6), ass. (0->32)
t22 = sin(pkin(5));
t23 = cos(pkin(5));
t24 = sin(qJ(3));
t25 = cos(qJ(3));
t12 = t25 * t22 + t24 * t23;
t34 = t25 * t23;
t37 = t24 * t22;
t11 = -t34 + t37;
t7 = t11 * qJD(1);
t39 = qJD(4) + t7;
t38 = t11 * qJDD(1);
t33 = mrSges(4,1) + mrSges(5,1);
t32 = -mrSges(4,2) + mrSges(5,3);
t31 = m(3) + m(4) + m(5);
t30 = qJD(1) * qJD(3);
t29 = t12 * qJDD(1) + t30 * t34;
t28 = qJD(1) * t37;
t27 = m(5) * pkin(3) + t33;
t26 = -m(5) * qJ(4) - t32;
t8 = t12 * qJD(1);
t21 = pkin(5) + qJ(3);
t20 = cos(t21);
t19 = sin(t21);
t10 = t11 * qJD(3);
t9 = t12 * qJD(3);
t6 = qJD(3) * qJ(4) + t8;
t5 = -qJD(3) * pkin(3) + t39;
t4 = -t12 * t30 - t38;
t3 = -qJD(3) * t28 + t29;
t2 = -qJDD(3) * pkin(3) + t8 * qJD(3) + qJDD(4) + t38;
t1 = qJDD(3) * qJ(4) + (qJD(4) - t28) * qJD(3) + t29;
t13 = [m(4) * (-t8 * t10 - t4 * t11 + t3 * t12 + t7 * t9) + m(5) * (t1 * t12 - t6 * t10 + t2 * t11 + t5 * t9) + (-t33 * t11 + t32 * t12) * qJDD(3) + (-t32 * t10 - t33 * t9) * qJD(3) + (-m(2) - t31) * g(2) + (m(2) + m(3) * (t22 ^ 2 + t23 ^ 2)) * qJDD(1); (-g(3) + qJDD(2)) * t31; t4 * mrSges(4,1) - t2 * mrSges(5,1) - t3 * mrSges(4,2) + t1 * mrSges(5,3) + (t26 * t19 - t27 * t20) * g(2) + (t27 * t19 + t26 * t20) * g(1) + (pkin(3) * mrSges(5,1) + qJ(4) * mrSges(5,3) + Ifges(5,2) + Ifges(4,3)) * qJDD(3) + (-t7 * mrSges(4,2) + t39 * mrSges(5,3) + t33 * t8) * qJD(3) + (-t2 * pkin(3) + t1 * qJ(4) + t39 * t6 - t5 * t8) * m(5); -qJD(3) ^ 2 * mrSges(5,3) - qJDD(3) * mrSges(5,1) + (-g(1) * t19 + g(2) * t20 - t6 * qJD(3) + t2) * m(5);];
tau  = t13;
