% Calculate vector of cutting forces with Newton-Euler
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
% f_new [3x3]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 18:18
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S2RR1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(2,1),zeros(3,1),zeros(1,1),zeros(3,1),zeros(3,3),zeros(3,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR1_invdynf_fixb_snew_vp2: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR1_invdynf_fixb_snew_vp2: qJD has to be [2x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [2 1]), ...
  'S2RR1_invdynf_fixb_snew_vp2: qJDD has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S2RR1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR1_invdynf_fixb_snew_vp2: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'S2RR1_invdynf_fixb_snew_vp2: m has to be [3x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [3,3]), ...
  'S2RR1_invdynf_fixb_snew_vp2: mrSges has to be [3x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [3 6]), ...
  'S2RR1_invdynf_fixb_snew_vp2: Ifges has to be [3x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:17:41
% EndTime: 2019-05-04 18:17:41
% DurationCPUTime: 0.08s
% Computational Cost: add. (196->38), mult. (399->58), div. (0->0), fcn. (186->4), ass. (0->23)
t18 = sin(qJ(1));
t20 = cos(qJ(1));
t28 = t18 * g(1) + t20 * g(3);
t17 = sin(qJ(2));
t27 = qJD(1) * t17;
t19 = cos(qJ(2));
t26 = qJD(1) * t19;
t25 = qJD(1) * qJD(2);
t24 = t20 * g(1) - t18 * g(3);
t10 = -t17 * qJDD(1) - t19 * t25;
t13 = -qJD(2) * mrSges(3,2) - mrSges(3,3) * t26;
t8 = -qJDD(1) * pkin(1) + t28;
t9 = (mrSges(3,1) * t19 - mrSges(3,2) * t17) * qJD(1);
t3 = m(3) * (t19 * g(2) - t17 * t8) - t10 * mrSges(3,3) + qJDD(2) * mrSges(3,1) + t9 * t27 + qJD(2) * t13;
t11 = -t19 * qJDD(1) + t17 * t25;
t12 = qJD(2) * mrSges(3,1) + mrSges(3,3) * t27;
t4 = m(3) * (t17 * g(2) + t19 * t8) + t11 * mrSges(3,3) - qJDD(2) * mrSges(3,2) - t9 * t26 - qJD(2) * t12;
t23 = -t17 * t4 - t19 * t3;
t21 = qJD(1) ^ 2;
t22 = -t11 * mrSges(3,1) - t12 * t27 + t13 * t26 + m(3) * (-t21 * pkin(1) + t24) + t10 * mrSges(3,2);
t2 = m(2) * t24 + qJDD(1) * mrSges(2,1) - t21 * mrSges(2,2) + t22;
t1 = m(2) * t28 - t21 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t17 * t3 + t19 * t4;
t5 = [-m(1) * g(1) - t18 * t1 - t20 * t2, t1, t4; (-m(1) - m(2)) * g(2) + t23, t2, t3; -m(1) * g(3) - t20 * t1 + t18 * t2, -m(2) * g(2) + t23, t22;];
f_new  = t5;
