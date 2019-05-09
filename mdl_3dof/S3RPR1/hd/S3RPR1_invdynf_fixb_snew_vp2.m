% Calculate vector of cutting forces with Newton-Euler
% S3RPR1
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% qJDD [3x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d3]';
% m_mdh [4x1]
%   mass of all robot links (including the base)
% mrSges [4x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% f_new [3x4]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 18:28
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S3RPR1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPR1_invdynf_fixb_snew_vp2: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RPR1_invdynf_fixb_snew_vp2: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3RPR1_invdynf_fixb_snew_vp2: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RPR1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RPR1_invdynf_fixb_snew_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'S3RPR1_invdynf_fixb_snew_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'S3RPR1_invdynf_fixb_snew_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'S3RPR1_invdynf_fixb_snew_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:28:43
% EndTime: 2019-05-04 18:28:44
% DurationCPUTime: 0.09s
% Computational Cost: add. (368->42), mult. (526->50), div. (0->0), fcn. (156->4), ass. (0->25)
t31 = -m(3) - m(4);
t30 = -pkin(1) - pkin(2);
t29 = mrSges(2,1) + mrSges(3,1);
t28 = -m(2) + t31;
t18 = sin(qJ(1));
t20 = cos(qJ(1));
t27 = t18 * g(1) - t20 * g(2);
t26 = -t20 * g(1) - t18 * g(2);
t17 = sin(qJ(3));
t19 = cos(qJ(3));
t21 = qJD(1) ^ 2;
t24 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t26;
t14 = -qJD(1) + qJD(3);
t12 = t14 ^ 2;
t13 = -qJDD(1) + qJDD(3);
t6 = t30 * t21 + t24;
t22 = -t21 * qJ(2) + qJDD(2) - t27;
t8 = t30 * qJDD(1) + t22;
t4 = m(4) * (-t17 * t6 + t19 * t8) + t13 * mrSges(4,1) - (t12 * mrSges(4,2));
t5 = m(4) * (t17 * t8 + t19 * t6) - t13 * mrSges(4,2) - t12 * mrSges(4,1);
t25 = -t17 * t4 + t19 * t5 + qJDD(1) * mrSges(3,3) + m(3) * (-t21 * pkin(1) + t24);
t23 = -m(3) * (-qJDD(1) * pkin(1) + t22) - t17 * t5 - t19 * t4;
t2 = m(2) * t27 + (-mrSges(2,2) + mrSges(3,3)) * t21 + t29 * qJDD(1) + t23;
t1 = m(2) * t26 - qJDD(1) * mrSges(2,2) - t29 * t21 + t25;
t3 = [-m(1) * g(1) + t20 * t1 - t18 * t2, t1, -t21 * mrSges(3,1) + t25, t5; -m(1) * g(2) + t18 * t1 + t20 * t2, t2, t31 * g(3), t4; (-m(1) + t28) * g(3), t28 * g(3), -qJDD(1) * mrSges(3,1) - t21 * mrSges(3,3) - t23, m(4) * g(3);];
f_new  = t3;
