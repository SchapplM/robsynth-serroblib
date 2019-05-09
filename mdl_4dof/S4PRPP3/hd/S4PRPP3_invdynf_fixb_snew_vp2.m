% Calculate vector of cutting forces with Newton-Euler
% S4PRPP3
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2]';
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
% f_new [3x5]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 18:56
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4PRPP3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(4,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP3_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPP3_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPP3_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPP3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S4PRPP3_invdynf_fixb_snew_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPP3_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPP3_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPP3_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:55:42
% EndTime: 2019-05-04 18:55:42
% DurationCPUTime: 0.12s
% Computational Cost: add. (285->47), mult. (370->44), div. (0->0), fcn. (80->2), ass. (0->24)
t37 = -m(3) - m(4);
t36 = -pkin(2) - pkin(3);
t22 = qJD(2) ^ 2;
t19 = -g(2) + qJDD(1);
t20 = sin(qJ(2));
t21 = cos(qJ(2));
t27 = -t21 * g(1) + t20 * t19;
t24 = qJDD(2) * qJ(3) + (2 * qJD(3) * qJD(2)) + t27;
t35 = qJDD(2) * mrSges(5,2) + m(5) * (t36 * t22 + t24);
t34 = -mrSges(4,1) - mrSges(5,1);
t33 = -mrSges(5,2) - mrSges(4,3);
t32 = m(2) - t37;
t26 = qJDD(2) * mrSges(4,3) + m(4) * (-t22 * pkin(2) + t24) + t35;
t30 = mrSges(3,1) - t34;
t3 = m(3) * t27 - qJDD(2) * mrSges(3,2) - t30 * t22 + t26;
t25 = t20 * g(1) + t21 * t19;
t23 = -t22 * qJ(3) + qJDD(3) - t25;
t7 = m(5) * (t36 * qJDD(2) + t23);
t29 = m(4) * (-qJDD(2) * pkin(2) + t23) + t7;
t5 = m(3) * t25 + (-mrSges(3,2) - t33) * t22 + t30 * qJDD(2) - t29;
t31 = m(2) * t19 + t20 * t3 + t21 * t5;
t28 = -t20 * t5 + t21 * t3;
t12 = m(5) * (g(3) + qJDD(4));
t1 = [(-m(1) - m(2)) * g(1) + t28, -m(2) * g(1) + t28, t3, t34 * t22 + t26, -t22 * mrSges(5,1) + t35; -m(1) * g(2) + t31, t32 * g(3) + t12, t5, -m(4) * g(3) - t12, -qJDD(2) * mrSges(5,1) - t22 * mrSges(5,2) + t7; -t12 + (-m(1) - t32) * g(3), t31, t37 * g(3) - t12, t34 * qJDD(2) + t33 * t22 + t29, t12;];
f_new  = t1;
