% Calculate vector of cutting forces with Newton-Euler
% S4RPRP2
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
%   pkin=[a2,a3,a4,d1,d3]';
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
% Datum: 2019-05-04 19:15
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4RPRP2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP2_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP2_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP2_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP2_invdynf_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP2_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP2_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP2_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:15:07
% EndTime: 2019-05-04 19:15:07
% DurationCPUTime: 0.14s
% Computational Cost: add. (735->57), mult. (981->59), div. (0->0), fcn. (276->4), ass. (0->31)
t43 = -m(3) - m(4);
t42 = -pkin(1) - pkin(2);
t29 = qJD(1) ^ 2;
t26 = sin(qJ(1));
t28 = cos(qJ(1));
t34 = -t28 * g(1) - t26 * g(2);
t32 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t34;
t11 = t29 * t42 + t32;
t36 = t26 * g(1) - t28 * g(2);
t30 = -t29 * qJ(2) + qJDD(2) - t36;
t13 = qJDD(1) * t42 + t30;
t25 = sin(qJ(3));
t27 = cos(qJ(3));
t41 = t27 * t11 + t25 * t13;
t20 = -qJDD(1) + qJDD(3);
t35 = -t25 * t11 + t27 * t13;
t40 = t20 * mrSges(5,1) + m(5) * (t20 * pkin(3) + t35);
t39 = mrSges(2,1) + mrSges(3,1);
t38 = (-mrSges(4,2) - mrSges(5,2));
t37 = -m(2) + t43;
t21 = -qJD(1) + qJD(3);
t19 = t21 ^ 2;
t4 = m(4) * t35 + t20 * mrSges(4,1) + (t19 * t38) + t40;
t7 = m(5) * (-t19 * pkin(3) + t41);
t5 = m(4) * t41 + t7 + t38 * t20 + (-mrSges(4,1) - mrSges(5,1)) * t19;
t33 = -t25 * t4 + t27 * t5 + m(3) * (-t29 * pkin(1) + t32) + qJDD(1) * mrSges(3,3);
t31 = -m(3) * (-qJDD(1) * pkin(1) + t30) - t25 * t5 - t27 * t4;
t16 = m(5) * (g(3) + qJDD(4));
t2 = m(2) * t36 + (-mrSges(2,2) + mrSges(3,3)) * t29 + t39 * qJDD(1) + t31;
t1 = m(2) * t34 - qJDD(1) * mrSges(2,2) - t29 * t39 + t33;
t3 = [-m(1) * g(1) + t28 * t1 - t26 * t2, t1, -t29 * mrSges(3,1) + t33, t5, -t19 * mrSges(5,1) - t20 * mrSges(5,2) + t7; -m(1) * g(2) + t26 * t1 + t28 * t2, t2, g(3) * t43 - t16, t4, -t19 * mrSges(5,2) + t40; -t16 + (-m(1) + t37) * g(3), g(3) * t37 - t16, -qJDD(1) * mrSges(3,1) - t29 * mrSges(3,3) - t31, m(4) * g(3) + t16, t16;];
f_new  = t3;
