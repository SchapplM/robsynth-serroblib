% Calculate vector of cutting forces with Newton-Euler
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
% f_new [3x5]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 18:43
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4PPRP1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP1_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRP1_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRP1_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRP1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPRP1_invdynf_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRP1_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRP1_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRP1_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:43:08
% EndTime: 2019-05-04 18:43:08
% DurationCPUTime: 0.09s
% Computational Cost: add. (374->38), mult. (539->46), div. (0->0), fcn. (276->4), ass. (0->25)
t22 = qJD(3) ^ 2;
t18 = sin(pkin(5));
t19 = cos(pkin(5));
t26 = t18 * g(1) - t19 * g(2);
t10 = qJDD(2) - t26;
t12 = -t19 * g(1) - t18 * g(2);
t20 = sin(qJ(3));
t21 = cos(qJ(3));
t25 = t21 * t10 - t20 * t12;
t33 = m(5) * (-qJDD(3) * pkin(3) - t22 * qJ(4) + qJDD(4) - t25);
t32 = t20 * t10 + t21 * t12;
t31 = qJDD(3) * mrSges(5,3) + m(5) * (-t22 * pkin(3) + qJDD(3) * qJ(4) + (2 * qJD(4) * qJD(3)) + t32);
t30 = mrSges(4,1) + mrSges(5,1);
t17 = g(3) - qJDD(1);
t13 = m(5) * t17;
t29 = m(4) * t17 + t13;
t4 = m(4) * t32 - qJDD(3) * mrSges(4,2) - t30 * t22 + t31;
t5 = m(4) * t25 - t33 + (-mrSges(4,2) + mrSges(5,3)) * t22 + t30 * qJDD(3);
t28 = m(3) * t12 - t20 * t5 + t21 * t4;
t27 = -m(3) * t17 - t29;
t24 = -m(2) * t17 + t27;
t23 = m(3) * t10 + t20 * t4 + t21 * t5;
t2 = m(2) * t12 + t28;
t1 = m(2) * t26 - t23;
t3 = [-m(1) * g(1) - t18 * t1 + t19 * t2, t2, t28, t4, -t22 * mrSges(5,1) + t31; -m(1) * g(2) + t19 * t1 + t18 * t2, t1, t27, t5, t13; -m(1) * g(3) + t24, t24, t23, t29, -qJDD(3) * mrSges(5,1) - t22 * mrSges(5,3) + t33;];
f_new  = t3;
