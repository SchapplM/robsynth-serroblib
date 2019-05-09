% Calculate vector of cutting forces with Newton-Euler
% S4PRPP2
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
%   pkin=[a2,a3,a4,d2,theta3]';
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
% Datum: 2019-05-04 18:54
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4PRPP2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP2_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPP2_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPP2_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPP2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP2_invdynf_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPP2_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPP2_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPP2_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:54:21
% EndTime: 2019-05-04 18:54:21
% DurationCPUTime: 0.12s
% Computational Cost: add. (599->47), mult. (803->52), div. (0->0), fcn. (324->4), ass. (0->26)
t27 = qJD(2) ^ 2;
t22 = -g(2) + qJDD(1);
t25 = sin(qJ(2));
t26 = cos(qJ(2));
t32 = t25 * g(1) + t26 * t22;
t12 = qJDD(2) * pkin(2) + t32;
t29 = -t26 * g(1) + t25 * t22;
t13 = -t27 * pkin(2) + t29;
t23 = sin(pkin(5));
t24 = cos(pkin(5));
t28 = t24 * t12 - t23 * t13;
t38 = m(5) * (-qJDD(2) * pkin(3) - t27 * qJ(4) + qJDD(4) - t28);
t37 = m(2) + m(3);
t34 = t23 * t12 + t24 * t13;
t36 = qJDD(2) * mrSges(5,3) + m(5) * (-t27 * pkin(3) + qJDD(2) * qJ(4) + (2 * qJD(4) * qJD(2)) + t34);
t35 = mrSges(4,1) + mrSges(5,1);
t21 = -g(3) + qJDD(3);
t16 = m(5) * t21;
t33 = m(4) * t21 + t16;
t6 = m(4) * t34 - qJDD(2) * mrSges(4,2) - t35 * t27 + t36;
t7 = m(4) * t28 - t38 + (-mrSges(4,2) + mrSges(5,3)) * t27 + t35 * qJDD(2);
t4 = m(3) * t32 + qJDD(2) * mrSges(3,1) - t27 * mrSges(3,2) + t23 * t6 + t24 * t7;
t5 = m(3) * t29 - t27 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t23 * t7 + t24 * t6;
t31 = m(2) * t22 + t25 * t5 + t26 * t4;
t30 = -t25 * t4 + t26 * t5;
t1 = [(-m(1) - m(2)) * g(1) + t30, -m(2) * g(1) + t30, t5, t6, -t27 * mrSges(5,1) + t36; -m(1) * g(2) + t31, t37 * g(3) - t33, t4, t7, t16; (-m(1) - t37) * g(3) + t33, t31, -m(3) * g(3) + t33, t33, -qJDD(2) * mrSges(5,1) - t27 * mrSges(5,3) + t38;];
f_new  = t1;
