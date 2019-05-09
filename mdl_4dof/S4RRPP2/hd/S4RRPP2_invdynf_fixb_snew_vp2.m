% Calculate vector of cutting forces with Newton-Euler
% S4RRPP2
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
%   pkin=[a2,a3,a4,d1,d2]';
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
% Datum: 2019-05-04 19:20
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4RRPP2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP2_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP2_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPP2_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP2_invdynf_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP2_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPP2_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPP2_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:20:13
% EndTime: 2019-05-04 19:20:14
% DurationCPUTime: 0.15s
% Computational Cost: add. (734->58), mult. (856->59), div. (0->0), fcn. (316->4), ass. (0->31)
t22 = qJDD(1) + qJDD(2);
t23 = (qJD(1) + qJD(2));
t27 = sin(qJ(1));
t29 = cos(qJ(1));
t38 = t27 * g(1) - t29 * g(2);
t13 = qJDD(1) * pkin(1) + t38;
t30 = qJD(1) ^ 2;
t34 = -t29 * g(1) - t27 * g(2);
t14 = -t30 * pkin(1) + t34;
t26 = sin(qJ(2));
t28 = cos(qJ(2));
t41 = t26 * t13 + t28 * t14;
t36 = t22 * qJ(3) + (2 * qJD(3) * t23) + t41;
t46 = t23 ^ 2;
t47 = t22 * mrSges(4,3) + m(4) * (-(pkin(2) * t46) + t36);
t45 = -m(3) - m(4);
t44 = -pkin(2) - pkin(3);
t43 = m(5) * (g(3) + qJDD(4));
t42 = mrSges(3,1) + mrSges(4,1);
t40 = -m(2) + t45;
t39 = -(t46 * mrSges(5,1)) + m(5) * ((t44 * t46) + t36);
t37 = t28 * t13 - t26 * t14;
t31 = -qJ(3) * t46 + qJDD(3) - t37;
t35 = (t46 * mrSges(5,2)) + t22 * mrSges(5,1) - m(5) * (t44 * t22 + t31);
t33 = t22 * mrSges(5,2) + t39;
t32 = -m(4) * (-t22 * pkin(2) + t31) + t35;
t4 = m(3) * t37 + t42 * t22 + (-mrSges(3,2) + mrSges(4,3)) * t46 + t32;
t3 = m(3) * t41 + (-mrSges(3,2) + mrSges(5,2)) * t22 - t42 * t46 + t39 + t47;
t2 = m(2) * t34 - t30 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t26 * t4 + t28 * t3;
t1 = m(2) * t38 + qJDD(1) * mrSges(2,1) - t30 * mrSges(2,2) + t26 * t3 + t28 * t4;
t5 = [-m(1) * g(1) - t27 * t1 + t29 * t2, t2, t3, -(mrSges(4,1) * t46) + t33 + t47, t33; -m(1) * g(2) + t29 * t1 + t27 * t2, t1, t4, -m(4) * g(3) - t43, -t35; -t43 + (-m(1) + t40) * g(3), t40 * g(3) - t43, t45 * g(3) - t43, -t22 * mrSges(4,1) - mrSges(4,3) * t46 - t32, t43;];
f_new  = t5;
