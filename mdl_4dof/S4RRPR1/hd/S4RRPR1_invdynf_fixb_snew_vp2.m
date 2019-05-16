% Calculate vector of cutting forces with Newton-Euler
% S4RRPR1
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2019-05-04 19:21
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4RRPR1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR1_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR1_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR1_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR1_invdynf_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR1_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR1_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR1_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:21:31
% EndTime: 2019-05-04 19:21:31
% DurationCPUTime: 0.24s
% Computational Cost: add. (3379->57), mult. (4559->72), div. (0->0), fcn. (2480->8), ass. (0->40)
t47 = -m(2) - m(3);
t28 = qJDD(1) + qJDD(2);
t35 = sin(qJ(1));
t38 = cos(qJ(1));
t43 = t35 * g(1) - t38 * g(2);
t19 = qJDD(1) * pkin(1) + t43;
t39 = qJD(1) ^ 2;
t40 = -t38 * g(1) - t35 * g(2);
t20 = -t39 * pkin(1) + t40;
t34 = sin(qJ(2));
t37 = cos(qJ(2));
t41 = t37 * t19 - t34 * t20;
t14 = t28 * pkin(2) + t41;
t29 = qJD(1) + qJD(2);
t27 = t29 ^ 2;
t45 = t34 * t19 + t37 * t20;
t15 = -t27 * pkin(2) + t45;
t31 = sin(pkin(7));
t32 = cos(pkin(7));
t46 = t31 * t14 + t32 * t15;
t30 = -g(3) + qJDD(3);
t24 = m(5) * t30;
t44 = m(4) * t30 + t24;
t42 = t32 * t14 - t31 * t15;
t36 = cos(qJ(4));
t33 = sin(qJ(4));
t23 = qJD(4) + t29;
t22 = qJDD(4) + t28;
t21 = t23 ^ 2;
t10 = -t27 * pkin(3) + t46;
t9 = t28 * pkin(3) + t42;
t8 = m(5) * (t36 * t10 + t33 * t9) - t22 * mrSges(5,2) - t21 * mrSges(5,1);
t7 = m(5) * (-t33 * t10 + t36 * t9) + t22 * mrSges(5,1) - t21 * mrSges(5,2);
t6 = m(4) * t46 - t27 * mrSges(4,1) - t28 * mrSges(4,2) - t33 * t7 + t36 * t8;
t5 = m(4) * t42 + t28 * mrSges(4,1) - t27 * mrSges(4,2) + t33 * t8 + t36 * t7;
t4 = m(3) * t45 - t27 * mrSges(3,1) - t28 * mrSges(3,2) - t31 * t5 + t32 * t6;
t3 = m(3) * t41 + t28 * mrSges(3,1) - t27 * mrSges(3,2) + t31 * t6 + t32 * t5;
t2 = m(2) * t40 - t39 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t34 * t3 + t37 * t4;
t1 = m(2) * t43 + qJDD(1) * mrSges(2,1) - t39 * mrSges(2,2) + t37 * t3 + t34 * t4;
t11 = [-m(1) * g(1) - t35 * t1 + t38 * t2, t2, t4, t6, t8; -m(1) * g(2) + t38 * t1 + t35 * t2, t1, t3, t5, t7; (-m(1) + t47) * g(3) + t44, t47 * g(3) + t44, -m(3) * g(3) + t44, t44, t24;];
f_new  = t11;
