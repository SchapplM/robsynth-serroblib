% Calculate vector of cutting forces with Newton-Euler
% S4RRRP1
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% Datum: 2019-05-04 19:24
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4RRRP1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP1_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP1_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP1_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP1_invdynf_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP1_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP1_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP1_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:24:05
% EndTime: 2019-05-04 19:24:06
% DurationCPUTime: 0.17s
% Computational Cost: add. (1631->56), mult. (2001->64), div. (0->0), fcn. (1008->6), ass. (0->36)
t44 = -m(3) - m(4);
t26 = qJDD(1) + qJDD(2);
t30 = sin(qJ(1));
t33 = cos(qJ(1));
t38 = t30 * g(1) - t33 * g(2);
t17 = qJDD(1) * pkin(1) + t38;
t34 = qJD(1) ^ 2;
t35 = -t33 * g(1) - t30 * g(2);
t18 = -t34 * pkin(1) + t35;
t29 = sin(qJ(2));
t32 = cos(qJ(2));
t36 = t32 * t17 - t29 * t18;
t12 = t26 * pkin(2) + t36;
t27 = qJD(1) + qJD(2);
t25 = t27 ^ 2;
t40 = t29 * t17 + t32 * t18;
t13 = -t25 * pkin(2) + t40;
t28 = sin(qJ(3));
t31 = cos(qJ(3));
t43 = t28 * t12 + t31 * t13;
t21 = qJDD(3) + t26;
t37 = t31 * t12 - t28 * t13;
t42 = t21 * mrSges(5,1) + m(5) * (t21 * pkin(3) + t37);
t41 = -mrSges(4,2) - mrSges(5,2);
t39 = -m(2) + t44;
t23 = m(5) * (-g(3) + qJDD(4));
t22 = qJD(3) + t27;
t20 = t22 ^ 2;
t8 = m(5) * (-t20 * pkin(3) + t43);
t6 = m(4) * t43 + t8 + t41 * t21 + (-mrSges(4,1) - mrSges(5,1)) * t20;
t5 = m(4) * t37 + t21 * mrSges(4,1) + t41 * t20 + t42;
t4 = m(3) * t40 - t25 * mrSges(3,1) - t26 * mrSges(3,2) - t28 * t5 + t31 * t6;
t3 = m(3) * t36 + t26 * mrSges(3,1) - t25 * mrSges(3,2) + t28 * t6 + t31 * t5;
t2 = m(2) * t35 - t34 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t29 * t3 + t32 * t4;
t1 = m(2) * t38 + qJDD(1) * mrSges(2,1) - t34 * mrSges(2,2) + t29 * t4 + t32 * t3;
t7 = [-m(1) * g(1) - t30 * t1 + t33 * t2, t2, t4, t6, -t20 * mrSges(5,1) - t21 * mrSges(5,2) + t8; -m(1) * g(2) + t33 * t1 + t30 * t2, t1, t3, t5, -t20 * mrSges(5,2) + t42; t23 + (-m(1) + t39) * g(3), t39 * g(3) + t23, t44 * g(3) + t23, -m(4) * g(3) + t23, t23;];
f_new  = t7;
