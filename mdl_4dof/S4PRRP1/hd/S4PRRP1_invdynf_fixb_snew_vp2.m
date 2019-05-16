% Calculate vector of cutting forces with Newton-Euler
% S4PRRP1
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
%   pkin=[a2,a3,a4,d2,d3,theta1]';
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
% Datum: 2019-05-04 19:02
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4PRRP1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP1_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP1_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP1_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP1_invdynf_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP1_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP1_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRP1_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:02:12
% EndTime: 2019-05-04 19:02:13
% DurationCPUTime: 0.16s
% Computational Cost: add. (1156->49), mult. (1649->61), div. (0->0), fcn. (1008->6), ass. (0->33)
t25 = (qJD(2) + qJD(3));
t23 = t25 ^ 2;
t24 = qJDD(2) + qJDD(3);
t27 = sin(pkin(6));
t28 = cos(pkin(6));
t16 = t27 * g(1) - t28 * g(2);
t17 = -t28 * g(1) - t27 * g(2);
t30 = sin(qJ(2));
t32 = cos(qJ(2));
t36 = t32 * t16 - t30 * t17;
t11 = qJDD(2) * pkin(2) + t36;
t33 = qJD(2) ^ 2;
t39 = t30 * t16 + t32 * t17;
t12 = -t33 * pkin(2) + t39;
t29 = sin(qJ(3));
t31 = cos(qJ(3));
t34 = t31 * t11 - t29 * t12;
t43 = m(5) * (-t24 * pkin(3) - t23 * qJ(4) + qJDD(4) - t34);
t42 = t29 * t11 + t31 * t12;
t41 = t24 * mrSges(5,3) + m(5) * (-t23 * pkin(3) + t24 * qJ(4) + (2 * qJD(4) * t25) + t42);
t40 = mrSges(4,1) + mrSges(5,1);
t26 = -g(3) + qJDD(1);
t19 = m(5) * t26;
t38 = m(4) * t26 + t19;
t37 = m(3) * t26 + t38;
t35 = m(2) * t26 + t37;
t6 = m(4) * t34 - t43 + t40 * t24 + (-mrSges(4,2) + mrSges(5,3)) * t23;
t5 = m(4) * t42 - t24 * mrSges(4,2) - t23 * t40 + t41;
t4 = m(3) * t39 - t33 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t29 * t6 + t31 * t5;
t3 = m(3) * t36 + qJDD(2) * mrSges(3,1) - t33 * mrSges(3,2) + t29 * t5 + t31 * t6;
t2 = m(2) * t17 - t30 * t3 + t32 * t4;
t1 = m(2) * t16 + t32 * t3 + t30 * t4;
t7 = [-m(1) * g(1) - t27 * t1 + t28 * t2, t2, t4, t5, -t23 * mrSges(5,1) + t41; -m(1) * g(2) + t28 * t1 + t27 * t2, t1, t3, t6, t19; -m(1) * g(3) + t35, t35, t37, t38, -t24 * mrSges(5,1) - t23 * mrSges(5,3) + t43;];
f_new  = t7;
