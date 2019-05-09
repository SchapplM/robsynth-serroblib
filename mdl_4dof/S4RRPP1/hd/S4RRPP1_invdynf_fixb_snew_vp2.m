% Calculate vector of cutting forces with Newton-Euler
% S4RRPP1
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
%   pkin=[a2,a3,a4,d1,d2,theta3]';
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
% Datum: 2019-05-04 19:19
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4RRPP1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP1_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP1_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPP1_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP1_invdynf_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP1_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPP1_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPP1_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:19:10
% EndTime: 2019-05-04 19:19:11
% DurationCPUTime: 0.17s
% Computational Cost: add. (1591->57), mult. (2081->67), div. (0->0), fcn. (1008->6), ass. (0->34)
t24 = (qJD(1) + qJD(2));
t22 = t24 ^ 2;
t23 = qJDD(1) + qJDD(2);
t29 = sin(qJ(1));
t31 = cos(qJ(1));
t36 = t29 * g(1) - t31 * g(2);
t16 = qJDD(1) * pkin(1) + t36;
t32 = qJD(1) ^ 2;
t34 = -t31 * g(1) - t29 * g(2);
t17 = -t32 * pkin(1) + t34;
t28 = sin(qJ(2));
t30 = cos(qJ(2));
t35 = t30 * t16 - t28 * t17;
t11 = t23 * pkin(2) + t35;
t38 = t28 * t16 + t30 * t17;
t12 = -t22 * pkin(2) + t38;
t26 = sin(pkin(6));
t27 = cos(pkin(6));
t33 = t27 * t11 - t26 * t12;
t43 = m(5) * (-t23 * pkin(3) - t22 * qJ(4) + qJDD(4) - t33);
t42 = -m(2) - m(3);
t41 = t26 * t11 + t27 * t12;
t40 = t23 * mrSges(5,3) + m(5) * (-t22 * pkin(3) + t23 * qJ(4) + (2 * qJD(4) * t24) + t41);
t39 = mrSges(4,1) + mrSges(5,1);
t25 = -g(3) + qJDD(3);
t19 = m(5) * t25;
t37 = m(4) * t25 + t19;
t6 = m(4) * t33 - t43 + t39 * t23 + (-mrSges(4,2) + mrSges(5,3)) * t22;
t5 = m(4) * t41 - t23 * mrSges(4,2) - t39 * t22 + t40;
t4 = m(3) * t38 - t22 * mrSges(3,1) - t23 * mrSges(3,2) - t26 * t6 + t27 * t5;
t3 = m(3) * t35 + t23 * mrSges(3,1) - t22 * mrSges(3,2) + t26 * t5 + t27 * t6;
t2 = m(2) * t34 - t32 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t28 * t3 + t30 * t4;
t1 = m(2) * t36 + qJDD(1) * mrSges(2,1) - t32 * mrSges(2,2) + t28 * t4 + t30 * t3;
t7 = [-m(1) * g(1) - t29 * t1 + t31 * t2, t2, t4, t5, -t22 * mrSges(5,1) + t40; -m(1) * g(2) + t31 * t1 + t29 * t2, t1, t3, t6, t19; (-m(1) + t42) * g(3) + t37, t42 * g(3) + t37, -m(3) * g(3) + t37, t37, -t23 * mrSges(5,1) - t22 * mrSges(5,3) + t43;];
f_new  = t7;
