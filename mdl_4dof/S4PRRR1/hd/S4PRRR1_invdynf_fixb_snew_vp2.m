% Calculate vector of cutting forces with Newton-Euler
% S4PRRR1
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
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2019-05-04 19:05
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4PRRR1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR1_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR1_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR1_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR1_invdynf_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR1_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR1_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR1_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:04:23
% EndTime: 2019-05-04 19:04:23
% DurationCPUTime: 0.21s
% Computational Cost: add. (2392->49), mult. (3575->66), div. (0->0), fcn. (2480->8), ass. (0->39)
t32 = sin(pkin(7));
t33 = cos(pkin(7));
t19 = t32 * g(1) - t33 * g(2);
t20 = -t33 * g(1) - t32 * g(2);
t36 = sin(qJ(2));
t39 = cos(qJ(2));
t42 = t39 * t19 - t36 * t20;
t14 = qJDD(2) * pkin(2) + t42;
t40 = qJD(2) ^ 2;
t46 = t36 * t19 + t39 * t20;
t15 = -t40 * pkin(2) + t46;
t35 = sin(qJ(3));
t38 = cos(qJ(3));
t47 = t35 * t14 + t38 * t15;
t31 = -g(3) + qJDD(1);
t24 = m(5) * t31;
t45 = m(4) * t31 + t24;
t30 = qJD(2) + qJD(3);
t29 = qJDD(2) + qJDD(3);
t44 = m(3) * t31 + t45;
t43 = t38 * t14 - t35 * t15;
t41 = m(2) * t31 + t44;
t37 = cos(qJ(4));
t34 = sin(qJ(4));
t28 = t30 ^ 2;
t23 = qJD(4) + t30;
t22 = qJDD(4) + t29;
t21 = t23 ^ 2;
t10 = -t28 * pkin(3) + t47;
t9 = t29 * pkin(3) + t43;
t8 = m(5) * (t37 * t10 + t34 * t9) - t22 * mrSges(5,2) - t21 * mrSges(5,1);
t7 = m(5) * (-t34 * t10 + t37 * t9) + t22 * mrSges(5,1) - t21 * mrSges(5,2);
t6 = m(4) * t47 - t28 * mrSges(4,1) - t29 * mrSges(4,2) - t34 * t7 + t37 * t8;
t5 = m(4) * t43 + t29 * mrSges(4,1) - t28 * mrSges(4,2) + t34 * t8 + t37 * t7;
t4 = m(3) * t46 - t40 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t35 * t5 + t38 * t6;
t3 = m(3) * t42 + qJDD(2) * mrSges(3,1) - t40 * mrSges(3,2) + t35 * t6 + t38 * t5;
t2 = m(2) * t20 - t36 * t3 + t39 * t4;
t1 = m(2) * t19 + t39 * t3 + t36 * t4;
t11 = [-m(1) * g(1) - t32 * t1 + t33 * t2, t2, t4, t6, t8; -m(1) * g(2) + t33 * t1 + t32 * t2, t1, t3, t5, t7; -m(1) * g(3) + t41, t41, t44, t45, t24;];
f_new  = t11;
