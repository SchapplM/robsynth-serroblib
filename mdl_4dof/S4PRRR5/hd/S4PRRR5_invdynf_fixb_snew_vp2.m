% Calculate vector of cutting forces with Newton-Euler
% S4PRRR5
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
% Datum: 2019-12-31 16:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4PRRR5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR5_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR5_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR5_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR5_invdynf_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR5_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR5_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR5_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:33:36
% EndTime: 2019-12-31 16:33:37
% DurationCPUTime: 0.26s
% Computational Cost: add. (1980->64), mult. (2520->90), div. (0->0), fcn. (1424->8), ass. (0->42)
t33 = qJDD(2) + qJDD(3);
t37 = sin(qJ(4));
t40 = cos(qJ(4));
t34 = qJD(2) + qJD(3);
t50 = qJD(4) * t34;
t20 = t37 * t33 + t40 * t50;
t21 = t40 * t33 - t37 * t50;
t55 = t34 * t37;
t23 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t55;
t54 = t34 * t40;
t24 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t54;
t32 = t34 ^ 2;
t36 = sin(pkin(7));
t51 = cos(pkin(7));
t27 = -t51 * g(1) - t36 * g(2);
t35 = -g(3) + qJDD(1);
t39 = sin(qJ(2));
t42 = cos(qJ(2));
t47 = -t39 * t27 + t42 * t35;
t17 = qJDD(2) * pkin(2) + t47;
t43 = qJD(2) ^ 2;
t52 = t42 * t27 + t39 * t35;
t18 = -t43 * pkin(2) + t52;
t38 = sin(qJ(3));
t41 = cos(qJ(3));
t46 = t41 * t17 - t38 * t18;
t56 = (t37 * t23 - t40 * t24) * t34 + m(5) * (-t33 * pkin(3) - t32 * pkin(6) - t46) - t21 * mrSges(5,1) + t20 * mrSges(5,2);
t53 = t38 * t17 + t41 * t18;
t14 = -t32 * pkin(3) + t33 * pkin(6) + t53;
t19 = (-mrSges(5,1) * t40 + mrSges(5,2) * t37) * t34;
t26 = t36 * g(1) - t51 * g(2);
t11 = m(5) * (-t37 * t14 - t40 * t26) - t20 * mrSges(5,3) + qJDD(4) * mrSges(5,1) - t19 * t55 + qJD(4) * t24;
t12 = m(5) * (t40 * t14 - t37 * t26) + t21 * mrSges(5,3) - qJDD(4) * mrSges(5,2) + t19 * t54 - qJD(4) * t23;
t6 = m(4) * t53 - t32 * mrSges(4,1) - t33 * mrSges(4,2) - t37 * t11 + t40 * t12;
t8 = m(4) * t46 + t33 * mrSges(4,1) - t32 * mrSges(4,2) - t56;
t4 = m(3) * t47 + qJDD(2) * mrSges(3,1) - t43 * mrSges(3,2) + t38 * t6 + t41 * t8;
t5 = m(3) * t52 - t43 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t38 * t8 + t41 * t6;
t49 = m(2) * t35 + t39 * t5 + t42 * t4;
t48 = m(4) * t26 - t40 * t11 - t37 * t12;
t7 = (m(2) + m(3)) * t26 + t48;
t1 = m(2) * t27 - t39 * t4 + t42 * t5;
t2 = [-m(1) * g(1) + t51 * t1 - t36 * t7, t1, t5, t6, t12; -m(1) * g(2) + t36 * t1 + t51 * t7, t7, t4, t8, t11; -m(1) * g(3) + t49, t49, -m(3) * t26 - t48, -t48, t56;];
f_new = t2;
