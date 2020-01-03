% Calculate vector of cutting forces with Newton-Euler
% S4RRRR1
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
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 17:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4RRRR1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR1_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR1_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR1_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR1_invdynf_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR1_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR1_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR1_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:22:09
% EndTime: 2019-12-31 17:22:10
% DurationCPUTime: 0.29s
% Computational Cost: add. (3237->72), mult. (3602->97), div. (0->0), fcn. (1830->8), ass. (0->47)
t32 = qJDD(1) + qJDD(2);
t28 = qJDD(3) + t32;
t34 = sin(qJ(4));
t38 = cos(qJ(4));
t33 = qJD(1) + qJD(2);
t29 = qJD(3) + t33;
t50 = qJD(4) * t29;
t18 = t34 * t28 + t38 * t50;
t19 = t38 * t28 - t34 * t50;
t54 = t29 * t34;
t23 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t54;
t53 = t29 * t38;
t24 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t53;
t27 = t29 ^ 2;
t37 = sin(qJ(1));
t41 = cos(qJ(1));
t48 = t37 * g(1) - t41 * g(2);
t25 = qJDD(1) * pkin(1) + t48;
t42 = qJD(1) ^ 2;
t46 = -t41 * g(1) - t37 * g(2);
t26 = -t42 * pkin(1) + t46;
t36 = sin(qJ(2));
t40 = cos(qJ(2));
t47 = t40 * t25 - t36 * t26;
t15 = t32 * pkin(2) + t47;
t31 = t33 ^ 2;
t51 = t36 * t25 + t40 * t26;
t16 = -t31 * pkin(2) + t51;
t35 = sin(qJ(3));
t39 = cos(qJ(3));
t45 = t39 * t15 - t35 * t16;
t57 = (t34 * t23 - t38 * t24) * t29 + m(5) * (-t28 * pkin(3) - t27 * pkin(7) - t45) - t19 * mrSges(5,1) + t18 * mrSges(5,2);
t56 = -m(3) - m(4);
t52 = t35 * t15 + t39 * t16;
t12 = -t27 * pkin(3) + t28 * pkin(7) + t52;
t17 = (-mrSges(5,1) * t38 + mrSges(5,2) * t34) * t29;
t10 = m(5) * (-t34 * g(3) + t38 * t12) + t19 * mrSges(5,3) - qJDD(4) * mrSges(5,2) + t17 * t53 - qJD(4) * t23;
t9 = m(5) * (-t38 * g(3) - t34 * t12) - t18 * mrSges(5,3) + qJDD(4) * mrSges(5,1) - t17 * t54 + qJD(4) * t24;
t55 = t34 * t10 + t38 * t9;
t49 = -m(2) + t56;
t6 = m(4) * t45 + t28 * mrSges(4,1) - t27 * mrSges(4,2) - t57;
t5 = m(4) * t52 - t27 * mrSges(4,1) - t28 * mrSges(4,2) + t38 * t10 - t34 * t9;
t4 = m(3) * t51 - t31 * mrSges(3,1) - t32 * mrSges(3,2) - t35 * t6 + t39 * t5;
t3 = m(3) * t47 + t32 * mrSges(3,1) - t31 * mrSges(3,2) + t35 * t5 + t39 * t6;
t2 = m(2) * t46 - t42 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t36 * t3 + t40 * t4;
t1 = m(2) * t48 + qJDD(1) * mrSges(2,1) - t42 * mrSges(2,2) + t40 * t3 + t36 * t4;
t7 = [-m(1) * g(1) - t37 * t1 + t41 * t2, t2, t4, t5, t10; -m(1) * g(2) + t41 * t1 + t37 * t2, t1, t3, t6, t9; (-m(1) + t49) * g(3) + t55, t49 * g(3) + t55, t56 * g(3) + t55, -m(4) * g(3) + t55, t57;];
f_new = t7;
