% Calculate vector of cutting forces with Newton-Euler
% S4RRPR3
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
% Datum: 2019-12-31 17:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4RRPR3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR3_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR3_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR3_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR3_invdynf_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR3_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR3_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR3_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:01:30
% EndTime: 2019-12-31 17:01:30
% DurationCPUTime: 0.27s
% Computational Cost: add. (2741->70), mult. (3602->96), div. (0->0), fcn. (1830->8), ass. (0->44)
t30 = qJDD(1) + qJDD(2);
t35 = sin(qJ(4));
t38 = cos(qJ(4));
t31 = qJD(1) + qJD(2);
t49 = qJD(4) * t31;
t18 = t35 * t30 + t38 * t49;
t19 = t38 * t30 - t35 * t49;
t53 = t31 * t35;
t25 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t53;
t52 = t31 * t38;
t26 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t52;
t29 = t31 ^ 2;
t37 = sin(qJ(1));
t40 = cos(qJ(1));
t47 = t37 * g(1) - t40 * g(2);
t23 = qJDD(1) * pkin(1) + t47;
t41 = qJD(1) ^ 2;
t45 = -t40 * g(1) - t37 * g(2);
t24 = -t41 * pkin(1) + t45;
t36 = sin(qJ(2));
t39 = cos(qJ(2));
t46 = t39 * t23 - t36 * t24;
t15 = t30 * pkin(2) + t46;
t50 = t36 * t23 + t39 * t24;
t16 = -t29 * pkin(2) + t50;
t33 = sin(pkin(7));
t34 = cos(pkin(7));
t44 = t34 * t15 - t33 * t16;
t55 = (t35 * t25 - t38 * t26) * t31 + m(5) * (-t30 * pkin(3) - t29 * pkin(6) - t44) - t19 * mrSges(5,1) + t18 * mrSges(5,2);
t54 = -m(2) - m(3);
t51 = t33 * t15 + t34 * t16;
t12 = -t29 * pkin(3) + t30 * pkin(6) + t51;
t17 = (-mrSges(5,1) * t38 + mrSges(5,2) * t35) * t31;
t32 = -g(3) + qJDD(3);
t10 = m(5) * (t38 * t12 + t35 * t32) + t19 * mrSges(5,3) - qJDD(4) * mrSges(5,2) + t17 * t52 - qJD(4) * t25;
t9 = m(5) * (-t35 * t12 + t38 * t32) - t18 * mrSges(5,3) + qJDD(4) * mrSges(5,1) - t17 * t53 + qJD(4) * t26;
t48 = m(4) * t32 + t35 * t10 + t38 * t9;
t6 = m(4) * t44 + t30 * mrSges(4,1) - t29 * mrSges(4,2) - t55;
t5 = m(4) * t51 - t29 * mrSges(4,1) - t30 * mrSges(4,2) + t38 * t10 - t35 * t9;
t4 = m(3) * t50 - t29 * mrSges(3,1) - t30 * mrSges(3,2) - t33 * t6 + t34 * t5;
t3 = m(3) * t46 + t30 * mrSges(3,1) - t29 * mrSges(3,2) + t33 * t5 + t34 * t6;
t2 = m(2) * t45 - t41 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t36 * t3 + t39 * t4;
t1 = m(2) * t47 + qJDD(1) * mrSges(2,1) - t41 * mrSges(2,2) + t39 * t3 + t36 * t4;
t7 = [-m(1) * g(1) - t37 * t1 + t40 * t2, t2, t4, t5, t10; -m(1) * g(2) + t40 * t1 + t37 * t2, t1, t3, t6, t9; (-m(1) + t54) * g(3) + t48, g(3) * t54 + t48, -m(3) * g(3) + t48, t48, t55;];
f_new = t7;
