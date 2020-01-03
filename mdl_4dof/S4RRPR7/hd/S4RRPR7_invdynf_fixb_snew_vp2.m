% Calculate vector of cutting forces with Newton-Euler
% S4RRPR7
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
% Datum: 2019-12-31 17:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4RRPR7_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR7_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR7_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR7_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR7_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR7_invdynf_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR7_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR7_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR7_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:05:59
% EndTime: 2019-12-31 17:06:00
% DurationCPUTime: 0.64s
% Computational Cost: add. (5140->125), mult. (11622->172), div. (0->0), fcn. (7329->8), ass. (0->65)
t54 = sin(pkin(7));
t55 = cos(pkin(7));
t57 = sin(qJ(2));
t60 = cos(qJ(2));
t41 = (t54 * t57 - t55 * t60) * qJD(1);
t76 = qJD(1) * qJD(2);
t47 = t57 * qJDD(1) + t60 * t76;
t48 = t60 * qJDD(1) - t57 * t76;
t78 = qJD(1) * t57;
t50 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t78;
t77 = qJD(1) * t60;
t51 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t77;
t63 = qJD(1) ^ 2;
t42 = (t54 * t60 + t55 * t57) * qJD(1);
t29 = t41 * pkin(3) - t42 * pkin(6);
t62 = qJD(2) ^ 2;
t58 = sin(qJ(1));
t61 = cos(qJ(1));
t72 = -t61 * g(1) - t58 * g(2);
t44 = -t63 * pkin(1) + qJDD(1) * pkin(5) + t72;
t79 = t57 * t44;
t80 = pkin(2) * t63;
t20 = qJDD(2) * pkin(2) - t47 * qJ(3) - t79 + (qJ(3) * t76 + t57 * t80 - g(3)) * t60;
t49 = qJD(2) * pkin(2) - qJ(3) * t78;
t53 = t60 ^ 2;
t73 = -t57 * g(3) + t60 * t44;
t21 = t48 * qJ(3) - qJD(2) * t49 - t53 * t80 + t73;
t82 = 2 * qJD(3);
t75 = t54 * t20 + t55 * t21 - t41 * t82;
t14 = -t62 * pkin(3) + qJDD(2) * pkin(6) - t41 * t29 + t75;
t32 = -t54 * t47 + t55 * t48;
t33 = t55 * t47 + t54 * t48;
t74 = t58 * g(1) - t61 * g(2);
t68 = -qJDD(1) * pkin(1) - t74;
t65 = -t48 * pkin(2) + qJDD(3) + t49 * t78 + (-qJ(3) * t53 - pkin(5)) * t63 + t68;
t15 = (qJD(2) * t41 - t33) * pkin(6) + (qJD(2) * t42 - t32) * pkin(3) + t65;
t56 = sin(qJ(4));
t59 = cos(qJ(4));
t34 = t59 * qJD(2) - t56 * t42;
t17 = t34 * qJD(4) + t56 * qJDD(2) + t59 * t33;
t35 = t56 * qJD(2) + t59 * t42;
t22 = -t34 * mrSges(5,1) + t35 * mrSges(5,2);
t40 = qJD(4) + t41;
t25 = -t40 * mrSges(5,2) + t34 * mrSges(5,3);
t31 = qJDD(4) - t32;
t11 = m(5) * (-t56 * t14 + t59 * t15) - t17 * mrSges(5,3) + t31 * mrSges(5,1) - t35 * t22 + t40 * t25;
t16 = -t35 * qJD(4) + t59 * qJDD(2) - t56 * t33;
t26 = t40 * mrSges(5,1) - t35 * mrSges(5,3);
t12 = m(5) * (t59 * t14 + t56 * t15) + t16 * mrSges(5,3) - t31 * mrSges(5,2) + t34 * t22 - t40 * t26;
t36 = -qJD(2) * mrSges(4,2) - t41 * mrSges(4,3);
t37 = qJD(2) * mrSges(4,1) - t42 * mrSges(4,3);
t67 = -m(4) * t65 + t32 * mrSges(4,1) - t33 * mrSges(4,2) - t59 * t11 - t56 * t12 - t41 * t36 - t42 * t37;
t83 = (t57 * t50 - t60 * t51) * qJD(1) + m(3) * (-t63 * pkin(5) + t68) - t48 * mrSges(3,1) + t47 * mrSges(3,2) - t67;
t46 = (-mrSges(3,1) * t60 + mrSges(3,2) * t57) * qJD(1);
t28 = t41 * mrSges(4,1) + t42 * mrSges(4,2);
t7 = m(4) * t75 - qJDD(2) * mrSges(4,2) + t32 * mrSges(4,3) - qJD(2) * t37 - t56 * t11 + t59 * t12 - t41 * t28;
t71 = -t55 * t20 + t54 * t21;
t66 = m(5) * (-qJDD(2) * pkin(3) - t62 * pkin(6) + (t82 + t29) * t42 + t71) - t16 * mrSges(5,1) + t17 * mrSges(5,2) - t34 * t25 + t35 * t26;
t8 = m(4) * (-0.2e1 * qJD(3) * t42 - t71) - t33 * mrSges(4,3) + qJDD(2) * mrSges(4,1) - t42 * t28 + qJD(2) * t36 - t66;
t4 = m(3) * (-t60 * g(3) - t79) - t47 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t46 * t78 + qJD(2) * t51 + t54 * t7 + t55 * t8;
t5 = m(3) * t73 - qJDD(2) * mrSges(3,2) + t48 * mrSges(3,3) - qJD(2) * t50 + t46 * t77 - t54 * t8 + t55 * t7;
t81 = t60 * t4 + t57 * t5;
t6 = m(2) * t74 + qJDD(1) * mrSges(2,1) - t63 * mrSges(2,2) - t83;
t1 = m(2) * t72 - t63 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t57 * t4 + t60 * t5;
t2 = [-m(1) * g(1) + t61 * t1 - t58 * t6, t1, t5, t7, t12; -m(1) * g(2) + t58 * t1 + t61 * t6, t6, t4, t8, t11; (-m(1) - m(2)) * g(3) + t81, -m(2) * g(3) + t81, t83, -t67, t66;];
f_new = t2;
