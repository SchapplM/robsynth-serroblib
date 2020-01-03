% Calculate vector of cutting forces with Newton-Euler
% S4RRRP7
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
% Datum: 2019-12-31 17:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4RRRP7_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP7_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP7_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP7_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP7_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP7_invdynf_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP7_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP7_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP7_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:20:10
% EndTime: 2019-12-31 17:20:11
% DurationCPUTime: 0.42s
% Computational Cost: add. (2763->119), mult. (5402->150), div. (0->0), fcn. (3071->6), ass. (0->58)
t57 = qJD(1) ^ 2;
t53 = sin(qJ(1));
t55 = cos(qJ(1));
t66 = t53 * g(1) - t55 * g(2);
t33 = -qJDD(1) * pkin(1) - t57 * pkin(5) - t66;
t52 = sin(qJ(2));
t54 = cos(qJ(2));
t69 = qJD(1) * qJD(2);
t63 = t54 * t69;
t41 = t52 * qJDD(1) + t63;
t64 = t52 * t69;
t42 = t54 * qJDD(1) - t64;
t71 = qJD(1) * t52;
t43 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t71;
t70 = t54 * qJD(1);
t44 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t70;
t51 = sin(qJ(3));
t78 = cos(qJ(3));
t38 = t51 * qJD(2) + t78 * t71;
t20 = t38 * qJD(3) - t78 * qJDD(2) + t51 * t41;
t46 = qJD(3) - t70;
t27 = mrSges(4,1) * t46 - mrSges(4,3) * t38;
t36 = qJDD(3) - t42;
t37 = -t78 * qJD(2) + t51 * t71;
t23 = pkin(3) * t37 - qJ(4) * t38;
t28 = -mrSges(5,1) * t46 + mrSges(5,2) * t38;
t45 = t46 ^ 2;
t14 = (-t41 - t63) * pkin(6) + (-t42 + t64) * pkin(2) + t33;
t40 = (-pkin(2) * t54 - pkin(6) * t52) * qJD(1);
t56 = qJD(2) ^ 2;
t62 = -t55 * g(1) - t53 * g(2);
t34 = -t57 * pkin(1) + qJDD(1) * pkin(5) + t62;
t65 = -t52 * g(3) + t54 * t34;
t17 = -t56 * pkin(2) + qJDD(2) * pkin(6) + t40 * t70 + t65;
t75 = t51 * t14 + t78 * t17;
t68 = t46 * t28 + t36 * mrSges(5,3) + m(5) * (-pkin(3) * t45 + qJ(4) * t36 + 0.2e1 * qJD(4) * t46 - t37 * t23 + t75);
t24 = mrSges(5,1) * t37 - mrSges(5,3) * t38;
t74 = -mrSges(4,1) * t37 - mrSges(4,2) * t38 - t24;
t76 = -mrSges(4,3) - mrSges(5,2);
t7 = m(4) * t75 - t36 * mrSges(4,2) + t76 * t20 - t46 * t27 + t74 * t37 + t68;
t21 = -t37 * qJD(3) + t51 * qJDD(2) + t78 * t41;
t26 = -mrSges(4,2) * t46 - mrSges(4,3) * t37;
t29 = -mrSges(5,2) * t37 + mrSges(5,3) * t46;
t60 = t78 * t14 - t51 * t17;
t79 = m(5) * (-t36 * pkin(3) - t45 * qJ(4) + t38 * t23 + qJDD(4) - t60);
t8 = m(4) * t60 - t79 + (t26 + t29) * t46 + t74 * t38 + (mrSges(4,1) + mrSges(5,1)) * t36 + t76 * t21;
t82 = m(3) * t33 - t42 * mrSges(3,1) + t41 * mrSges(3,2) + (t43 * t52 - t44 * t54) * qJD(1) + t51 * t7 + t78 * t8;
t72 = -t54 * g(3) - t52 * t34;
t16 = -qJDD(2) * pkin(2) - t56 * pkin(6) + t40 * t71 - t72;
t67 = m(5) * (-0.2e1 * qJD(4) * t38 + (t37 * t46 - t21) * qJ(4) + (t38 * t46 + t20) * pkin(3) + t16) + t20 * mrSges(5,1) + t37 * t29;
t81 = m(4) * t16 + t20 * mrSges(4,1) + (mrSges(4,2) - mrSges(5,3)) * t21 + t37 * t26 + (t27 - t28) * t38 + t67;
t39 = (-mrSges(3,1) * t54 + mrSges(3,2) * t52) * qJD(1);
t4 = m(3) * t65 - qJDD(2) * mrSges(3,2) + t42 * mrSges(3,3) - qJD(2) * t43 + t39 * t70 - t51 * t8 + t78 * t7;
t6 = m(3) * t72 + qJDD(2) * mrSges(3,1) - t41 * mrSges(3,3) + qJD(2) * t44 - t39 * t71 - t81;
t80 = t52 * t4 + t54 * t6;
t2 = m(2) * t66 + qJDD(1) * mrSges(2,1) - t57 * mrSges(2,2) - t82;
t1 = m(2) * t62 - t57 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t54 * t4 - t52 * t6;
t3 = [-m(1) * g(1) + t55 * t1 - t53 * t2, t1, t4, t7, -t20 * mrSges(5,2) - t37 * t24 + t68; -m(1) * g(2) + t53 * t1 + t55 * t2, t2, t6, t8, -t21 * mrSges(5,3) - t38 * t28 + t67; (-m(1) - m(2)) * g(3) + t80, -m(2) * g(3) + t80, t82, t81, -t36 * mrSges(5,1) + t21 * mrSges(5,2) + t38 * t24 - t46 * t29 + t79;];
f_new = t3;
