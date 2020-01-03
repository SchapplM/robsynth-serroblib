% Calculate vector of cutting forces with Newton-Euler
% S4RRPR10
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
%   pkin=[a2,a3,a4,d1,d2,d4]';
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
% Datum: 2019-12-31 17:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4RRPR10_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR10_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR10_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR10_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR10_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR10_invdynf_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR10_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR10_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR10_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:11:11
% EndTime: 2019-12-31 17:11:12
% DurationCPUTime: 0.47s
% Computational Cost: add. (2093->123), mult. (4286->154), div. (0->0), fcn. (2115->6), ass. (0->61)
t53 = sin(qJ(2));
t56 = cos(qJ(2));
t73 = qJD(1) * qJD(2);
t69 = t56 * t73;
t38 = t53 * qJDD(1) + t69;
t70 = t53 * t73;
t39 = t56 * qJDD(1) - t70;
t75 = qJD(1) * t56;
t41 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t75;
t54 = sin(qJ(1));
t57 = cos(qJ(1));
t72 = t54 * g(1) - t57 * g(2);
t66 = -qJDD(1) * pkin(1) - t72;
t42 = -mrSges(4,1) * t75 - qJD(2) * mrSges(4,3);
t52 = sin(qJ(4));
t55 = cos(qJ(4));
t74 = t53 * qJD(1);
t84 = -2 * qJD(3);
t60 = pkin(2) * t70 + t74 * t84 + (-t38 - t69) * qJ(3) + t66;
t44 = pkin(3) * t74 - qJD(2) * pkin(6);
t51 = t56 ^ 2;
t59 = qJD(1) ^ 2;
t10 = -t44 * t74 + (-pkin(3) * t51 - pkin(5)) * t59 + (-pkin(2) - pkin(6)) * t39 + t60;
t35 = (-pkin(2) * t56 - qJ(3) * t53) * qJD(1);
t58 = qJD(2) ^ 2;
t68 = -t57 * g(1) - t54 * g(2);
t29 = -t59 * pkin(1) + qJDD(1) * pkin(5) + t68;
t78 = -t56 * g(3) - t53 * t29;
t17 = -qJDD(2) * pkin(2) - t58 * qJ(3) + t35 * t74 + qJDD(3) - t78;
t13 = (-t53 * t56 * t59 - qJDD(2)) * pkin(6) + (t38 - t69) * pkin(3) + t17;
t33 = -t52 * qJD(2) - t55 * t75;
t21 = t33 * qJD(4) + t55 * qJDD(2) - t52 * t39;
t34 = t55 * qJD(2) - t52 * t75;
t22 = -t33 * mrSges(5,1) + t34 * mrSges(5,2);
t46 = qJD(4) + t74;
t23 = -t46 * mrSges(5,2) + t33 * mrSges(5,3);
t32 = qJDD(4) + t38;
t8 = m(5) * (-t52 * t10 + t55 * t13) - t21 * mrSges(5,3) + t32 * mrSges(5,1) - t34 * t22 + t46 * t23;
t82 = t59 * pkin(5);
t20 = -t34 * qJD(4) - t52 * qJDD(2) - t55 * t39;
t24 = t46 * mrSges(5,1) - t34 * mrSges(5,3);
t9 = m(5) * (t55 * t10 + t52 * t13) + t20 * mrSges(5,3) - t32 * mrSges(5,2) + t33 * t22 - t46 * t24;
t67 = t52 * t8 - t55 * t9 - m(4) * (-t39 * pkin(2) + t60 - t82) - t42 * t75 + t38 * mrSges(4,3);
t43 = mrSges(4,1) * t74 + qJD(2) * mrSges(4,2);
t76 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t74 - t43;
t80 = mrSges(3,1) - mrSges(4,2);
t87 = (-t56 * t41 + t76 * t53) * qJD(1) - t80 * t39 + m(3) * (t66 - t82) + t38 * mrSges(3,2) - t67;
t65 = -m(4) * t17 - t52 * t9 - t55 * t8;
t36 = (mrSges(4,2) * t56 - mrSges(4,3) * t53) * qJD(1);
t77 = t36 + (-mrSges(3,1) * t56 + mrSges(3,2) * t53) * qJD(1);
t79 = mrSges(3,3) + mrSges(4,1);
t4 = m(3) * t78 - t79 * t38 + t80 * qJDD(2) + (t41 - t42) * qJD(2) - t77 * t74 + t65;
t71 = -t53 * g(3) + t56 * t29;
t61 = -t58 * pkin(2) + qJDD(2) * qJ(3) + t35 * t75 + t71;
t64 = -t20 * mrSges(5,1) - t33 * t23 + m(5) * (-t51 * t59 * pkin(6) + t39 * pkin(3) + ((2 * qJD(3)) + t44) * qJD(2) + t61) + t21 * mrSges(5,2) + t34 * t24;
t63 = -m(4) * (qJD(2) * t84 - t61) + t64;
t6 = m(3) * t71 + t79 * t39 + (-mrSges(3,2) + mrSges(4,3)) * qJDD(2) - t76 * qJD(2) + t77 * t75 + t63;
t83 = t56 * t4 + t53 * t6;
t2 = m(2) * t72 + qJDD(1) * mrSges(2,1) - t59 * mrSges(2,2) - t87;
t1 = m(2) * t68 - t59 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t53 * t4 + t56 * t6;
t3 = [-m(1) * g(1) + t57 * t1 - t54 * t2, t1, t6, t39 * mrSges(4,2) - t43 * t74 - t67, t9; -m(1) * g(2) + t54 * t1 + t57 * t2, t2, t4, -t39 * mrSges(4,1) - qJDD(2) * mrSges(4,3) - qJD(2) * t43 - t36 * t75 - t63, t8; (-m(1) - m(2)) * g(3) + t83, -m(2) * g(3) + t83, t87, t38 * mrSges(4,1) + qJDD(2) * mrSges(4,2) + qJD(2) * t42 + t36 * t74 - t65, t64;];
f_new = t3;
