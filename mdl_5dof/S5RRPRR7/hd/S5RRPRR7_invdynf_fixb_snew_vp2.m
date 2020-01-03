% Calculate vector of cutting forces with Newton-Euler
% S5RRPRR7
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% f_new [3x6]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRPRR7_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR7_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR7_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR7_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR7_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR7_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR7_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR7_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR7_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:15:30
% EndTime: 2019-12-31 20:15:31
% DurationCPUTime: 0.51s
% Computational Cost: add. (6250->113), mult. (7591->145), div. (0->0), fcn. (4106->8), ass. (0->64)
t50 = qJDD(1) + qJDD(2);
t52 = (qJD(1) + qJD(2));
t57 = sin(qJ(1));
t61 = cos(qJ(1));
t72 = t57 * g(1) - t61 * g(2);
t39 = qJDD(1) * pkin(1) + t72;
t62 = qJD(1) ^ 2;
t68 = -t61 * g(1) - t57 * g(2);
t40 = -t62 * pkin(1) + t68;
t56 = sin(qJ(2));
t60 = cos(qJ(2));
t80 = t56 * t39 + t60 * t40;
t69 = t50 * qJ(3) + (2 * qJD(3) * t52) + t80;
t86 = -m(3) - m(4);
t85 = -pkin(2) - pkin(7);
t55 = sin(qJ(4));
t84 = t52 * t55;
t59 = cos(qJ(4));
t83 = t52 * t59;
t82 = (mrSges(3,1) - mrSges(4,2));
t81 = -mrSges(3,2) + mrSges(4,3);
t48 = t52 ^ 2;
t70 = t60 * t39 - t56 * t40;
t66 = -t48 * qJ(3) + qJDD(3) - t70;
t21 = t85 * t50 + t66;
t79 = t55 * g(3) + t59 * t21;
t77 = qJD(4) * t52;
t76 = -m(2) + t86;
t74 = t55 * t77;
t33 = (mrSges(5,1) * t55 + mrSges(5,2) * t59) * t52;
t35 = t59 * t50 - t74;
t41 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t84;
t54 = sin(qJ(5));
t58 = cos(qJ(5));
t11 = (-t35 - t74) * pkin(8) + (-t48 * t55 * t59 + qJDD(4)) * pkin(4) + t79;
t34 = -t55 * t50 - t59 * t77;
t43 = qJD(4) * pkin(4) - pkin(8) * t83;
t53 = t55 ^ 2;
t71 = -t59 * g(3) + t55 * t21;
t13 = -t53 * t48 * pkin(4) + t34 * pkin(8) - qJD(4) * t43 + t71;
t28 = (-t54 * t59 - t55 * t58) * t52;
t16 = t28 * qJD(5) + t54 * t34 + t58 * t35;
t29 = (-t54 * t55 + t58 * t59) * t52;
t25 = -t28 * mrSges(6,1) + t29 * mrSges(6,2);
t51 = qJD(4) + qJD(5);
t26 = -t51 * mrSges(6,2) + t28 * mrSges(6,3);
t49 = qJDD(4) + qJDD(5);
t8 = m(6) * (t58 * t11 - t54 * t13) - t16 * mrSges(6,3) + t49 * mrSges(6,1) - t29 * t25 + t51 * t26;
t15 = -t29 * qJD(5) + t58 * t34 - t54 * t35;
t27 = t51 * mrSges(6,1) - t29 * mrSges(6,3);
t9 = m(6) * (t54 * t11 + t58 * t13) + t15 * mrSges(6,3) - t49 * mrSges(6,2) + t28 * t25 - t51 * t27;
t5 = m(5) * t79 + qJDD(4) * mrSges(5,1) - t35 * mrSges(5,3) + qJD(4) * t41 - t33 * t83 + t54 * t9 + t58 * t8;
t42 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t83;
t6 = m(5) * t71 - qJDD(4) * mrSges(5,2) + t34 * mrSges(5,3) - qJD(4) * t42 - t33 * t84 - t54 * t8 + t58 * t9;
t73 = -t55 * t5 + t59 * t6;
t67 = -m(4) * (-t50 * pkin(2) + t66) - t59 * t5 - t55 * t6;
t65 = -t15 * mrSges(6,1) - t28 * t26 + m(6) * (t43 * t83 - t34 * pkin(4) + (-pkin(8) * t53 + t85) * t48 + t69) + t16 * mrSges(6,2) + t29 * t27;
t64 = -t34 * mrSges(5,1) + m(5) * (t85 * t48 + t69) + t41 * t84 + t42 * t83 + t35 * mrSges(5,2) + t65;
t63 = -m(4) * (t48 * pkin(2) - t69) + t64;
t7 = m(3) * t80 - (t82 * t48) + t81 * t50 + t63;
t3 = m(3) * t70 + t81 * t48 + t82 * t50 + t67;
t2 = m(2) * t68 - t62 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t56 * t3 + t60 * t7;
t1 = m(2) * t72 + qJDD(1) * mrSges(2,1) - t62 * mrSges(2,2) + t60 * t3 + t56 * t7;
t4 = [-m(1) * g(1) - t57 * t1 + t61 * t2, t2, t7, -m(4) * g(3) + t73, t6, t9; -m(1) * g(2) + t61 * t1 + t57 * t2, t1, t3, -(t48 * mrSges(4,2)) - t50 * mrSges(4,3) - t63, t5, t8; (-m(1) + t76) * g(3) + t73, t76 * g(3) + t73, t86 * g(3) + t73, t50 * mrSges(4,2) - t48 * mrSges(4,3) - t67, t64, t65;];
f_new = t4;
