% Calculate vector of cutting forces with Newton-Euler
% S6RRPRRR14V3
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% f_new [3x7]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-07 03:56
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPRRR14V3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(1,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14V3_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14V3_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR14V3_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR14V3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S6RRPRRR14V3_invdynf_fixb_snew_vp2: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR14V3_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR14V3_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR14V3_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 03:48:22
% EndTime: 2019-05-07 03:48:26
% DurationCPUTime: 1.31s
% Computational Cost: add. (14744->163), mult. (27885->209), div. (0->0), fcn. (20524->10), ass. (0->81)
t77 = sin(qJ(2));
t82 = cos(qJ(2));
t95 = qJD(1) * qJD(2);
t93 = t82 * t95;
t61 = t77 * qJDD(1) + t93;
t62 = t82 * qJDD(1) - t77 * t95;
t97 = qJD(1) * t77;
t64 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t97;
t65 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t97;
t78 = sin(qJ(1));
t83 = cos(qJ(1));
t68 = t78 * g(1) - t83 * g(2);
t39 = -0.2e1 * qJD(3) * t97 + (-t61 - t93) * qJ(3) - t68;
t84 = qJD(1) ^ 2;
t69 = -t83 * g(1) - t78 * g(2);
t94 = -t77 * g(3) + t82 * t69;
t43 = 0.2e1 * qJD(3) * qJD(2) + (-t77 * t82 * t84 + qJDD(2)) * qJ(3) + t94;
t76 = sin(qJ(4));
t81 = cos(qJ(4));
t28 = -t81 * t39 + t76 * t43;
t56 = t81 * qJD(2) - t76 * t97;
t41 = t56 * qJD(4) + t76 * qJDD(2) + t81 * t61;
t57 = t76 * qJD(2) + t81 * t97;
t44 = -t56 * mrSges(5,1) + t57 * mrSges(5,2);
t96 = t82 * qJD(1);
t70 = qJD(4) - t96;
t50 = -t70 * mrSges(5,2) + t56 * mrSges(5,3);
t55 = qJDD(4) - t62;
t75 = sin(qJ(5));
t80 = cos(qJ(5));
t48 = -t75 * t57 + t80 * t70;
t25 = t48 * qJD(5) + t80 * t41 + t75 * t55;
t49 = t80 * t57 + t75 * t70;
t53 = qJD(5) - t56;
t74 = sin(qJ(6));
t79 = cos(qJ(6));
t31 = -t74 * t49 + t79 * t53;
t40 = -t57 * qJD(4) + t81 * qJDD(2) - t76 * t61;
t38 = qJDD(5) - t40;
t19 = t31 * qJD(6) + t79 * t25 + t74 * t38;
t29 = t76 * t39 + t81 * t43;
t91 = -t82 * g(3) - t77 * t69;
t47 = qJDD(3) + (-t77 ^ 2 * t84 - qJD(2) ^ 2) * qJ(3) - t91;
t21 = t80 * t29 + t75 * t47;
t32 = t79 * t49 + t74 * t53;
t22 = -t31 * mrSges(7,1) + t32 * mrSges(7,2);
t24 = -t49 * qJD(5) - t75 * t41 + t80 * t55;
t23 = qJDD(6) - t24;
t46 = qJD(6) - t48;
t26 = -t46 * mrSges(7,2) + t31 * mrSges(7,3);
t16 = m(7) * (-t74 * t21 + t79 * t28) - t19 * mrSges(7,3) + t23 * mrSges(7,1) - t32 * t22 + t46 * t26;
t18 = -t32 * qJD(6) - t74 * t25 + t79 * t38;
t27 = t46 * mrSges(7,1) - t32 * mrSges(7,3);
t17 = m(7) * (t79 * t21 + t74 * t28) + t18 * mrSges(7,3) - t23 * mrSges(7,2) + t31 * t22 - t46 * t27;
t33 = -t53 * mrSges(6,2) + t48 * mrSges(6,3);
t34 = t53 * mrSges(6,1) - t49 * mrSges(6,3);
t85 = t24 * mrSges(6,1) - t25 * mrSges(6,2) - t79 * t16 - t74 * t17 + t48 * t33 - t49 * t34;
t11 = t55 * mrSges(5,1) - t41 * mrSges(5,3) - t57 * t44 + t70 * t50 + (-m(5) - m(6)) * t28 + t85;
t30 = -t48 * mrSges(6,1) + t49 * mrSges(6,2);
t13 = m(6) * t21 - t38 * mrSges(6,2) + t24 * mrSges(6,3) - t74 * t16 + t79 * t17 + t48 * t30 - t53 * t34;
t20 = t75 * t29 - t80 * t47;
t87 = t18 * mrSges(7,1) - t19 * mrSges(7,2) + t31 * t26 - t32 * t27;
t15 = t38 * mrSges(6,1) - t25 * mrSges(6,3) - t49 * t30 + t53 * t33 + (-m(6) - m(7)) * t20 + t87;
t51 = t70 * mrSges(5,1) - t57 * mrSges(5,3);
t9 = m(5) * t29 - t55 * mrSges(5,2) + t40 * mrSges(5,3) + t80 * t13 - t75 * t15 + t56 * t44 - t70 * t51;
t92 = m(4) * t39 - t62 * mrSges(4,1) + t81 * t11 + t76 * t9;
t67 = mrSges(4,2) * t96 + qJD(2) * mrSges(4,3);
t98 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t96 + t67;
t105 = -(-t98 * t82 + (t64 - t65) * t77) * qJD(1) + t62 * mrSges(3,1) - (mrSges(3,2) - mrSges(4,3)) * t61 - t92;
t100 = mrSges(3,3) + mrSges(4,2);
t60 = (-mrSges(3,1) * t82 + mrSges(3,2) * t77) * qJD(1);
t59 = (-mrSges(4,1) * t82 - mrSges(4,3) * t77) * qJD(1);
t89 = m(4) * t43 + qJDD(2) * mrSges(4,3) + qJD(2) * t65 - t76 * t11 + t59 * t96 + t81 * t9;
t4 = m(3) * t94 - qJDD(2) * mrSges(3,2) - qJD(2) * t64 + t100 * t62 + t60 * t96 + t89;
t88 = -m(5) * t47 + t40 * mrSges(5,1) - t41 * mrSges(5,2) - t75 * t13 - t80 * t15 + t56 * t50 - t57 * t51;
t86 = -m(4) * t47 + t88;
t6 = m(3) * t91 - t100 * t61 + (mrSges(3,1) + mrSges(4,1)) * qJDD(2) + t98 * qJD(2) + (-t59 - t60) * t97 + t86;
t102 = t77 * t4 + t82 * t6;
t2 = qJDD(1) * mrSges(2,1) - t84 * mrSges(2,2) + (m(2) + m(3)) * t68 + t105;
t1 = m(2) * t69 - t84 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t82 * t4 - t77 * t6;
t3 = [-m(1) * g(1) + t83 * t1 - t78 * t2, t1, t4, t62 * mrSges(4,2) + t89, t9, t13, t17; -m(1) * g(2) + t78 * t1 + t83 * t2, t2, t6, -t61 * mrSges(4,3) + (-t65 * t77 - t67 * t82) * qJD(1) + t92, t11, t15, t16; (-m(1) - m(2)) * g(3) + t102, -m(2) * g(3) + t102, -m(3) * t68 - t105, -qJDD(2) * mrSges(4,1) + t61 * mrSges(4,2) - qJD(2) * t67 + t59 * t97 - t86, -t88, m(6) * t28 - t85, m(7) * t20 - t87;];
f_new  = t3;
