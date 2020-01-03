% Calculate vector of cutting forces with Newton-Euler
% S5RRPPR7
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRPPR7_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR7_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR7_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR7_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR7_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR7_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR7_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR7_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR7_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:34:59
% EndTime: 2019-12-31 19:35:02
% DurationCPUTime: 1.03s
% Computational Cost: add. (8078->166), mult. (18595->209), div. (0->0), fcn. (11732->8), ass. (0->82)
t79 = cos(qJ(2));
t100 = qJD(1) * t79;
t76 = sin(qJ(2));
t101 = qJD(1) * t76;
t102 = cos(pkin(8));
t74 = sin(pkin(8));
t59 = -t102 * t100 + t74 * t101;
t60 = (t102 * t76 + t74 * t79) * qJD(1);
t37 = t59 * pkin(3) - t60 * qJ(4);
t119 = (2 * qJD(3)) + t37;
t98 = qJD(1) * qJD(2);
t66 = t76 * qJDD(1) + t79 * t98;
t67 = t79 * qJDD(1) - t76 * t98;
t69 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t101;
t70 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t100;
t82 = qJD(1) ^ 2;
t49 = t59 * mrSges(5,1) - qJD(2) * mrSges(5,3);
t103 = -qJD(2) * mrSges(4,2) - t59 * mrSges(4,3) - t49;
t107 = mrSges(4,1) - mrSges(5,2);
t43 = -t102 * t67 + t74 * t66;
t44 = t102 * t66 + t74 * t67;
t48 = qJD(2) * mrSges(4,1) - t60 * mrSges(4,3);
t68 = qJD(2) * pkin(2) - qJ(3) * t101;
t73 = t79 ^ 2;
t77 = sin(qJ(1));
t80 = cos(qJ(1));
t97 = t77 * g(1) - t80 * g(2);
t92 = -qJDD(1) * pkin(1) - t97;
t85 = -t67 * pkin(2) + qJDD(3) + t68 * t101 + (-qJ(3) * t73 - pkin(6)) * t82 + t92;
t81 = qJD(2) ^ 2;
t95 = -t80 * g(1) - t77 * g(2);
t63 = -t82 * pkin(1) + qJDD(1) * pkin(6) + t95;
t108 = t76 * t63;
t111 = pkin(2) * t82;
t27 = qJDD(2) * pkin(2) - t66 * qJ(3) - t108 + (qJ(3) * t98 + t76 * t111 - g(3)) * t79;
t96 = -t76 * g(3) + t79 * t63;
t28 = t67 * qJ(3) - qJD(2) * t68 - t73 * t111 + t96;
t93 = t102 * t27 - t74 * t28;
t17 = -qJDD(2) * pkin(3) - t81 * qJ(4) + t119 * t60 + qJDD(4) - t93;
t99 = qJD(2) * t59;
t12 = (t59 * t60 - qJDD(2)) * pkin(7) + (t44 + t99) * pkin(4) + t17;
t51 = t60 * pkin(4) - qJD(2) * pkin(7);
t58 = t59 ^ 2;
t114 = -2 * qJD(4);
t83 = (-t44 + t99) * qJ(4) + t85 + (qJD(2) * pkin(3) + t114) * t60;
t15 = -t58 * pkin(4) - t60 * t51 + (pkin(3) + pkin(7)) * t43 + t83;
t75 = sin(qJ(5));
t78 = cos(qJ(5));
t45 = -t75 * qJD(2) + t78 * t59;
t23 = t45 * qJD(5) + t78 * qJDD(2) + t75 * t43;
t46 = t78 * qJD(2) + t75 * t59;
t29 = -t45 * mrSges(6,1) + t46 * mrSges(6,2);
t57 = qJD(5) + t60;
t33 = -t57 * mrSges(6,2) + t45 * mrSges(6,3);
t42 = qJDD(5) + t44;
t10 = m(6) * (t78 * t12 - t75 * t15) - t23 * mrSges(6,3) + t42 * mrSges(6,1) - t46 * t29 + t57 * t33;
t22 = -t46 * qJD(5) - t75 * qJDD(2) + t78 * t43;
t34 = t57 * mrSges(6,1) - t46 * mrSges(6,3);
t11 = m(6) * (t75 * t12 + t78 * t15) + t22 * mrSges(6,3) - t42 * mrSges(6,2) + t45 * t29 - t57 * t34;
t50 = t60 * mrSges(5,1) + qJD(2) * mrSges(5,2);
t91 = t75 * t10 - t78 * t11 - m(5) * (t43 * pkin(3) + t83) + t60 * t50 + t44 * mrSges(5,3);
t84 = m(4) * t85 + t44 * mrSges(4,2) + t103 * t59 + t107 * t43 + t60 * t48 - t91;
t118 = (t76 * t69 - t79 * t70) * qJD(1) - t67 * mrSges(3,1) + t66 * mrSges(3,2) + m(3) * (-t82 * pkin(6) + t92) + t84;
t116 = -2 * qJD(3);
t65 = (-mrSges(3,1) * t79 + mrSges(3,2) * t76) * qJD(1);
t39 = -t59 * mrSges(5,2) - t60 * mrSges(5,3);
t104 = -t59 * mrSges(4,1) - t60 * mrSges(4,2) - t39;
t106 = -mrSges(4,3) - mrSges(5,1);
t90 = -m(5) * t17 - t78 * t10 - t75 * t11;
t7 = m(4) * t93 + (m(4) * t116 + t104) * t60 + t106 * t44 + t107 * qJDD(2) + t103 * qJD(2) + t90;
t105 = t102 * t28 + t74 * t27;
t54 = t59 * t116;
t89 = t81 * pkin(3) - qJDD(2) * qJ(4) - t105;
t88 = -t22 * mrSges(6,1) - t45 * t33 + m(6) * (-t43 * pkin(4) - t58 * pkin(7) - t59 * t37 + t54 + ((2 * qJD(4)) + t51) * qJD(2) - t89) + t46 * t34 + t23 * mrSges(6,2);
t86 = -m(5) * (qJD(2) * t114 + t119 * t59 + t89) + t88;
t8 = m(4) * (t54 + t105) + t104 * t59 + t106 * t43 + (-mrSges(4,2) + mrSges(5,3)) * qJDD(2) + (-t48 + t50) * qJD(2) + t86;
t4 = m(3) * (-t79 * g(3) - t108) - t66 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t65 * t101 + qJD(2) * t70 + t74 * t8 + t102 * t7;
t5 = m(3) * t96 - qJDD(2) * mrSges(3,2) + t67 * mrSges(3,3) - qJD(2) * t69 + t65 * t100 + t102 * t8 - t74 * t7;
t113 = t79 * t4 + t76 * t5;
t6 = m(2) * t97 + qJDD(1) * mrSges(2,1) - t82 * mrSges(2,2) - t118;
t1 = m(2) * t95 - t82 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t76 * t4 + t79 * t5;
t2 = [-m(1) * g(1) + t80 * t1 - t77 * t6, t1, t5, t8, -t43 * mrSges(5,2) - t59 * t49 - t91, t11; -m(1) * g(2) + t77 * t1 + t80 * t6, t6, t4, t7, t43 * mrSges(5,1) - qJDD(2) * mrSges(5,3) - qJD(2) * t50 + t59 * t39 - t86, t10; (-m(1) - m(2)) * g(3) + t113, -m(2) * g(3) + t113, t118, t84, t44 * mrSges(5,1) + qJDD(2) * mrSges(5,2) + qJD(2) * t49 + t60 * t39 - t90, t88;];
f_new = t2;
