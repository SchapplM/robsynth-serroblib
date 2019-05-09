% Calculate vector of cutting forces with Newton-Euler
% S6PPRRRR3
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
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
% Datum: 2019-05-04 21:20
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PPRRRR3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR3_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRRR3_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRR3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_invdynf_fixb_snew_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRR3_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRR3_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRRR3_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 21:12:27
% EndTime: 2019-05-04 21:12:34
% DurationCPUTime: 6.84s
% Computational Cost: add. (143270->152), mult. (272623->221), div. (0->0), fcn. (231037->18), ass. (0->99)
t78 = sin(pkin(13));
t83 = cos(pkin(13));
t71 = -t83 * g(1) - t78 * g(2);
t77 = sin(pkin(14));
t82 = cos(pkin(14));
t70 = t78 * g(1) - t83 * g(2);
t76 = -g(3) + qJDD(1);
t81 = sin(pkin(6));
t86 = cos(pkin(6));
t99 = t70 * t86 + t76 * t81;
t48 = -t77 * t71 + t99 * t82;
t58 = -t81 * t70 + t86 * t76 + qJDD(2);
t80 = sin(pkin(7));
t85 = cos(pkin(7));
t123 = t48 * t85 + t58 * t80;
t49 = t82 * t71 + t99 * t77;
t90 = sin(qJ(3));
t94 = cos(qJ(3));
t104 = t123 * t94 - t90 * t49;
t79 = sin(pkin(8));
t119 = pkin(10) * t79;
t95 = qJD(3) ^ 2;
t30 = qJDD(3) * pkin(3) + t95 * t119 + t104;
t37 = -t80 * t48 + t85 * t58;
t84 = cos(pkin(8));
t122 = t30 * t84 + t37 * t79;
t111 = qJD(3) * t79;
t93 = cos(qJ(4));
t106 = t93 * t111;
t108 = t123 * t90 + t94 * t49;
t31 = -t95 * pkin(3) + qJDD(3) * t119 + t108;
t89 = sin(qJ(4));
t109 = t122 * t89 + t93 * t31;
t65 = (-pkin(4) * t93 - pkin(11) * t89) * t111;
t74 = t84 * qJD(3) + qJD(4);
t72 = t74 ^ 2;
t73 = t84 * qJDD(3) + qJDD(4);
t24 = -t72 * pkin(4) + t73 * pkin(11) + t65 * t106 + t109;
t36 = t84 * t37;
t110 = qJD(3) * qJD(4);
t66 = (qJDD(3) * t89 + t93 * t110) * t79;
t67 = (-qJDD(3) * t93 + t89 * t110) * t79;
t26 = t67 * pkin(4) - t66 * pkin(11) + t36 + (-t30 + (pkin(4) * t89 - pkin(11) * t93) * t74 * qJD(3)) * t79;
t88 = sin(qJ(5));
t92 = cos(qJ(5));
t112 = t92 * t24 + t88 * t26;
t107 = t89 * t111;
t59 = -t88 * t107 + t92 * t74;
t60 = t92 * t107 + t88 * t74;
t44 = -t59 * pkin(5) - t60 * pkin(12);
t61 = qJDD(5) + t67;
t69 = qJD(5) - t106;
t68 = t69 ^ 2;
t20 = -t68 * pkin(5) + t61 * pkin(12) + t59 * t44 + t112;
t120 = t122 * t93 - t89 * t31;
t23 = -t73 * pkin(4) - t72 * pkin(11) + t65 * t107 - t120;
t41 = -t60 * qJD(5) - t88 * t66 + t92 * t73;
t42 = t59 * qJD(5) + t92 * t66 + t88 * t73;
t21 = (-t59 * t69 - t42) * pkin(12) + (t60 * t69 - t41) * pkin(5) + t23;
t87 = sin(qJ(6));
t91 = cos(qJ(6));
t50 = -t87 * t60 + t91 * t69;
t33 = t50 * qJD(6) + t91 * t42 + t87 * t61;
t51 = t91 * t60 + t87 * t69;
t34 = -t50 * mrSges(7,1) + t51 * mrSges(7,2);
t57 = qJD(6) - t59;
t38 = -t57 * mrSges(7,2) + t50 * mrSges(7,3);
t40 = qJDD(6) - t41;
t17 = m(7) * (-t87 * t20 + t91 * t21) - t33 * mrSges(7,3) + t40 * mrSges(7,1) - t51 * t34 + t57 * t38;
t32 = -t51 * qJD(6) - t87 * t42 + t91 * t61;
t39 = t57 * mrSges(7,1) - t51 * mrSges(7,3);
t18 = m(7) * (t91 * t20 + t87 * t21) + t32 * mrSges(7,3) - t40 * mrSges(7,2) + t50 * t34 - t57 * t39;
t43 = -t59 * mrSges(6,1) + t60 * mrSges(6,2);
t53 = t69 * mrSges(6,1) - t60 * mrSges(6,3);
t15 = m(6) * t112 - t61 * mrSges(6,2) + t41 * mrSges(6,3) - t87 * t17 + t91 * t18 + t59 * t43 - t69 * t53;
t101 = -t88 * t24 + t92 * t26;
t52 = -t69 * mrSges(6,2) + t59 * mrSges(6,3);
t97 = m(7) * (-t61 * pkin(5) - t68 * pkin(12) + t60 * t44 - t101) - t32 * mrSges(7,1) + t33 * mrSges(7,2) - t50 * t38 + t51 * t39;
t16 = m(6) * t101 + t61 * mrSges(6,1) - t42 * mrSges(6,3) - t60 * t43 + t69 * t52 - t97;
t62 = t74 * mrSges(5,1) - mrSges(5,3) * t107;
t64 = (-mrSges(5,1) * t93 + mrSges(5,2) * t89) * t111;
t12 = m(5) * t109 - t73 * mrSges(5,2) - t67 * mrSges(5,3) + t64 * t106 + t92 * t15 - t88 * t16 - t74 * t62;
t63 = -t74 * mrSges(5,2) + mrSges(5,3) * t106;
t96 = m(6) * t23 - t41 * mrSges(6,1) + t42 * mrSges(6,2) + t91 * t17 + t87 * t18 - t59 * t52 + t60 * t53;
t14 = m(5) * t120 + t73 * mrSges(5,1) - t66 * mrSges(5,3) - t64 * t107 + t74 * t63 - t96;
t102 = t89 * t12 + t93 * t14;
t13 = m(5) * (-t79 * t30 + t36) + t66 * mrSges(5,2) + t67 * mrSges(5,1) + t88 * t15 + t92 * t16 + (t62 * t89 - t63 * t93) * t111;
t10 = m(4) * t37 + t102 * t79 + t84 * t13;
t11 = m(4) * t108 - t95 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t93 * t12 - t89 * t14;
t9 = m(4) * t104 + qJDD(3) * mrSges(4,1) - t95 * mrSges(4,2) + t102 * t84 - t79 * t13;
t103 = t11 * t90 + t9 * t94;
t4 = m(3) * t48 - t80 * t10 + t103 * t85;
t8 = m(3) * t49 + t94 * t11 - t90 * t9;
t121 = t4 * t82 + t77 * t8;
t6 = m(3) * t58 + t85 * t10 + t103 * t80;
t105 = m(2) * t76 + t121 * t81 + t86 * t6;
t2 = m(2) * t71 - t77 * t4 + t82 * t8;
t1 = m(2) * t70 + t121 * t86 - t81 * t6;
t3 = [-m(1) * g(1) - t78 * t1 + t83 * t2, t2, t8, t11, t12, t15, t18; -m(1) * g(2) + t83 * t1 + t78 * t2, t1, t4, t9, t14, t16, t17; -m(1) * g(3) + t105, t105, t6, t10, t13, t96, t97;];
f_new  = t3;
