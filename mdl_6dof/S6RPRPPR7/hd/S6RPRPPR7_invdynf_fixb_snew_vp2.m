% Calculate vector of cutting forces with Newton-Euler
% S6RPRPPR7
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4]';
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
% Datum: 2019-05-05 17:17
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRPPR7_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR7_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR7_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR7_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR7_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR7_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR7_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR7_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR7_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:15:07
% EndTime: 2019-05-05 17:15:10
% DurationCPUTime: 1.20s
% Computational Cost: add. (10934->180), mult. (23996->218), div. (0->0), fcn. (14726->8), ass. (0->89)
t86 = sin(qJ(1));
t89 = cos(qJ(1));
t106 = -t89 * g(1) - t86 * g(2);
t103 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t106;
t91 = qJD(1) ^ 2;
t88 = cos(qJ(3));
t115 = qJD(1) * t88;
t85 = sin(qJ(3));
t116 = qJD(1) * t85;
t128 = -pkin(1) - pkin(7);
t113 = qJD(1) * qJD(3);
t72 = -t85 * qJDD(1) - t88 * t113;
t108 = t85 * t113;
t73 = t88 * qJDD(1) - t108;
t74 = -(qJD(3) * mrSges(4,2)) - mrSges(4,3) * t116;
t76 = (qJD(3) * mrSges(4,1)) - mrSges(4,3) * t115;
t117 = sin(pkin(9));
t83 = cos(pkin(9));
t64 = (t117 * t88 + t83 * t85) * qJD(1);
t114 = qJD(3) * t64;
t110 = t86 * g(1) - t89 * g(2);
t99 = -t91 * qJ(2) + qJDD(2) - t110;
t55 = t128 * qJDD(1) + t99;
t118 = t85 * g(3) + t88 * t55;
t26 = (-t73 - t108) * qJ(4) + (-t85 * t88 * t91 + qJDD(3)) * pkin(3) + t118;
t109 = -t88 * g(3) + t85 * t55;
t75 = qJD(3) * pkin(3) - qJ(4) * t115;
t82 = t85 ^ 2;
t27 = -t82 * t91 * pkin(3) + t72 * qJ(4) - qJD(3) * t75 + t109;
t107 = -t117 * t27 + t83 * t26;
t65 = t83 * t115 - t117 * t116;
t36 = t64 * pkin(4) - t65 * qJ(5);
t133 = (2 * qJD(4)) + t36;
t90 = qJD(3) ^ 2;
t16 = -qJDD(3) * pkin(4) - t90 * qJ(5) + t133 * t65 + qJDD(5) - t107;
t43 = t117 * t72 + t83 * t73;
t11 = (t64 * t65 - qJDD(3)) * pkin(8) + (t43 + t114) * pkin(5) + t16;
t42 = t117 * t73 - t83 * t72;
t54 = t65 * pkin(5) - qJD(3) * pkin(8);
t63 = t64 ^ 2;
t130 = -2 * qJD(5);
t94 = -t72 * pkin(3) + qJDD(4) + t75 * t115 + (-qJ(4) * t82 + t128) * t91 + t103;
t92 = (-t43 + t114) * qJ(5) + t94 + (qJD(3) * pkin(4) + t130) * t65;
t14 = (pkin(4) + pkin(8)) * t42 - t63 * pkin(5) - t65 * t54 + t92;
t84 = sin(qJ(6));
t87 = cos(qJ(6));
t45 = t87 * qJD(3) + t84 * t64;
t21 = -t45 * qJD(6) - t84 * qJDD(3) + t87 * t42;
t44 = -t84 * qJD(3) + t87 * t64;
t30 = -t44 * mrSges(7,1) + t45 * mrSges(7,2);
t61 = qJD(6) + t65;
t33 = t61 * mrSges(7,1) - t45 * mrSges(7,3);
t41 = qJDD(6) + t43;
t10 = m(7) * (t84 * t11 + t87 * t14) + t21 * mrSges(7,3) - t41 * mrSges(7,2) + t44 * t30 - t61 * t33;
t53 = t65 * mrSges(6,1) + (qJD(3) * mrSges(6,2));
t22 = t44 * qJD(6) + t87 * qJDD(3) + t84 * t42;
t32 = -t61 * mrSges(7,2) + t44 * mrSges(7,3);
t9 = m(7) * (t87 * t11 - t84 * t14) - t22 * mrSges(7,3) + t41 * mrSges(7,1) - t45 * t30 + t61 * t32;
t104 = -t87 * t10 + t84 * t9 - m(6) * (t42 * pkin(4) + t92) + t65 * t53 + t43 * mrSges(6,3);
t52 = t64 * mrSges(6,1) - qJD(3) * mrSges(6,3);
t119 = -(qJD(3) * mrSges(5,2)) - t64 * mrSges(5,3) - t52;
t124 = mrSges(5,1) - mrSges(6,2);
t51 = (qJD(3) * mrSges(5,1)) - t65 * mrSges(5,3);
t95 = m(5) * t94 + t43 * mrSges(5,2) + t119 * t64 + t124 * t42 + t65 * t51 - t104;
t93 = m(4) * (t128 * t91 + t103) + t74 * t116 + t76 * t115 + t73 * mrSges(4,2) + t95 - t72 * mrSges(4,1);
t134 = m(3) * ((t91 * pkin(1)) - t103) - t93;
t132 = -2 * qJD(4);
t129 = -m(2) - m(3);
t125 = (mrSges(2,1) - mrSges(3,2));
t123 = -mrSges(2,2) + mrSges(3,3);
t122 = -mrSges(5,3) - mrSges(6,1);
t121 = t117 * t26 + t83 * t27;
t38 = -t64 * mrSges(6,2) - t65 * mrSges(6,3);
t120 = -t64 * mrSges(5,1) - t65 * mrSges(5,2) - t38;
t101 = -m(6) * t16 - t84 * t10 - t87 * t9;
t6 = m(5) * t107 + (m(5) * t132 + t120) * t65 + t122 * t43 + t124 * qJDD(3) + t119 * qJD(3) + t101;
t59 = t64 * t132;
t100 = t90 * pkin(4) - qJDD(3) * qJ(5) - t121;
t98 = -t21 * mrSges(7,1) - t44 * t32 + m(7) * (-t42 * pkin(5) - t63 * pkin(8) - t64 * t36 + t59 + ((2 * qJD(5)) + t54) * qJD(3) - t100) + t45 * t33 + t22 * mrSges(7,2);
t96 = -m(6) * ((qJD(3) * t130) + t133 * t64 + t100) + t98;
t7 = m(5) * (t59 + t121) + t120 * t64 + t122 * t42 + (-mrSges(5,2) + mrSges(6,3)) * qJDD(3) + (-t51 + t53) * qJD(3) + t96;
t71 = (mrSges(4,1) * t85 + mrSges(4,2) * t88) * qJD(1);
t3 = m(4) * t118 + qJDD(3) * mrSges(4,1) - t73 * mrSges(4,3) + qJD(3) * t74 - t71 * t115 + t117 * t7 + t83 * t6;
t4 = m(4) * t109 - qJDD(3) * mrSges(4,2) + t72 * mrSges(4,3) - qJD(3) * t76 - t71 * t116 - t117 * t6 + t83 * t7;
t111 = -t85 * t3 + t88 * t4;
t102 = -m(3) * (-qJDD(1) * pkin(1) + t99) - t88 * t3 - t85 * t4;
t5 = m(2) * t106 + t123 * qJDD(1) - (t125 * t91) - t134;
t1 = m(2) * t110 + t125 * qJDD(1) + t123 * t91 + t102;
t2 = [-m(1) * g(1) - t86 * t1 + t89 * t5, t5, -m(3) * g(3) + t111, t4, t7, -t42 * mrSges(6,2) - t64 * t52 - t104, t10; -m(1) * g(2) + t89 * t1 + t86 * t5, t1, -(t91 * mrSges(3,2)) - qJDD(1) * mrSges(3,3) + t134, t3, t6, t42 * mrSges(6,1) - qJDD(3) * mrSges(6,3) - qJD(3) * t53 + t64 * t38 - t96, t9; (-m(1) + t129) * g(3) + t111, t129 * g(3) + t111, qJDD(1) * mrSges(3,2) - t91 * mrSges(3,3) - t102, t93, t95, t43 * mrSges(6,1) + qJDD(3) * mrSges(6,2) + qJD(3) * t52 + t65 * t38 - t101, t98;];
f_new  = t2;
