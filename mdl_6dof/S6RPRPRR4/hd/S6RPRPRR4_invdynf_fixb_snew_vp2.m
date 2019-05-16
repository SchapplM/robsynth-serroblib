% Calculate vector of cutting forces with Newton-Euler
% S6RPRPRR4
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
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
% Datum: 2019-05-05 18:46
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRPRR4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR4_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR4_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR4_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR4_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR4_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR4_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR4_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 18:42:43
% EndTime: 2019-05-05 18:42:47
% DurationCPUTime: 1.64s
% Computational Cost: add. (18515->179), mult. (36744->223), div. (0->0), fcn. (21150->10), ass. (0->94)
t137 = -2 * qJD(4);
t93 = sin(qJ(1));
t97 = cos(qJ(1));
t116 = t93 * g(1) - t97 * g(2);
t62 = qJDD(1) * pkin(1) + t116;
t111 = -t97 * g(1) - t93 * g(2);
t99 = qJD(1) ^ 2;
t66 = -t99 * pkin(1) + t111;
t88 = sin(pkin(10));
t89 = cos(pkin(10));
t112 = t89 * t62 - t88 * t66;
t107 = -qJDD(1) * pkin(2) - t112;
t119 = qJD(1) * qJD(3);
t96 = cos(qJ(3));
t114 = t96 * t119;
t92 = sin(qJ(3));
t115 = t92 * t119;
t67 = t92 * qJDD(1) + t114;
t81 = t92 * qJD(1);
t102 = pkin(3) * t115 + t81 * t137 + (-t67 - t114) * qJ(4) + t107;
t120 = qJD(1) * t96;
t130 = t99 * pkin(7);
t68 = t96 * qJDD(1) - t115;
t72 = -mrSges(5,1) * t120 - qJD(3) * mrSges(5,3);
t132 = -pkin(3) - pkin(8);
t74 = pkin(4) * t81 - qJD(3) * pkin(8);
t86 = t96 ^ 2;
t19 = -t74 * t81 + (-pkin(4) * t86 - pkin(7)) * t99 + t132 * t68 + t102;
t123 = t88 * t62 + t89 * t66;
t39 = -t99 * pkin(2) + qJDD(1) * pkin(7) + t123;
t36 = t92 * t39;
t63 = (-pkin(3) * t96 - qJ(4) * t92) * qJD(1);
t98 = qJD(3) ^ 2;
t110 = -t98 * qJ(4) + t63 * t81 + qJDD(4) + t36;
t131 = pkin(8) * t99;
t87 = -g(3) + qJDD(2);
t27 = t67 * pkin(4) + t132 * qJDD(3) + (-pkin(4) * t119 - t131 * t92 - t87) * t96 + t110;
t91 = sin(qJ(5));
t95 = cos(qJ(5));
t113 = -t91 * t19 + t95 * t27;
t60 = -t91 * qJD(3) - t120 * t95;
t43 = t60 * qJD(5) + t95 * qJDD(3) - t91 * t68;
t59 = qJDD(5) + t67;
t61 = t95 * qJD(3) - t120 * t91;
t77 = t81 + qJD(5);
t14 = (t60 * t77 - t43) * pkin(9) + (t60 * t61 + t59) * pkin(5) + t113;
t125 = t95 * t19 + t91 * t27;
t42 = -t61 * qJD(5) - t91 * qJDD(3) - t95 * t68;
t49 = t77 * pkin(5) - t61 * pkin(9);
t58 = t60 ^ 2;
t15 = -t58 * pkin(5) + t42 * pkin(9) - t77 * t49 + t125;
t90 = sin(qJ(6));
t94 = cos(qJ(6));
t44 = t94 * t60 - t90 * t61;
t22 = t44 * qJD(6) + t90 * t42 + t94 * t43;
t45 = t90 * t60 + t94 * t61;
t33 = -t44 * mrSges(7,1) + t45 * mrSges(7,2);
t75 = qJD(6) + t77;
t34 = -t75 * mrSges(7,2) + t44 * mrSges(7,3);
t55 = qJDD(6) + t59;
t12 = m(7) * (t94 * t14 - t90 * t15) - t22 * mrSges(7,3) + t55 * mrSges(7,1) - t45 * t33 + t75 * t34;
t21 = -t45 * qJD(6) + t94 * t42 - t90 * t43;
t35 = t75 * mrSges(7,1) - t45 * mrSges(7,3);
t13 = m(7) * (t90 * t14 + t94 * t15) + t21 * mrSges(7,3) - t55 * mrSges(7,2) + t44 * t33 - t75 * t35;
t46 = -t60 * mrSges(6,1) + t61 * mrSges(6,2);
t47 = -t77 * mrSges(6,2) + t60 * mrSges(6,3);
t8 = m(6) * t113 + t59 * mrSges(6,1) - t43 * mrSges(6,3) + t94 * t12 + t90 * t13 - t61 * t46 + t77 * t47;
t48 = t77 * mrSges(6,1) - t61 * mrSges(6,3);
t9 = m(6) * t125 - t59 * mrSges(6,2) + t42 * mrSges(6,3) - t90 * t12 + t94 * t13 + t60 * t46 - t77 * t48;
t109 = t91 * t8 - t95 * t9 - m(5) * (-t68 * pkin(3) + t102 - t130) - t72 * t120 + t67 * mrSges(5,3);
t73 = mrSges(5,1) * t81 + qJD(3) * mrSges(5,2);
t121 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t81 - t73;
t127 = mrSges(4,1) - mrSges(5,2);
t71 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t120;
t136 = (t121 * t92 - t96 * t71) * qJD(1) - t127 * t68 + m(4) * (t107 - t130) + t67 * mrSges(4,2) - t109;
t124 = t96 * t39 + t92 * t87;
t135 = t98 * pkin(3) - qJDD(3) * qJ(4) + qJD(3) * t137 - t63 * t120 - t124;
t128 = t96 * t87;
t126 = mrSges(4,3) + mrSges(5,1);
t64 = (mrSges(5,2) * t96 - mrSges(5,3) * t92) * qJD(1);
t122 = t64 + (-mrSges(4,1) * t96 + mrSges(4,2) * t92) * qJD(1);
t103 = t68 * pkin(4) + qJD(3) * t74 - t131 * t86 - t135;
t105 = -t21 * mrSges(7,1) - t44 * t34 + m(7) * (-t42 * pkin(5) - t58 * pkin(9) + t61 * t49 + t103) + t22 * mrSges(7,2) + t45 * t35;
t101 = m(6) * t103 - t42 * mrSges(6,1) + t43 * mrSges(6,2) - t60 * t47 + t61 * t48 + t105;
t100 = -m(5) * t135 + t101;
t11 = t122 * t120 + t100 + t126 * t68 + (-mrSges(4,2) + mrSges(5,3)) * qJDD(3) - t121 * qJD(3) + m(4) * t124;
t106 = -m(5) * (-qJDD(3) * pkin(3) + t110 - t128) - t95 * t8 - t91 * t9;
t6 = m(4) * (-t36 + t128) - t126 * t67 + t127 * qJDD(3) + (t71 - t72) * qJD(3) - t122 * t81 + t106;
t117 = m(3) * t87 + t92 * t11 + t96 * t6;
t4 = m(3) * t112 + qJDD(1) * mrSges(3,1) - t99 * mrSges(3,2) - t136;
t3 = m(3) * t123 - t99 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t96 * t11 - t92 * t6;
t2 = m(2) * t111 - t99 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t89 * t3 - t88 * t4;
t1 = m(2) * t116 + qJDD(1) * mrSges(2,1) - t99 * mrSges(2,2) + t88 * t3 + t89 * t4;
t5 = [-m(1) * g(1) - t93 * t1 + t97 * t2, t2, t3, t11, t68 * mrSges(5,2) - t73 * t81 - t109, t9, t13; -m(1) * g(2) + t97 * t1 + t93 * t2, t1, t4, t6, -t68 * mrSges(5,1) - qJDD(3) * mrSges(5,3) - qJD(3) * t73 - t120 * t64 - t100, t8, t12; (-m(1) - m(2)) * g(3) + t117, -m(2) * g(3) + t117, t117, t136, t67 * mrSges(5,1) + qJDD(3) * mrSges(5,2) + qJD(3) * t72 + t64 * t81 - t106, t101, t105;];
f_new  = t5;
