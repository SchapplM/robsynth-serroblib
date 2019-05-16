% Calculate vector of cutting forces with Newton-Euler
% S6RPRPRR12
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6]';
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
% Datum: 2019-05-05 20:47
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRPRR12_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR12_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR12_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR12_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR12_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRR12_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR12_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR12_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR12_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 20:44:42
% EndTime: 2019-05-05 20:44:46
% DurationCPUTime: 1.26s
% Computational Cost: add. (12586->180), mult. (25137->218), div. (0->0), fcn. (13950->8), ass. (0->92)
t86 = sin(qJ(1));
t90 = cos(qJ(1));
t109 = -t90 * g(1) - t86 * g(2);
t105 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t109;
t92 = qJD(1) ^ 2;
t131 = (-pkin(1) - pkin(7));
t116 = t131 * t92;
t85 = sin(qJ(3));
t120 = qJD(1) * t85;
t124 = mrSges(4,1) - mrSges(5,2);
t119 = qJD(1) * qJD(3);
t89 = cos(qJ(3));
t112 = t89 * t119;
t62 = t85 * qJDD(1) + t112;
t75 = t85 * t119;
t63 = t89 * qJDD(1) - t75;
t65 = -(qJD(3) * mrSges(4,2)) - mrSges(4,3) * t120;
t76 = t89 * qJD(1);
t66 = (qJD(3) * mrSges(4,1)) - mrSges(4,3) * t76;
t130 = pkin(3) + pkin(8);
t69 = pkin(4) * t76 - (qJD(3) * pkin(8));
t82 = t85 ^ 2;
t135 = -2 * qJD(4);
t98 = pkin(3) * t112 + t76 * t135 + t105 + (-t63 + t75) * qJ(4);
t19 = -t69 * t76 + t130 * t62 + (-pkin(4) * t82 + t131) * t92 + t98;
t114 = t86 * g(1) - t90 * g(2);
t102 = -t92 * qJ(2) + qJDD(2) - t114;
t46 = t131 * qJDD(1) + t102;
t126 = t89 * t46;
t59 = (pkin(3) * t85 - qJ(4) * t89) * qJD(1);
t91 = qJD(3) ^ 2;
t101 = -t91 * qJ(4) + t59 * t76 + qJDD(4) - t126;
t128 = pkin(8) * t92;
t24 = t63 * pkin(4) - t130 * qJDD(3) + (pkin(4) * t119 + t89 * t128 - g(3)) * t85 + t101;
t84 = sin(qJ(5));
t88 = cos(qJ(5));
t111 = -t84 * t19 + t88 * t24;
t57 = -t84 * qJD(3) + t88 * t120;
t36 = t57 * qJD(5) + t88 * qJDD(3) + t84 * t62;
t56 = qJDD(5) + t63;
t58 = t88 * qJD(3) + t84 * t120;
t73 = t76 + qJD(5);
t11 = (t57 * t73 - t36) * pkin(9) + (t57 * t58 + t56) * pkin(5) + t111;
t121 = t88 * t19 + t84 * t24;
t35 = -t58 * qJD(5) - t84 * qJDD(3) + t88 * t62;
t44 = t73 * pkin(5) - t58 * pkin(9);
t55 = t57 ^ 2;
t12 = -t55 * pkin(5) + t35 * pkin(9) - t73 * t44 + t121;
t83 = sin(qJ(6));
t87 = cos(qJ(6));
t38 = t83 * t57 + t87 * t58;
t16 = -t38 * qJD(6) + t87 * t35 - t83 * t36;
t37 = t87 * t57 - t83 * t58;
t27 = -t37 * mrSges(7,1) + t38 * mrSges(7,2);
t70 = qJD(6) + t73;
t32 = t70 * mrSges(7,1) - t38 * mrSges(7,3);
t50 = qJDD(6) + t56;
t10 = m(7) * (t83 * t11 + t87 * t12) + t16 * mrSges(7,3) - t50 * mrSges(7,2) + t37 * t27 - t70 * t32;
t39 = -t57 * mrSges(6,1) + t58 * mrSges(6,2);
t40 = -t73 * mrSges(6,2) + t57 * mrSges(6,3);
t17 = t37 * qJD(6) + t83 * t35 + t87 * t36;
t31 = -t70 * mrSges(7,2) + t37 * mrSges(7,3);
t9 = m(7) * (t87 * t11 - t83 * t12) - t17 * mrSges(7,3) + t50 * mrSges(7,1) - t38 * t27 + t70 * t31;
t5 = m(6) * t111 + t56 * mrSges(6,1) - t36 * mrSges(6,3) + t83 * t10 - t58 * t39 + t73 * t40 + t87 * t9;
t41 = t73 * mrSges(6,1) - t58 * mrSges(6,3);
t6 = m(6) * t121 - t56 * mrSges(6,2) + t35 * mrSges(6,3) + t87 * t10 + t57 * t39 - t73 * t41 - t83 * t9;
t67 = mrSges(5,1) * t120 - (qJD(3) * mrSges(5,3));
t68 = mrSges(5,1) * t76 + (qJD(3) * mrSges(5,2));
t97 = -(t85 * t67 + t89 * t68) * qJD(1) - t84 * t5 + t88 * t6 + m(5) * (t62 * pkin(3) + t116 + t98) - t63 * mrSges(5,3);
t93 = t124 * t62 + m(4) * (t116 + t105) + t65 * t120 + t66 * t76 + t63 * mrSges(4,2) + t97;
t136 = m(3) * ((t92 * pkin(1)) - t105) - t93;
t113 = -t89 * g(3) + t85 * t46;
t134 = t91 * pkin(3) - qJDD(3) * qJ(4) + (qJD(3) * t135) + t59 * t120 - t113;
t132 = -m(2) - m(3);
t127 = t85 * g(3);
t125 = (mrSges(2,1) - mrSges(3,2));
t123 = -mrSges(2,2) + mrSges(3,3);
t122 = -mrSges(4,3) - mrSges(5,1);
t103 = -m(5) * (-qJDD(3) * pkin(3) + t101 - t127) - t88 * t5 - t84 * t6;
t60 = (-mrSges(5,2) * t85 - mrSges(5,3) * t89) * qJD(1);
t110 = qJD(1) * (-t60 - (mrSges(4,1) * t85 + mrSges(4,2) * t89) * qJD(1));
t3 = m(4) * (t126 + t127) + t122 * t63 + t124 * qJDD(3) + (t65 - t67) * qJD(3) + t89 * t110 + t103;
t95 = -t62 * pkin(4) + qJD(3) * t69 - t82 * t128 - t134;
t100 = -t16 * mrSges(7,1) - t37 * t31 + m(7) * (-t35 * pkin(5) - t55 * pkin(9) + t58 * t44 + t95) + t17 * mrSges(7,2) + t38 * t32;
t96 = m(6) * t95 - t35 * mrSges(6,1) + t36 * mrSges(6,2) - t57 * t40 + t58 * t41 + t100;
t94 = -m(5) * t134 + t96;
t8 = t94 + t122 * t62 + (-mrSges(4,2) + mrSges(5,3)) * qJDD(3) + (-t66 + t68) * qJD(3) + m(4) * t113 + t85 * t110;
t115 = -t85 * t3 + t89 * t8;
t104 = -m(3) * (-qJDD(1) * pkin(1) + t102) - t89 * t3 - t85 * t8;
t2 = m(2) * t109 + t123 * qJDD(1) - (t125 * t92) - t136;
t1 = m(2) * t114 + t125 * qJDD(1) + t123 * t92 + t104;
t4 = [-m(1) * g(1) - t86 * t1 + t90 * t2, t2, -m(3) * g(3) + t115, t8, -t62 * mrSges(5,2) + t97, t6, t10; -m(1) * g(2) + t90 * t1 + t86 * t2, t1, -(t92 * mrSges(3,2)) - qJDD(1) * mrSges(3,3) + t136, t3, t62 * mrSges(5,1) - qJDD(3) * mrSges(5,3) - qJD(3) * t68 + t60 * t120 - t94, t5, t9; (-m(1) + t132) * g(3) + t115, t132 * g(3) + t115, qJDD(1) * mrSges(3,2) - t92 * mrSges(3,3) - t104, t93, t63 * mrSges(5,1) + qJDD(3) * mrSges(5,2) + qJD(3) * t67 + t60 * t76 - t103, t96, t100;];
f_new  = t4;
