% Calculate vector of cutting forces with Newton-Euler
% S6RRRPRP9
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
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
% Datum: 2019-05-07 08:48
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRPRP9_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP9_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP9_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP9_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP9_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP9_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP9_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP9_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP9_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 08:38:53
% EndTime: 2019-05-07 08:38:58
% DurationCPUTime: 1.61s
% Computational Cost: add. (18640->198), mult. (36403->239), div. (0->0), fcn. (23588->8), ass. (0->94)
t106 = sin(qJ(3));
t107 = sin(qJ(2));
t109 = cos(qJ(2));
t146 = cos(qJ(3));
t105 = sin(qJ(5));
t145 = cos(qJ(5));
t134 = qJD(1) * t107;
t85 = -t146 * qJD(2) + t106 * t134;
t133 = qJD(1) * t109;
t97 = qJD(3) - t133;
t144 = t85 * t97;
t132 = qJD(1) * qJD(2);
t128 = t109 * t132;
t129 = t107 * t132;
t112 = qJD(1) ^ 2;
t108 = sin(qJ(1));
t110 = cos(qJ(1));
t127 = t108 * g(1) - t110 * g(2);
t78 = -qJDD(1) * pkin(1) - t112 * pkin(7) - t127;
t89 = qJDD(1) * t107 + t128;
t90 = t109 * qJDD(1) - t129;
t40 = (-t89 - t128) * pkin(8) + (-t90 + t129) * pkin(2) + t78;
t111 = qJD(2) ^ 2;
t124 = -g(1) * t110 - g(2) * t108;
t79 = -pkin(1) * t112 + qJDD(1) * pkin(7) + t124;
t130 = -g(3) * t107 + t109 * t79;
t88 = (-pkin(2) * t109 - pkin(8) * t107) * qJD(1);
t45 = -pkin(2) * t111 + qJDD(2) * pkin(8) + t88 * t133 + t130;
t125 = -t106 * t45 + t146 * t40;
t86 = t106 * qJD(2) + t146 * t134;
t63 = pkin(3) * t85 - qJ(4) * t86;
t84 = qJDD(3) - t90;
t96 = t97 ^ 2;
t25 = -t84 * pkin(3) - t96 * qJ(4) + t86 * t63 + qJDD(4) - t125;
t59 = -t85 * qJD(3) + t106 * qJDD(2) + t146 * t89;
t17 = (-t59 - t144) * pkin(9) + (t85 * t86 - t84) * pkin(4) + t25;
t138 = t106 * t40 + t146 * t45;
t150 = 2 * qJD(4);
t120 = -pkin(3) * t96 + t84 * qJ(4) + t97 * t150 - t85 * t63 + t138;
t58 = qJD(3) * t86 - t146 * qJDD(2) + t106 * t89;
t70 = -pkin(4) * t97 - pkin(9) * t86;
t83 = t85 ^ 2;
t20 = -pkin(4) * t83 + t58 * pkin(9) + t70 * t97 + t120;
t140 = t105 * t17 + t145 * t20;
t60 = t105 * t86 - t145 * t85;
t61 = t105 * t85 + t145 * t86;
t36 = pkin(5) * t60 - qJ(6) * t61;
t95 = qJD(5) - t97;
t49 = -mrSges(7,1) * t95 + t61 * mrSges(7,2);
t82 = qJDD(5) - t84;
t94 = t95 ^ 2;
t131 = m(7) * (-pkin(5) * t94 + qJ(6) * t82 + 0.2e1 * qJD(6) * t95 - t60 * t36 + t140) + t95 * t49 + t82 * mrSges(7,3);
t37 = mrSges(7,1) * t60 - mrSges(7,3) * t61;
t139 = -mrSges(6,1) * t60 - mrSges(6,2) * t61 - t37;
t141 = -mrSges(6,3) - mrSges(7,2);
t30 = t61 * qJD(5) + t105 * t59 - t145 * t58;
t48 = mrSges(6,1) * t95 - t61 * mrSges(6,3);
t10 = m(6) * t140 - t82 * mrSges(6,2) + t139 * t60 + t141 * t30 - t95 * t48 + t131;
t121 = -t105 * t20 + t145 * t17;
t148 = m(7) * (-t82 * pkin(5) - t94 * qJ(6) + t61 * t36 + qJDD(6) - t121);
t31 = -t60 * qJD(5) + t105 * t58 + t145 * t59;
t46 = -t60 * mrSges(7,2) + mrSges(7,3) * t95;
t47 = -mrSges(6,2) * t95 - t60 * mrSges(6,3);
t11 = m(6) * t121 - t148 + (t47 + t46) * t95 + (mrSges(6,1) + mrSges(7,1)) * t82 + t139 * t61 + t141 * t31;
t68 = -mrSges(5,1) * t97 + mrSges(5,2) * t86;
t122 = m(5) * t120 + t84 * mrSges(5,3) + t145 * t10 - t105 * t11 + t97 * t68;
t64 = mrSges(5,1) * t85 - mrSges(5,3) * t86;
t137 = -mrSges(4,1) * t85 - mrSges(4,2) * t86 - t64;
t142 = -mrSges(4,3) - mrSges(5,2);
t67 = mrSges(4,1) * t97 - mrSges(4,3) * t86;
t5 = m(4) * t138 - t84 * mrSges(4,2) + t137 * t85 + t142 * t58 - t97 * t67 + t122;
t118 = m(5) * t25 + t105 * t10 + t145 * t11;
t66 = -mrSges(4,2) * t97 - mrSges(4,3) * t85;
t69 = -mrSges(5,2) * t85 + mrSges(5,3) * t97;
t6 = m(4) * t125 + (t66 + t69) * t97 + t137 * t86 + (mrSges(4,1) + mrSges(5,1)) * t84 + t142 * t59 - t118;
t91 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t134;
t92 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t133;
t152 = m(3) * t78 - t90 * mrSges(3,1) + t89 * mrSges(3,2) + (t107 * t91 - t109 * t92) * qJD(1) + t106 * t5 + t146 * t6;
t135 = -t109 * g(3) - t107 * t79;
t44 = -qJDD(2) * pkin(2) - t111 * pkin(8) + t88 * t134 - t135;
t119 = t58 * pkin(3) + t44 + (t144 - t59) * qJ(4);
t147 = pkin(3) * t97;
t114 = -t58 * pkin(4) - t83 * pkin(9) - t119 + (-t147 + t150 + t70) * t86;
t126 = m(7) * (t114 + (t60 * t95 - t31) * qJ(6) + (t61 * t95 + t30) * pkin(5) - 0.2e1 * qJD(6) * t61) - t31 * mrSges(7,3) + t30 * mrSges(7,1) - t61 * t49 + t60 * t46;
t117 = m(6) * t114 + t30 * mrSges(6,1) + t31 * mrSges(6,2) + t60 * t47 + t61 * t48 + t126;
t116 = m(5) * ((-(2 * qJD(4)) + t147) * t86 + t119) + t58 * mrSges(5,1) + t85 * t69 - t117;
t151 = m(4) * t44 + t58 * mrSges(4,1) + (t67 - t68) * t86 + (mrSges(4,2) - mrSges(5,3)) * t59 + t85 * t66 + t116;
t87 = (-mrSges(3,1) * t109 + mrSges(3,2) * t107) * qJD(1);
t4 = m(3) * t130 - qJDD(2) * mrSges(3,2) + t90 * mrSges(3,3) - qJD(2) * t91 - t106 * t6 + t87 * t133 + t146 * t5;
t8 = m(3) * t135 + qJDD(2) * mrSges(3,1) - t89 * mrSges(3,3) + qJD(2) * t92 - t87 * t134 - t151;
t149 = t107 * t4 + t109 * t8;
t2 = m(2) * t127 + qJDD(1) * mrSges(2,1) - t112 * mrSges(2,2) - t152;
t1 = m(2) * t124 - t112 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t107 * t8 + t109 * t4;
t3 = [-m(1) * g(1) + t1 * t110 - t108 * t2, t1, t4, t5, -t58 * mrSges(5,2) - t85 * t64 + t122, t10, -t30 * mrSges(7,2) - t60 * t37 + t131; -m(1) * g(2) + t1 * t108 + t110 * t2, t2, t8, t6, -t59 * mrSges(5,3) - t86 * t68 + t116, t11, t126; (-m(1) - m(2)) * g(3) + t149, -m(2) * g(3) + t149, t152, t151, -t84 * mrSges(5,1) + t59 * mrSges(5,2) + t86 * t64 - t97 * t69 + t118, t117, -t82 * mrSges(7,1) + t31 * mrSges(7,2) + t61 * t37 - t95 * t46 + t148;];
f_new  = t3;
