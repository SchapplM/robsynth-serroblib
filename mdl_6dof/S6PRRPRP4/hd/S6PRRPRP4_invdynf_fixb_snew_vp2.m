% Calculate vector of cutting forces with Newton-Euler
% S6PRRPRP4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
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
% Datum: 2019-05-05 04:05
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRRPRP4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP4_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRP4_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP4_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP4_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRP4_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 04:00:36
% EndTime: 2019-05-05 04:00:40
% DurationCPUTime: 1.24s
% Computational Cost: add. (11784->175), mult. (23093->214), div. (0->0), fcn. (13756->10), ass. (0->89)
t144 = -2 * qJD(4);
t94 = cos(qJ(3));
t119 = qJD(2) * t94;
t87 = sin(pkin(6));
t92 = sin(qJ(2));
t131 = t87 * t92;
t86 = sin(pkin(10));
t88 = cos(pkin(10));
t69 = g(1) * t86 - g(2) * t88;
t89 = cos(pkin(6));
t132 = t69 * t89;
t70 = -g(1) * t88 - g(2) * t86;
t85 = -g(3) + qJDD(1);
t95 = cos(qJ(2));
t113 = t85 * t131 + t92 * t132 + t95 * t70;
t97 = qJD(2) ^ 2;
t32 = -pkin(2) * t97 + qJDD(2) * pkin(8) + t113;
t50 = -t69 * t87 + t85 * t89;
t91 = sin(qJ(3));
t124 = t94 * t32 + t91 * t50;
t63 = (-pkin(3) * t94 - qJ(4) * t91) * qJD(2);
t96 = qJD(3) ^ 2;
t141 = pkin(3) * t96 - qJDD(3) * qJ(4) + (qJD(3) * t144) - t63 * t119 - t124;
t135 = pkin(9) * t97;
t118 = qJD(2) * qJD(3);
t110 = t91 * t118;
t67 = qJDD(2) * t94 - t110;
t120 = qJD(2) * t91;
t76 = pkin(4) * t120 - qJD(3) * pkin(9);
t84 = t94 ^ 2;
t100 = pkin(4) * t67 + qJD(3) * t76 - t84 * t135 - t141;
t90 = sin(qJ(5));
t93 = cos(qJ(5));
t62 = qJD(3) * t93 - t90 * t119;
t39 = -qJD(5) * t62 - qJDD(3) * t90 - t67 * t93;
t61 = -qJD(3) * t90 - t93 * t119;
t40 = qJD(5) * t61 + qJDD(3) * t93 - t67 * t90;
t78 = qJD(5) + t120;
t47 = pkin(5) * t78 - qJ(6) * t62;
t48 = mrSges(7,1) * t78 - mrSges(7,3) * t62;
t57 = t61 ^ 2;
t114 = m(7) * (-pkin(5) * t39 - qJ(6) * t57 + t47 * t62 + qJDD(6) + t100) + t40 * mrSges(7,2) + t62 * t48;
t45 = -mrSges(7,2) * t78 + mrSges(7,3) * t61;
t46 = -mrSges(6,2) * t78 + mrSges(6,3) * t61;
t49 = mrSges(6,1) * t78 - mrSges(6,3) * t62;
t99 = m(6) * t100 + t40 * mrSges(6,2) + (-t46 - t45) * t61 - (mrSges(6,1) + mrSges(7,1)) * t39 + t62 * t49 + t114;
t143 = -m(5) * t141 + t99;
t140 = (t85 * t87 + t132) * t95 - t92 * t70;
t102 = -qJDD(2) * pkin(2) - t140;
t29 = t91 * t32;
t106 = -t96 * qJ(4) + t63 * t120 + qJDD(4) + t29;
t137 = -pkin(3) - pkin(9);
t111 = t94 * t118;
t66 = qJDD(2) * t91 + t111;
t22 = t66 * pkin(4) + t137 * qJDD(3) + (-pkin(4) * t118 - t91 * t135 - t50) * t94 + t106;
t98 = pkin(3) * t110 + t120 * t144 + (-t66 - t111) * qJ(4) + t102;
t24 = -t76 * t120 + (-pkin(4) * t84 - pkin(8)) * t97 + t137 * t67 + t98;
t109 = t93 * t22 - t90 * t24;
t58 = qJDD(5) + t66;
t116 = m(7) * (-0.2e1 * qJD(6) * t62 + (t61 * t78 - t40) * qJ(6) + (t61 * t62 + t58) * pkin(5) + t109) + t78 * t45 + t58 * mrSges(7,1);
t42 = -mrSges(7,1) * t61 + mrSges(7,2) * t62;
t43 = -mrSges(6,1) * t61 + mrSges(6,2) * t62;
t11 = m(6) * t109 + t58 * mrSges(6,1) + t78 * t46 + (-t43 - t42) * t62 + (-mrSges(6,3) - mrSges(7,3)) * t40 + t116;
t125 = t90 * t22 + t93 * t24;
t115 = m(7) * (-pkin(5) * t57 + qJ(6) * t39 + 0.2e1 * qJD(6) * t61 - t47 * t78 + t125) + t61 * t42 + t39 * mrSges(7,3);
t13 = m(6) * t125 + t39 * mrSges(6,3) + t61 * t43 + (-t49 - t48) * t78 + (-mrSges(6,2) - mrSges(7,2)) * t58 + t115;
t133 = t97 * pkin(8);
t73 = -mrSges(5,1) * t119 - (qJD(3) * mrSges(5,3));
t104 = t90 * t11 - t93 * t13 - m(5) * (-t67 * pkin(3) - t133 + t98) - t73 * t119 + t66 * mrSges(5,3);
t74 = mrSges(5,1) * t120 + (qJD(3) * mrSges(5,2));
t121 = (qJD(3) * mrSges(4,1)) - mrSges(4,3) * t120 - t74;
t128 = mrSges(4,1) - mrSges(5,2);
t72 = -(qJD(3) * mrSges(4,2)) + mrSges(4,3) * t119;
t142 = (t121 * t91 - t94 * t72) * qJD(2) - t128 * t67 + m(4) * (t102 - t133) + t66 * mrSges(4,2) - t104;
t8 = m(3) * t140 + qJDD(2) * mrSges(3,1) - t97 * mrSges(3,2) - t142;
t134 = t8 * t95;
t130 = t94 * t50;
t126 = mrSges(4,3) + mrSges(5,1);
t64 = (mrSges(5,2) * t94 - mrSges(5,3) * t91) * qJD(2);
t122 = t64 + (-mrSges(4,1) * t94 + mrSges(4,2) * t91) * qJD(2);
t10 = m(4) * t124 + t126 * t67 + (-mrSges(4,2) + mrSges(5,3)) * qJDD(3) - t121 * qJD(3) + t122 * t119 + t143;
t103 = -m(5) * (-qJDD(3) * pkin(3) + t106 - t130) - t93 * t11 - t90 * t13;
t9 = m(4) * (-t29 + t130) - t126 * t66 + t128 * qJDD(3) + (t72 - t73) * qJD(3) - t122 * t120 + t103;
t4 = m(3) * t113 - t97 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t94 * t10 - t91 * t9;
t6 = m(3) * t50 + t10 * t91 + t9 * t94;
t112 = m(2) * t85 + t4 * t131 + t87 * t134 + t89 * t6;
t2 = m(2) * t70 + t4 * t95 - t8 * t92;
t1 = m(2) * t69 - t87 * t6 + (t4 * t92 + t134) * t89;
t3 = [-m(1) * g(1) - t1 * t86 + t2 * t88, t2, t4, t10, t67 * mrSges(5,2) - t74 * t120 - t104, t13, -t58 * mrSges(7,2) - t78 * t48 + t115; -m(1) * g(2) + t1 * t88 + t2 * t86, t1, t8, t9, -t67 * mrSges(5,1) - qJDD(3) * mrSges(5,3) - qJD(3) * t74 - t64 * t119 - t143, t11, -t40 * mrSges(7,3) - t62 * t42 + t116; -m(1) * g(3) + t112, t112, t6, t142, t66 * mrSges(5,1) + qJDD(3) * mrSges(5,2) + qJD(3) * t73 + t64 * t120 - t103, t99, -t39 * mrSges(7,1) - t61 * t45 + t114;];
f_new  = t3;
