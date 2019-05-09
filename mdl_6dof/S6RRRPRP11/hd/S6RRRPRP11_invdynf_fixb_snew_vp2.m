% Calculate vector of cutting forces with Newton-Euler
% S6RRRPRP11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
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
% Datum: 2019-05-07 09:21
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRPRP11_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP11_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP11_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP11_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP11_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP11_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP11_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP11_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP11_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 09:07:40
% EndTime: 2019-05-07 09:07:47
% DurationCPUTime: 2.64s
% Computational Cost: add. (32610->209), mult. (70044->265), div. (0->0), fcn. (53177->10), ass. (0->101)
t105 = sin(qJ(3));
t151 = cos(qJ(3));
t102 = sin(pkin(6));
t106 = sin(qJ(2));
t109 = cos(qJ(2));
t133 = qJD(1) * t109;
t103 = cos(pkin(6));
t111 = qJD(1) ^ 2;
t107 = sin(qJ(1));
t110 = cos(qJ(1));
t125 = t107 * g(1) - t110 * g(2);
t150 = pkin(8) * t102;
t85 = qJDD(1) * pkin(1) + t111 * t150 + t125;
t138 = t103 * t85;
t122 = -t110 * g(1) - t107 * g(2);
t86 = -t111 * pkin(1) + qJDD(1) * t150 + t122;
t139 = t106 * t138 + t109 * t86;
t134 = qJD(1) * t102;
t88 = (-pkin(2) * t109 - pkin(9) * t106) * t134;
t99 = t103 * qJD(1) + qJD(2);
t97 = t99 ^ 2;
t98 = t103 * qJDD(1) + qJDD(2);
t41 = -t97 * pkin(2) + t98 * pkin(9) + (-g(3) * t106 + t88 * t133) * t102 + t139;
t149 = t103 * g(3);
t131 = qJD(1) * qJD(2);
t89 = (qJDD(1) * t106 + t109 * t131) * t102;
t90 = (-qJDD(1) * t109 + t106 * t131) * t102;
t42 = t90 * pkin(2) - t89 * pkin(9) - t149 + (-t85 + (pkin(2) * t106 - pkin(9) * t109) * t99 * qJD(1)) * t102;
t143 = t105 * t42 + t151 * t41;
t127 = t106 * t134;
t77 = t105 * t127 - t151 * t99;
t78 = t105 * t99 + t151 * t127;
t59 = t77 * pkin(3) - t78 * qJ(4);
t82 = qJDD(3) + t90;
t126 = t102 * t133;
t95 = -qJD(3) + t126;
t94 = t95 ^ 2;
t154 = t94 * pkin(3) - t82 * qJ(4) + 0.2e1 * qJD(4) * t95 + t77 * t59 - t143;
t57 = t78 * qJD(3) + t105 * t89 - t151 * t98;
t70 = t78 * pkin(4) + t95 * pkin(10);
t76 = t77 ^ 2;
t115 = -t57 * pkin(4) - t76 * pkin(10) - t95 * t70 - t154;
t104 = sin(qJ(5));
t108 = cos(qJ(5));
t65 = t104 * t77 - t108 * t95;
t32 = -t65 * qJD(5) - t104 * t82 + t108 * t57;
t64 = t104 * t95 + t108 * t77;
t33 = t64 * qJD(5) + t104 * t57 + t108 * t82;
t75 = qJD(5) + t78;
t50 = pkin(5) * t75 - t65 * qJ(6);
t51 = mrSges(7,1) * t75 - t65 * mrSges(7,3);
t63 = t64 ^ 2;
t128 = m(7) * (-t32 * pkin(5) - t63 * qJ(6) + t65 * t50 + qJDD(6) + t115) + t33 * mrSges(7,2) + t65 * t51;
t48 = -mrSges(7,2) * t75 + t64 * mrSges(7,3);
t49 = -mrSges(6,2) * t75 + t64 * mrSges(6,3);
t52 = mrSges(6,1) * t75 - t65 * mrSges(6,3);
t114 = m(6) * t115 + t33 * mrSges(6,2) + (-t49 - t48) * t64 - (mrSges(6,1) + mrSges(7,1)) * t32 + t65 * t52 + t128;
t155 = -m(5) * t154 + t114;
t148 = t77 * t95;
t123 = -t105 * t41 + t151 * t42;
t26 = -t82 * pkin(3) - t94 * qJ(4) + t78 * t59 + qJDD(4) - t123;
t58 = -t77 * qJD(3) + t105 * t98 + t151 * t89;
t20 = (t77 * t78 - t82) * pkin(10) + (t58 - t148) * pkin(4) + t26;
t135 = t102 * t109;
t120 = -g(3) * t135 - t106 * t86 + t109 * t138;
t40 = -t98 * pkin(2) - t97 * pkin(9) + t88 * t127 - t120;
t112 = (-t58 - t148) * qJ(4) + t40 + (-t95 * pkin(3) - 0.2e1 * qJD(4)) * t78;
t24 = -t76 * pkin(4) - t78 * t70 + (pkin(3) + pkin(10)) * t57 + t112;
t124 = -t104 * t24 + t108 * t20;
t55 = qJDD(5) + t58;
t130 = m(7) * (-0.2e1 * qJD(6) * t65 + (t64 * t75 - t33) * qJ(6) + (t64 * t65 + t55) * pkin(5) + t124) + t75 * t48 + t55 * mrSges(7,1);
t44 = -mrSges(7,1) * t64 + mrSges(7,2) * t65;
t45 = -mrSges(6,1) * t64 + mrSges(6,2) * t65;
t11 = m(6) * t124 + t55 * mrSges(6,1) + t75 * t49 + (-t45 - t44) * t65 + (-mrSges(6,3) - mrSges(7,3)) * t33 + t130;
t144 = t104 * t20 + t108 * t24;
t129 = m(7) * (-t63 * pkin(5) + t32 * qJ(6) + 0.2e1 * qJD(6) * t64 - t50 * t75 + t144) + t32 * mrSges(7,3) + t64 * t44;
t13 = m(6) * t144 + t32 * mrSges(6,3) + t64 * t45 + (-t52 - t51) * t75 + (-mrSges(6,2) - mrSges(7,2)) * t55 + t129;
t69 = t78 * mrSges(5,1) - t95 * mrSges(5,2);
t118 = t104 * t11 - t108 * t13 - m(5) * (t57 * pkin(3) + t112) + t58 * mrSges(5,3) + t78 * t69;
t68 = t77 * mrSges(5,1) + t95 * mrSges(5,3);
t140 = t95 * mrSges(4,2) - t77 * mrSges(4,3) - t68;
t147 = mrSges(4,1) - mrSges(5,2);
t67 = -t95 * mrSges(4,1) - t78 * mrSges(4,3);
t153 = m(4) * t40 + t58 * mrSges(4,2) + t140 * t77 + t147 * t57 + t78 * t67 - t118;
t145 = -mrSges(4,3) - mrSges(5,1);
t61 = -t77 * mrSges(5,2) - t78 * mrSges(5,3);
t141 = -t77 * mrSges(4,1) - t78 * mrSges(4,2) - t61;
t136 = t102 * t106;
t10 = m(4) * t143 + (t67 - t69) * t95 + (-mrSges(4,2) + mrSges(5,3)) * t82 + t141 * t77 + t145 * t57 + t155;
t83 = t99 * mrSges(3,1) - mrSges(3,3) * t127;
t87 = (-mrSges(3,1) * t109 + mrSges(3,2) * t106) * t134;
t116 = -m(5) * t26 - t104 * t13 - t108 * t11;
t9 = m(4) * t123 - t140 * t95 + t141 * t78 + t145 * t58 + t147 * t82 + t116;
t4 = m(3) * (-g(3) * t136 + t139) - t90 * mrSges(3,3) - t98 * mrSges(3,2) + t87 * t126 - t99 * t83 + t151 * t10 - t105 * t9;
t84 = -t99 * mrSges(3,2) + mrSges(3,3) * t126;
t6 = m(3) * (-t102 * t85 - t149) + t89 * mrSges(3,2) + t90 * mrSges(3,1) + t105 * t10 + t151 * t9 + (t106 * t83 - t109 * t84) * t134;
t8 = m(3) * t120 + t98 * mrSges(3,1) - t89 * mrSges(3,3) - t87 * t127 + t99 * t84 - t153;
t132 = t103 * t6 + t8 * t135 + t4 * t136;
t2 = m(2) * t122 - t111 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t106 * t8 + t109 * t4;
t1 = m(2) * t125 + qJDD(1) * mrSges(2,1) - t111 * mrSges(2,2) - t102 * t6 + (t106 * t4 + t109 * t8) * t103;
t3 = [-m(1) * g(1) - t107 * t1 + t110 * t2, t2, t4, t10, -t57 * mrSges(5,2) - t77 * t68 - t118, t13, -t55 * mrSges(7,2) - t75 * t51 + t129; -m(1) * g(2) + t110 * t1 + t107 * t2, t1, t8, t9, t57 * mrSges(5,1) - t82 * mrSges(5,3) + t77 * t61 + t95 * t69 - t155, t11, -t33 * mrSges(7,3) - t65 * t44 + t130; (-m(1) - m(2)) * g(3) + t132, -m(2) * g(3) + t132, t6, t153, t58 * mrSges(5,1) + t82 * mrSges(5,2) + t78 * t61 - t95 * t68 - t116, t114, -t32 * mrSges(7,1) - t64 * t48 + t128;];
f_new  = t3;
