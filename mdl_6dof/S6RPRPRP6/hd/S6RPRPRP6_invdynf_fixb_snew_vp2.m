% Calculate vector of cutting forces with Newton-Euler
% S6RPRPRP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
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
% Datum: 2019-05-05 17:55
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRPRP6_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP6_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP6_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP6_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP6_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP6_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP6_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP6_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP6_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:52:28
% EndTime: 2019-05-05 17:52:32
% DurationCPUTime: 1.47s
% Computational Cost: add. (13851->188), mult. (33223->223), div. (0->0), fcn. (23162->8), ass. (0->90)
t149 = -2 * qJD(4);
t95 = sin(qJ(1));
t97 = cos(qJ(1));
t119 = t95 * g(1) - t97 * g(2);
t114 = qJDD(2) - t119;
t92 = cos(pkin(9));
t89 = t92 ^ 2;
t91 = sin(pkin(9));
t129 = t91 ^ 2 + t89;
t99 = qJD(1) ^ 2;
t101 = (-pkin(2) * t92 - pkin(1)) * qJDD(1) + (-t129 * pkin(7) - qJ(2)) * t99 + t114;
t139 = cos(qJ(3));
t94 = sin(qJ(3));
t107 = t139 * t91 + t92 * t94;
t79 = t107 * qJD(1);
t126 = t79 * qJD(3);
t144 = -t92 * t139 + t91 * t94;
t78 = t144 * qJD(1);
t127 = t78 * qJD(3);
t62 = t107 * qJDD(1) - t127;
t100 = pkin(3) * t126 + t79 * t149 + (-t62 + t127) * qJ(4) + t101;
t61 = t144 * qJDD(1) + t126;
t71 = t79 * pkin(4) - qJD(3) * pkin(8);
t77 = t78 ^ 2;
t17 = -t77 * pkin(4) - t79 * t71 + (pkin(3) + pkin(8)) * t61 + t100;
t125 = qJD(1) * qJD(2);
t118 = -t92 * g(3) - 0.2e1 * t91 * t125;
t128 = pkin(7) * qJDD(1);
t140 = pkin(2) * t99;
t115 = -t97 * g(1) - t95 * g(2);
t80 = -t99 * pkin(1) + qJDD(1) * qJ(2) + t115;
t48 = (t92 * t140 - t128 - t80) * t91 + t118;
t116 = -t91 * g(3) + (0.2e1 * t125 + t80) * t92;
t49 = t92 * t128 - t89 * t140 + t116;
t113 = t139 * t48 - t94 * t49;
t54 = t78 * pkin(3) - t79 * qJ(4);
t98 = qJD(3) ^ 2;
t26 = -qJDD(3) * pkin(3) - t98 * qJ(4) + t79 * t54 + qJDD(4) - t113;
t20 = (t78 * t79 - qJDD(3)) * pkin(8) + (t62 + t127) * pkin(4) + t26;
t93 = sin(qJ(5));
t96 = cos(qJ(5));
t134 = t96 * t17 + t93 * t20;
t65 = t96 * qJD(3) + t93 * t78;
t33 = -t65 * qJD(5) - t93 * qJDD(3) + t96 * t61;
t64 = -t93 * qJD(3) + t96 * t78;
t36 = -t64 * mrSges(7,1) + t65 * mrSges(7,2);
t75 = qJD(5) + t79;
t43 = t75 * pkin(5) - t65 * qJ(6);
t63 = t64 ^ 2;
t122 = m(7) * (-t63 * pkin(5) + t33 * qJ(6) + 0.2e1 * qJD(6) * t64 - t75 * t43 + t134) + t64 * t36 + t33 * mrSges(7,3);
t37 = -t64 * mrSges(6,1) + t65 * mrSges(6,2);
t44 = t75 * mrSges(7,1) - t65 * mrSges(7,3);
t45 = t75 * mrSges(6,1) - t65 * mrSges(6,3);
t59 = qJDD(5) + t62;
t10 = m(6) * t134 + t33 * mrSges(6,3) + t64 * t37 + (-t45 - t44) * t75 + (-mrSges(6,2) - mrSges(7,2)) * t59 + t122;
t70 = t79 * mrSges(5,1) + qJD(3) * mrSges(5,2);
t117 = -t93 * t17 + t96 * t20;
t34 = t64 * qJD(5) + t96 * qJDD(3) + t93 * t61;
t41 = -t75 * mrSges(7,2) + t64 * mrSges(7,3);
t123 = m(7) * (-0.2e1 * qJD(6) * t65 + (t64 * t75 - t34) * qJ(6) + (t64 * t65 + t59) * pkin(5) + t117) + t75 * t41 + t59 * mrSges(7,1);
t42 = -t75 * mrSges(6,2) + t64 * mrSges(6,3);
t8 = m(6) * t117 + t59 * mrSges(6,1) + t75 * t42 + (-t37 - t36) * t65 + (-mrSges(6,3) - mrSges(7,3)) * t34 + t123;
t109 = -t96 * t10 + t93 * t8 - m(5) * (t61 * pkin(3) + t100) + t79 * t70 + t62 * mrSges(5,3);
t69 = t78 * mrSges(5,1) - qJD(3) * mrSges(5,3);
t130 = -qJD(3) * mrSges(4,2) - t78 * mrSges(4,3) - t69;
t137 = mrSges(4,1) - mrSges(5,2);
t68 = qJD(3) * mrSges(4,1) - t79 * mrSges(4,3);
t102 = m(4) * t101 + t62 * mrSges(4,2) + t130 * t78 + t137 * t61 + t79 * t68 - t109;
t148 = -t102 - m(3) * (-qJDD(1) * pkin(1) - t99 * qJ(2) + t114);
t132 = t139 * t49 + t94 * t48;
t145 = t98 * pkin(3) - qJDD(3) * qJ(4) + qJD(3) * t149 + t78 * t54 - t132;
t104 = -t61 * pkin(4) - t77 * pkin(8) + qJD(3) * t71 - t145;
t121 = m(7) * (-t33 * pkin(5) - t63 * qJ(6) + t65 * t43 + qJDD(6) + t104) + t34 * mrSges(7,2) + t65 * t44;
t103 = m(6) * t104 + t34 * mrSges(6,2) + (-t42 - t41) * t64 - (mrSges(6,1) + mrSges(7,1)) * t33 + t65 * t45 + t121;
t147 = -m(5) * t145 + t103;
t146 = t129 * mrSges(3,3);
t56 = -t78 * mrSges(5,2) - t79 * mrSges(5,3);
t131 = -t78 * mrSges(4,1) - t79 * mrSges(4,2) - t56;
t135 = -mrSges(4,3) - mrSges(5,1);
t11 = m(4) * t132 + t131 * t78 + t135 * t61 + (-mrSges(4,2) + mrSges(5,3)) * qJDD(3) + (-t68 + t70) * qJD(3) + t147;
t112 = -t92 * mrSges(3,1) + t91 * mrSges(3,2);
t110 = qJDD(1) * mrSges(3,3) + t99 * t112;
t106 = -m(5) * t26 - t93 * t10 - t96 * t8;
t7 = m(4) * t113 + t130 * qJD(3) + t137 * qJDD(3) + t131 * t79 + t135 * t62 + t106;
t4 = m(3) * t118 + t94 * t11 + t139 * t7 + (-m(3) * t80 - t110) * t91;
t5 = m(3) * t116 + t139 * t11 + t110 * t92 - t94 * t7;
t143 = t92 * t4 + t91 * t5;
t6 = m(2) * t119 + (-mrSges(2,2) + t146) * t99 + (mrSges(2,1) - t112) * qJDD(1) + t148;
t1 = m(2) * t115 - t99 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t91 * t4 + t92 * t5;
t2 = [-m(1) * g(1) + t97 * t1 - t95 * t6, t1, t5, t11, -t61 * mrSges(5,2) - t78 * t69 - t109, t10, -t59 * mrSges(7,2) - t75 * t44 + t122; -m(1) * g(2) + t95 * t1 + t97 * t6, t6, t4, t7, t61 * mrSges(5,1) - qJDD(3) * mrSges(5,3) - qJD(3) * t70 + t78 * t56 - t147, t8, -t34 * mrSges(7,3) - t65 * t36 + t123; (-m(1) - m(2)) * g(3) + t143, -m(2) * g(3) + t143, t112 * qJDD(1) - t99 * t146 - t148, t102, t62 * mrSges(5,1) + qJDD(3) * mrSges(5,2) + qJD(3) * t69 + t79 * t56 - t106, t103, -t33 * mrSges(7,1) - t64 * t41 + t121;];
f_new  = t2;
