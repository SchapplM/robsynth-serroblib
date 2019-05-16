% Calculate vector of cutting forces with Newton-Euler
% S6RPRRPP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
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
% Datum: 2019-05-05 21:40
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRRPP5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP5_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP5_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP5_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP5_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP5_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP5_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP5_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:36:15
% EndTime: 2019-05-05 21:36:19
% DurationCPUTime: 1.62s
% Computational Cost: add. (15844->192), mult. (37159->223), div. (0->0), fcn. (26578->8), ass. (0->92)
t142 = cos(qJ(3));
t95 = sin(pkin(9));
t96 = cos(pkin(9));
t98 = sin(qJ(3));
t150 = t96 * t142 - t95 * t98;
t102 = qJD(1) ^ 2;
t94 = t96 ^ 2;
t130 = t95 ^ 2 + t94;
t149 = t130 * mrSges(3,3);
t141 = cos(qJ(4));
t101 = qJD(3) ^ 2;
t126 = qJD(1) * qJD(2);
t121 = -t96 * g(3) - 0.2e1 * t95 * t126;
t129 = pkin(7) * qJDD(1);
t140 = pkin(2) * t102;
t100 = cos(qJ(1));
t99 = sin(qJ(1));
t115 = -t100 * g(1) - t99 * g(2);
t84 = -t102 * pkin(1) + qJDD(1) * qJ(2) + t115;
t58 = (t96 * t140 - t129 - t84) * t95 + t121;
t119 = -t95 * g(3) + (0.2e1 * t126 + t84) * t96;
t59 = t96 * t129 - t94 * t140 + t119;
t131 = t142 * t59 + t98 * t58;
t82 = t150 * qJD(1);
t111 = t142 * t95 + t96 * t98;
t83 = t111 * qJD(1);
t68 = -t82 * pkin(3) - t83 * pkin(8);
t25 = -t101 * pkin(3) + qJDD(3) * pkin(8) + t82 * t68 + t131;
t122 = t99 * g(1) - t100 * g(2);
t117 = qJDD(2) - t122;
t103 = (-pkin(2) * t96 - pkin(1)) * qJDD(1) + (-t130 * pkin(7) - qJ(2)) * t102 + t117;
t127 = t83 * qJD(3);
t128 = t82 * qJD(3);
t70 = t150 * qJDD(1) - t127;
t71 = t111 * qJDD(1) + t128;
t27 = (-t71 - t128) * pkin(8) + (-t70 + t127) * pkin(3) + t103;
t97 = sin(qJ(4));
t135 = t141 * t25 + t97 * t27;
t144 = 2 * qJD(5);
t73 = -t141 * qJD(3) + t97 * t83;
t74 = t97 * qJD(3) + t141 * t83;
t42 = t73 * pkin(4) - t74 * qJ(5);
t67 = qJDD(4) - t70;
t80 = qJD(4) - t82;
t79 = t80 ^ 2;
t110 = -t79 * pkin(4) + t67 * qJ(5) + t80 * t144 - t73 * t42 + t135;
t55 = -t80 * mrSges(6,1) + t74 * mrSges(6,2);
t148 = m(6) * t110 + t67 * mrSges(6,3) + t80 * t55;
t132 = t142 * t58 - t98 * t59;
t107 = qJDD(3) * pkin(3) + t101 * pkin(8) - t83 * t68 + t132;
t139 = t73 * t80;
t39 = -t73 * qJD(4) + t97 * qJDD(3) + t141 * t71;
t147 = (-t39 + t139) * qJ(5) - t107;
t38 = t74 * qJD(4) - t141 * qJDD(3) + t97 * t71;
t50 = t80 * mrSges(7,2) + t73 * mrSges(7,3);
t52 = -t80 * pkin(5) - t74 * qJ(6);
t53 = -t80 * mrSges(7,1) - t74 * mrSges(7,3);
t72 = t73 ^ 2;
t118 = m(7) * (-t72 * qJ(6) + qJDD(6) + (-pkin(4) - pkin(5)) * t38 + (-pkin(4) * t80 + t144 + t52) * t74 - t147) + t39 * mrSges(7,2) - t38 * mrSges(7,1) + t74 * t53 - t73 * t50;
t145 = -0.2e1 * t74;
t49 = -t73 * mrSges(6,2) + t80 * mrSges(6,3);
t109 = m(6) * (qJD(5) * t145 + (t74 * t80 + t38) * pkin(4) + t147) + t38 * mrSges(6,1) + t73 * t49 - t118;
t51 = -t80 * mrSges(5,2) - t73 * mrSges(5,3);
t54 = t80 * mrSges(5,1) - t74 * mrSges(5,3);
t146 = -m(5) * t107 + t38 * mrSges(5,1) + (t54 - t55) * t74 + (mrSges(5,2) - mrSges(6,3)) * t39 + t73 * t51 + t109;
t114 = -t96 * mrSges(3,1) + t95 * mrSges(3,2);
t112 = qJDD(1) * mrSges(3,3) + t102 * t114;
t116 = t141 * t27 - t97 * t25;
t19 = -t67 * pkin(4) - t79 * qJ(5) + t74 * t42 + qJDD(5) - t116;
t124 = t80 * t50 + t67 * mrSges(7,1) - m(7) * (qJD(6) * t145 + (-t39 - t139) * qJ(6) + (t73 * t74 - t67) * pkin(5) + t19);
t113 = m(6) * t19 - t124;
t43 = t73 * mrSges(6,1) - t74 * mrSges(6,3);
t134 = -t73 * mrSges(5,1) - t74 * mrSges(5,2) - t43;
t136 = -mrSges(5,3) - mrSges(6,2);
t44 = -t73 * mrSges(7,1) + t74 * mrSges(7,2);
t10 = m(5) * t116 + (t51 + t49) * t80 + (mrSges(5,1) + mrSges(6,1)) * t67 + (t44 + t134) * t74 + (mrSges(7,3) + t136) * t39 - t113;
t125 = m(7) * (-t72 * pkin(5) + t38 * qJ(6) + 0.2e1 * qJD(6) * t73 + t80 * t52 + t110) + t73 * t44 + t38 * mrSges(7,3);
t12 = m(5) * t135 + (-t54 + t53) * t80 + t134 * t73 + (-mrSges(5,2) + mrSges(7,2)) * t67 + t136 * t38 + t125 + t148;
t65 = -t82 * mrSges(4,1) + t83 * mrSges(4,2);
t76 = qJD(3) * mrSges(4,1) - t83 * mrSges(4,3);
t7 = m(4) * t131 - qJDD(3) * mrSges(4,2) + t70 * mrSges(4,3) - qJD(3) * t76 - t97 * t10 + t141 * t12 + t82 * t65;
t75 = -qJD(3) * mrSges(4,2) + t82 * mrSges(4,3);
t8 = m(4) * t132 + qJDD(3) * mrSges(4,1) - t71 * mrSges(4,3) + qJD(3) * t75 - t83 * t65 - t146;
t4 = m(3) * t121 + t98 * t7 + t142 * t8 + (-m(3) * t84 - t112) * t95;
t5 = m(3) * t119 + t112 * t96 + t142 * t7 - t98 * t8;
t143 = t96 * t4 + t95 * t5;
t108 = t67 * mrSges(7,2) + t80 * t53 + t125;
t106 = -m(4) * t103 + t70 * mrSges(4,1) - t71 * mrSges(4,2) - t141 * t10 - t97 * t12 + t82 * t75 - t83 * t76;
t105 = m(3) * (-qJDD(1) * pkin(1) - t102 * qJ(2) + t117) - t106;
t6 = m(2) * t122 + (-mrSges(2,2) + t149) * t102 + (mrSges(2,1) - t114) * qJDD(1) - t105;
t1 = m(2) * t115 - t102 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t95 * t4 + t96 * t5;
t2 = [-m(1) * g(1) + t100 * t1 - t99 * t6, t1, t5, t7, t12, -t38 * mrSges(6,2) - t73 * t43 + t108 + t148, t108; -m(1) * g(2) + t99 * t1 + t100 * t6, t6, t4, t8, t10, -t39 * mrSges(6,3) - t74 * t55 + t109, -t39 * mrSges(7,3) - t74 * t44 - t124; (-m(1) - m(2)) * g(3) + t143, -m(2) * g(3) + t143, t114 * qJDD(1) - t102 * t149 + t105, -t106, t146, -t67 * mrSges(6,1) - t80 * t49 + (t43 - t44) * t74 + (mrSges(6,2) - mrSges(7,3)) * t39 + t113, t118;];
f_new  = t2;
