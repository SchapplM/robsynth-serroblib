% Calculate vector of cutting forces with Newton-Euler
% S6RPRRPP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta5]';
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
% Datum: 2019-05-05 21:45
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRRPP6_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP6_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP6_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP6_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP6_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP6_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP6_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP6_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP6_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:41:24
% EndTime: 2019-05-05 21:41:27
% DurationCPUTime: 1.42s
% Computational Cost: add. (17140->174), mult. (34129->214), div. (0->0), fcn. (21291->8), ass. (0->86)
t83 = sin(qJ(1));
t86 = cos(qJ(1));
t101 = -t86 * g(1) - t83 * g(2);
t124 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t101;
t123 = -2 * qJD(5);
t122 = -m(2) - m(3);
t121 = (-pkin(1) - pkin(7));
t114 = cos(pkin(9));
t85 = cos(qJ(3));
t112 = qJD(1) * t85;
t81 = sin(qJ(4));
t84 = cos(qJ(4));
t67 = qJD(3) * t84 - t81 * t112;
t68 = qJD(3) * t81 + t84 * t112;
t80 = sin(pkin(9));
t47 = -t114 * t67 + t68 * t80;
t48 = t114 * t68 + t80 * t67;
t28 = pkin(5) * t47 - qJ(6) * t48;
t111 = qJD(1) * qJD(3);
t103 = t85 * t111;
t82 = sin(qJ(3));
t71 = -qJDD(1) * t82 - t103;
t66 = qJDD(4) - t71;
t113 = qJD(1) * t82;
t76 = qJD(4) + t113;
t75 = t76 ^ 2;
t104 = t82 * t111;
t72 = qJDD(1) * t85 - t104;
t88 = qJD(1) ^ 2;
t92 = (t121 * t88) - t124;
t33 = (-t72 + t104) * pkin(8) + (-t71 + t103) * pkin(3) + t92;
t107 = t83 * g(1) - t86 * g(2);
t96 = -t88 * qJ(2) + qJDD(2) - t107;
t56 = t121 * qJDD(1) + t96;
t106 = -g(3) * t85 + t82 * t56;
t70 = (pkin(3) * t82 - pkin(8) * t85) * qJD(1);
t87 = qJD(3) ^ 2;
t37 = -pkin(3) * t87 + qJDD(3) * pkin(8) - t70 * t113 + t106;
t102 = t84 * t33 - t81 * t37;
t46 = qJD(4) * t67 + qJDD(3) * t81 + t72 * t84;
t17 = (t67 * t76 - t46) * qJ(5) + (t67 * t68 + t66) * pkin(4) + t102;
t115 = t81 * t33 + t84 * t37;
t45 = -qJD(4) * t68 + qJDD(3) * t84 - t72 * t81;
t51 = pkin(4) * t76 - qJ(5) * t68;
t65 = t67 ^ 2;
t19 = -pkin(4) * t65 + qJ(5) * t45 - t51 * t76 + t115;
t98 = t114 * t17 - t80 * t19;
t120 = m(7) * (-t66 * pkin(5) - t75 * qJ(6) + qJDD(6) + ((2 * qJD(5)) + t28) * t48 - t98);
t119 = (mrSges(2,1) - mrSges(3,2));
t118 = -mrSges(2,2) + mrSges(3,3);
t117 = -mrSges(6,3) - mrSges(7,2);
t29 = mrSges(7,1) * t47 - mrSges(7,3) * t48;
t116 = -mrSges(6,1) * t47 - mrSges(6,2) * t48 - t29;
t108 = t114 * t19 + t47 * t123 + t80 * t17;
t41 = -mrSges(7,1) * t76 + mrSges(7,2) * t48;
t109 = m(7) * (-pkin(5) * t75 + qJ(6) * t66 + 0.2e1 * qJD(6) * t76 - t28 * t47 + t108) + t76 * t41 + t66 * mrSges(7,3);
t100 = t82 * g(3) + t85 * t56;
t69 = (mrSges(4,1) * t82 + mrSges(4,2) * t85) * qJD(1);
t73 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t113;
t36 = -qJDD(3) * pkin(3) - t87 * pkin(8) + t70 * t112 - t100;
t50 = -mrSges(5,2) * t76 + mrSges(5,3) * t67;
t52 = mrSges(5,1) * t76 - mrSges(5,3) * t68;
t24 = -t114 * t45 + t80 * t46;
t25 = t114 * t46 + t80 * t45;
t38 = -mrSges(6,2) * t76 - mrSges(6,3) * t47;
t40 = mrSges(6,1) * t76 - mrSges(6,3) * t48;
t90 = -t45 * pkin(4) - t65 * qJ(5) + t68 * t51 + qJDD(5) + t36;
t39 = -mrSges(7,2) * t47 + mrSges(7,3) * t76;
t95 = t25 * mrSges(7,3) + t48 * t41 - m(7) * (-0.2e1 * qJD(6) * t48 + (t47 * t76 - t25) * qJ(6) + (t48 * t76 + t24) * pkin(5) + t90) - t24 * mrSges(7,1) - t47 * t39;
t91 = m(6) * t90 + t24 * mrSges(6,1) + t25 * mrSges(6,2) + t47 * t38 + t48 * t40 - t95;
t89 = m(5) * t36 - t45 * mrSges(5,1) + t46 * mrSges(5,2) - t67 * t50 + t68 * t52 + t91;
t11 = m(4) * t100 + qJDD(3) * mrSges(4,1) - t72 * mrSges(4,3) + qJD(3) * t73 - t69 * t112 - t89;
t10 = m(6) * t98 - t120 + (t38 + t39) * t76 + (mrSges(6,1) + mrSges(7,1)) * t66 + (m(6) * t123 + t116) * t48 + t117 * t25;
t49 = -mrSges(5,1) * t67 + mrSges(5,2) * t68;
t9 = m(6) * t108 - t66 * mrSges(6,2) + t116 * t47 + t117 * t24 - t76 * t40 + t109;
t7 = m(5) * t102 + t66 * mrSges(5,1) - t46 * mrSges(5,3) + t114 * t10 - t68 * t49 + t76 * t50 + t80 * t9;
t74 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t112;
t8 = m(5) * t115 - t66 * mrSges(5,2) + t45 * mrSges(5,3) - t80 * t10 + t114 * t9 + t67 * t49 - t76 * t52;
t4 = m(4) * t106 - qJDD(3) * mrSges(4,2) + t71 * mrSges(4,3) - qJD(3) * t74 - t69 * t113 - t81 * t7 + t84 * t8;
t105 = -t82 * t11 + t85 * t4;
t97 = -m(3) * (-qJDD(1) * pkin(1) + t96) - t85 * t11 - t82 * t4;
t94 = m(4) * t92 - t71 * mrSges(4,1) + t72 * mrSges(4,2) + t74 * t112 + t73 * t113 + t84 * t7 + t81 * t8;
t93 = -m(3) * (t88 * pkin(1) + t124) + t94;
t2 = m(2) * t101 + t118 * qJDD(1) - (t119 * t88) + t93;
t1 = m(2) * t107 + t119 * qJDD(1) + t118 * t88 + t97;
t3 = [-m(1) * g(1) - t1 * t83 + t2 * t86, t2, -m(3) * g(3) + t105, t4, t8, t9, -t24 * mrSges(7,2) - t47 * t29 + t109; -m(1) * g(2) + t1 * t86 + t2 * t83, t1, -(t88 * mrSges(3,2)) - qJDD(1) * mrSges(3,3) - t93, t11, t7, t10, -t95; (-m(1) + t122) * g(3) + t105, t122 * g(3) + t105, qJDD(1) * mrSges(3,2) - t88 * mrSges(3,3) - t97, t94, t89, t91, -t66 * mrSges(7,1) + t25 * mrSges(7,2) + t48 * t29 - t76 * t39 + t120;];
f_new  = t3;
