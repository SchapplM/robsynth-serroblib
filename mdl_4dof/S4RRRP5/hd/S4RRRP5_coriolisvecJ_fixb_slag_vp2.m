% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRRP5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP5_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP5_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP5_coriolisvecJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP5_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:16:38
% EndTime: 2019-12-31 17:16:42
% DurationCPUTime: 1.61s
% Computational Cost: add. (1236->212), mult. (3357->287), div. (0->0), fcn. (1998->4), ass. (0->103)
t135 = Ifges(4,1) + Ifges(5,1);
t138 = Ifges(5,4) + Ifges(4,5);
t139 = mrSges(4,1) + mrSges(5,1);
t137 = -Ifges(4,6) + Ifges(5,6);
t90 = cos(qJ(2));
t104 = qJD(1) * t90;
t88 = sin(qJ(2));
t105 = qJD(1) * t88;
t87 = sin(qJ(3));
t89 = cos(qJ(3));
t62 = -t104 * t89 + t105 * t87;
t116 = t62 * Ifges(5,5);
t59 = Ifges(4,4) * t62;
t69 = t87 * t90 + t89 * t88;
t63 = t69 * qJD(1);
t86 = qJD(2) + qJD(3);
t136 = t135 * t63 + t138 * t86 + t116 - t59;
t110 = -Ifges(4,4) + Ifges(5,5);
t126 = -pkin(6) - pkin(5);
t78 = t126 * t90;
t72 = qJD(1) * t78;
t113 = t72 * t87;
t77 = t126 * t88;
t71 = qJD(1) * t77;
t66 = qJD(2) * pkin(2) + t71;
t42 = t66 * t89 + t113;
t134 = -t42 + qJD(4);
t132 = -t62 / 0.2e1;
t131 = t62 / 0.2e1;
t129 = t63 / 0.2e1;
t128 = -t86 / 0.2e1;
t84 = -pkin(2) * t90 - pkin(1);
t76 = qJD(1) * t84;
t125 = m(4) * t76;
t124 = pkin(1) * mrSges(3,1);
t123 = pkin(1) * mrSges(3,2);
t112 = t89 * t72;
t43 = t66 * t87 - t112;
t98 = qJD(2) * t126;
t94 = qJD(1) * t98;
t67 = t88 * t94;
t92 = t90 * t94;
t7 = qJD(3) * t43 + t67 * t87 - t89 * t92;
t93 = t89 * t77 + t78 * t87;
t122 = t93 * t7;
t121 = t69 * t7;
t120 = mrSges(5,2) * t63;
t119 = mrSges(4,3) * t62;
t118 = Ifges(3,4) * t88;
t21 = -pkin(3) * t86 + t134;
t117 = t21 * t62;
t115 = t63 * mrSges(4,3);
t114 = t63 * Ifges(4,4);
t111 = mrSges(5,2) + mrSges(4,3);
t50 = -mrSges(4,2) * t86 - t119;
t53 = -t62 * mrSges(5,2) + mrSges(5,3) * t86;
t109 = t50 + t53;
t108 = t139 * t86 - t115 - t120;
t107 = Ifges(3,5) * qJD(2);
t106 = Ifges(3,6) * qJD(2);
t103 = qJD(2) * mrSges(3,1);
t102 = qJD(2) * mrSges(3,2);
t101 = qJD(2) * t88;
t100 = qJD(3) * t89;
t99 = pkin(2) * t101;
t97 = t107 / 0.2e1;
t96 = -t106 / 0.2e1;
t95 = t90 * t98;
t36 = pkin(3) * t63 + qJ(4) * t62;
t49 = t77 * t87 - t78 * t89;
t68 = t87 * t88 - t89 * t90;
t6 = qJD(3) * t113 + t66 * t100 + t89 * t67 + t87 * t92;
t47 = t86 * t69;
t46 = t86 * t68;
t20 = t62 * pkin(3) - t63 * qJ(4) + t76;
t58 = Ifges(5,5) * t63;
t25 = t86 * Ifges(5,6) + t62 * Ifges(5,3) + t58;
t26 = -t62 * Ifges(4,2) + t86 * Ifges(4,6) + t114;
t29 = qJ(4) * t86 + t43;
t34 = t46 * qJD(1);
t35 = t47 * qJD(1);
t4 = qJD(4) * t86 + t6;
t91 = -t6 * mrSges(4,2) + t4 * mrSges(5,3) + t29 * t120 + t26 * t129 - t20 * (mrSges(5,1) * t63 + mrSges(5,3) * t62) - t76 * (mrSges(4,1) * t63 - mrSges(4,2) * t62) + (Ifges(5,3) * t63 - t116) * t132 - t42 * t119 - t139 * t7 + t137 * t35 - t138 * t34 + (t137 * t63 - t138 * t62) * t128 + (-Ifges(4,2) * t63 + t136 - t59) * t131 - (-t135 * t62 - t114 + t25 + t58) * t63 / 0.2e1;
t85 = Ifges(3,4) * t104;
t83 = -pkin(2) * t89 - pkin(3);
t81 = pkin(2) * t87 + qJ(4);
t79 = pkin(2) * t100 + qJD(4);
t75 = mrSges(3,3) * t104 - t102;
t74 = -mrSges(3,3) * t105 + t103;
t73 = t88 * t98;
t61 = Ifges(3,1) * t105 + t107 + t85;
t60 = t106 + (Ifges(3,2) * t90 + t118) * qJD(1);
t45 = t71 * t89 + t113;
t44 = t71 * t87 - t112;
t41 = t68 * pkin(3) - t69 * qJ(4) + t84;
t40 = mrSges(4,1) * t62 + mrSges(4,2) * t63;
t39 = mrSges(5,1) * t62 - mrSges(5,3) * t63;
t22 = pkin(2) * t105 + t36;
t9 = qJD(3) * t49 + t73 * t87 - t89 * t95;
t8 = qJD(3) * t93 + t89 * t73 + t87 * t95;
t2 = pkin(3) * t47 + qJ(4) * t46 - qJD(4) * t69 + t99;
t1 = pkin(3) * t35 + qJ(4) * t34 + qJD(1) * t99 - qJD(4) * t63;
t3 = [t1 * (t68 * mrSges(5,1) - t69 * mrSges(5,3)) + t2 * t39 - t108 * t9 + t109 * t8 + m(4) * (-t42 * t9 + t43 * t8 + t49 * t6 - t122) + m(5) * (t1 * t41 + t2 * t20 + t21 * t9 + t29 * t8 + t4 * t49 - t122) + (-t6 * t68 + t121) * mrSges(4,3) + (-t4 * t68 + t121) * mrSges(5,2) + (t61 / 0.2e1 - pkin(5) * t74 + t97 + (-0.2e1 * t123 + 0.3e1 / 0.2e1 * Ifges(3,4) * t90) * qJD(1)) * t90 * qJD(2) + (mrSges(4,1) * t84 + mrSges(5,1) * t41 + t110 * t69 + (Ifges(4,2) + Ifges(5,3)) * t68 - t111 * t49) * t35 - (mrSges(4,2) * t84 - mrSges(5,3) * t41 + t110 * t68 - t111 * t93 + t135 * t69) * t34 + (-t60 / 0.2e1 - pkin(5) * t75 + t96 + (-0.2e1 * t124 - 0.3e1 / 0.2e1 * t118 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t90) * qJD(1) + (0.2e1 * t125 + t40 + qJD(1) * (mrSges(4,1) * t68 + mrSges(4,2) * t69)) * pkin(2)) * t101 - (t76 * mrSges(4,2) + t21 * mrSges(5,2) - t42 * mrSges(4,3) - t20 * mrSges(5,3) + Ifges(4,4) * t132 + Ifges(5,5) * t131) * t46 + (t76 * mrSges(4,1) + Ifges(5,3) * t131 - Ifges(4,2) * t132 + t20 * mrSges(5,1) + t25 / 0.2e1 - t26 / 0.2e1 - t43 * mrSges(4,3) - t29 * mrSges(5,2)) * t47 - t136 * t46 / 0.2e1 + (t110 * t47 - t135 * t46) * t129 + (t137 * t47 - t138 * t46) * t86 / 0.2e1; -t109 * t45 + t108 * t44 + m(5) * (t29 * t79 + t4 * t81 + t7 * t83) - m(4) * (-t42 * t44 + t43 * t45) + ((t97 - t85 / 0.2e1 - t61 / 0.2e1 + qJD(1) * t123 + (t74 - t103) * pkin(5)) * t90 + (t96 + t60 / 0.2e1 + (t124 + t118 / 0.2e1 + (Ifges(3,2) / 0.2e1 - Ifges(3,1) / 0.2e1) * t90) * qJD(1) + (t75 + t102) * pkin(5) + (-t40 - t125) * pkin(2)) * t88) * qJD(1) + (m(4) * (t6 * t87 - t7 * t89) + (t34 * t89 - t35 * t87) * mrSges(4,3) + ((m(4) * t43 + t50) * t89 + (-m(4) * t42 + m(5) * t21 - t108) * t87) * qJD(3)) * pkin(2) - m(5) * (t20 * t22 + t21 * t44 + t29 * t45) + (-t34 * t83 - t35 * t81 + t117) * mrSges(5,2) + t91 + t79 * t53 + t43 * t115 - t22 * t39; -t109 * t42 + (t108 + t115) * t43 + t91 + (pkin(3) * t34 - qJ(4) * t35 + t117) * mrSges(5,2) + qJD(4) * t53 - t36 * t39 + (-pkin(3) * t7 + qJ(4) * t4 + t134 * t29 - t20 * t36 - t21 * t43) * m(5); -t34 * mrSges(5,2) + t63 * t39 - t86 * t53 + 0.2e1 * (t7 / 0.2e1 + t20 * t129 + t29 * t128) * m(5);];
tauc = t3(:);
