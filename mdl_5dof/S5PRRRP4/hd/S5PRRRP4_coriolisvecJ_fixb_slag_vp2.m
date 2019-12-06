% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRRP4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP4_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP4_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP4_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP4_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:45:37
% EndTime: 2019-12-05 16:45:41
% DurationCPUTime: 1.24s
% Computational Cost: add. (1310->199), mult. (2683->279), div. (0->0), fcn. (1522->6), ass. (0->97)
t150 = mrSges(5,1) + mrSges(6,1);
t78 = cos(qJ(4));
t118 = qJD(4) * t78;
t79 = cos(qJ(3));
t120 = qJD(3) * t79;
t76 = sin(qJ(3));
t121 = qJD(3) * t76;
t80 = cos(qJ(2));
t122 = qJD(1) * t80;
t65 = qJD(2) * pkin(2) + t122;
t77 = sin(qJ(2));
t58 = t76 * t77 - t79 * t80;
t89 = t58 * qJD(2);
t14 = t65 * t120 + (-t77 * t121 - t89) * qJD(1);
t123 = qJD(1) * t77;
t38 = t79 * t123 + t65 * t76;
t72 = qJD(2) + qJD(3);
t27 = pkin(7) * t72 + t38;
t75 = sin(qJ(4));
t7 = t27 * t118 + t14 * t75;
t138 = t7 * t75;
t11 = t78 * t14;
t129 = t27 * t75;
t4 = t11 + (qJD(5) - t129) * qJD(4);
t101 = t4 * t78 + t138;
t18 = -qJD(4) * pkin(4) + qJD(5) + t129;
t149 = t18 * t118 + t101;
t126 = t72 * t78;
t113 = mrSges(6,2) * t126;
t63 = qJD(4) * mrSges(6,3) + t113;
t124 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t126 + t63;
t127 = t72 * t75;
t125 = (-mrSges(6,2) - mrSges(5,3)) * t127 + t150 * qJD(4);
t148 = t124 * t75 + t125 * t78;
t147 = -mrSges(5,1) * t78 + mrSges(5,2) * t75 - mrSges(4,1);
t59 = t76 * t80 + t77 * t79;
t90 = t59 * qJD(2);
t15 = t65 * t121 + (t77 * t120 + t90) * qJD(1);
t92 = (mrSges(5,1) * t75 + mrSges(5,2) * t78) * qJD(4);
t44 = t72 * t92;
t146 = m(5) * t15 + t44;
t119 = qJD(4) * t75;
t6 = -t27 * t119 + t11;
t140 = t6 * t78;
t100 = t138 + t140;
t115 = qJD(4) * qJ(5);
t19 = t27 * t78 + t115;
t116 = t19 * qJD(4);
t145 = m(6) * (-t75 * t116 + t149) + m(5) * t100;
t95 = t18 * t75 + t19 * t78;
t110 = (t75 ^ 2 + t78 ^ 2) * t27;
t142 = pkin(2) * t79;
t139 = t7 * mrSges(5,3);
t137 = Ifges(5,4) * t75;
t135 = Ifges(6,5) * t75;
t134 = Ifges(6,5) * t78;
t133 = t15 * t58;
t131 = t19 * mrSges(6,2);
t117 = qJD(5) * t75;
t111 = t127 / 0.2e1;
t108 = t119 / 0.2e1;
t106 = t147 * t72;
t68 = t76 * t123;
t37 = t65 * t79 - t68;
t105 = t125 * qJD(4);
t104 = t72 * t108;
t98 = -mrSges(6,1) * t78 - mrSges(6,3) * t75;
t50 = t98 * t72;
t102 = -t106 - t50;
t97 = pkin(4) * t75 - qJ(5) * t78;
t96 = t18 * t78 - t19 * t75;
t64 = -pkin(4) * t78 - qJ(5) * t75 - pkin(3);
t94 = (Ifges(6,1) * t75 - t134) * t72;
t93 = (Ifges(5,2) * t78 + t137) * t72;
t91 = (mrSges(6,1) * t75 - mrSges(6,3) * t78) * qJD(4);
t49 = pkin(4) * t119 - t78 * t115 - t117;
t84 = -t72 * mrSges(4,2) + t124 * t78 - t125 * t75;
t13 = t64 * t72 - t37;
t26 = -pkin(3) * t72 - t37;
t66 = Ifges(6,5) * t127;
t39 = Ifges(6,6) * qJD(4) - Ifges(6,3) * t126 + t66;
t40 = Ifges(5,6) * qJD(4) + t93;
t41 = Ifges(6,4) * qJD(4) + t94;
t67 = Ifges(5,4) * t126;
t42 = Ifges(5,1) * t127 + Ifges(5,5) * qJD(4) + t67;
t8 = (t97 * qJD(4) - t117) * t72 + t15;
t82 = t8 * t98 - t72 * (Ifges(6,3) * t75 + t134) * t118 + (-Ifges(6,3) * t78 + t135) * t104 + t13 * t91 + t26 * t92 + t39 * t108 + mrSges(5,3) * t140 + t147 * t15 + 0.2e1 * (t135 - t137 + (Ifges(5,1) + Ifges(6,1)) * t78) * t104 + ((Ifges(6,4) + Ifges(5,5)) * t78 + (-Ifges(5,6) + Ifges(6,6)) * t75) * qJD(4) ^ 2 / 0.2e1 - (t93 + t40) * t119 / 0.2e1 + t149 * mrSges(6,2) + ((0.3e1 * Ifges(5,4) * t78 + (Ifges(5,1) - 0.2e1 * Ifges(5,2)) * t75) * t72 + t94 + t42 + t41) * t118 / 0.2e1;
t69 = pkin(2) * t76 + pkin(7);
t57 = t64 - t142;
t56 = t79 * t122 - t68;
t55 = t59 * qJD(1);
t52 = t97 * t72;
t43 = t72 * t91;
t32 = pkin(2) * t121 + t49;
t21 = t59 * qJD(3) + t90;
t20 = -t58 * qJD(3) - t89;
t1 = [(-mrSges(3,1) * t77 - mrSges(3,2) * t80) * qJD(2) ^ 2 + (t43 + t44) * t58 - t102 * t21 - t148 * t59 * qJD(4) + t84 * t20 + m(4) * (t14 * t59 + t20 * t38 - t21 * t37 + t133) + m(5) * (t100 * t59 + t20 * t110 + t21 * t26 + t133) + m(6) * (t13 * t21 + t58 * t8 + t95 * t20 + (t96 * qJD(4) + t101) * t59); t82 + (m(5) * t110 * t120 + (m(4) * t14 + (-m(4) * t37 + m(5) * t26 + t106) * qJD(3)) * t76 + (-m(4) * t15 + (m(4) * t38 + m(6) * t95 + t84) * qJD(3)) * t79) * pkin(2) - m(6) * (t13 * t55 + t95 * t56) - m(5) * (t56 * t110 + t26 * t55) + m(6) * (t13 * t32 + t57 * t8) + t69 * t145 + (t139 + t125 * t56 + (-t124 * t69 - t131) * qJD(4)) * t75 + (-t69 * t105 - t124 * t56) * t78 + t57 * t43 + t32 * t50 + t102 * t55 - m(4) * (-t37 * t55 + t38 * t56) + (t56 * t72 - t14) * mrSges(4,2) + t146 * (-pkin(3) - t142); t82 + pkin(7) * t145 + (t139 + t125 * t37 + (-t124 * pkin(7) - t131) * qJD(4)) * t75 + (-pkin(7) * t105 - t124 * t37) * t78 + t102 * t38 + (t37 * t72 - t14) * mrSges(4,2) - m(5) * (t37 * t110 + t26 * t38) + t64 * t43 + t49 * t50 - t146 * pkin(3) + (-t95 * t37 + t8 * t64 + (-t38 + t49) * t13) * m(6); -t6 * mrSges(5,2) + t4 * mrSges(6,3) + qJD(5) * t63 - t52 * t50 - t150 * t7 + t148 * t27 + ((-t26 * mrSges(5,2) + t13 * mrSges(6,3) - t41 / 0.2e1 - t42 / 0.2e1 - t67 / 0.2e1 - t18 * mrSges(6,2) + Ifges(6,5) * t126 / 0.2e1) * t78 + (-t26 * mrSges(5,1) - t13 * mrSges(6,1) - t66 / 0.2e1 - t39 / 0.2e1 + t40 / 0.2e1 + t131 + Ifges(5,4) * t111 + (-Ifges(6,1) / 0.2e1 + Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1 - Ifges(5,1) / 0.2e1) * t126) * t75 + ((Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1 - pkin(4) * mrSges(6,2)) * t78 + (Ifges(6,6) / 0.2e1 - Ifges(5,6) / 0.2e1 - qJ(5) * mrSges(6,2)) * t75) * qJD(4)) * t72 + (-t7 * pkin(4) + t4 * qJ(5) + t19 * qJD(5) - t13 * t52 - t96 * t27) * m(6); t50 * t127 + (-t63 + t113) * qJD(4) + 0.2e1 * (t7 / 0.2e1 + t13 * t111 - t116 / 0.2e1) * m(6);];
tauc = t1(:);
