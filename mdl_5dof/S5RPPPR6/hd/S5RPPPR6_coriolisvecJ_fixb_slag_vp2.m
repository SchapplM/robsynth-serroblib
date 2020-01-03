% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
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
% Datum: 2019-12-31 17:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPPR6_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR6_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR6_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR6_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR6_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR6_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR6_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:47:37
% EndTime: 2019-12-31 17:47:41
% DurationCPUTime: 1.56s
% Computational Cost: add. (1115->181), mult. (2845->287), div. (0->0), fcn. (1802->6), ass. (0->102)
t71 = sin(pkin(7));
t68 = t71 ^ 2;
t73 = cos(pkin(7));
t69 = t73 ^ 2;
t134 = (mrSges(4,1) + mrSges(3,3)) * (t68 + t69);
t131 = pkin(3) + qJ(2);
t74 = sin(qJ(5));
t111 = t73 * t74;
t75 = cos(qJ(5));
t113 = t71 * t75;
t70 = sin(pkin(8));
t44 = t70 * t111 + t113;
t36 = t44 * qJD(1);
t100 = t71 * qJD(1);
t110 = t73 * t75;
t94 = t70 * t110;
t40 = qJD(1) * t94 - t74 * t100;
t84 = t70 * t73 * mrSges(5,3) + t71 * mrSges(5,1);
t107 = -mrSges(6,1) * t36 - mrSges(6,2) * t40 - t84 * qJD(1);
t90 = -qJ(3) * t71 - pkin(1);
t46 = (-pkin(2) - qJ(4)) * t73 + t90;
t26 = qJD(1) * t46 + qJD(2);
t96 = qJ(2) * qJD(1);
t58 = t71 * t96 + qJD(3);
t49 = pkin(3) * t100 + t58;
t72 = cos(pkin(8));
t14 = -t26 * t70 + t49 * t72;
t11 = -pkin(4) * t100 - t14;
t130 = -m(5) * t14 + m(6) * t11 + t107;
t98 = t73 * qJD(1);
t55 = t72 * t98 + qJD(5);
t19 = -mrSges(6,2) * t55 + mrSges(6,3) * t36;
t20 = mrSges(6,1) * t55 + mrSges(6,3) * t40;
t15 = t72 * t26 + t70 * t49;
t12 = pkin(6) * t100 + t15;
t50 = pkin(3) * t98 + t73 * t96 + qJD(4);
t81 = (pkin(4) * t72 + pkin(6) * t70) * t73;
t21 = qJD(1) * t81 + t50;
t5 = -t12 * t74 + t21 * t75;
t6 = t12 * t75 + t21 * t74;
t129 = m(6) * (-t5 * t74 + t6 * t75) + t75 * t19 - t74 * t20;
t112 = t72 * t73;
t128 = t112 / 0.2e1;
t126 = qJD(2) * t69;
t99 = t71 * qJD(3);
t54 = -qJD(4) * t73 - t99;
t52 = t54 * qJD(1);
t95 = qJD(1) * qJD(2);
t92 = t71 * t95;
t23 = t52 * t72 + t70 * t92;
t91 = t73 * t95;
t1 = t5 * qJD(5) + t23 * t75 + t74 * t91;
t2 = -t6 * qJD(5) - t23 * t74 + t75 * t91;
t122 = t2 * mrSges(6,1) - t1 * mrSges(6,2);
t120 = -t40 / 0.2e1;
t117 = Ifges(6,4) * t40;
t22 = t52 * t70 - t72 * t92;
t116 = t22 * t70;
t115 = t22 * t72;
t114 = t71 * t74;
t37 = t44 * qJD(5);
t31 = qJD(1) * t37;
t82 = t94 - t114;
t38 = t82 * qJD(5);
t32 = qJD(1) * t38;
t106 = Ifges(6,5) * t31 + Ifges(6,6) * t32;
t56 = t131 * t71;
t105 = t72 * t46 + t70 * t56;
t104 = t131 * t73;
t76 = qJD(1) ^ 2;
t102 = qJ(2) * t76;
t101 = qJD(2) * t71;
t97 = t73 * qJD(2);
t89 = qJ(2) * t95;
t18 = pkin(6) * t71 + t105;
t25 = t81 + t104;
t8 = t18 * t75 + t25 * t74;
t7 = -t18 * t74 + t25 * t75;
t86 = -t46 * t70 + t56 * t72;
t85 = -pkin(2) * t73 + t90;
t83 = -t71 * mrSges(5,2) - mrSges(5,3) * t112;
t80 = qJD(1) * (mrSges(5,1) * t72 - mrSges(5,2) * t70);
t53 = (mrSges(4,2) * t73 - mrSges(4,3) * t71) * qJD(1);
t79 = t1 * t75 - t2 * t74 + (-t5 * t75 - t6 * t74) * qJD(5);
t78 = (-t19 * t74 - t20 * t75) * qJD(5) + (t31 * t74 + t32 * t75) * mrSges(6,3);
t48 = t83 * qJD(1);
t77 = t130 * t70 + (m(5) * t15 + t129 + t48) * t72;
t63 = t69 * t102;
t60 = t68 * t89;
t43 = t85 * qJD(1) + qJD(2);
t41 = t73 * t80;
t39 = (t70 * t113 + t111) * qJD(1);
t35 = (-t70 * t114 + t110) * qJD(1);
t33 = Ifges(6,4) * t36;
t30 = t70 * t101 + t54 * t72;
t17 = -pkin(4) * t71 - t86;
t13 = -mrSges(6,1) * t32 + mrSges(6,2) * t31;
t10 = -Ifges(6,1) * t40 + Ifges(6,5) * t55 + t33;
t9 = Ifges(6,2) * t36 + Ifges(6,6) * t55 - t117;
t4 = -t8 * qJD(5) - t30 * t74 + t75 * t97;
t3 = t7 * qJD(5) + t30 * t75 + t74 * t97;
t16 = [t130 * (-t72 * t101 + t54 * t70) + m(5) * (t23 * t105 + t15 * t30 + (qJD(1) * t104 + t50) * t97) + m(6) * (t1 * t8 + t2 * t7 + t3 * t6 + t4 * t5) + 0.2e1 * t95 * t134 + (-m(5) * t86 + m(6) * t17 - t44 * mrSges(6,1) - mrSges(6,2) * t82 - t84) * t22 + (t1 * t44 + t2 * t82 - t5 * t37 + t6 * t38) * mrSges(6,3) + (t8 * mrSges(6,3) - Ifges(6,4) * t82 + Ifges(6,2) * t44 + Ifges(6,6) * t128) * t32 + (-t7 * mrSges(6,3) - Ifges(6,1) * t82 + Ifges(6,4) * t44 + Ifges(6,5) * t128) * t31 - 0.2e1 * t53 * t99 + 0.2e1 * m(3) * (t69 * t89 + t60) + t23 * t83 + t55 * (Ifges(6,5) * t37 + Ifges(6,6) * t38) / 0.2e1 + t30 * t48 + t37 * t10 / 0.2e1 + t36 * (Ifges(6,4) * t37 + Ifges(6,2) * t38) / 0.2e1 + t11 * (-t38 * mrSges(6,1) + t37 * mrSges(6,2)) + t38 * t9 / 0.2e1 + t17 * t13 + t3 * t19 + t4 * t20 + (t106 / 0.2e1 + t122) * t112 + (Ifges(6,1) * t37 + Ifges(6,4) * t38) * t120 + t41 * t97 + t80 * t126 + m(4) * (t60 + (t58 * qJD(2) - t43 * qJD(3)) * t71 + (0.2e1 * qJ(2) * t126 - t85 * t99) * qJD(1)); t70 * t13 - t39 * t19 - t35 * t20 + (-t73 * t41 + (-m(4) * qJD(3) - t48 * t70) * t71) * qJD(1) + (t107 * t100 + t78) * t72 - m(3) * (t68 * t102 + t63) - m(4) * (t58 * t100 + t63) - t76 * t134 + (-t5 * t35 - t6 * t39 + t116 + (t100 * t11 + t79) * t72) * m(6) + (-(t50 * t73 + (t14 * t72 + t15 * t70) * t71) * qJD(1) + t23 * t72 + t116) * m(5); -t72 * t13 + t78 * t70 + m(5) * (t23 * t70 - t115) + m(6) * (t70 * t79 - t115) + (t53 + (qJD(2) + t43) * m(4) + t77) * t100; m(6) * (t1 * t74 + t2 * t75) + (-t31 * t75 + t32 * t74) * mrSges(6,3) + (m(5) * qJD(2) + t77) * t98 + t129 * qJD(5); -t11 * (-mrSges(6,1) * t40 + mrSges(6,2) * t36) + t40 * (Ifges(6,1) * t36 + t117) / 0.2e1 + t9 * t120 - t55 * (Ifges(6,5) * t36 + Ifges(6,6) * t40) / 0.2e1 - t5 * t19 + t6 * t20 + (t36 * t5 - t40 * t6) * mrSges(6,3) + t106 - (Ifges(6,2) * t40 + t10 + t33) * t36 / 0.2e1 + t122;];
tauc = t16(:);
