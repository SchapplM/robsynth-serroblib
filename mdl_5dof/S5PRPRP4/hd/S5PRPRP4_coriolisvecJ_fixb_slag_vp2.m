% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRPRP4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP4_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP4_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP4_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP4_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:35:04
% EndTime: 2019-12-05 15:35:08
% DurationCPUTime: 1.01s
% Computational Cost: add. (730->175), mult. (1868->253), div. (0->0), fcn. (1087->6), ass. (0->83)
t52 = cos(qJ(4));
t69 = qJD(2) * t52;
t98 = -t69 / 0.2e1;
t50 = sin(qJ(4));
t97 = t52 * mrSges(5,1) - t50 * mrSges(5,2) + mrSges(4,1);
t96 = mrSges(5,1) + mrSges(6,1);
t90 = mrSges(6,2) + mrSges(5,3);
t95 = -Ifges(5,1) / 0.2e1;
t94 = Ifges(6,6) / 0.2e1;
t93 = Ifges(5,4) * t98;
t92 = -qJD(2) / 0.2e1;
t91 = -qJD(4) / 0.2e1;
t53 = cos(qJ(2));
t71 = qJD(1) * t53;
t39 = qJD(2) * pkin(2) + t71;
t48 = sin(pkin(8));
t49 = cos(pkin(8));
t51 = sin(qJ(2));
t72 = qJD(1) * t51;
t14 = t48 * t39 + t49 * t72;
t11 = qJD(2) * pkin(6) + t14;
t81 = t11 * t50;
t8 = qJD(3) * t52 - t81;
t73 = -t8 + qJD(5);
t88 = 2 * m(5);
t87 = 2 * m(6);
t40 = t48 * t72;
t20 = t49 * t71 - t40;
t86 = -t20 / 0.2e1;
t41 = pkin(2) * t48 + pkin(6);
t85 = t41 / 0.2e1;
t84 = pkin(2) * t49;
t29 = t48 * t51 - t49 * t53;
t21 = t29 * qJD(2);
t16 = qJD(1) * t21;
t9 = qJD(3) * t50 + t11 * t52;
t3 = qJD(4) * t9 - t16 * t50;
t83 = t3 * t50;
t82 = t3 * t52;
t30 = t48 * t53 + t49 * t51;
t19 = t30 * qJD(2);
t15 = qJD(1) * t19;
t80 = t15 * t29;
t79 = t21 * t50;
t78 = t21 * t52;
t67 = qJD(4) * t52;
t77 = qJD(3) * t67 - t52 * t16;
t70 = qJD(2) * t50;
t76 = t96 * qJD(4) - t90 * t70;
t37 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t69;
t38 = mrSges(6,2) * t69 + qJD(4) * mrSges(6,3);
t75 = t37 + t38;
t68 = qJD(4) * t50;
t66 = qJD(2) * qJD(4);
t65 = t3 * t85;
t64 = Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1;
t63 = -0.3e1 / 0.2e1 * Ifges(6,5) + 0.3e1 / 0.2e1 * Ifges(5,4);
t62 = t94 - Ifges(5,6) / 0.2e1;
t61 = t70 / 0.2e1;
t13 = t39 * t49 - t40;
t31 = (-t52 * mrSges(6,1) - t50 * mrSges(6,3)) * qJD(2);
t60 = t97 * qJD(2) - t31;
t59 = pkin(4) * t50 - qJ(5) * t52;
t58 = -pkin(4) * t52 - qJ(5) * t50 - pkin(3);
t18 = t30 * qJD(1);
t57 = -t50 * t76 + t75 * t52;
t17 = qJD(4) * t59 - qJD(5) * t50;
t10 = -qJD(2) * pkin(3) - t13;
t44 = Ifges(6,5) * t70;
t6 = qJD(4) * qJ(5) + t9;
t7 = qJD(2) * t58 - t13;
t56 = t10 * mrSges(5,1) + t7 * mrSges(6,1) + qJD(4) * t94 + Ifges(6,3) * t98 + t44 / 0.2e1 + Ifges(5,6) * t91 + (Ifges(5,4) * t50 + Ifges(5,2) * t52) * t92 - t6 * mrSges(6,2) - t9 * mrSges(5,3);
t5 = -qJD(4) * pkin(4) + t73;
t55 = t7 * mrSges(6,3) + t8 * mrSges(5,3) + (Ifges(6,1) * t50 - Ifges(6,5) * t52) * t92 + t70 * t95 + t93 - t10 * mrSges(5,2) - t5 * mrSges(6,2) + (Ifges(6,4) + Ifges(5,5)) * t91;
t42 = -pkin(3) - t84;
t33 = t59 * qJD(2);
t28 = t58 - t84;
t27 = (mrSges(5,1) * t50 + mrSges(5,2) * t52) * t66;
t26 = (mrSges(6,1) * t50 - mrSges(6,3) * t52) * t66;
t4 = (t18 + t17) * qJD(2);
t2 = -t11 * t68 + t77;
t1 = (qJD(5) - t81) * qJD(4) + t77;
t12 = [(-mrSges(3,1) * t51 - mrSges(3,2) * t53) * qJD(2) ^ 2 + (t26 + t27) * t29 - t60 * t19 - (-qJD(2) * mrSges(4,2) + t57) * t21 + m(4) * (-t13 * t19 - t14 * t21 + t80) + m(5) * (t10 * t19 - t9 * t78 + t8 * t79 + t80) + m(6) * (t19 * t7 + t29 * t4 - t5 * t79 - t6 * t78) + ((-t50 * t75 - t52 * t76) * qJD(4) - m(4) * t16 + m(5) * (t2 * t52 - t8 * t67 - t9 * t68 + t83) + m(6) * (t1 * t52 + t5 * t67 - t6 * t68 + t83)) * t30; t17 * t31 + t28 * t26 + t42 * t27 + (t20 * qJD(2) + t16) * mrSges(4,2) + m(6) * (t17 * t7 + t28 * t4) + (-t4 * mrSges(6,1) + t1 * mrSges(6,2) + t2 * mrSges(5,3) - t75 * t20 + (t1 * t85 + t6 * t86) * t87 + (t2 * t85 + t86 * t9) * t88 + (t64 * qJD(4) + t63 * t69 + (-m(5) * t8 + m(6) * t5 - t76) * t41 - t55) * qJD(4)) * t52 + (-t4 * mrSges(6,3) + t90 * t3 + t76 * t20 + (t5 * t86 + t65) * t87 + (t65 + t20 * t8 / 0.2e1) * t88 + (t62 * qJD(4) + (-m(5) * t9 - m(6) * t6 - t75) * t41 + (-t63 * t50 + (-0.3e1 / 0.2e1 * Ifges(6,3) - 0.3e1 / 0.2e1 * Ifges(5,2) + 0.3e1 / 0.2e1 * Ifges(5,1) + 0.3e1 / 0.2e1 * Ifges(6,1)) * t52) * qJD(2) + t56) * qJD(4)) * t50 + (t42 * t88 / 0.2e1 - t97) * t15 + (t60 - m(6) * t7 - t10 * t88 / 0.2e1) * t18 + (t13 * t18 - t14 * t20 + (-t15 * t49 - t16 * t48) * pkin(2)) * m(4); m(5) * (t2 * t50 - t82) + m(6) * (t1 * t50 - t82) + (m(5) * (-t50 * t8 + t52 * t9) + m(6) * (t5 * t50 + t52 * t6) + t57 + t90 * qJD(2) * (-t50 ^ 2 - t52 ^ 2)) * qJD(4); -t2 * mrSges(5,2) + t1 * mrSges(6,3) - t33 * t31 - t8 * t37 + t76 * t9 + t73 * t38 - t96 * t3 + ((t93 + Ifges(6,5) * t69 / 0.2e1 + (-pkin(4) * mrSges(6,2) + t64) * qJD(4) + t55) * t52 + (-t44 / 0.2e1 + Ifges(5,4) * t61 + (-qJ(5) * mrSges(6,2) + t62) * qJD(4) + (-Ifges(6,1) / 0.2e1 + Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1 + t95) * t69 - t56) * t50) * qJD(2) + (-t3 * pkin(4) + t1 * qJ(5) - t7 * t33 - t5 * t9 + t73 * t6) * m(6); -qJD(4) * t38 + (mrSges(6,2) * t67 + t31 * t50) * qJD(2) + (t3 / 0.2e1 + t7 * t61 + t6 * t91) * t87;];
tauc = t12(:);
