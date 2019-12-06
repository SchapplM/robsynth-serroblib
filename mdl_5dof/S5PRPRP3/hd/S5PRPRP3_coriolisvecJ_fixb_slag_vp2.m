% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRPRP3
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
% Datum: 2019-12-05 15:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRPRP3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP3_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP3_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP3_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP3_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:32:40
% EndTime: 2019-12-05 15:32:44
% DurationCPUTime: 1.14s
% Computational Cost: add. (734->179), mult. (1896->240), div. (0->0), fcn. (1106->6), ass. (0->82)
t102 = Ifges(5,4) + Ifges(6,4);
t60 = cos(qJ(2));
t82 = qJD(1) * t60;
t43 = qJD(2) * pkin(2) + t82;
t55 = sin(pkin(8));
t58 = sin(qJ(2));
t83 = qJD(1) * t58;
t44 = t55 * t83;
t56 = cos(pkin(8));
t17 = t43 * t56 - t44;
t59 = cos(qJ(4));
t71 = -pkin(4) * t59 - pkin(3);
t10 = qJD(2) * t71 + qJD(5) - t17;
t12 = -qJD(2) * pkin(3) - t17;
t57 = sin(qJ(4));
t35 = (-t59 * mrSges(6,1) + t57 * mrSges(6,2)) * qJD(2);
t18 = t55 * t43 + t56 * t83;
t13 = qJD(2) * pkin(6) + t18;
t68 = qJ(5) * qJD(2) + t13;
t79 = qJD(3) * t57;
t7 = t59 * t68 + t79;
t9 = t13 * t59 + t79;
t92 = m(6) * t10;
t96 = qJD(4) / 0.2e1;
t100 = -(t35 + t92) * pkin(4) - t12 * mrSges(5,1) - t10 * mrSges(6,1) + t9 * mrSges(5,3) + t7 * mrSges(6,3) + ((Ifges(5,2) + Ifges(6,2)) * t59 + t102 * t57) * qJD(2) / 0.2e1 + (Ifges(6,6) + Ifges(5,6)) * t96;
t99 = t59 * mrSges(5,1) - t57 * mrSges(5,2) + mrSges(4,1);
t80 = qJD(2) * t59;
t94 = -t102 * t80 / 0.2e1;
t47 = pkin(2) * t55 + pkin(6);
t93 = m(5) * t47;
t91 = pkin(2) * t56;
t34 = t55 * t60 + t56 * t58;
t22 = t34 * qJD(2);
t19 = qJD(1) * t22;
t33 = t55 * t58 - t56 * t60;
t90 = t19 * t33;
t24 = t33 * qJD(2);
t89 = t24 * t57;
t88 = t24 * t59;
t81 = qJD(2) * t57;
t39 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t81;
t40 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t81;
t87 = -t39 - t40;
t41 = -qJD(4) * mrSges(6,2) + mrSges(6,3) * t80;
t42 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t80;
t86 = t41 + t42;
t75 = qJD(2) * qJD(4);
t29 = (mrSges(6,1) * t57 + mrSges(6,2) * t59) * t75;
t84 = qJ(5) + t47;
t78 = qJD(4) * t57;
t77 = qJD(4) * t59;
t76 = qJD(5) * t59;
t52 = t59 * qJD(3);
t74 = -0.3e1 / 0.2e1 * Ifges(6,4) - 0.3e1 / 0.2e1 * Ifges(5,4);
t73 = Ifges(6,5) / 0.2e1 + Ifges(5,5) / 0.2e1;
t8 = -t13 * t57 + t52;
t72 = -m(5) * t8 - t40;
t70 = mrSges(5,3) + t93;
t69 = qJD(4) * t84;
t67 = t99 * qJD(2) - t35;
t66 = (-Ifges(6,6) / 0.2e1 - Ifges(5,6) / 0.2e1) * qJD(4);
t20 = qJD(1) * t24;
t3 = qJD(4) * t52 - t13 * t78 - t59 * t20;
t6 = -t57 * t68 + t52;
t5 = qJD(4) * pkin(4) + t6;
t65 = m(6) * (-t5 * t57 + t59 * t7);
t21 = t34 * qJD(1);
t64 = t57 * t87 + t59 * t86;
t63 = t12 * mrSges(5,2) + t10 * mrSges(6,2) - t8 * mrSges(5,3) - t5 * mrSges(6,3) - t94 + (Ifges(5,1) + Ifges(6,1)) * t81 / 0.2e1 + (Ifges(5,5) + Ifges(6,5)) * t96;
t48 = -pkin(3) - t91;
t38 = t71 - t91;
t32 = t84 * t59;
t31 = t84 * t57;
t30 = (mrSges(5,1) * t57 + mrSges(5,2) * t59) * t75;
t23 = t56 * t82 - t44;
t15 = -qJD(5) * t57 - t59 * t69;
t14 = -t57 * t69 + t76;
t11 = (pkin(4) * t78 + t21) * qJD(2);
t4 = -qJD(4) * t9 + t20 * t57;
t2 = (-qJD(2) * qJD(5) + t20) * t57 - t7 * qJD(4);
t1 = (-qJ(5) * t78 + t76) * qJD(2) + t3;
t16 = [(-mrSges(3,1) * t58 - mrSges(3,2) * t60) * qJD(2) ^ 2 + (t29 + t30) * t33 - t67 * t22 - (-qJD(2) * mrSges(4,2) + t64) * t24 + m(4) * (-t17 * t22 - t18 * t24 + t90) + m(5) * (t12 * t22 + t8 * t89 - t9 * t88 + t90) + m(6) * (t10 * t22 + t11 * t33 + t5 * t89 - t7 * t88) + ((-t57 * t86 + t59 * t87) * qJD(4) - m(4) * t20 + m(5) * (t3 * t59 - t4 * t57 - t8 * t77 - t9 * t78) + m(6) * (t1 * t59 - t2 * t57 - t5 * t77 - t7 * t78)) * t34; t14 * t41 + t15 * t39 + t38 * t29 + t48 * t30 + (t23 * qJD(2) + t20) * mrSges(4,2) + m(6) * (t1 * t32 + t11 * t38 + t14 * t7 + t15 * t5 - t2 * t31) + (t11 * mrSges(6,2) - t2 * mrSges(6,3) - t70 * t4 + (m(6) * t5 + t39 - t72) * t23 + (-t47 * t42 + t66 - t9 * t93 + (-t32 * mrSges(6,3) + t57 * t74) * qJD(2) - t100) * qJD(4)) * t57 + (-t11 * mrSges(6,1) + t1 * mrSges(6,3) + t70 * t3 + (-m(5) * t9 - m(6) * t7 - t86) * t23 + (t73 * qJD(4) + t72 * t47 + (t31 * mrSges(6,3) - t74 * t59 + (0.3e1 / 0.2e1 * Ifges(6,1) + 0.3e1 / 0.2e1 * Ifges(5,1) - 0.3e1 / 0.2e1 * Ifges(5,2) - 0.3e1 / 0.2e1 * Ifges(6,2)) * t57) * qJD(2) + t63) * qJD(4)) * t59 + (m(5) * t48 - t99) * t19 + (-t18 * t23 + (-t19 * t56 - t20 * t55) * pkin(2)) * m(4) + (m(4) * t17 - m(5) * t12 + t67 - t92) * t21; m(5) * (t3 * t57 + t4 * t59) + m(6) * (t1 * t57 + t2 * t59) + (m(5) * (-t57 * t8 + t59 * t9) + t65 + t64 + (mrSges(5,3) + mrSges(6,3)) * qJD(2) * (-t57 ^ 2 - t59 ^ 2)) * qJD(4); t4 * mrSges(5,1) - t3 * mrSges(5,2) - t1 * mrSges(6,2) + t9 * t40 - t6 * t41 - t8 * t42 + (m(6) * pkin(4) + mrSges(6,1)) * t2 + (t39 - m(6) * (-t5 + t6)) * t7 + (((-pkin(4) * mrSges(6,3) + t73) * qJD(4) - t63 + t94) * t59 + ((Ifges(6,4) / 0.2e1 + Ifges(5,4) / 0.2e1) * t81 + t66 + (Ifges(6,2) / 0.2e1 + Ifges(5,2) / 0.2e1 - Ifges(6,1) / 0.2e1 - Ifges(5,1) / 0.2e1) * t80 + t100) * t57) * qJD(2); m(6) * t11 + (t57 * t39 - t59 * t41 - t65) * qJD(2) + t29;];
tauc = t16(:);
