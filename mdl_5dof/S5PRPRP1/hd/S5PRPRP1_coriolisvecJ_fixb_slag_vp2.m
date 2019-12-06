% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRPRP1
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
% Datum: 2019-12-05 15:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRPRP1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP1_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP1_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP1_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP1_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:28:19
% EndTime: 2019-12-05 15:28:23
% DurationCPUTime: 1.10s
% Computational Cost: add. (841->170), mult. (2270->221), div. (0->0), fcn. (1488->4), ass. (0->80)
t99 = Ifges(5,1) + Ifges(6,1);
t98 = Ifges(6,4) + Ifges(5,5);
t55 = sin(pkin(8));
t56 = cos(pkin(8));
t65 = (t55 ^ 2 + t56 ^ 2) * qJD(2);
t101 = mrSges(4,3) * t65;
t100 = mrSges(6,1) + mrSges(5,1);
t78 = mrSges(6,2) + mrSges(5,3);
t77 = -Ifges(5,4) + Ifges(6,5);
t97 = Ifges(6,6) - Ifges(5,6);
t57 = sin(qJ(4));
t86 = cos(qJ(4));
t67 = t86 * t56;
t64 = qJD(2) * t67;
t71 = qJD(2) * t55;
t35 = t57 * t71 - t64;
t34 = Ifges(5,4) * t35;
t41 = t86 * t55 + t57 * t56;
t36 = t41 * qJD(2);
t84 = Ifges(6,5) * t35;
t96 = t98 * qJD(4) + t99 * t36 - t34 + t84;
t76 = pkin(6) + qJ(3);
t45 = t76 * t55;
t52 = t56 * qJD(1);
t31 = -qJD(2) * t45 + t52;
t68 = t86 * t31;
t70 = qJ(3) * qJD(2);
t43 = t55 * qJD(1) + t56 * t70;
t32 = pkin(6) * qJD(2) * t56 + t43;
t81 = t57 * t32;
t11 = t68 - t81;
t95 = -t11 + qJD(5);
t94 = -t35 / 0.2e1;
t93 = t35 / 0.2e1;
t91 = t36 / 0.2e1;
t12 = t57 * t31 + t86 * t32;
t59 = t41 * qJD(3);
t3 = qJD(2) * t59 + t12 * qJD(4);
t46 = t76 * t56;
t61 = -t86 * t45 - t57 * t46;
t89 = t61 * t3;
t60 = -t57 * t55 + t67;
t88 = t3 * t60;
t87 = t3 * t41;
t85 = Ifges(5,4) * t36;
t83 = t35 * mrSges(5,3);
t82 = t36 * mrSges(5,3);
t80 = -qJD(4) / 0.2e1;
t79 = qJD(4) / 0.2e1;
t75 = qJD(3) * t64 + qJD(4) * t68;
t26 = -t35 * mrSges(6,2) + qJD(4) * mrSges(6,3);
t74 = -qJD(4) * mrSges(5,2) + t26 - t83;
t73 = -t36 * mrSges(6,2) + t100 * qJD(4) - t82;
t69 = -pkin(3) * t56 - pkin(2);
t66 = qJD(3) * t71;
t62 = -(-t55 * t70 + t52) * t55 + t43 * t56;
t21 = -t57 * t45 + t86 * t46;
t44 = t69 * qJD(2) + qJD(3);
t38 = t41 * qJD(4);
t37 = t60 * qJD(4);
t33 = Ifges(6,5) * t36;
t30 = qJD(2) * t38;
t29 = qJD(2) * t37;
t28 = t30 * mrSges(6,1);
t27 = t29 * mrSges(5,2);
t19 = -pkin(4) * t60 - qJ(5) * t41 + t69;
t18 = mrSges(6,1) * t35 - mrSges(6,3) * t36;
t17 = pkin(4) * t36 + qJ(5) * t35;
t14 = -Ifges(5,2) * t35 + Ifges(5,6) * qJD(4) + t85;
t13 = Ifges(6,6) * qJD(4) + Ifges(6,3) * t35 + t33;
t10 = t21 * qJD(4) + t59;
t9 = t60 * qJD(3) + t61 * qJD(4);
t8 = qJD(4) * qJ(5) + t12;
t7 = pkin(4) * t35 - qJ(5) * t36 + t44;
t6 = -qJD(4) * pkin(4) + t95;
t5 = pkin(4) * t38 - qJ(5) * t37 - qJD(5) * t41;
t4 = pkin(4) * t30 - qJ(5) * t29 - qJD(5) * t36;
t2 = (-qJD(4) * t32 - t66) * t57 + t75;
t1 = -t57 * t66 + (qJD(5) - t81) * qJD(4) + t75;
t15 = [-t73 * t38 + t74 * t37 + m(5) * (-t11 * t38 + t12 * t37 + t2 * t41 - t88) + m(6) * (t1 * t41 + t37 * t8 + t38 * t6 - t88) + t78 * (-t29 * t60 - t41 * t30); t19 * t28 + t4 * (-mrSges(6,1) * t60 - t41 * mrSges(6,3)) + t5 * t18 + (m(4) * (qJ(3) * t65 + t62) + 0.2e1 * t101) * qJD(3) + t69 * t27 - t73 * t10 + t74 * t9 + (-mrSges(6,3) * t19 + t99 * t41 - t60 * t77 - t61 * t78) * t29 + (t69 * mrSges(5,1) + t77 * t41 - (Ifges(5,2) + Ifges(6,3)) * t60 - t78 * t21) * t30 + (t1 * t60 + t87) * mrSges(6,2) + (t2 * t60 + t87) * mrSges(5,3) + m(6) * (t1 * t21 + t10 * t6 + t19 * t4 + t5 * t7 + t8 * t9 - t89) + m(5) * (-t10 * t11 + t12 * t9 + t2 * t21 - t89) + (Ifges(6,3) * t93 - Ifges(5,2) * t94 + t44 * mrSges(5,1) + t7 * mrSges(6,1) + t13 / 0.2e1 - t14 / 0.2e1 - t8 * mrSges(6,2) - t12 * mrSges(5,3) + t77 * t91 + t97 * t79) * t38 + (t96 / 0.2e1 + t44 * mrSges(5,2) + t6 * mrSges(6,2) - t11 * mrSges(5,3) - t7 * mrSges(6,3) + Ifges(5,4) * t94 + Ifges(6,5) * t93 + t98 * t79 + t99 * t91) * t37; t30 * mrSges(5,1) - t29 * mrSges(6,3) + t27 + t28 + t73 * t36 + t74 * t35 - m(5) * (-t11 * t36 - t12 * t35) + (-m(4) * t62 - t101) * qJD(2) + (t8 * t35 - t6 * t36 + t4) * m(6); -t44 * (mrSges(5,1) * t36 - mrSges(5,2) * t35) + (Ifges(6,3) * t36 - t84) * t94 - t7 * (t36 * mrSges(6,1) + t35 * mrSges(6,3)) + t14 * t91 - t17 * t18 + qJD(5) * t26 + t1 * mrSges(6,3) - t2 * mrSges(5,2) + t97 * t30 - t100 * t3 + t98 * t29 + (t73 + t82) * t12 + (-t74 - t83) * t11 + (-pkin(4) * t29 - qJ(5) * t30 + t35 * t6 + t36 * t8) * mrSges(6,2) + (-t98 * t35 + t97 * t36) * t80 + (-t3 * pkin(4) + t1 * qJ(5) - t6 * t12 - t7 * t17 + t95 * t8) * m(6) + (-Ifges(5,2) * t36 - t34 + t96) * t93 - (-t99 * t35 + t13 + t33 - t85) * t36 / 0.2e1; t29 * mrSges(6,2) - qJD(4) * t26 + t36 * t18 + 0.2e1 * (t3 / 0.2e1 + t8 * t80 + t7 * t91) * m(6);];
tauc = t15(:);
