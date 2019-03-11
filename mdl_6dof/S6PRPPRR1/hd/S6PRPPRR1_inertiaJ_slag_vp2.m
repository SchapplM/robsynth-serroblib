% Calculate joint inertia matrix for
% S6PRPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
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
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPPRR1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_inertiaJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR1_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPPRR1_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPPRR1_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:14:11
% EndTime: 2019-03-08 19:14:12
% DurationCPUTime: 0.48s
% Computational Cost: add. (804->168), mult. (1785->251), div. (0->0), fcn. (1921->12), ass. (0->74)
t56 = cos(pkin(12));
t53 = sin(pkin(12));
t60 = sin(qJ(5));
t77 = t60 * t53;
t84 = cos(qJ(5));
t31 = -t84 * t56 + t77;
t71 = t84 * t53;
t33 = t60 * t56 + t71;
t18 = t31 * mrSges(6,1) + t33 * mrSges(6,2);
t35 = -t56 * mrSges(5,1) + t53 * mrSges(5,2);
t96 = t18 + t35;
t59 = sin(qJ(6));
t80 = t33 * t59;
t16 = -t31 * mrSges(7,2) - mrSges(7,3) * t80;
t62 = cos(qJ(6));
t79 = t33 * t62;
t17 = t31 * mrSges(7,1) - mrSges(7,3) * t79;
t95 = t62 * t16 - t59 * t17;
t36 = -t62 * mrSges(7,1) + t59 * mrSges(7,2);
t94 = -m(7) * pkin(5) - mrSges(6,1) + t36;
t54 = sin(pkin(11));
t55 = sin(pkin(6));
t57 = cos(pkin(11));
t61 = sin(qJ(2));
t63 = cos(qJ(2));
t24 = (t54 * t63 + t57 * t61) * t55;
t58 = cos(pkin(6));
t19 = -t24 * t53 + t58 * t56;
t20 = t24 * t56 + t58 * t53;
t5 = -t84 * t19 + t60 * t20;
t93 = t5 ^ 2;
t41 = t54 * pkin(2) + qJ(4);
t85 = pkin(8) + t41;
t27 = t85 * t56;
t11 = t60 * t27 + t85 * t71;
t92 = t11 ^ 2;
t22 = (t54 * t61 - t57 * t63) * t55;
t21 = t22 ^ 2;
t91 = t31 ^ 2;
t49 = t56 ^ 2;
t90 = 0.2e1 * t11;
t89 = t62 / 0.2e1;
t88 = pkin(5) * t31;
t87 = t11 * t5;
t86 = t31 * t5;
t83 = Ifges(7,4) * t59;
t82 = Ifges(7,4) * t62;
t81 = t31 * t11;
t75 = Ifges(7,5) * t79 + Ifges(7,3) * t31;
t74 = Ifges(7,5) * t59 + Ifges(7,6) * t62;
t73 = t53 ^ 2 + t49;
t72 = t59 ^ 2 + t62 ^ 2;
t43 = -t57 * pkin(2) - pkin(3);
t70 = t72 * t33;
t7 = t60 * t19 + t84 * t20;
t1 = t22 * t62 - t59 * t7;
t2 = t22 * t59 + t62 * t7;
t69 = -t1 * t59 + t2 * t62;
t34 = -t56 * pkin(4) + t43;
t10 = -t33 * pkin(9) + t34 + t88;
t13 = t84 * t27 - t85 * t77;
t3 = t62 * t10 - t59 * t13;
t4 = t59 * t10 + t62 * t13;
t68 = -t3 * t59 + t4 * t62;
t67 = mrSges(7,1) * t59 + mrSges(7,2) * t62;
t66 = -t19 * t53 + t20 * t56;
t50 = t58 ^ 2;
t38 = Ifges(7,1) * t59 + t82;
t37 = Ifges(7,2) * t62 + t83;
t30 = t33 ^ 2;
t14 = t67 * t33;
t9 = Ifges(7,5) * t31 + (Ifges(7,1) * t62 - t83) * t33;
t8 = Ifges(7,6) * t31 + (-Ifges(7,2) * t59 + t82) * t33;
t6 = [m(2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t93) + m(6) * (t7 ^ 2 + t21 + t93) + m(5) * (t19 ^ 2 + t20 ^ 2 + t21) + m(4) * (t24 ^ 2 + t21 + t50) + m(3) * (t50 + (t61 ^ 2 + t63 ^ 2) * t55 ^ 2); -t24 * mrSges(4,2) + t1 * t17 + t5 * t14 + t2 * t16 + (t63 * mrSges(3,1) - t61 * mrSges(3,2)) * t55 + (-t7 * t31 + t5 * t33) * mrSges(6,3) + t66 * mrSges(5,3) + (-mrSges(4,1) + t96) * t22 + m(7) * (t3 * t1 + t4 * t2 + t87) + m(6) * (t13 * t7 + t34 * t22 + t87) + m(5) * (t43 * t22 + t66 * t41) + m(4) * (-t22 * t57 + t24 * t54) * pkin(2); Ifges(5,2) * t49 + t14 * t90 + 0.2e1 * t4 * t16 + 0.2e1 * t3 * t17 + 0.2e1 * t34 * t18 + 0.2e1 * t43 * t35 + Ifges(3,3) + Ifges(4,3) + (Ifges(5,1) * t53 + 0.2e1 * Ifges(5,4) * t56) * t53 + (-0.2e1 * t13 * mrSges(6,3) + Ifges(6,2) * t31 + t75) * t31 + (mrSges(6,3) * t90 + Ifges(6,1) * t33 - t59 * t8 + t62 * t9 + (-Ifges(7,6) * t59 - (2 * Ifges(6,4))) * t31) * t33 + m(7) * (t3 ^ 2 + t4 ^ 2 + t92) + m(6) * (t13 ^ 2 + t34 ^ 2 + t92) + m(5) * (t73 * t41 ^ 2 + t43 ^ 2) + m(4) * (t54 ^ 2 + t57 ^ 2) * pkin(2) ^ 2 + 0.2e1 * (t57 * mrSges(4,1) - t54 * mrSges(4,2)) * pkin(2) + 0.2e1 * t73 * t41 * mrSges(5,3); m(4) * t58 + m(7) * (t69 * t33 + t86) + m(6) * (t33 * t7 + t86) + m(5) * (t56 * t19 + t53 * t20); t31 * t14 + t95 * t33 + m(7) * (t68 * t33 + t81) + m(6) * (t33 * t13 + t81); m(4) + m(5) * t73 + m(6) * (t30 + t91) + m(7) * (t72 * t30 + t91); m(7) * (t62 * t1 + t59 * t2) + 0.2e1 * (m(6) / 0.2e1 + m(5) / 0.2e1) * t22; t59 * t16 + t62 * t17 + m(7) * (t62 * t3 + t59 * t4) + m(6) * t34 + m(5) * t43 + t96; 0; m(7) * t72 + m(5) + m(6); -t7 * mrSges(6,2) + (m(7) * pkin(9) + mrSges(7,3)) * t69 + t94 * t5; t59 * t9 / 0.2e1 + t8 * t89 - pkin(5) * t14 - t13 * mrSges(6,2) + t68 * mrSges(7,3) + (-t59 * t37 / 0.2e1 + t38 * t89 + Ifges(6,5)) * t33 + (t74 / 0.2e1 - Ifges(6,6)) * t31 + t94 * t11 + (m(7) * t68 + t95) * pkin(9); t31 * t36 + m(7) * (pkin(9) * t70 - t88) + mrSges(7,3) * t70 - t18; 0; Ifges(6,3) - 0.2e1 * pkin(5) * t36 + m(7) * (t72 * pkin(9) ^ 2 + pkin(5) ^ 2) + t59 * t38 + t62 * t37 + 0.2e1 * t72 * pkin(9) * mrSges(7,3); t1 * mrSges(7,1) - t2 * mrSges(7,2); t3 * mrSges(7,1) - t4 * mrSges(7,2) - Ifges(7,6) * t80 + t75; -t14; -t36; -t67 * pkin(9) + t74; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t6(1) t6(2) t6(4) t6(7) t6(11) t6(16); t6(2) t6(3) t6(5) t6(8) t6(12) t6(17); t6(4) t6(5) t6(6) t6(9) t6(13) t6(18); t6(7) t6(8) t6(9) t6(10) t6(14) t6(19); t6(11) t6(12) t6(13) t6(14) t6(15) t6(20); t6(16) t6(17) t6(18) t6(19) t6(20) t6(21);];
Mq  = res;
