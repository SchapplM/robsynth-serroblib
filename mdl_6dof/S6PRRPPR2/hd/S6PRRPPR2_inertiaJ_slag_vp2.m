% Calculate joint inertia matrix for
% S6PRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
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
% Datum: 2019-03-08 21:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPPR2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR2_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR2_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPPR2_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:04:10
% EndTime: 2019-03-08 21:04:11
% DurationCPUTime: 0.64s
% Computational Cost: add. (826->194), mult. (1714->265), div. (0->0), fcn. (1775->10), ass. (0->81)
t99 = m(5) + m(6);
t59 = cos(pkin(6));
t61 = sin(qJ(3));
t64 = cos(qJ(3));
t57 = sin(pkin(6));
t62 = sin(qJ(2));
t86 = t57 * t62;
t28 = t59 * t64 - t61 * t86;
t29 = t59 * t61 + t64 * t86;
t56 = sin(pkin(11));
t58 = cos(pkin(11));
t13 = t56 * t28 + t58 * t29;
t10 = t13 ^ 2;
t95 = t56 * pkin(3);
t44 = qJ(5) + t95;
t97 = t44 ^ 2;
t55 = t64 ^ 2;
t96 = m(5) * pkin(3);
t94 = t58 * pkin(3);
t60 = sin(qJ(6));
t63 = cos(qJ(6));
t79 = t60 ^ 2 + t63 ^ 2;
t36 = m(7) * t79;
t93 = m(6) + t36;
t92 = Ifges(7,4) * t60;
t91 = Ifges(7,4) * t63;
t90 = Ifges(7,6) * t60;
t33 = t56 * t61 - t58 * t64;
t89 = t33 * t60;
t88 = t33 * t63;
t87 = t44 * t13;
t65 = cos(qJ(2));
t85 = t57 * t65;
t84 = mrSges(6,1) + mrSges(5,3);
t83 = -mrSges(6,2) + mrSges(5,1);
t82 = -qJ(4) - pkin(8);
t38 = t60 * mrSges(7,1) + t63 * mrSges(7,2);
t81 = t38 + mrSges(6,3);
t80 = t61 ^ 2 + t55;
t34 = t56 * t64 + t58 * t61;
t77 = Ifges(7,5) * t89 + Ifges(7,6) * t88 + Ifges(7,3) * t34;
t39 = t82 * t64;
t75 = t82 * t61;
t21 = -t56 * t39 - t58 * t75;
t23 = -t58 * t39 + t56 * t75;
t76 = t21 ^ 2 + t23 ^ 2;
t49 = -t64 * pkin(3) - pkin(2);
t48 = -pkin(4) - t94;
t74 = t79 * mrSges(7,3);
t68 = -t34 * qJ(5) + t49;
t7 = (pkin(4) + pkin(9)) * t33 + t68;
t8 = t34 * pkin(5) + t21;
t1 = -t60 * t7 + t63 * t8;
t2 = t60 * t8 + t63 * t7;
t73 = t63 * t1 + t60 * t2;
t11 = -t58 * t28 + t56 * t29;
t3 = t63 * t11 + t60 * t85;
t4 = t60 * t11 - t63 * t85;
t72 = t63 * t3 + t60 * t4;
t71 = t63 * mrSges(7,1) - t60 * mrSges(7,2);
t70 = t21 * t11 + t23 * t13;
t69 = -t28 * t61 + t29 * t64;
t51 = t57 ^ 2;
t50 = Ifges(7,5) * t63;
t46 = t51 * t65 ^ 2;
t43 = -pkin(9) + t48;
t41 = Ifges(7,1) * t63 - t92;
t40 = -Ifges(7,2) * t60 + t91;
t37 = -t64 * mrSges(4,1) + t61 * mrSges(4,2);
t31 = t34 * mrSges(6,3);
t30 = t34 * mrSges(5,2);
t20 = -t33 * mrSges(6,2) - t31;
t19 = t33 * mrSges(5,1) + t30;
t18 = -t34 * mrSges(7,2) + mrSges(7,3) * t88;
t17 = t34 * mrSges(7,1) - mrSges(7,3) * t89;
t16 = t33 * pkin(4) + t68;
t15 = t71 * t33;
t9 = -t33 * pkin(5) + t23;
t6 = Ifges(7,5) * t34 + (Ifges(7,1) * t60 + t91) * t33;
t5 = Ifges(7,6) * t34 + (Ifges(7,2) * t63 + t92) * t33;
t12 = [m(2) + m(7) * (t3 ^ 2 + t4 ^ 2 + t10) + m(4) * (t28 ^ 2 + t29 ^ 2 + t46) + m(3) * (t51 * t62 ^ 2 + t59 ^ 2 + t46) + t99 * (t11 ^ 2 + t10 + t46); t3 * t17 + t4 * t18 + t84 * t34 * t11 + t69 * mrSges(4,3) + (-t84 * t33 - t15) * t13 + (-t62 * mrSges(3,2) + (mrSges(3,1) - t19 - t20 - t37) * t65) * t57 + m(6) * (-t16 * t85 + t70) + m(7) * (t1 * t3 + t9 * t13 + t2 * t4) + m(5) * (-t49 * t85 + t70) + m(4) * (pkin(2) * t85 + t69 * pkin(8)); Ifges(4,2) * t55 - 0.2e1 * pkin(2) * t37 + 0.2e1 * t1 * t17 - 0.2e1 * t9 * t15 + 0.2e1 * t16 * t20 + 0.2e1 * t2 * t18 + 0.2e1 * t49 * t19 + Ifges(3,3) + (Ifges(4,1) * t61 + 0.2e1 * Ifges(4,4) * t64) * t61 + 0.2e1 * t80 * pkin(8) * mrSges(4,3) + m(4) * (t80 * pkin(8) ^ 2 + pkin(2) ^ 2) + m(5) * (t49 ^ 2 + t76) + m(6) * (t16 ^ 2 + t76) + m(7) * (t1 ^ 2 + t2 ^ 2 + t9 ^ 2) + ((Ifges(5,1) + Ifges(6,2)) * t34 + 0.2e1 * t84 * t21 + t77) * t34 + (t63 * t5 + t60 * t6 + (Ifges(5,2) + Ifges(6,3)) * t33 + 0.2e1 * (-Ifges(5,4) - Ifges(6,6)) * t34 - 0.2e1 * t84 * t23) * t33; t28 * mrSges(4,1) - t29 * mrSges(4,2) - t83 * t11 - t72 * mrSges(7,3) + (-mrSges(5,2) + t81) * t13 + m(6) * (t48 * t11 + t87) + m(7) * (t72 * t43 + t87) + (-t11 * t58 + t13 * t56) * t96; Ifges(4,5) * t61 + Ifges(4,6) * t64 - t44 * t15 + t9 * t38 + (-mrSges(5,2) + mrSges(6,3)) * t23 - t83 * t21 + (-t61 * mrSges(4,1) - t64 * mrSges(4,2)) * pkin(8) + (t43 * t17 - t1 * mrSges(7,3) + t6 / 0.2e1) * t63 + (t43 * t18 - t2 * mrSges(7,3) - t5 / 0.2e1) * t60 + m(7) * (t73 * t43 + t44 * t9) + m(6) * (t48 * t21 + t44 * t23) + (-t21 * t58 + t23 * t56) * t96 + (t48 * mrSges(6,1) - t90 / 0.2e1 + t50 / 0.2e1 + Ifges(5,5) - Ifges(6,4) - mrSges(5,3) * t94) * t34 + (t60 * t41 / 0.2e1 + t63 * t40 / 0.2e1 - t44 * mrSges(6,1) + Ifges(6,5) - Ifges(5,6) - mrSges(5,3) * t95) * t33; 0.2e1 * t48 * mrSges(6,2) - t60 * t40 + t63 * t41 + Ifges(6,1) + Ifges(4,3) + Ifges(5,3) + m(7) * (t79 * t43 ^ 2 + t97) + m(6) * (t48 ^ 2 + t97) + m(5) * (t56 ^ 2 + t58 ^ 2) * pkin(3) ^ 2 + 0.2e1 * t81 * t44 + 0.2e1 * (t58 * mrSges(5,1) - t56 * mrSges(5,2)) * pkin(3) - 0.2e1 * t43 * t74; m(7) * (-t60 * t3 + t63 * t4) - t99 * t85; -t60 * t17 + t63 * t18 + t30 - t31 + t83 * t33 + m(7) * (-t60 * t1 + t63 * t2) + m(6) * t16 + m(5) * t49; 0; m(5) + t93; m(6) * t11 + m(7) * t72; m(6) * t21 + m(7) * t73 + t34 * mrSges(6,1) + t63 * t17 + t60 * t18; m(6) * t48 + t43 * t36 + mrSges(6,2) - t74; 0; t93; t3 * mrSges(7,1) - t4 * mrSges(7,2); t1 * mrSges(7,1) - t2 * mrSges(7,2) + t77; t71 * t43 + t50 - t90; -t38; t71; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t12(1) t12(2) t12(4) t12(7) t12(11) t12(16); t12(2) t12(3) t12(5) t12(8) t12(12) t12(17); t12(4) t12(5) t12(6) t12(9) t12(13) t12(18); t12(7) t12(8) t12(9) t12(10) t12(14) t12(19); t12(11) t12(12) t12(13) t12(14) t12(15) t12(20); t12(16) t12(17) t12(18) t12(19) t12(20) t12(21);];
Mq  = res;
