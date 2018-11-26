% Calculate joint inertia matrix for
% S6RRPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 16:48
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRR1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR1_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR1_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR1_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR1_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR1_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:48:27
% EndTime: 2018-11-23 16:48:27
% DurationCPUTime: 0.73s
% Computational Cost: add. (1502->218), mult. (2713->293), div. (0->0), fcn. (2919->8), ass. (0->81)
t57 = sin(qJ(6));
t60 = cos(qJ(6));
t80 = t57 ^ 2 + t60 ^ 2;
t74 = t80 * mrSges(7,3);
t39 = -mrSges(7,1) * t60 + mrSges(7,2) * t57;
t81 = -t39 + mrSges(6,1);
t55 = sin(pkin(10));
t56 = cos(pkin(10));
t59 = sin(qJ(2));
t62 = cos(qJ(2));
t32 = t55 * t59 - t56 * t62;
t33 = t55 * t62 + t56 * t59;
t58 = sin(qJ(5));
t61 = cos(qJ(5));
t18 = -t32 * t61 + t33 * t58;
t19 = t32 * t58 + t33 * t61;
t90 = t19 * t57;
t10 = -mrSges(7,2) * t18 - mrSges(7,3) * t90;
t89 = t19 * t60;
t11 = mrSges(7,1) * t18 - mrSges(7,3) * t89;
t104 = t60 * t10 - t57 * t11;
t103 = m(7) * pkin(5) + t81;
t83 = -qJ(3) - pkin(7);
t75 = t83 * t59;
t76 = t83 * t62;
t23 = t55 * t75 - t56 * t76;
t14 = pkin(8) * t32 + t23;
t21 = -t55 * t76 - t56 * t75;
t66 = -t33 * pkin(8) + t21;
t6 = t14 * t58 - t61 * t66;
t102 = t6 ^ 2;
t101 = 0.2e1 * t6;
t47 = -pkin(2) * t62 - pkin(1);
t69 = qJ(4) * t33 - t47;
t12 = (-pkin(3) - pkin(4)) * t32 + t69;
t100 = 0.2e1 * t12;
t99 = -0.2e1 * t39;
t98 = 0.2e1 * t47;
t97 = t60 / 0.2e1;
t96 = pkin(2) * t55;
t95 = pkin(2) * t56;
t94 = t6 * t61;
t8 = t61 * t14 + t58 * t66;
t93 = t8 * mrSges(6,2);
t92 = Ifges(7,4) * t57;
t91 = Ifges(7,4) * t60;
t46 = -pkin(3) - t95;
t43 = -pkin(4) + t46;
t44 = qJ(4) + t96;
t27 = t43 * t61 - t44 * t58;
t88 = t27 * mrSges(6,1);
t28 = t43 * t58 + t44 * t61;
t87 = t28 * mrSges(6,2);
t84 = mrSges(5,2) + mrSges(4,3);
t82 = Ifges(7,5) * t89 + Ifges(7,3) * t18;
t79 = t59 ^ 2 + t62 ^ 2;
t40 = Ifges(7,5) * t57 + Ifges(7,6) * t60;
t78 = t40 / 0.2e1 - Ifges(6,6);
t77 = t21 ^ 2 + t23 ^ 2;
t26 = -pkin(9) + t28;
t73 = t80 * t26;
t72 = t80 * t58;
t3 = pkin(5) * t18 - pkin(9) * t19 + t12;
t1 = t3 * t60 - t57 * t8;
t2 = t3 * t57 + t60 * t8;
t71 = -t1 * t57 + t2 * t60;
t70 = mrSges(7,1) * t57 + mrSges(7,2) * t60;
t41 = Ifges(7,2) * t60 + t92;
t42 = Ifges(7,1) * t57 + t91;
t68 = t60 * t41 + t57 * t42 + Ifges(6,3);
t67 = t42 * t97 - t57 * t41 / 0.2e1 + Ifges(6,5);
t53 = t61 ^ 2;
t50 = t58 ^ 2;
t31 = t33 * mrSges(4,2);
t30 = t32 * mrSges(5,1);
t25 = pkin(5) - t27;
t16 = pkin(3) * t32 - t69;
t9 = t70 * t19;
t5 = Ifges(7,5) * t18 + (Ifges(7,1) * t60 - t92) * t19;
t4 = Ifges(7,6) * t18 + (-Ifges(7,2) * t57 + t91) * t19;
t7 = [-0.2e1 * pkin(1) * (-t62 * mrSges(3,1) + t59 * mrSges(3,2)) + t59 * (Ifges(3,1) * t59 + Ifges(3,4) * t62) + t62 * (Ifges(3,4) * t59 + Ifges(3,2) * t62) + t31 * t98 + 0.2e1 * t16 * t30 + t9 * t101 + 0.2e1 * t2 * t10 + 0.2e1 * t1 * t11 + Ifges(2,3) + 0.2e1 * t79 * pkin(7) * mrSges(3,3) + (mrSges(6,1) * t100 - 0.2e1 * mrSges(6,3) * t8 + Ifges(6,2) * t18 + t82) * t18 + (mrSges(6,2) * t100 + mrSges(6,3) * t101 + Ifges(6,1) * t19 - t57 * t4 + t60 * t5 + (-Ifges(7,6) * t57 - (2 * Ifges(6,4))) * t18) * t19 + (-0.2e1 * t16 * mrSges(5,3) + (Ifges(5,1) + Ifges(4,1)) * t33 + 0.2e1 * t84 * t21) * t33 + (mrSges(4,1) * t98 + (Ifges(5,3) + Ifges(4,2)) * t32 + 0.2e1 * (-Ifges(4,4) + Ifges(5,5)) * t33 - 0.2e1 * t84 * t23) * t32 + m(3) * (pkin(7) ^ 2 * t79 + pkin(1) ^ 2) + m(4) * (t47 ^ 2 + t77) + m(5) * (t16 ^ 2 + t77) + m(7) * (t1 ^ 2 + t2 ^ 2 + t102) + m(6) * (t12 ^ 2 + t8 ^ 2 + t102); t93 + Ifges(3,5) * t59 + Ifges(3,6) * t62 + t25 * t9 + t81 * t6 + (-mrSges(4,2) + mrSges(5,3)) * t23 + (-mrSges(5,1) - mrSges(4,1)) * t21 + (-mrSges(3,1) * t59 - mrSges(3,2) * t62) * pkin(7) + (-t4 / 0.2e1 + t26 * t10 - t2 * mrSges(7,3)) * t60 + (-t5 / 0.2e1 + t1 * mrSges(7,3) - t26 * t11) * t57 + (-t28 * mrSges(6,3) - t78) * t18 + (mrSges(5,2) * t46 - mrSges(4,3) * t95 + Ifges(5,4) + Ifges(4,5)) * t33 + (-mrSges(5,2) * t44 - mrSges(4,3) * t96 - Ifges(4,6) + Ifges(5,6)) * t32 + (-t27 * mrSges(6,3) - t67) * t19 + m(7) * (t25 * t6 + t26 * t71) + m(5) * (t21 * t46 + t23 * t44) + m(6) * (-t27 * t6 + t28 * t8) + m(4) * (-t21 * t56 + t23 * t55) * pkin(2); -0.2e1 * t46 * mrSges(5,1) - 0.2e1 * t88 + 0.2e1 * t87 + 0.2e1 * t44 * mrSges(5,3) + t25 * t99 + Ifges(5,2) + Ifges(3,3) + Ifges(4,3) + m(7) * (t26 ^ 2 * t80 + t25 ^ 2) + m(6) * (t27 ^ 2 + t28 ^ 2) + m(5) * (t44 ^ 2 + t46 ^ 2) + m(4) * (t55 ^ 2 + t56 ^ 2) * pkin(2) ^ 2 + t68 + 0.2e1 * (mrSges(4,1) * t56 - mrSges(4,2) * t55) * pkin(2) - 0.2e1 * t26 * t74; t32 * mrSges(4,1) - t18 * mrSges(6,1) - t19 * mrSges(6,2) - t33 * mrSges(5,3) - t57 * t10 - t60 * t11 + t30 + t31 + m(7) * (-t1 * t60 - t2 * t57) - m(6) * t12 + m(5) * t16 + m(4) * t47; 0; m(7) * t80 + m(4) + m(5) + m(6); t33 * mrSges(5,2) + (-t19 * mrSges(6,3) - t9) * t61 + (-t18 * mrSges(6,3) + t104) * t58 + m(7) * (t58 * t71 - t94) + m(6) * (t58 * t8 - t94) + m(5) * t21; -mrSges(5,1) - t81 * t61 + (mrSges(6,2) - t74) * t58 + m(7) * (-t25 * t61 + t26 * t72) + m(6) * (t27 * t61 + t28 * t58) + m(5) * t46; 0; m(5) + m(6) * (t50 + t53) + m(7) * (t50 * t80 + t53); t4 * t97 + t57 * t5 / 0.2e1 - pkin(5) * t9 - t93 + t78 * t18 + t71 * mrSges(7,3) + t67 * t19 - t103 * t6 + (m(7) * t71 + t104) * pkin(9); m(7) * (-pkin(5) * t25 + pkin(9) * t73) - t87 + t88 + (pkin(5) + t25) * t39 + (-pkin(9) * t80 + t73) * mrSges(7,3) - t68; 0; -t58 * mrSges(6,2) + (m(7) * pkin(9) + mrSges(7,3)) * t72 + t103 * t61; pkin(5) * t99 + m(7) * (pkin(9) ^ 2 * t80 + pkin(5) ^ 2) + 0.2e1 * pkin(9) * t74 + t68; mrSges(7,1) * t1 - mrSges(7,2) * t2 - Ifges(7,6) * t90 + t82; -t26 * t70 - t40; t39; -t70 * t58; -pkin(9) * t70 + t40; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t7(1) t7(2) t7(4) t7(7) t7(11) t7(16); t7(2) t7(3) t7(5) t7(8) t7(12) t7(17); t7(4) t7(5) t7(6) t7(9) t7(13) t7(18); t7(7) t7(8) t7(9) t7(10) t7(14) t7(19); t7(11) t7(12) t7(13) t7(14) t7(15) t7(20); t7(16) t7(17) t7(18) t7(19) t7(20) t7(21);];
Mq  = res;
