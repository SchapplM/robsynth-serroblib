% Calculate joint inertia matrix for
% S6RPRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4,theta5]';
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
% Datum: 2019-03-09 02:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPPR1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR1_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPPR1_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR1_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR1_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR1_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:38:01
% EndTime: 2019-03-09 02:38:03
% DurationCPUTime: 0.72s
% Computational Cost: add. (1339->212), mult. (2502->317), div. (0->0), fcn. (2644->10), ass. (0->80)
t77 = cos(qJ(3));
t103 = t77 ^ 2;
t70 = sin(pkin(9));
t60 = pkin(1) * t70 + pkin(7);
t86 = qJ(4) + t60;
t43 = t86 * t77;
t69 = sin(pkin(10));
t72 = cos(pkin(10));
t75 = sin(qJ(3));
t83 = t86 * t75;
t25 = t43 * t69 + t72 * t83;
t102 = t25 ^ 2;
t46 = t69 * t75 - t72 * t77;
t44 = t46 ^ 2;
t101 = 0.2e1 * t25;
t73 = cos(pkin(9));
t62 = -pkin(1) * t73 - pkin(2);
t52 = -pkin(3) * t77 + t62;
t100 = 0.2e1 * t52;
t71 = cos(pkin(11));
t98 = t71 / 0.2e1;
t58 = pkin(3) * t69 + qJ(5);
t97 = pkin(8) + t58;
t49 = t69 * t77 + t72 * t75;
t91 = t49 * t71;
t68 = sin(pkin(11));
t92 = t49 * t68;
t28 = mrSges(6,1) * t92 + mrSges(6,2) * t91;
t74 = sin(qJ(6));
t76 = cos(qJ(6));
t50 = t68 * t76 + t71 * t74;
t21 = t50 * t49;
t48 = -t68 * t74 + t71 * t76;
t22 = t48 * t49;
t9 = t21 * mrSges(7,1) + t22 * mrSges(7,2);
t96 = t28 + t9;
t95 = Ifges(6,4) * t68;
t94 = Ifges(6,4) * t71;
t93 = t25 * t46;
t20 = pkin(4) * t46 - qJ(5) * t49 + t52;
t27 = t72 * t43 - t69 * t83;
t8 = t68 * t20 + t71 * t27;
t90 = Ifges(7,5) * t50 + Ifges(7,6) * t48;
t53 = -t71 * mrSges(6,1) + t68 * mrSges(6,2);
t89 = t53 - mrSges(5,1);
t88 = t68 ^ 2 + t71 ^ 2;
t87 = t75 ^ 2 + t103;
t85 = Ifges(7,5) * t22 - Ifges(7,6) * t21 + Ifges(7,3) * t46;
t61 = -pkin(3) * t72 - pkin(4);
t84 = t88 * t58;
t7 = t71 * t20 - t27 * t68;
t82 = -t7 * t68 + t8 * t71;
t81 = -t77 * mrSges(4,1) + t75 * mrSges(4,2);
t31 = -t48 * mrSges(7,1) + t50 * mrSges(7,2);
t29 = -mrSges(6,2) * t46 - mrSges(6,3) * t92;
t30 = mrSges(6,1) * t46 - mrSges(6,3) * t91;
t80 = t71 * t29 - t68 * t30;
t55 = Ifges(6,1) * t68 + t94;
t54 = Ifges(6,2) * t71 + t95;
t51 = -pkin(5) * t71 + t61;
t45 = t49 ^ 2;
t39 = t49 * mrSges(5,2);
t38 = t97 * t71;
t37 = t97 * t68;
t33 = Ifges(7,1) * t50 + Ifges(7,4) * t48;
t32 = Ifges(7,4) * t50 + Ifges(7,2) * t48;
t24 = -t37 * t74 + t38 * t76;
t23 = -t37 * t76 - t38 * t74;
t14 = Ifges(6,5) * t46 + (Ifges(6,1) * t71 - t95) * t49;
t13 = Ifges(6,6) * t46 + (-Ifges(6,2) * t68 + t94) * t49;
t12 = pkin(5) * t92 + t25;
t11 = mrSges(7,1) * t46 - mrSges(7,3) * t22;
t10 = -mrSges(7,2) * t46 - mrSges(7,3) * t21;
t6 = Ifges(7,1) * t22 - Ifges(7,4) * t21 + Ifges(7,5) * t46;
t5 = Ifges(7,4) * t22 - Ifges(7,2) * t21 + Ifges(7,6) * t46;
t4 = -pkin(8) * t92 + t8;
t3 = pkin(5) * t46 - pkin(8) * t91 + t7;
t2 = t3 * t74 + t4 * t76;
t1 = t3 * t76 - t4 * t74;
t15 = [0.2e1 * t62 * t81 + Ifges(4,2) * t103 + t39 * t100 + t28 * t101 + 0.2e1 * t8 * t29 + 0.2e1 * t7 * t30 - t21 * t5 + t22 * t6 + 0.2e1 * t2 * t10 + 0.2e1 * t1 * t11 + 0.2e1 * t12 * t9 + Ifges(2,3) + Ifges(3,3) + (mrSges(5,1) * t100 - 0.2e1 * t27 * mrSges(5,3) + (Ifges(6,3) + Ifges(5,2)) * t46 + t85) * t46 + (mrSges(5,3) * t101 + Ifges(5,1) * t49 - t68 * t13 + t71 * t14 + (Ifges(6,5) * t71 - Ifges(6,6) * t68 - (2 * Ifges(5,4))) * t46) * t49 + m(5) * (t27 ^ 2 + t52 ^ 2 + t102) + m(4) * (t60 ^ 2 * t87 + t62 ^ 2) + m(6) * (t7 ^ 2 + t8 ^ 2 + t102) + m(7) * (t1 ^ 2 + t12 ^ 2 + t2 ^ 2) + m(3) * (t70 ^ 2 + t73 ^ 2) * pkin(1) ^ 2 + (Ifges(4,1) * t75 + 0.2e1 * Ifges(4,4) * t77) * t75 + 0.2e1 * (mrSges(3,1) * t73 - mrSges(3,2) * t70) * pkin(1) + 0.2e1 * t87 * t60 * mrSges(4,3); t22 * t10 - t21 * t11 + t80 * t49 + t96 * t46 + m(7) * (-t1 * t21 + t12 * t46 + t2 * t22) + m(5) * (t27 * t49 + t93) + m(6) * (t49 * t82 + t93); m(3) + m(7) * (t21 ^ 2 + t22 ^ 2 + t44) + m(5) * (t45 + t44) + m(6) * (t45 * t88 + t44) + m(4) * t87; m(7) * (t1 * t23 + t12 * t51 + t2 * t24) + (-t75 * mrSges(4,1) - t77 * mrSges(4,2)) * t60 + (-t68 * t54 / 0.2e1 + t55 * t98 + Ifges(5,5)) * t49 + t89 * t25 + (-t1 * t50 + t2 * t48) * mrSges(7,3) + t80 * t58 + t82 * mrSges(6,3) + m(6) * (t25 * t61 + t58 * t82) + (m(5) * (-t25 * t72 + t27 * t69) + (-t69 * t46 - t72 * t49) * mrSges(5,3)) * pkin(3) + (Ifges(6,5) * t68 + Ifges(6,6) * t71 + t90) * t46 / 0.2e1 + Ifges(4,5) * t75 + Ifges(4,6) * t77 + t68 * t14 / 0.2e1 + t50 * t6 / 0.2e1 + t51 * t9 + t61 * t28 - Ifges(5,6) * t46 + t48 * t5 / 0.2e1 + t23 * t11 + t24 * t10 - t27 * mrSges(5,2) + t12 * t31 - t21 * t32 / 0.2e1 + t22 * t33 / 0.2e1 + t13 * t98; -t39 + (t21 * t50 + t22 * t48) * mrSges(7,3) + t88 * t49 * mrSges(6,3) + (t31 + t89) * t46 + m(7) * (-t21 * t23 + t22 * t24 + t46 * t51) + m(6) * (t46 * t61 + t49 * t84) + m(5) * (-t46 * t72 + t49 * t69) * pkin(3) - t81; 0.2e1 * t51 * t31 + t48 * t32 + t50 * t33 + 0.2e1 * t61 * t53 + t71 * t54 + t68 * t55 + Ifges(4,3) + Ifges(5,3) + m(7) * (t23 ^ 2 + t24 ^ 2 + t51 ^ 2) + m(6) * (t58 ^ 2 * t88 + t61 ^ 2) + m(5) * (t69 ^ 2 + t72 ^ 2) * pkin(3) ^ 2 + 0.2e1 * (mrSges(5,1) * t72 - mrSges(5,2) * t69) * pkin(3) + 0.2e1 * (-t23 * t50 + t24 * t48) * mrSges(7,3) + 0.2e1 * mrSges(6,3) * t84; t46 * mrSges(5,1) + t50 * t10 + t48 * t11 + t68 * t29 + t71 * t30 + t39 + m(7) * (t1 * t48 + t2 * t50) + m(6) * (t68 * t8 + t7 * t71) + m(5) * t52; m(7) * (-t21 * t48 + t22 * t50); m(7) * (t23 * t48 + t24 * t50); m(5) + m(6) * t88 + m(7) * (t48 ^ 2 + t50 ^ 2); m(6) * t25 + m(7) * t12 + t96; 0.2e1 * (m(7) / 0.2e1 + m(6) / 0.2e1) * t46; m(6) * t61 + m(7) * t51 + t31 + t53; 0; m(6) + m(7); mrSges(7,1) * t1 - mrSges(7,2) * t2 + t85; -t9; mrSges(7,1) * t23 - t24 * mrSges(7,2) + t90; -t31; 0; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t15(1) t15(2) t15(4) t15(7) t15(11) t15(16); t15(2) t15(3) t15(5) t15(8) t15(12) t15(17); t15(4) t15(5) t15(6) t15(9) t15(13) t15(18); t15(7) t15(8) t15(9) t15(10) t15(14) t15(19); t15(11) t15(12) t15(13) t15(14) t15(15) t15(20); t15(16) t15(17) t15(18) t15(19) t15(20) t15(21);];
Mq  = res;
