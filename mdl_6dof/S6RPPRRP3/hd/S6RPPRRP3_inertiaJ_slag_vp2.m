% Calculate joint inertia matrix for
% S6RPPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2]';
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
% Datum: 2019-03-09 02:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRP3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP3_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP3_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP3_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP3_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP3_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:02:44
% EndTime: 2019-03-09 02:02:46
% DurationCPUTime: 0.58s
% Computational Cost: add. (543->179), mult. (1038->237), div. (0->0), fcn. (731->6), ass. (0->81)
t107 = Ifges(7,2) + Ifges(6,3);
t58 = sin(qJ(5));
t60 = cos(qJ(5));
t80 = t58 ^ 2 + t60 ^ 2;
t57 = cos(pkin(9));
t44 = -t57 * pkin(1) - pkin(2);
t33 = -pkin(7) + t44;
t106 = -0.2e1 * t33;
t105 = m(6) + m(7);
t104 = (mrSges(7,2) + mrSges(6,3)) * t80;
t56 = sin(pkin(9));
t42 = t56 * pkin(1) + qJ(3);
t103 = t42 ^ 2;
t102 = 0.2e1 * t42;
t59 = sin(qJ(4));
t100 = t59 * pkin(4);
t61 = cos(qJ(4));
t99 = t61 * pkin(8);
t53 = t59 ^ 2;
t55 = t61 ^ 2;
t79 = t55 + t53;
t98 = m(5) * t79 + m(4);
t97 = Ifges(6,4) * t58;
t96 = Ifges(6,4) * t60;
t95 = Ifges(7,5) * t58;
t94 = Ifges(7,5) * t60;
t93 = Ifges(6,6) * t59;
t92 = Ifges(7,6) * t59;
t91 = t33 * t59;
t90 = t58 * t59;
t89 = t58 * t61;
t14 = t42 - t99 + t100;
t88 = t60 * t14;
t86 = t61 * mrSges(7,2);
t85 = t61 * mrSges(6,3);
t4 = t58 * t14 + t60 * t91;
t27 = -t60 * mrSges(6,1) + t58 * mrSges(6,2);
t84 = t27 - mrSges(5,1);
t83 = t80 * pkin(8) * t59;
t82 = t80 * t99;
t81 = t80 * pkin(8) ^ 2;
t26 = -t60 * mrSges(7,1) - t58 * mrSges(7,3);
t77 = -t26 - t84;
t69 = t58 * mrSges(7,1) - t60 * mrSges(7,3);
t15 = t69 * t61;
t68 = -pkin(5) * t58 + qJ(6) * t60;
t5 = (-t33 - t68) * t61;
t76 = m(7) * t5 + t15;
t75 = t79 * mrSges(5,3);
t21 = -t59 * mrSges(7,1) + t60 * t86;
t73 = Ifges(7,6) * t89 + (Ifges(7,4) + Ifges(6,5)) * t60 * t61 + t107 * t59;
t1 = t59 * qJ(6) + t4;
t2 = -t88 + (t33 * t58 - pkin(5)) * t59;
t72 = t1 * t60 + t2 * t58;
t3 = -t33 * t90 + t88;
t71 = -t3 * t58 + t4 * t60;
t70 = t58 * mrSges(6,1) + t60 * mrSges(6,2);
t19 = -t59 * mrSges(6,2) - t58 * t85;
t20 = t59 * mrSges(6,1) - t60 * t85;
t22 = t59 * mrSges(7,3) - t58 * t86;
t67 = (t19 + t22) * t60 + (-t20 + t21) * t58;
t66 = -mrSges(5,2) + t104;
t65 = m(7) * t72 + t67;
t64 = m(7) * t68 - t69 - t70;
t49 = Ifges(7,4) * t58;
t48 = Ifges(6,5) * t58;
t46 = Ifges(6,6) * t60;
t32 = t33 ^ 2;
t31 = Ifges(6,1) * t58 + t96;
t30 = Ifges(7,1) * t58 - t94;
t29 = Ifges(6,2) * t60 + t97;
t28 = -Ifges(7,3) * t60 + t95;
t25 = t55 * t33;
t24 = t55 * t32;
t23 = -t60 * pkin(5) - t58 * qJ(6) - pkin(4);
t16 = t70 * t61;
t10 = Ifges(6,5) * t59 + (Ifges(6,1) * t60 - t97) * t61;
t9 = Ifges(7,4) * t59 + (Ifges(7,1) * t60 + t95) * t61;
t8 = t93 + (-Ifges(6,2) * t58 + t96) * t61;
t7 = t92 + (Ifges(7,3) * t58 + t94) * t61;
t6 = [0.2e1 * t44 * mrSges(4,2) + mrSges(4,3) * t102 + 0.2e1 * t1 * t22 + 0.2e1 * t5 * t15 + 0.2e1 * t4 * t19 + 0.2e1 * t2 * t21 + 0.2e1 * t3 * t20 + Ifges(4,1) + Ifges(2,3) + Ifges(3,3) + (mrSges(5,1) * t102 + Ifges(5,2) * t59 + t73) * t59 + (mrSges(5,2) * t102 + Ifges(5,1) * t61 - 0.2e1 * Ifges(5,4) * t59 + t16 * t106 + (t10 + t9) * t60 + (t7 - t8 - t93) * t58) * t61 + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) + m(6) * (t3 ^ 2 + t4 ^ 2 + t24) + m(4) * (t44 ^ 2 + t103) + m(5) * (t53 * t32 + t103 + t24) + m(3) * (t56 ^ 2 + t57 ^ 2) * pkin(1) ^ 2 + 0.2e1 * (t57 * mrSges(3,1) - t56 * mrSges(3,2)) * pkin(1) + t75 * t106; (t16 + t76) * t59 + (m(6) * (t71 - t91) + t65) * t61; m(3) + (t80 * t55 + t53) * t105 + t98; mrSges(4,2) + (-t15 - t16) * t61 - t75 + t67 * t59 + m(7) * (-t61 * t5 + t59 * t72) + m(6) * (t59 * t71 + t25) + m(5) * (t53 * t33 + t25) + m(4) * t44; (-0.1e1 + t80) * t61 * t59 * t105; (t80 * t53 + t55) * t105 + t98; -pkin(4) * t16 + t5 * t26 + t76 * t23 + (-t33 * mrSges(5,2) + t49 / 0.2e1 + t48 / 0.2e1 + t46 / 0.2e1 - Ifges(5,6)) * t59 + (t1 * mrSges(7,2) + t4 * mrSges(6,3) - t92 / 0.2e1 - t7 / 0.2e1 + t8 / 0.2e1) * t60 + (t9 / 0.2e1 + t10 / 0.2e1 + t2 * mrSges(7,2) - t3 * mrSges(6,3)) * t58 + (m(6) * t71 + t65) * pkin(8) + (Ifges(5,5) + (t30 / 0.2e1 + t31 / 0.2e1) * t60 + (t28 / 0.2e1 - t29 / 0.2e1) * t58 + (m(6) * pkin(4) - t84) * t33) * t61; -t77 * t59 + m(6) * (t82 - t100) + m(7) * (t23 * t59 + t82) + t66 * t61; t77 * t61 + m(7) * (-t23 * t61 + t83) + m(6) * (pkin(4) * t61 + t83) + t66 * t59; -0.2e1 * pkin(4) * t27 + 0.2e1 * t23 * t26 + Ifges(5,3) + (-t28 + t29) * t60 + (t30 + t31) * t58 + m(7) * (t23 ^ 2 + t81) + m(6) * (pkin(4) ^ 2 + t81) + 0.2e1 * pkin(8) * t104; -Ifges(6,6) * t89 + m(7) * (-pkin(5) * t2 + qJ(6) * t1) + qJ(6) * t22 + t1 * mrSges(7,3) - t4 * mrSges(6,2) + t3 * mrSges(6,1) - t2 * mrSges(7,1) - pkin(5) * t21 + t73; t64 * t61; t64 * t59; mrSges(7,2) * t68 - Ifges(7,6) * t60 + pkin(8) * t64 + t46 + t48 + t49; 0.2e1 * pkin(5) * mrSges(7,1) + 0.2e1 * qJ(6) * mrSges(7,3) + m(7) * (pkin(5) ^ 2 + qJ(6) ^ 2) + t107; m(7) * t2 + t21; m(7) * t89; m(7) * t90; (m(7) * pkin(8) + mrSges(7,2)) * t58; -m(7) * pkin(5) - mrSges(7,1); m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t6(1) t6(2) t6(4) t6(7) t6(11) t6(16); t6(2) t6(3) t6(5) t6(8) t6(12) t6(17); t6(4) t6(5) t6(6) t6(9) t6(13) t6(18); t6(7) t6(8) t6(9) t6(10) t6(14) t6(19); t6(11) t6(12) t6(13) t6(14) t6(15) t6(20); t6(16) t6(17) t6(18) t6(19) t6(20) t6(21);];
Mq  = res;
