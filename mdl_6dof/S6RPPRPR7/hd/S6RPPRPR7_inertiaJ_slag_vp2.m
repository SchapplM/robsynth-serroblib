% Calculate joint inertia matrix for
% S6RPPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 01:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRPR7_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR7_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR7_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR7_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR7_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR7_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:52:54
% EndTime: 2019-03-09 01:52:56
% DurationCPUTime: 0.65s
% Computational Cost: add. (1228->202), mult. (2226->292), div. (0->0), fcn. (2389->8), ass. (0->80)
t74 = -pkin(1) - qJ(3);
t102 = -pkin(7) + t74;
t71 = sin(pkin(9));
t52 = t102 * t71;
t76 = sin(qJ(4));
t100 = cos(qJ(4));
t73 = cos(pkin(9));
t85 = t100 * t73;
t29 = -t102 * t85 + t76 * t52;
t108 = t29 ^ 2;
t48 = t76 * t71 - t85;
t107 = t48 ^ 2;
t68 = t73 ^ 2;
t106 = -2 * mrSges(5,3);
t70 = sin(pkin(10));
t103 = t70 / 0.2e1;
t72 = cos(pkin(10));
t95 = t48 * t72;
t96 = t48 * t70;
t24 = -mrSges(6,1) * t96 - mrSges(6,2) * t95;
t75 = sin(qJ(6));
t77 = cos(qJ(6));
t50 = t77 * t70 + t75 * t72;
t21 = t50 * t48;
t80 = t75 * t70 - t77 * t72;
t23 = t80 * t48;
t7 = -t21 * mrSges(7,1) + t23 * mrSges(7,2);
t101 = t24 + t7;
t99 = Ifges(6,4) * t70;
t98 = Ifges(6,4) * t72;
t97 = t48 * t29;
t94 = t76 * t73;
t93 = pkin(8) + qJ(5);
t49 = t100 * t71 + t94;
t58 = t71 * pkin(3) + qJ(2);
t25 = t49 * pkin(4) + t48 * qJ(5) + t58;
t31 = t100 * t52 + t102 * t94;
t9 = t70 * t25 + t72 * t31;
t92 = Ifges(7,5) * t50 - Ifges(7,6) * t80;
t54 = -t72 * mrSges(6,1) + t70 * mrSges(6,2);
t91 = t54 - mrSges(5,1);
t90 = t71 * mrSges(4,1) + t73 * mrSges(4,2);
t89 = t70 ^ 2 + t72 ^ 2;
t88 = t71 ^ 2 + t68;
t87 = Ifges(7,5) * t23 + Ifges(7,6) * t21 + Ifges(7,3) * t49;
t86 = m(4) * t88;
t84 = t88 * mrSges(4,3);
t8 = t72 * t25 - t70 * t31;
t83 = qJ(5) * t89;
t82 = -t70 * t8 + t72 * t9;
t32 = mrSges(7,1) * t80 + t50 * mrSges(7,2);
t27 = -t49 * mrSges(6,2) + mrSges(6,3) * t96;
t28 = t49 * mrSges(6,1) + mrSges(6,3) * t95;
t81 = t72 * t27 - t70 * t28;
t79 = qJ(2) ^ 2;
t60 = -t72 * pkin(5) - pkin(4);
t57 = Ifges(6,1) * t70 + t98;
t56 = Ifges(6,2) * t72 + t99;
t55 = t93 * t72;
t53 = t93 * t70;
t45 = t49 ^ 2;
t40 = t48 * mrSges(5,2);
t36 = -t75 * t53 + t77 * t55;
t35 = -t77 * t53 - t75 * t55;
t34 = Ifges(7,1) * t50 - Ifges(7,4) * t80;
t33 = Ifges(7,4) * t50 - Ifges(7,2) * t80;
t22 = t80 * t49;
t20 = t50 * t49;
t14 = -pkin(5) * t96 + t29;
t13 = Ifges(6,5) * t49 + (-Ifges(6,1) * t72 + t99) * t48;
t12 = Ifges(6,6) * t49 + (Ifges(6,2) * t70 - t98) * t48;
t11 = t49 * mrSges(7,1) - t23 * mrSges(7,3);
t10 = -t49 * mrSges(7,2) + t21 * mrSges(7,3);
t6 = pkin(8) * t96 + t9;
t5 = Ifges(7,1) * t23 + Ifges(7,4) * t21 + Ifges(7,5) * t49;
t4 = Ifges(7,4) * t23 + Ifges(7,2) * t21 + Ifges(7,6) * t49;
t3 = t49 * pkin(5) + pkin(8) * t95 + t8;
t2 = t75 * t3 + t77 * t6;
t1 = t77 * t3 - t75 * t6;
t15 = [Ifges(4,1) * t68 - 0.2e1 * t58 * t40 + 0.2e1 * t29 * t24 + t21 * t4 + t23 * t5 + 0.2e1 * t9 * t27 + 0.2e1 * t8 * t28 + 0.2e1 * t2 * t10 + 0.2e1 * t1 * t11 + 0.2e1 * t14 * t7 + Ifges(3,1) + Ifges(2,3) - (2 * pkin(1) * mrSges(3,2)) + (-0.2e1 * Ifges(4,4) * t73 + Ifges(4,2) * t71) * t71 - 0.2e1 * t74 * t84 + (0.2e1 * t58 * mrSges(5,1) + t31 * t106 + (Ifges(6,3) + Ifges(5,2)) * t49 + t87) * t49 + (t29 * t106 + Ifges(5,1) * t48 + t70 * t12 - t72 * t13 + (-Ifges(6,5) * t72 + Ifges(6,6) * t70 + (2 * Ifges(5,4))) * t49) * t48 + m(4) * (t88 * t74 ^ 2 + t79) + m(3) * ((pkin(1) ^ 2) + t79) + m(5) * (t31 ^ 2 + t58 ^ 2 + t108) + m(6) * (t8 ^ 2 + t9 ^ 2 + t108) + m(7) * (t1 ^ 2 + t14 ^ 2 + t2 ^ 2) + 0.2e1 * (t90 + mrSges(3,3)) * qJ(2); -m(3) * pkin(1) - t107 * mrSges(5,3) - t22 * t10 - t20 * t11 + mrSges(3,2) + t101 * t48 - t84 + (-mrSges(5,3) * t49 + t81) * t49 + m(7) * (-t20 * t1 + t48 * t14 - t22 * t2) + m(6) * (t82 * t49 + t97) + m(5) * (t49 * t31 + t97) + t74 * t86; m(3) + m(7) * (t20 ^ 2 + t22 ^ 2 + t107) + m(6) * (t89 * t45 + t107) + m(5) * (t45 + t107) + t86; m(4) * qJ(2) + t49 * mrSges(5,1) + t50 * t10 - t80 * t11 + t70 * t27 + t72 * t28 - t40 + m(7) * (-t1 * t80 + t50 * t2) + m(6) * (t70 * t9 + t72 * t8) + m(5) * t58 + t90; m(7) * (t20 * t80 - t50 * t22); m(4) + m(5) + m(6) * t89 + m(7) * (t50 ^ 2 + t80 ^ 2); t13 * t103 + t72 * t12 / 0.2e1 + t60 * t7 - Ifges(5,6) * t49 + t50 * t5 / 0.2e1 - t80 * t4 / 0.2e1 - t31 * mrSges(5,2) + t14 * t32 + t21 * t33 / 0.2e1 + t23 * t34 / 0.2e1 + t35 * t11 + t36 * t10 - pkin(4) * t24 + t91 * t29 + t81 * qJ(5) + (-t1 * t50 - t2 * t80) * mrSges(7,3) + t82 * mrSges(6,3) + m(6) * (-pkin(4) * t29 + t82 * qJ(5)) + m(7) * (t35 * t1 + t60 * t14 + t36 * t2) + (-t72 * t57 / 0.2e1 + t56 * t103 - Ifges(5,5)) * t48 + (Ifges(6,5) * t70 + Ifges(6,6) * t72 + t92) * t49 / 0.2e1; (t20 * t50 + t22 * t80) * mrSges(7,3) + (t89 * mrSges(6,3) - mrSges(5,2)) * t49 + (t32 + t91) * t48 + m(7) * (-t20 * t35 - t22 * t36 + t48 * t60) + m(6) * (-pkin(4) * t48 + t49 * t83); m(7) * (-t35 * t80 + t36 * t50); -0.2e1 * pkin(4) * t54 + 0.2e1 * t60 * t32 - t80 * t33 + t50 * t34 + t72 * t56 + t70 * t57 + Ifges(5,3) + m(7) * (t35 ^ 2 + t36 ^ 2 + t60 ^ 2) + m(6) * (t89 * qJ(5) ^ 2 + pkin(4) ^ 2) + 0.2e1 * (-t35 * t50 - t36 * t80) * mrSges(7,3) + 0.2e1 * mrSges(6,3) * t83; m(6) * t29 + m(7) * t14 + t101; 0.2e1 * (m(7) / 0.2e1 + m(6) / 0.2e1) * t48; 0; -m(6) * pkin(4) + m(7) * t60 + t32 + t54; m(6) + m(7); mrSges(7,1) * t1 - mrSges(7,2) * t2 + t87; -mrSges(7,1) * t20 + mrSges(7,2) * t22; -t32; t35 * mrSges(7,1) - t36 * mrSges(7,2) + t92; 0; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t15(1) t15(2) t15(4) t15(7) t15(11) t15(16); t15(2) t15(3) t15(5) t15(8) t15(12) t15(17); t15(4) t15(5) t15(6) t15(9) t15(13) t15(18); t15(7) t15(8) t15(9) t15(10) t15(14) t15(19); t15(11) t15(12) t15(13) t15(14) t15(15) t15(20); t15(16) t15(17) t15(18) t15(19) t15(20) t15(21);];
Mq  = res;
