% Calculate joint inertia matrix for
% S6RPRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
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
% Datum: 2019-03-09 02:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPPR2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR2_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR2_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR2_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR2_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR2_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:41:16
% EndTime: 2019-03-09 02:41:18
% DurationCPUTime: 0.59s
% Computational Cost: add. (789->167), mult. (1423->227), div. (0->0), fcn. (1360->8), ass. (0->73)
t55 = cos(qJ(3));
t92 = t55 ^ 2;
t91 = m(6) + m(5);
t48 = sin(pkin(10));
t50 = cos(pkin(10));
t53 = sin(qJ(3));
t25 = t48 * t53 - t50 * t55;
t90 = 0.2e1 * t25;
t27 = t48 * t55 + t50 * t53;
t51 = cos(pkin(9));
t42 = -t51 * pkin(1) - pkin(2);
t30 = -t55 * pkin(3) + t42;
t58 = -t27 * qJ(5) + t30;
t8 = t25 * pkin(4) + t58;
t88 = -0.2e1 * t8;
t23 = t25 ^ 2;
t24 = t27 ^ 2;
t84 = t48 * pkin(3);
t37 = qJ(5) + t84;
t87 = t37 ^ 2;
t86 = 0.2e1 * t30;
t85 = m(5) * pkin(3);
t83 = t50 * pkin(3);
t52 = sin(qJ(6));
t54 = cos(qJ(6));
t71 = t52 ^ 2 + t54 ^ 2;
t29 = m(7) * t71;
t82 = m(6) + t29;
t81 = Ifges(7,4) * t52;
t80 = Ifges(7,4) * t54;
t79 = Ifges(7,6) * t52;
t78 = t25 * t52;
t77 = t25 * t54;
t76 = t37 * t27;
t75 = mrSges(6,1) + mrSges(5,3);
t74 = -mrSges(6,2) + mrSges(5,1);
t19 = t27 * mrSges(5,2);
t20 = t27 * mrSges(6,3);
t73 = t19 - t20;
t72 = t53 ^ 2 + t92;
t49 = sin(pkin(9));
t40 = t49 * pkin(1) + pkin(7);
t70 = qJ(4) + t40;
t69 = Ifges(7,5) * t78 + Ifges(7,6) * t77 + Ifges(7,3) * t27;
t22 = t70 * t55;
t65 = t70 * t53;
t11 = t50 * t22 - t48 * t65;
t9 = t48 * t22 + t50 * t65;
t68 = t11 ^ 2 + t9 ^ 2;
t41 = -pkin(4) - t83;
t67 = t71 * mrSges(7,3);
t35 = -pkin(8) + t41;
t66 = t71 * t35;
t3 = (pkin(4) + pkin(8)) * t25 + t58;
t4 = t27 * pkin(5) + t9;
t1 = -t52 * t3 + t54 * t4;
t2 = t54 * t3 + t52 * t4;
t64 = t54 * t1 + t52 * t2;
t62 = mrSges(6,2) - t67;
t61 = -t55 * mrSges(4,1) + t53 * mrSges(4,2);
t60 = t54 * mrSges(7,1) - t52 * mrSges(7,2);
t14 = t27 * mrSges(7,1) - mrSges(7,3) * t78;
t15 = -t27 * mrSges(7,2) + mrSges(7,3) * t77;
t59 = t54 * t14 + t52 * t15;
t43 = Ifges(7,5) * t54;
t33 = Ifges(7,1) * t54 - t81;
t32 = -Ifges(7,2) * t52 + t80;
t31 = t52 * mrSges(7,1) + t54 * mrSges(7,2);
t13 = t60 * t25;
t7 = Ifges(7,5) * t27 + (Ifges(7,1) * t52 + t80) * t25;
t6 = Ifges(7,6) * t27 + (Ifges(7,2) * t54 + t81) * t25;
t5 = -t25 * pkin(5) + t11;
t10 = [0.2e1 * t42 * t61 + Ifges(4,2) * t92 + t19 * t86 + t20 * t88 - 0.2e1 * t5 * t13 + 0.2e1 * t1 * t14 + 0.2e1 * t2 * t15 + Ifges(2,3) + Ifges(3,3) + 0.2e1 * (t51 * mrSges(3,1) - t49 * mrSges(3,2)) * pkin(1) + 0.2e1 * t72 * t40 * mrSges(4,3) + (mrSges(5,1) * t86 + mrSges(6,2) * t88 + t52 * t7 + t54 * t6 + (Ifges(5,2) + Ifges(6,3)) * t25 - 0.2e1 * t75 * t11) * t25 + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) + m(6) * (t8 ^ 2 + t68) + m(5) * (t30 ^ 2 + t68) + m(4) * (t72 * t40 ^ 2 + t42 ^ 2) + m(3) * (t49 ^ 2 + t51 ^ 2) * pkin(1) ^ 2 + ((Ifges(5,1) + Ifges(6,2)) * t27 + t69 + 0.2e1 * t75 * t9 + (-Ifges(5,4) - Ifges(6,6)) * t90) * t27 + (Ifges(4,1) * t53 + 0.2e1 * Ifges(4,4) * t55) * t53; -t27 * t13 + t59 * t25 + m(7) * (t64 * t25 + t5 * t27) + t91 * (t11 * t27 + t9 * t25); m(3) + m(7) * (t71 * t23 + t24) + m(4) * t72 + t91 * (t23 + t24); Ifges(4,5) * t53 + Ifges(4,6) * t55 - t37 * t13 + t5 * t31 - t74 * t9 + (-t53 * mrSges(4,1) - t55 * mrSges(4,2)) * t40 + (-mrSges(5,2) + mrSges(6,3)) * t11 + (-t1 * mrSges(7,3) + t35 * t14 + t7 / 0.2e1) * t54 + (-t2 * mrSges(7,3) + t35 * t15 - t6 / 0.2e1) * t52 + m(7) * (t64 * t35 + t37 * t5) + m(6) * (t37 * t11 + t41 * t9) + (t11 * t48 - t50 * t9) * t85 + (-t79 / 0.2e1 + t43 / 0.2e1 + Ifges(5,5) - Ifges(6,4) - mrSges(5,3) * t83 + t41 * mrSges(6,1)) * t27 + (t52 * t33 / 0.2e1 + t54 * t32 / 0.2e1 - Ifges(5,6) + Ifges(6,5) - mrSges(5,3) * t84 - t37 * mrSges(6,1)) * t25; t27 * t31 + (-mrSges(5,1) + t62) * t25 + m(6) * (t41 * t25 + t76) + m(7) * (t25 * t66 + t76) + (-t25 * t50 + t27 * t48) * t85 - t61 - t73; 0.2e1 * t41 * mrSges(6,2) - t52 * t32 + t54 * t33 + Ifges(6,1) + Ifges(4,3) + Ifges(5,3) + m(7) * (t71 * t35 ^ 2 + t87) + m(6) * (t41 ^ 2 + t87) + m(5) * (t48 ^ 2 + t50 ^ 2) * pkin(3) ^ 2 + 0.2e1 * (t31 + mrSges(6,3)) * t37 + 0.2e1 * (t50 * mrSges(5,1) - t48 * mrSges(5,2)) * pkin(3) - 0.2e1 * t35 * t67; -t52 * t14 + t54 * t15 + t74 * t25 + m(7) * (-t52 * t1 + t54 * t2) + m(6) * t8 + m(5) * t30 + t73; 0; 0; m(5) + t82; m(6) * t9 + m(7) * t64 + t27 * mrSges(6,1) + t59; (m(6) / 0.2e1 + t29 / 0.2e1) * t90; m(6) * t41 + m(7) * t66 + t62; 0; t82; t1 * mrSges(7,1) - t2 * mrSges(7,2) + t69; t13; t60 * t35 + t43 - t79; -t31; t60; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t10(1) t10(2) t10(4) t10(7) t10(11) t10(16); t10(2) t10(3) t10(5) t10(8) t10(12) t10(17); t10(4) t10(5) t10(6) t10(9) t10(13) t10(18); t10(7) t10(8) t10(9) t10(10) t10(14) t10(19); t10(11) t10(12) t10(13) t10(14) t10(15) t10(20); t10(16) t10(17) t10(18) t10(19) t10(20) t10(21);];
Mq  = res;
