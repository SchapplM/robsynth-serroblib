% Calculate joint inertia matrix for
% S6RPPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta5]';
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
% Datum: 2019-03-09 01:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRPR5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR5_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR5_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR5_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR5_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR5_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:48:37
% EndTime: 2019-03-09 01:48:38
% DurationCPUTime: 0.53s
% Computational Cost: add. (655->192), mult. (1200->265), div. (0->0), fcn. (1058->6), ass. (0->75)
t60 = (qJ(2) - pkin(7));
t92 = -2 * t60;
t58 = sin(pkin(9));
t59 = cos(pkin(9));
t63 = sin(qJ(6));
t84 = cos(qJ(6));
t34 = t84 * t58 + t63 * t59;
t68 = -t63 * t58 + t84 * t59;
t8 = -mrSges(7,1) * t68 + t34 * mrSges(7,2);
t61 = (pkin(1) + qJ(3));
t91 = t61 ^ 2;
t90 = 2 * t61;
t89 = m(6) * pkin(4);
t88 = t58 / 0.2e1;
t87 = t59 / 0.2e1;
t86 = m(6) + m(7);
t65 = cos(qJ(4));
t79 = t59 * t65;
t80 = t58 * t65;
t26 = mrSges(6,1) * t80 + mrSges(6,2) * t79;
t23 = t34 * t65;
t25 = t68 * t65;
t5 = t23 * mrSges(7,1) + t25 * mrSges(7,2);
t85 = -t26 - t5;
t83 = Ifges(6,4) * t58;
t82 = Ifges(6,4) * t59;
t64 = sin(qJ(4));
t78 = t60 * t64;
t40 = -t59 * mrSges(6,1) + t58 * mrSges(6,2);
t77 = mrSges(5,1) - t40;
t76 = pkin(8) + qJ(5);
t38 = t64 * pkin(4) - t65 * qJ(5) + t61;
t16 = t58 * t38 + t59 * t78;
t75 = t58 ^ 2 + t59 ^ 2;
t56 = t64 ^ 2;
t57 = t65 ^ 2;
t74 = t57 + t56;
t73 = Ifges(7,5) * t25 - Ifges(7,6) * t23 + Ifges(7,3) * t64;
t72 = t74 * mrSges(5,3);
t71 = qJ(5) * t75;
t28 = t59 * t38;
t15 = -t58 * t78 + t28;
t70 = -t15 * t58 + t16 * t59;
t36 = -t64 * mrSges(6,2) - mrSges(6,3) * t80;
t37 = t64 * mrSges(6,1) - mrSges(6,3) * t79;
t69 = t59 * t36 - t58 * t37;
t67 = (qJ(2) ^ 2);
t55 = t60 ^ 2;
t50 = t57 * t60;
t49 = t57 * t55;
t48 = -t59 * pkin(5) - pkin(4);
t43 = Ifges(6,1) * t58 + t82;
t42 = Ifges(6,2) * t59 + t83;
t41 = t76 * t59;
t39 = t76 * t58;
t32 = (pkin(5) * t58 - t60) * t65;
t31 = Ifges(7,5) * t34;
t30 = Ifges(7,6) * t68;
t24 = t68 * t64;
t22 = t34 * t64;
t21 = Ifges(6,5) * t64 + (Ifges(6,1) * t59 - t83) * t65;
t20 = Ifges(6,6) * t64 + (-Ifges(6,2) * t58 + t82) * t65;
t14 = t64 * mrSges(7,1) - t25 * mrSges(7,3);
t13 = -t64 * mrSges(7,2) - t23 * mrSges(7,3);
t12 = -t63 * t39 + t84 * t41;
t11 = -t84 * t39 - t63 * t41;
t10 = Ifges(7,1) * t34 + Ifges(7,4) * t68;
t9 = Ifges(7,4) * t34 + Ifges(7,2) * t68;
t7 = -pkin(8) * t80 + t16;
t6 = -pkin(8) * t79 + t28 + (-t58 * t60 + pkin(5)) * t64;
t4 = Ifges(7,1) * t25 - Ifges(7,4) * t23 + Ifges(7,5) * t64;
t3 = Ifges(7,4) * t25 - Ifges(7,2) * t23 + Ifges(7,6) * t64;
t2 = t63 * t6 + t84 * t7;
t1 = t84 * t6 - t63 * t7;
t17 = [-(2 * pkin(1) * mrSges(3,2)) + (mrSges(4,3) * t90) + 0.2e1 * t1 * t14 + 0.2e1 * t2 * t13 + 0.2e1 * t15 * t37 + 0.2e1 * t16 * t36 - t23 * t3 + t25 * t4 + 0.2e1 * t32 * t5 + Ifges(3,1) + Ifges(4,1) + Ifges(2,3) + (mrSges(5,1) * t90 + (Ifges(6,3) + Ifges(5,2)) * t64 + t73) * t64 + ((mrSges(5,2) * t90) + Ifges(5,1) * t65 - t58 * t20 + t59 * t21 + t26 * t92 + (Ifges(6,5) * t59 - Ifges(6,6) * t58 - (2 * Ifges(5,4))) * t64) * t65 + m(4) * (t67 + t91) + (m(3) * (pkin(1) ^ 2 + t67)) + m(5) * (t56 * t55 + t49 + t91) + m(6) * (t15 ^ 2 + t16 ^ 2 + t49) + m(7) * (t1 ^ 2 + t2 ^ 2 + t32 ^ 2) + (2 * (mrSges(3,3) + mrSges(4,2)) * qJ(2)) + t72 * t92; -(m(3) * pkin(1)) - t64 * mrSges(5,1) - t65 * mrSges(5,2) - t34 * t13 - t68 * t14 - t58 * t36 - t59 * t37 + mrSges(3,2) - mrSges(4,3) + m(7) * (-t1 * t68 - t34 * t2) + m(6) * (-t59 * t15 - t58 * t16) + (-m(5) / 0.2e1 - m(4) / 0.2e1) * t90; m(3) + m(4) + m(5) + m(6) * t75 + m(7) * (t34 ^ 2 + t68 ^ 2); m(4) * qJ(2) + t24 * t13 - t22 * t14 + mrSges(4,2) + t85 * t65 + t69 * t64 - t72 + m(7) * (-t22 * t1 + t24 * t2 - t65 * t32) + m(6) * (t70 * t64 + t50) + m(5) * (t56 * t60 + t50); m(7) * (t22 * t68 - t24 * t34); m(4) + m(5) * t74 + m(6) * (t75 * t56 + t57) + m(7) * (t22 ^ 2 + t24 ^ 2 + t57); m(7) * (t11 * t1 + t12 * t2 + t48 * t32) + t21 * t88 + t20 * t87 + t48 * t5 + t68 * t3 / 0.2e1 + t34 * t4 / 0.2e1 + t25 * t10 / 0.2e1 - pkin(4) * t26 + t32 * t8 + t12 * t13 + t11 * t14 - t23 * t9 / 0.2e1 + (-(t60 * mrSges(5,2)) + Ifges(6,5) * t88 + Ifges(6,6) * t87 + t31 / 0.2e1 + t30 / 0.2e1 - Ifges(5,6)) * t64 + (-t1 * t34 + t2 * t68) * mrSges(7,3) + t70 * mrSges(6,3) + (m(6) * t70 + t69) * qJ(5) + (t43 * t87 - t58 * t42 / 0.2e1 + Ifges(5,5) + (t77 + t89) * t60) * t65; m(7) * (-t11 * t68 - t12 * t34); (t22 * t34 + t24 * t68) * mrSges(7,3) + (-t8 + t77) * t65 + (t75 * mrSges(6,3) - mrSges(5,2)) * t64 + m(6) * (pkin(4) * t65 + t64 * t71) + m(7) * (-t11 * t22 + t12 * t24 - t48 * t65); -0.2e1 * pkin(4) * t40 + t34 * t10 + t68 * t9 + t59 * t42 + t58 * t43 + 0.2e1 * t48 * t8 + Ifges(5,3) + m(7) * (t11 ^ 2 + t12 ^ 2 + t48 ^ 2) + m(6) * (t75 * qJ(5) ^ 2 + pkin(4) ^ 2) + 0.2e1 * (-t11 * t34 + t12 * t68) * mrSges(7,3) + 0.2e1 * mrSges(6,3) * t71; -m(6) * t65 * t60 + m(7) * t32 - t85; 0; -t86 * t65; m(7) * t48 + t40 + t8 - t89; t86; t1 * mrSges(7,1) - t2 * mrSges(7,2) + t73; t8; -t22 * mrSges(7,1) - t24 * mrSges(7,2); t11 * mrSges(7,1) - t12 * mrSges(7,2) + t30 + t31; 0; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t17(1) t17(2) t17(4) t17(7) t17(11) t17(16); t17(2) t17(3) t17(5) t17(8) t17(12) t17(17); t17(4) t17(5) t17(6) t17(9) t17(13) t17(18); t17(7) t17(8) t17(9) t17(10) t17(14) t17(19); t17(11) t17(12) t17(13) t17(14) t17(15) t17(20); t17(16) t17(17) t17(18) t17(19) t17(20) t17(21);];
Mq  = res;
