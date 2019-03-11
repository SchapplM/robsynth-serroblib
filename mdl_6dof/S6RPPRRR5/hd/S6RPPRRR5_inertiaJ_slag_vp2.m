% Calculate joint inertia matrix for
% S6RPPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
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
% Datum: 2019-03-09 02:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRR5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR5_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR5_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR5_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR5_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR5_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:28:16
% EndTime: 2019-03-09 02:28:18
% DurationCPUTime: 0.55s
% Computational Cost: add. (704->159), mult. (1213->219), div. (0->0), fcn. (1092->6), ass. (0->68)
t55 = sin(qJ(6));
t58 = cos(qJ(6));
t80 = t55 ^ 2 + t58 ^ 2;
t104 = mrSges(7,3) * t80;
t37 = -t58 * mrSges(7,1) + t55 * mrSges(7,2);
t103 = -mrSges(6,1) + t37;
t60 = cos(qJ(4));
t102 = t60 ^ 2;
t57 = sin(qJ(4));
t101 = t57 * mrSges(5,1) + t60 * mrSges(5,2) + mrSges(4,3);
t56 = sin(qJ(5));
t59 = cos(qJ(5));
t33 = t56 * t57 - t59 * t60;
t34 = t56 * t60 + t59 * t57;
t85 = t55 * mrSges(7,3);
t11 = -t34 * mrSges(7,2) + t33 * t85;
t86 = t33 * t58;
t12 = t34 * mrSges(7,1) + mrSges(7,3) * t86;
t100 = t58 * t11 - t55 * t12;
t69 = -mrSges(7,1) * t55 - mrSges(7,2) * t58;
t10 = t69 * t33;
t52 = qJ(2) - pkin(7);
t92 = -pkin(8) + t52;
t36 = t92 * t57;
t76 = t92 * t60;
t14 = t56 * t36 - t59 * t76;
t99 = m(7) * t14 + t10;
t16 = t59 * t36 + t56 * t76;
t53 = pkin(1) + qJ(3);
t40 = t57 * pkin(4) + t53;
t9 = t34 * pkin(5) + t33 * pkin(9) + t40;
t2 = -t55 * t16 + t58 * t9;
t3 = t58 * t16 + t55 * t9;
t93 = t3 * t58;
t72 = -t2 * t55 + t93;
t98 = m(7) * t72 + t100;
t97 = t14 ^ 2;
t30 = t33 ^ 2;
t96 = t53 ^ 2;
t95 = -2 * mrSges(6,3);
t94 = 0.2e1 * t53;
t90 = Ifges(7,4) * t55;
t89 = Ifges(7,4) * t58;
t88 = t33 * t14;
t87 = t33 * t55;
t82 = -Ifges(7,5) * t86 + Ifges(7,3) * t34;
t81 = Ifges(7,5) * t55 + Ifges(7,6) * t58;
t79 = t57 ^ 2 + t102;
t38 = Ifges(7,2) * t58 + t90;
t39 = Ifges(7,1) * t55 + t89;
t78 = t58 * t38 + t55 * t39 + Ifges(6,3);
t77 = m(5) * t79;
t75 = t80 * pkin(9);
t43 = t56 * pkin(4) + pkin(9);
t74 = t80 * t43;
t73 = t79 * mrSges(5,3);
t71 = t60 * mrSges(5,1) - t57 * mrSges(5,2);
t68 = t33 * t59 - t34 * t56;
t67 = 0.2e1 * t104;
t66 = t103 * t33 + (-mrSges(6,2) + t104) * t34;
t65 = (t59 * mrSges(6,1) - t56 * mrSges(6,2)) * pkin(4);
t6 = Ifges(7,6) * t34 + (Ifges(7,2) * t55 - t89) * t33;
t7 = Ifges(7,5) * t34 + (-Ifges(7,1) * t58 + t90) * t33;
t64 = -t16 * mrSges(6,2) + mrSges(7,3) * t93 - t2 * t85 + t38 * t87 / 0.2e1 - t39 * t86 / 0.2e1 - Ifges(6,5) * t33 + t55 * t7 / 0.2e1 + t58 * t6 / 0.2e1 + (t81 / 0.2e1 - Ifges(6,6)) * t34 + t103 * t14;
t61 = (qJ(2) ^ 2);
t44 = -t59 * pkin(4) - pkin(5);
t29 = t34 ^ 2;
t1 = [Ifges(5,1) * t102 + 0.2e1 * t2 * t12 + 0.2e1 * t14 * t10 + 0.2e1 * t3 * t11 + Ifges(3,1) + Ifges(4,1) + Ifges(2,3) - (2 * pkin(1) * mrSges(3,2)) + (0.2e1 * t40 * mrSges(6,1) + Ifges(6,2) * t34 + t16 * t95 + t82) * t34 + (-0.2e1 * t40 * mrSges(6,2) + t14 * t95 + Ifges(6,1) * t33 + t55 * t6 - t58 * t7 + (Ifges(7,6) * t55 + (2 * Ifges(6,4))) * t34) * t33 + (m(3) * (pkin(1) ^ 2 + t61)) + m(4) * (t61 + t96) + m(5) * (t52 ^ 2 * t79 + t96) + m(6) * (t16 ^ 2 + t40 ^ 2 + t97) + m(7) * (t2 ^ 2 + t3 ^ 2 + t97) + (-0.2e1 * Ifges(5,4) * t60 + Ifges(5,2) * t57) * t57 + t101 * t94 + (2 * (mrSges(3,3) + mrSges(4,2)) * qJ(2)) - 0.2e1 * t52 * t73; -(m(3) * pkin(1)) - t34 * mrSges(6,1) + t33 * mrSges(6,2) - t55 * t11 - t58 * t12 + mrSges(3,2) + m(7) * (-t58 * t2 - t55 * t3) - m(6) * t40 + (-m(5) / 0.2e1 - m(4) / 0.2e1) * t94 - t101; m(7) * t80 + m(3) + m(4) + m(5) + m(6); m(4) * qJ(2) - t30 * mrSges(6,3) + t33 * t10 + mrSges(4,2) - t73 + (-mrSges(6,3) * t34 + t100) * t34 + m(7) * (t34 * t72 + t88) + m(6) * (t34 * t16 + t88) + t52 * t77; 0; m(4) + t77 + m(6) * (t29 + t30) + m(7) * (t29 * t80 + t30); t64 + (m(6) * (-t14 * t59 + t16 * t56) + t68 * mrSges(6,3)) * pkin(4) + t71 * t52 - Ifges(5,6) * t57 + Ifges(5,5) * t60 + t99 * t44 + t98 * t43; 0; m(7) * (t44 * t33 + t34 * t74) - m(6) * t68 * pkin(4) + t66 + t71; 0.2e1 * t44 * t37 + Ifges(5,3) + 0.2e1 * t65 + t43 * t67 + m(7) * (t43 ^ 2 * t80 + t44 ^ 2) + m(6) * (t56 ^ 2 + t59 ^ 2) * pkin(4) ^ 2 + t78; -t99 * pkin(5) + t98 * pkin(9) + t64; 0; m(7) * (-pkin(5) * t33 + t34 * t75) + t66; m(7) * (-pkin(5) * t44 + pkin(9) * t74) + (t44 - pkin(5)) * t37 + t65 + (t74 + t75) * mrSges(7,3) + t78; -0.2e1 * pkin(5) * t37 + m(7) * (pkin(9) ^ 2 * t80 + pkin(5) ^ 2) + pkin(9) * t67 + t78; t2 * mrSges(7,1) - t3 * mrSges(7,2) + Ifges(7,6) * t87 + t82; t37; t69 * t34; t43 * t69 + t81; pkin(9) * t69 + t81; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
