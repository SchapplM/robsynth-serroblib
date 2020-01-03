% Calculate joint inertia matrix for
% S5RRPRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP9_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP9_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP9_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP9_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP9_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP9_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:05:43
% EndTime: 2019-12-31 20:05:45
% DurationCPUTime: 0.67s
% Computational Cost: add. (735->202), mult. (1511->272), div. (0->0), fcn. (1409->6), ass. (0->77)
t103 = 2 * pkin(6);
t102 = -m(6) * pkin(4) - mrSges(6,1);
t78 = sin(pkin(8));
t101 = -t78 / 0.2e1;
t79 = cos(pkin(8));
t100 = t79 / 0.2e1;
t82 = cos(qJ(2));
t99 = pkin(6) * t82;
t81 = sin(qJ(2));
t73 = t81 * pkin(6);
t98 = cos(qJ(4));
t97 = Ifges(4,4) * t78;
t96 = Ifges(4,4) * t79;
t95 = t78 * t81;
t94 = t79 * t81;
t92 = -Ifges(5,3) - Ifges(6,2);
t91 = pkin(7) + qJ(3);
t59 = -t82 * pkin(2) - t81 * qJ(3) - pkin(1);
t52 = t79 * t59;
t14 = -pkin(7) * t94 + t52 + (-pkin(6) * t78 - pkin(3)) * t82;
t32 = t78 * t59 + t79 * t99;
t22 = -pkin(7) * t95 + t32;
t80 = sin(qJ(4));
t4 = t80 * t14 + t98 * t22;
t85 = -t80 * t78 + t98 * t79;
t43 = t85 * t81;
t30 = t82 * mrSges(6,1) + t43 * mrSges(6,2);
t44 = mrSges(4,1) * t95 + mrSges(4,2) * t94;
t58 = pkin(3) * t95 + t73;
t90 = t78 ^ 2 + t79 ^ 2;
t61 = t91 * t79;
t88 = t91 * t78;
t24 = t80 * t61 + t98 * t88;
t26 = t98 * t61 - t80 * t88;
t89 = t24 ^ 2 + t26 ^ 2;
t69 = -t79 * pkin(3) - pkin(2);
t60 = -t79 * mrSges(4,1) + t78 * mrSges(4,2);
t55 = t98 * t78 + t80 * t79;
t42 = t55 * t81;
t11 = t42 * mrSges(5,1) + t43 * mrSges(5,2);
t17 = -mrSges(5,1) * t85 + t55 * mrSges(5,2);
t10 = t42 * mrSges(6,1) - t43 * mrSges(6,3);
t16 = -mrSges(6,1) * t85 - t55 * mrSges(6,3);
t87 = (-Ifges(6,4) - Ifges(5,5)) * t43 + (Ifges(5,6) - Ifges(6,6)) * t42;
t31 = -t78 * t99 + t52;
t86 = -t31 * t78 + t32 * t79;
t3 = t98 * t14 - t80 * t22;
t84 = pkin(6) ^ 2;
t77 = t82 ^ 2;
t76 = t81 ^ 2;
t72 = t76 * t84;
t63 = Ifges(4,1) * t78 + t96;
t62 = Ifges(4,2) * t79 + t97;
t57 = -t82 * mrSges(4,1) - mrSges(4,3) * t94;
t56 = t82 * mrSges(4,2) - mrSges(4,3) * t95;
t50 = Ifges(6,4) * t55;
t49 = Ifges(5,5) * t55;
t48 = Ifges(5,6) * t85;
t47 = Ifges(6,6) * t85;
t41 = -Ifges(4,5) * t82 + (Ifges(4,1) * t79 - t97) * t81;
t40 = -Ifges(4,6) * t82 + (-Ifges(4,2) * t78 + t96) * t81;
t29 = -t82 * mrSges(5,1) - t43 * mrSges(5,3);
t28 = t82 * mrSges(5,2) - t42 * mrSges(5,3);
t27 = -t42 * mrSges(6,2) - t82 * mrSges(6,3);
t21 = Ifges(5,1) * t55 + Ifges(5,4) * t85;
t20 = Ifges(6,1) * t55 - Ifges(6,5) * t85;
t19 = Ifges(5,4) * t55 + Ifges(5,2) * t85;
t18 = Ifges(6,5) * t55 - Ifges(6,3) * t85;
t13 = -pkin(4) * t85 - t55 * qJ(5) + t69;
t9 = Ifges(5,1) * t43 - Ifges(5,4) * t42 - Ifges(5,5) * t82;
t8 = Ifges(6,1) * t43 - Ifges(6,4) * t82 + Ifges(6,5) * t42;
t7 = Ifges(5,4) * t43 - Ifges(5,2) * t42 - Ifges(5,6) * t82;
t6 = Ifges(6,5) * t43 - Ifges(6,6) * t82 + Ifges(6,3) * t42;
t5 = t42 * pkin(4) - t43 * qJ(5) + t58;
t2 = t82 * pkin(4) - t3;
t1 = -t82 * qJ(5) + t4;
t12 = [0.2e1 * t1 * t27 + 0.2e1 * t5 * t10 + 0.2e1 * t58 * t11 + 0.2e1 * t2 * t30 + 0.2e1 * t4 * t28 + 0.2e1 * t3 * t29 + 0.2e1 * t31 * t57 + 0.2e1 * t32 * t56 + Ifges(2,3) + (t8 + t9) * t43 + (t6 - t7) * t42 + (t76 + t77) * mrSges(3,3) * t103 + (-0.2e1 * pkin(1) * mrSges(3,2) + Ifges(3,1) * t81 + t44 * t103 - t78 * t40 + t79 * t41) * t81 + m(6) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) + m(5) * (t3 ^ 2 + t4 ^ 2 + t58 ^ 2) + m(4) * (t31 ^ 2 + t32 ^ 2 + t72) + m(3) * (pkin(1) ^ 2 + t77 * t84 + t72) + (0.2e1 * pkin(1) * mrSges(3,1) + (Ifges(4,3) + Ifges(3,2) - t92) * t82 + (-Ifges(4,5) * t79 + Ifges(4,6) * t78 + (2 * Ifges(3,4))) * t81 + t87) * t82; t78 * t41 / 0.2e1 + t40 * t100 + t69 * t11 + t58 * t17 - pkin(2) * t44 + t13 * t10 + t5 * t16 + (Ifges(4,5) * t101 - Ifges(4,6) * t79 / 0.2e1 - t49 / 0.2e1 - t48 / 0.2e1 - t50 / 0.2e1 + t47 / 0.2e1 + Ifges(3,6) - pkin(6) * mrSges(3,2)) * t82 + (t20 / 0.2e1 + t21 / 0.2e1) * t43 + (t18 / 0.2e1 - t19 / 0.2e1) * t42 + (t27 + t28) * t26 + (-t29 + t30) * t24 + (t79 * t56 - t78 * t57) * qJ(3) + t86 * mrSges(4,3) + (Ifges(3,5) + t62 * t101 + t63 * t100 + (t60 - mrSges(3,1)) * pkin(6)) * t81 + m(6) * (t26 * t1 + t13 * t5 + t24 * t2) + m(5) * (-t24 * t3 + t26 * t4 + t69 * t58) + m(4) * (-pkin(2) * t73 + t86 * qJ(3)) + (t9 / 0.2e1 + t8 / 0.2e1 + t2 * mrSges(6,2) - t3 * mrSges(5,3)) * t55 - (t6 / 0.2e1 - t7 / 0.2e1 - t1 * mrSges(6,2) - t4 * mrSges(5,3)) * t85; -0.2e1 * pkin(2) * t60 + 0.2e1 * t13 * t16 + 0.2e1 * t69 * t17 + t79 * t62 + t78 * t63 + Ifges(3,3) + m(6) * (t13 ^ 2 + t89) + m(5) * (t69 ^ 2 + t89) + m(4) * (t90 * qJ(3) ^ 2 + pkin(2) ^ 2) + (t20 + t21) * t55 - (t18 - t19) * t85 + 0.2e1 * t90 * qJ(3) * mrSges(4,3) + 0.2e1 * (t24 * t55 + t26 * t85) * (mrSges(6,2) + mrSges(5,3)); m(4) * t73 + m(5) * t58 + m(6) * t5 + t10 + t11 + t44; -m(4) * pkin(2) + m(5) * t69 + m(6) * t13 + t16 + t17 + t60; m(4) + m(5) + m(6); -pkin(4) * t30 + m(6) * (-pkin(4) * t2 + qJ(5) * t1) + qJ(5) * t27 + t1 * mrSges(6,3) - t4 * mrSges(5,2) + t3 * mrSges(5,1) - t2 * mrSges(6,1) + t92 * t82 - t87; t50 - t47 + t49 + t48 + (-pkin(4) * t55 + qJ(5) * t85) * mrSges(6,2) + (m(6) * qJ(5) - mrSges(5,2) + mrSges(6,3)) * t26 + (-mrSges(5,1) + t102) * t24; 0; 0.2e1 * pkin(4) * mrSges(6,1) + 0.2e1 * qJ(5) * mrSges(6,3) + m(6) * (pkin(4) ^ 2 + qJ(5) ^ 2) - t92; m(6) * t2 + t30; m(6) * t24 + t55 * mrSges(6,2); 0; t102; m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t12(1), t12(2), t12(4), t12(7), t12(11); t12(2), t12(3), t12(5), t12(8), t12(12); t12(4), t12(5), t12(6), t12(9), t12(13); t12(7), t12(8), t12(9), t12(10), t12(14); t12(11), t12(12), t12(13), t12(14), t12(15);];
Mq = res;
