% Calculate joint inertia matrix for
% S6RRPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5]';
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
% Datum: 2019-03-09 08:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRP3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP3_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPPRP3_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP3_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP3_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRP3_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:33:36
% EndTime: 2019-03-09 08:33:38
% DurationCPUTime: 0.71s
% Computational Cost: add. (606->214), mult. (988->260), div. (0->0), fcn. (702->4), ass. (0->74)
t94 = Ifges(6,6) + Ifges(7,6);
t93 = Ifges(6,3) + Ifges(7,3);
t58 = sin(qJ(2));
t60 = cos(qJ(2));
t92 = t58 ^ 2 + t60 ^ 2;
t91 = 2 * mrSges(7,3);
t40 = t60 * qJ(4);
t29 = t60 * pkin(7) - t40;
t90 = t29 ^ 2;
t89 = -2 * mrSges(5,3);
t25 = -t60 * pkin(2) - t58 * qJ(3) - pkin(1);
t88 = -0.2e1 * t25;
t87 = m(7) * pkin(5);
t61 = -pkin(2) - pkin(3);
t57 = sin(qJ(5));
t59 = cos(qJ(5));
t35 = t57 ^ 2 + t59 ^ 2;
t18 = m(6) * t35;
t86 = t59 * pkin(5);
t28 = (pkin(7) - qJ(4)) * t58;
t15 = t60 * pkin(3) - t25;
t7 = t58 * pkin(4) + t60 * pkin(8) + t15;
t4 = t59 * t28 + t57 * t7;
t85 = Ifges(6,4) * t57;
t84 = Ifges(6,4) * t59;
t83 = Ifges(7,4) * t57;
t82 = Ifges(7,4) * t59;
t81 = t57 * t60;
t80 = t59 * t60;
t79 = -mrSges(6,2) - mrSges(7,2);
t78 = -mrSges(5,3) + mrSges(4,2);
t77 = -Ifges(6,5) - Ifges(7,5);
t56 = qJ(3) + pkin(4);
t19 = -t58 * mrSges(7,2) + mrSges(7,3) * t81;
t20 = -t58 * mrSges(6,2) + mrSges(6,3) * t81;
t76 = t19 + t20;
t21 = t58 * mrSges(7,1) + mrSges(7,3) * t80;
t22 = t58 * mrSges(6,1) + mrSges(6,3) * t80;
t75 = t21 + t22;
t74 = t92 * pkin(7) ^ 2;
t73 = qJ(6) * t60;
t51 = -pkin(8) + t61;
t72 = qJ(6) - t51;
t71 = m(7) * t35 + m(5) + t18;
t70 = -m(4) * pkin(2) - mrSges(4,1);
t69 = mrSges(7,1) + t87;
t3 = -t57 * t28 + t59 * t7;
t68 = t77 * t59;
t67 = t35 * mrSges(6,3);
t42 = t59 * mrSges(7,1);
t26 = -t57 * mrSges(7,2) + t42;
t66 = t93 * t58 + t94 * t81;
t64 = -t57 * t3 + t59 * t4;
t13 = (-t57 * mrSges(7,1) - t59 * mrSges(7,2)) * t60;
t62 = qJ(3) ^ 2;
t43 = t59 * mrSges(6,1);
t41 = t58 * mrSges(5,1);
t36 = t56 + t86;
t34 = -Ifges(6,1) * t57 - t84;
t33 = -Ifges(7,1) * t57 - t82;
t32 = -Ifges(6,2) * t59 - t85;
t31 = -Ifges(7,2) * t59 - t83;
t27 = -t57 * mrSges(6,2) + t43;
t24 = t72 * t59;
t23 = t72 * t57;
t14 = (-t57 * mrSges(6,1) - t59 * mrSges(6,2)) * t60;
t12 = -t40 + (-pkin(5) * t57 + pkin(7)) * t60;
t11 = Ifges(6,5) * t58 + (-Ifges(6,1) * t59 + t85) * t60;
t10 = Ifges(7,5) * t58 + (-Ifges(7,1) * t59 + t83) * t60;
t9 = Ifges(6,6) * t58 + (Ifges(6,2) * t57 - t84) * t60;
t8 = Ifges(7,6) * t58 + (Ifges(7,2) * t57 - t82) * t60;
t2 = t57 * t73 + t4;
t1 = t58 * pkin(5) + t59 * t73 + t3;
t5 = [0.2e1 * t1 * t21 + 0.2e1 * t12 * t13 + 0.2e1 * t29 * t14 + 0.2e1 * t15 * t41 + 0.2e1 * t2 * t19 + 0.2e1 * t4 * t20 + 0.2e1 * t3 * t22 + Ifges(2,3) + m(4) * (t25 ^ 2 + t74) + m(3) * (pkin(1) ^ 2 + t74) + m(6) * (t3 ^ 2 + t4 ^ 2 + t90) + m(5) * (t15 ^ 2 + t28 ^ 2 + t90) + m(7) * (t1 ^ 2 + t12 ^ 2 + t2 ^ 2) + (-0.2e1 * pkin(1) * mrSges(3,2) + mrSges(4,3) * t88 + t28 * t89 + (Ifges(5,2) + Ifges(4,1) + Ifges(3,1)) * t58 + t66) * t58 + (0.2e1 * pkin(1) * mrSges(3,1) + mrSges(4,1) * t88 - 0.2e1 * t15 * mrSges(5,2) + t29 * t89 + (Ifges(5,1) + Ifges(4,3) + Ifges(3,2)) * t60 + (-t10 - t11) * t59 + (t8 + t9) * t57 + ((2 * Ifges(3,4)) + (2 * Ifges(5,4)) - (2 * Ifges(4,5)) + t68) * t58) * t60 + 0.2e1 * (mrSges(4,2) + mrSges(3,3)) * pkin(7) * t92; t28 * mrSges(5,2) + t12 * t26 + t36 * t13 + t56 * t14 - t24 * t19 + t23 * t21 + (t27 + mrSges(5,1)) * t29 + (t51 * t20 - t2 * mrSges(7,3) - t4 * mrSges(6,3) - t8 / 0.2e1 - t9 / 0.2e1) * t59 + (-t51 * t22 + t3 * mrSges(6,3) + t1 * mrSges(7,3) - t10 / 0.2e1 - t11 / 0.2e1) * t57 + m(5) * (qJ(3) * t29 + t61 * t28) + m(6) * (t56 * t29 + t64 * t51) + m(7) * (t23 * t1 + t36 * t12 - t24 * t2) + (-pkin(2) * mrSges(4,2) - t61 * mrSges(5,3) + Ifges(4,4) + Ifges(3,5) + Ifges(5,6) + (-Ifges(6,6) / 0.2e1 - Ifges(7,6) / 0.2e1) * t59 + (-Ifges(6,5) / 0.2e1 - Ifges(7,5) / 0.2e1) * t57 + (-mrSges(3,1) + t70) * pkin(7)) * t58 + (Ifges(3,6) + Ifges(5,5) - Ifges(4,6) + (-t33 / 0.2e1 - t34 / 0.2e1) * t59 + (t31 / 0.2e1 + t32 / 0.2e1) * t57 + t78 * qJ(3) + (m(4) * qJ(3) - mrSges(3,2) + mrSges(4,3)) * pkin(7)) * t60; 0.2e1 * pkin(2) * mrSges(4,1) + 0.2e1 * t61 * mrSges(5,2) + 0.2e1 * t36 * t26 + 0.2e1 * t56 * t27 + Ifges(4,2) + Ifges(3,3) + Ifges(5,3) + (t24 * t91 - t31 - t32) * t59 + (t23 * t91 - t33 - t34) * t57 + m(7) * (t23 ^ 2 + t24 ^ 2 + t36 ^ 2) + m(6) * (t35 * t51 ^ 2 + t56 ^ 2) + m(5) * (t61 ^ 2 + t62) + m(4) * (pkin(2) ^ 2 + t62) + 0.2e1 * (mrSges(5,1) + mrSges(4,3)) * qJ(3) - 0.2e1 * t51 * t67; t76 * t59 - t75 * t57 + (m(4) * pkin(7) + t78) * t58 + m(7) * (-t57 * t1 + t59 * t2) + m(6) * t64 + m(5) * t28; mrSges(5,2) - t35 * mrSges(7,3) - t67 + m(7) * (-t57 * t23 - t59 * t24) + m(5) * t61 + t51 * t18 + t70; m(4) + t71; -t60 * mrSges(5,2) + t41 + t75 * t59 + t76 * t57 + m(7) * (t59 * t1 + t57 * t2) + m(6) * (t59 * t3 + t57 * t4) + m(5) * t15; m(7) * (t59 * t23 - t57 * t24); 0; t71; t3 * mrSges(6,1) + t1 * mrSges(7,1) - t4 * mrSges(6,2) - t2 * mrSges(7,2) + t60 * t68 + (m(7) * t1 + t21) * pkin(5) + t66; t24 * mrSges(7,2) + t69 * t23 + (-mrSges(6,2) * t51 - t94) * t59 + (-mrSges(6,1) * t51 + mrSges(7,3) * pkin(5) + t77) * t57; t79 * t59 + (-mrSges(6,1) - t69) * t57; m(7) * t86 + t79 * t57 + t42 + t43; (0.2e1 * mrSges(7,1) + t87) * pkin(5) + t93; m(7) * t12 + t13; m(7) * t36 + t26; 0; 0; 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t5(1) t5(2) t5(4) t5(7) t5(11) t5(16); t5(2) t5(3) t5(5) t5(8) t5(12) t5(17); t5(4) t5(5) t5(6) t5(9) t5(13) t5(18); t5(7) t5(8) t5(9) t5(10) t5(14) t5(19); t5(11) t5(12) t5(13) t5(14) t5(15) t5(20); t5(16) t5(17) t5(18) t5(19) t5(20) t5(21);];
Mq  = res;
