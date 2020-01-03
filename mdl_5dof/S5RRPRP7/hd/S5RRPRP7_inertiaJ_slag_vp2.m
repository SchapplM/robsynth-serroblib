% Calculate joint inertia matrix for
% S5RRPRP7
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
% Datum: 2019-12-31 20:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP7_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP7_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP7_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP7_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP7_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP7_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:59:56
% EndTime: 2019-12-31 19:59:58
% DurationCPUTime: 0.66s
% Computational Cost: add. (686->175), mult. (1317->243), div. (0->0), fcn. (1273->6), ass. (0->70)
t92 = Ifges(6,2) + Ifges(5,3);
t57 = sin(qJ(4));
t59 = cos(qJ(4));
t91 = t57 ^ 2 + t59 ^ 2;
t90 = 0.2e1 * t91;
t89 = m(4) * pkin(2);
t60 = cos(qJ(2));
t75 = -qJ(3) - pkin(6);
t35 = t75 * t60;
t55 = sin(pkin(8));
t56 = cos(pkin(8));
t58 = sin(qJ(2));
t69 = t75 * t58;
t19 = -t55 * t35 - t56 * t69;
t88 = t19 ^ 2;
t87 = 0.2e1 * t19;
t47 = -t60 * pkin(2) - pkin(1);
t86 = 0.2e1 * t47;
t85 = t55 * pkin(2);
t84 = t56 * pkin(2);
t83 = Ifges(5,4) * t57;
t82 = Ifges(5,4) * t59;
t81 = Ifges(6,5) * t57;
t80 = Ifges(6,5) * t59;
t30 = t55 * t58 - t56 * t60;
t79 = Ifges(5,6) * t30;
t78 = Ifges(6,6) * t30;
t31 = t55 * t60 + t56 * t58;
t77 = t31 * t57;
t76 = t31 * t59;
t14 = -t30 * mrSges(5,2) - mrSges(5,3) * t77;
t17 = -mrSges(6,2) * t77 + t30 * mrSges(6,3);
t74 = t14 + t17;
t15 = t30 * mrSges(5,1) - mrSges(5,3) * t76;
t16 = -t30 * mrSges(6,1) + mrSges(6,2) * t76;
t73 = t15 - t16;
t13 = t30 * pkin(3) - t31 * pkin(7) + t47;
t21 = -t56 * t35 + t55 * t69;
t4 = t57 * t13 + t59 * t21;
t45 = pkin(7) + t85;
t72 = t91 * t45 ^ 2;
t70 = t58 ^ 2 + t60 ^ 2;
t46 = -pkin(3) - t84;
t67 = Ifges(6,6) * t77 + (Ifges(6,4) + Ifges(5,5)) * t76 + t92 * t30;
t34 = -t59 * mrSges(5,1) + t57 * mrSges(5,2);
t66 = t57 * mrSges(5,1) + t59 * mrSges(5,2);
t33 = -t59 * mrSges(6,1) - t57 * mrSges(6,3);
t65 = t57 * mrSges(6,1) - t59 * mrSges(6,3);
t64 = t59 * pkin(4) + t57 * qJ(5);
t63 = pkin(4) * t57 - qJ(5) * t59;
t3 = t59 * t13 - t57 * t21;
t50 = Ifges(6,4) * t57;
t49 = Ifges(5,5) * t57;
t48 = Ifges(5,6) * t59;
t39 = Ifges(5,1) * t57 + t82;
t38 = Ifges(6,1) * t57 - t80;
t37 = Ifges(5,2) * t59 + t83;
t36 = -Ifges(6,3) * t59 + t81;
t29 = t46 - t64;
t26 = t31 * mrSges(4,2);
t12 = t66 * t31;
t11 = t65 * t31;
t9 = Ifges(5,5) * t30 + (Ifges(5,1) * t59 - t83) * t31;
t8 = Ifges(6,4) * t30 + (Ifges(6,1) * t59 + t81) * t31;
t7 = t79 + (-Ifges(5,2) * t57 + t82) * t31;
t6 = t78 + (Ifges(6,3) * t57 + t80) * t31;
t5 = t63 * t31 + t19;
t2 = -t30 * pkin(4) - t3;
t1 = t30 * qJ(5) + t4;
t10 = [t58 * (Ifges(3,1) * t58 + Ifges(3,4) * t60) + t60 * (Ifges(3,4) * t58 + Ifges(3,2) * t60) - 0.2e1 * pkin(1) * (-t60 * mrSges(3,1) + t58 * mrSges(3,2)) + t26 * t86 + t12 * t87 + 0.2e1 * t5 * t11 + 0.2e1 * t4 * t14 + 0.2e1 * t3 * t15 + 0.2e1 * t2 * t16 + 0.2e1 * t1 * t17 + Ifges(2,3) + 0.2e1 * t70 * pkin(6) * mrSges(3,3) + (mrSges(4,1) * t86 - 0.2e1 * t21 * mrSges(4,3) + Ifges(4,2) * t30 + t67) * t30 + m(6) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) + m(5) * (t3 ^ 2 + t4 ^ 2 + t88) + m(4) * (t21 ^ 2 + t47 ^ 2 + t88) + m(3) * (t70 * pkin(6) ^ 2 + pkin(1) ^ 2) + (mrSges(4,3) * t87 + Ifges(4,1) * t31 - 0.2e1 * Ifges(4,4) * t30 + (t8 + t9) * t59 + (t6 - t7 - t79) * t57) * t31; -t21 * mrSges(4,2) + Ifges(3,5) * t58 + Ifges(3,6) * t60 + t29 * t11 + t46 * t12 + t5 * t33 + (t34 - mrSges(4,1)) * t19 + (-t58 * mrSges(3,1) - t60 * mrSges(3,2)) * pkin(6) + (t50 / 0.2e1 + t49 / 0.2e1 + t48 / 0.2e1 - Ifges(4,6) - mrSges(4,3) * t85) * t30 + (-t78 / 0.2e1 - t6 / 0.2e1 + t7 / 0.2e1 + t1 * mrSges(6,2) + t4 * mrSges(5,3) + t74 * t45) * t59 + (t8 / 0.2e1 + t9 / 0.2e1 + t2 * mrSges(6,2) - t3 * mrSges(5,3) - t73 * t45) * t57 + m(6) * (t29 * t5 + (t1 * t59 + t2 * t57) * t45) + m(5) * (t46 * t19 + (-t3 * t57 + t4 * t59) * t45) + (-t19 * t56 + t21 * t55) * t89 + (-mrSges(4,3) * t84 + Ifges(4,5) + (t38 / 0.2e1 + t39 / 0.2e1) * t59 + (t36 / 0.2e1 - t37 / 0.2e1) * t57) * t31; 0.2e1 * t29 * t33 + 0.2e1 * t46 * t34 + Ifges(3,3) + Ifges(4,3) + (-t36 + t37) * t59 + (t38 + t39) * t57 + m(6) * (t29 ^ 2 + t72) + m(5) * (t46 ^ 2 + t72) + (mrSges(6,2) + mrSges(5,3)) * t45 * t90 + (0.2e1 * t56 * mrSges(4,1) - 0.2e1 * t55 * mrSges(4,2) + (t55 ^ 2 + t56 ^ 2) * t89) * pkin(2); t30 * mrSges(4,1) + t26 + t73 * t59 + t74 * t57 + m(6) * (t57 * t1 - t59 * t2) + m(5) * (t59 * t3 + t57 * t4) + m(4) * t47; 0; m(4) + (m(5) / 0.2e1 + m(6) / 0.2e1) * t90; -Ifges(5,6) * t77 - t4 * mrSges(5,2) - pkin(4) * t16 + qJ(5) * t17 + m(6) * (-pkin(4) * t2 + qJ(5) * t1) + t1 * mrSges(6,3) - t2 * mrSges(6,1) + t3 * mrSges(5,1) + t67; -Ifges(6,6) * t59 + t48 + t49 + t50 - t63 * mrSges(6,2) + (-m(6) * t63 - t65 - t66) * t45; m(6) * t64 - t33 - t34; 0.2e1 * pkin(4) * mrSges(6,1) + 0.2e1 * qJ(5) * mrSges(6,3) + m(6) * (pkin(4) ^ 2 + qJ(5) ^ 2) + t92; m(6) * t2 + t16; (m(6) * t45 + mrSges(6,2)) * t57; -m(6) * t59; -m(6) * pkin(4) - mrSges(6,1); m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t10(1), t10(2), t10(4), t10(7), t10(11); t10(2), t10(3), t10(5), t10(8), t10(12); t10(4), t10(5), t10(6), t10(9), t10(13); t10(7), t10(8), t10(9), t10(10), t10(14); t10(11), t10(12), t10(13), t10(14), t10(15);];
Mq = res;
