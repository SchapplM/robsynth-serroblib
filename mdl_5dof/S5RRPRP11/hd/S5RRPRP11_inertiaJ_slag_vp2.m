% Calculate joint inertia matrix for
% S5RRPRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% Datum: 2019-12-31 20:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP11_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP11_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP11_inertiaJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP11_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP11_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP11_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:12:24
% EndTime: 2019-12-31 20:12:25
% DurationCPUTime: 0.52s
% Computational Cost: add. (418->164), mult. (763->211), div. (0->0), fcn. (524->4), ass. (0->63)
t83 = pkin(3) + pkin(6);
t82 = Ifges(6,2) + Ifges(5,3);
t49 = sin(qJ(2));
t51 = cos(qJ(2));
t81 = t49 ^ 2 + t51 ^ 2;
t48 = sin(qJ(4));
t50 = cos(qJ(4));
t68 = -t48 ^ 2 - t50 ^ 2;
t80 = (mrSges(6,2) + mrSges(5,3)) * t68;
t52 = -pkin(2) - pkin(7);
t79 = Ifges(5,4) * t48;
t78 = Ifges(5,4) * t50;
t77 = Ifges(6,5) * t48;
t76 = Ifges(6,5) * t50;
t75 = Ifges(5,6) * t48;
t74 = t48 * t51;
t73 = t50 * t51;
t65 = -t49 * qJ(3) - pkin(1);
t13 = t51 * t52 + t65;
t28 = t83 * t49;
t4 = t50 * t13 + t48 * t28;
t16 = t49 * mrSges(5,1) + mrSges(5,3) * t74;
t17 = -t49 * mrSges(6,1) - mrSges(6,2) * t74;
t72 = t16 - t17;
t18 = -t49 * mrSges(5,2) - mrSges(5,3) * t73;
t19 = -mrSges(6,2) * t73 + t49 * mrSges(6,3);
t71 = t18 + t19;
t70 = t68 * t52 ^ 2;
t69 = t81 * pkin(6) ^ 2;
t29 = t83 * t51;
t67 = Ifges(6,6) * t73 + t82 * t49;
t66 = -m(4) * pkin(2) + mrSges(4,2);
t1 = t49 * qJ(5) + t4;
t3 = -t48 * t13 + t50 * t28;
t2 = -t49 * pkin(4) - t3;
t62 = t48 * t1 - t50 * t2;
t61 = t50 * t3 + t48 * t4;
t60 = t50 * mrSges(5,1) - t48 * mrSges(5,2);
t59 = t50 * mrSges(6,1) + t48 * mrSges(6,3);
t58 = pkin(4) * t50 + qJ(5) * t48;
t57 = -0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * t68;
t56 = -Ifges(5,6) * t50 + (-Ifges(6,4) - Ifges(5,5)) * t48;
t55 = m(6) * t58 + t59 + t60;
t53 = qJ(3) ^ 2;
t38 = Ifges(6,4) * t50;
t37 = Ifges(5,5) * t50;
t35 = Ifges(6,6) * t48;
t27 = Ifges(5,1) * t50 - t79;
t26 = Ifges(6,1) * t50 + t77;
t25 = -Ifges(5,2) * t48 + t78;
t24 = Ifges(6,3) * t48 + t76;
t23 = t48 * mrSges(5,1) + t50 * mrSges(5,2);
t22 = t48 * mrSges(6,1) - t50 * mrSges(6,3);
t21 = -t51 * pkin(2) + t65;
t20 = t48 * pkin(4) - t50 * qJ(5) + qJ(3);
t12 = t60 * t51;
t11 = t59 * t51;
t9 = Ifges(5,5) * t49 + (-Ifges(5,1) * t48 - t78) * t51;
t8 = Ifges(6,4) * t49 + (-Ifges(6,1) * t48 + t76) * t51;
t7 = Ifges(5,6) * t49 + (-Ifges(5,2) * t50 - t79) * t51;
t6 = Ifges(6,6) * t49 + (Ifges(6,3) * t50 - t77) * t51;
t5 = t51 * t58 + t29;
t10 = [0.2e1 * t1 * t19 + 0.2e1 * t5 * t11 + 0.2e1 * t29 * t12 + 0.2e1 * t3 * t16 + 0.2e1 * t2 * t17 + 0.2e1 * t4 * t18 + Ifges(2,3) + m(6) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) + m(5) * (t29 ^ 2 + t3 ^ 2 + t4 ^ 2) + m(3) * (pkin(1) ^ 2 + t69) + m(4) * (t21 ^ 2 + t69) + (-0.2e1 * pkin(1) * mrSges(3,2) - 0.2e1 * t21 * mrSges(4,3) + (Ifges(3,1) + Ifges(4,2)) * t49 + t67) * t49 + (0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * t21 * mrSges(4,2) + (Ifges(3,2) + Ifges(4,3)) * t51 + (t6 - t7) * t50 + (-t8 - t9) * t48 + ((2 * Ifges(3,4)) + (2 * Ifges(4,6)) + t56) * t49) * t51 + 0.2e1 * (mrSges(4,1) + mrSges(3,3)) * pkin(6) * t81; qJ(3) * t12 + t20 * t11 + t5 * t22 + t29 * t23 + (t8 / 0.2e1 + t9 / 0.2e1 + t2 * mrSges(6,2) - t3 * mrSges(5,3) + t72 * t52) * t50 + (t6 / 0.2e1 - t7 / 0.2e1 - t1 * mrSges(6,2) - t4 * mrSges(5,3) + t71 * t52) * t48 + m(6) * (t20 * t5 + t52 * t62) + m(5) * (qJ(3) * t29 + t52 * t61) + (-t75 / 0.2e1 + t37 / 0.2e1 + t38 / 0.2e1 + t35 / 0.2e1 - Ifges(4,4) + Ifges(3,5) - pkin(2) * mrSges(4,1) + (-mrSges(3,1) + t66) * pkin(6)) * t49 + (qJ(3) * mrSges(4,1) - Ifges(4,5) + Ifges(3,6) + (t24 / 0.2e1 - t25 / 0.2e1) * t50 + (-t26 / 0.2e1 - t27 / 0.2e1) * t48 + (m(4) * qJ(3) - mrSges(3,2) + mrSges(4,3)) * pkin(6)) * t51; -0.2e1 * pkin(2) * mrSges(4,2) + 0.2e1 * t20 * t22 + Ifges(4,1) + Ifges(3,3) + (t26 + t27) * t50 + (t24 - t25) * t48 + 0.2e1 * (t23 + mrSges(4,3)) * qJ(3) + m(6) * (t20 ^ 2 - t70) + m(5) * (t53 - t70) + m(4) * (pkin(2) ^ 2 + t53) + 0.2e1 * t52 * t80; t72 * t50 + (m(4) * pkin(6) + mrSges(4,1)) * t49 + t71 * t48 + m(6) * t62 + m(5) * t61; t52 * t57 + t66 + t80; m(4) + t57; qJ(5) * t19 + m(6) * (-pkin(4) * t2 + qJ(5) * t1) + t1 * mrSges(6,3) - t4 * mrSges(5,2) + t3 * mrSges(5,1) - t2 * mrSges(6,1) - pkin(4) * t17 + t56 * t51 + t67; -mrSges(6,2) * t58 + t52 * t55 + t35 + t37 + t38 - t75; t55; 0.2e1 * pkin(4) * mrSges(6,1) + 0.2e1 * qJ(5) * mrSges(6,3) + m(6) * (pkin(4) ^ 2 + qJ(5) ^ 2) + t82; m(6) * t2 + t17; (-m(6) * t52 + mrSges(6,2)) * t50; -m(6) * t50; -m(6) * pkin(4) - mrSges(6,1); m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t10(1), t10(2), t10(4), t10(7), t10(11); t10(2), t10(3), t10(5), t10(8), t10(12); t10(4), t10(5), t10(6), t10(9), t10(13); t10(7), t10(8), t10(9), t10(10), t10(14); t10(11), t10(12), t10(13), t10(14), t10(15);];
Mq = res;
