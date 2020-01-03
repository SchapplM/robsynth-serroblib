% Calculate joint inertia matrix for
% S5RPRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP11_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP11_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP11_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP11_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP11_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP11_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:53:28
% EndTime: 2019-12-31 18:53:29
% DurationCPUTime: 0.58s
% Computational Cost: add. (638->162), mult. (1237->222), div. (0->0), fcn. (1212->6), ass. (0->67)
t88 = Ifges(6,2) + Ifges(5,3);
t55 = sin(qJ(4));
t57 = cos(qJ(4));
t87 = t55 ^ 2 + t57 ^ 2;
t86 = 0.2e1 * t87;
t54 = cos(pkin(8));
t73 = pkin(6) + qJ(2);
t33 = t73 * t54;
t56 = sin(qJ(3));
t53 = sin(pkin(8));
t67 = t73 * t53;
t82 = cos(qJ(3));
t19 = t56 * t33 + t82 * t67;
t85 = t19 ^ 2;
t50 = t54 ^ 2;
t84 = 0.2e1 * t19;
t42 = -t54 * pkin(2) - pkin(1);
t83 = 0.2e1 * t42;
t81 = Ifges(5,4) * t55;
t80 = Ifges(5,4) * t57;
t79 = Ifges(6,5) * t55;
t78 = Ifges(6,5) * t57;
t30 = t56 * t53 - t54 * t82;
t77 = Ifges(5,6) * t30;
t76 = Ifges(6,6) * t30;
t31 = t53 * t82 + t56 * t54;
t75 = t31 * t55;
t74 = t31 * t57;
t14 = -t30 * mrSges(5,2) - mrSges(5,3) * t75;
t17 = -mrSges(6,2) * t75 + t30 * mrSges(6,3);
t72 = t14 + t17;
t15 = t30 * mrSges(5,1) - mrSges(5,3) * t74;
t16 = -t30 * mrSges(6,1) + mrSges(6,2) * t74;
t71 = t15 - t16;
t13 = t30 * pkin(3) - t31 * pkin(7) + t42;
t21 = t33 * t82 - t56 * t67;
t4 = t55 * t13 + t57 * t21;
t70 = t87 * pkin(7) ^ 2;
t69 = t53 ^ 2 + t50;
t66 = -t54 * mrSges(3,1) + t53 * mrSges(3,2);
t64 = Ifges(6,6) * t75 + (Ifges(6,4) + Ifges(5,5)) * t74 + t88 * t30;
t35 = -t57 * mrSges(5,1) + t55 * mrSges(5,2);
t63 = t55 * mrSges(5,1) + t57 * mrSges(5,2);
t34 = -t57 * mrSges(6,1) - t55 * mrSges(6,3);
t62 = t55 * mrSges(6,1) - t57 * mrSges(6,3);
t61 = t57 * pkin(4) + t55 * qJ(5);
t60 = pkin(4) * t55 - qJ(5) * t57;
t3 = t57 * t13 - t55 * t21;
t46 = Ifges(6,4) * t55;
t45 = Ifges(5,5) * t55;
t44 = Ifges(5,6) * t57;
t39 = Ifges(5,1) * t55 + t80;
t38 = Ifges(6,1) * t55 - t78;
t37 = Ifges(5,2) * t57 + t81;
t36 = -Ifges(6,3) * t57 + t79;
t32 = -pkin(3) - t61;
t26 = t31 * mrSges(4,2);
t12 = t63 * t31;
t11 = t62 * t31;
t9 = Ifges(5,5) * t30 + (Ifges(5,1) * t57 - t81) * t31;
t8 = Ifges(6,4) * t30 + (Ifges(6,1) * t57 + t79) * t31;
t7 = t77 + (-Ifges(5,2) * t55 + t80) * t31;
t6 = t76 + (Ifges(6,3) * t55 + t78) * t31;
t5 = t31 * t60 + t19;
t2 = -t30 * pkin(4) - t3;
t1 = t30 * qJ(5) + t4;
t10 = [-0.2e1 * pkin(1) * t66 + Ifges(3,2) * t50 + t26 * t83 + 0.2e1 * t5 * t11 + 0.2e1 * t4 * t14 + 0.2e1 * t3 * t15 + 0.2e1 * t2 * t16 + 0.2e1 * t1 * t17 + t12 * t84 + Ifges(2,3) + (Ifges(3,1) * t53 + 0.2e1 * Ifges(3,4) * t54) * t53 + 0.2e1 * t69 * qJ(2) * mrSges(3,3) + (mrSges(4,1) * t83 - 0.2e1 * t21 * mrSges(4,3) + Ifges(4,2) * t30 + t64) * t30 + m(6) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) + m(5) * (t3 ^ 2 + t4 ^ 2 + t85) + m(4) * (t21 ^ 2 + t42 ^ 2 + t85) + m(3) * (qJ(2) ^ 2 * t69 + pkin(1) ^ 2) + (mrSges(4,3) * t84 + Ifges(4,1) * t31 - 0.2e1 * Ifges(4,4) * t30 + (t8 + t9) * t57 + (t6 - t7 - t77) * t55) * t31; -m(3) * pkin(1) + t30 * mrSges(4,1) + t26 + t71 * t57 + t72 * t55 + m(6) * (t55 * t1 - t57 * t2) + m(5) * (t57 * t3 + t55 * t4) + m(4) * t42 + t66; m(3) + m(4) + (m(5) / 0.2e1 + m(6) / 0.2e1) * t86; -t21 * mrSges(4,2) - pkin(3) * t12 + t32 * t11 + t5 * t34 + (t46 / 0.2e1 + t45 / 0.2e1 + t44 / 0.2e1 - Ifges(4,6)) * t30 + (t35 - mrSges(4,1)) * t19 + (-t76 / 0.2e1 - t6 / 0.2e1 + t7 / 0.2e1 + t1 * mrSges(6,2) + t4 * mrSges(5,3) + t72 * pkin(7)) * t57 + (t8 / 0.2e1 + t9 / 0.2e1 + t2 * mrSges(6,2) - t3 * mrSges(5,3) - t71 * pkin(7)) * t55 + m(5) * (-pkin(3) * t19 + (-t3 * t55 + t4 * t57) * pkin(7)) + m(6) * (t32 * t5 + (t1 * t57 + t2 * t55) * pkin(7)) + (Ifges(4,5) + (t38 / 0.2e1 + t39 / 0.2e1) * t57 + (t36 / 0.2e1 - t37 / 0.2e1) * t55) * t31; 0; -0.2e1 * pkin(3) * t35 + 0.2e1 * t32 * t34 + Ifges(4,3) + (-t36 + t37) * t57 + (t39 + t38) * t55 + m(6) * (t32 ^ 2 + t70) + m(5) * (pkin(3) ^ 2 + t70) + (mrSges(6,2) + mrSges(5,3)) * pkin(7) * t86; -Ifges(5,6) * t75 - pkin(4) * t16 + m(6) * (-pkin(4) * t2 + qJ(5) * t1) + qJ(5) * t17 + t1 * mrSges(6,3) - t4 * mrSges(5,2) + t3 * mrSges(5,1) - t2 * mrSges(6,1) + t64; m(6) * t61 - t34 - t35; -Ifges(6,6) * t57 + t44 + t45 + t46 - t60 * mrSges(6,2) + (-m(6) * t60 - t62 - t63) * pkin(7); 0.2e1 * pkin(4) * mrSges(6,1) + 0.2e1 * qJ(5) * mrSges(6,3) + m(6) * (pkin(4) ^ 2 + qJ(5) ^ 2) + t88; m(6) * t2 + t16; -m(6) * t57; (m(6) * pkin(7) + mrSges(6,2)) * t55; -m(6) * pkin(4) - mrSges(6,1); m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t10(1), t10(2), t10(4), t10(7), t10(11); t10(2), t10(3), t10(5), t10(8), t10(12); t10(4), t10(5), t10(6), t10(9), t10(13); t10(7), t10(8), t10(9), t10(10), t10(14); t10(11), t10(12), t10(13), t10(14), t10(15);];
Mq = res;
