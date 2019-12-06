% Calculate joint inertia matrix for
% S5PRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRP7_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP7_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP7_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP7_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP7_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP7_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:54:08
% EndTime: 2019-12-05 16:54:11
% DurationCPUTime: 0.68s
% Computational Cost: add. (396->171), mult. (914->234), div. (0->0), fcn. (786->8), ass. (0->72)
t93 = Ifges(5,5) + Ifges(6,5);
t73 = Ifges(5,6) + Ifges(6,6);
t92 = 2 * pkin(7);
t91 = -2 * mrSges(6,3);
t55 = sin(qJ(4));
t58 = cos(qJ(4));
t90 = t93 * t55 + t73 * t58;
t27 = -t58 * mrSges(5,1) + t55 * mrSges(5,2);
t89 = -m(5) * pkin(3) - mrSges(4,1) + t27;
t26 = -t58 * mrSges(6,1) + t55 * mrSges(6,2);
t41 = -t58 * pkin(4) - pkin(3);
t88 = m(6) * t41 + t26;
t54 = cos(pkin(5));
t56 = sin(qJ(3));
t59 = cos(qJ(3));
t53 = sin(pkin(5));
t57 = sin(qJ(2));
t77 = t53 * t57;
t12 = -t54 * t59 + t56 * t77;
t87 = t12 ^ 2;
t86 = m(6) * pkin(4);
t84 = pkin(7) * t59;
t83 = Ifges(5,4) * t55;
t82 = Ifges(5,4) * t58;
t81 = Ifges(6,4) * t55;
t80 = Ifges(6,4) * t58;
t79 = t12 * t56;
t14 = t54 * t56 + t59 * t77;
t78 = t14 * t59;
t60 = cos(qJ(2));
t76 = t53 * t60;
t75 = t55 * t56;
t74 = t56 * t58;
t72 = Ifges(5,3) + Ifges(6,3);
t71 = -qJ(5) - pkin(8);
t15 = mrSges(6,1) * t75 + mrSges(6,2) * t74;
t69 = t93 * t74;
t24 = -t59 * pkin(3) - t56 * pkin(8) - pkin(2);
t7 = t55 * t24 + t58 * t84;
t66 = t55 ^ 2 + t58 ^ 2;
t65 = qJ(5) * t56;
t63 = mrSges(5,1) * t55 + mrSges(5,2) * t58;
t62 = pkin(7) ^ 2;
t52 = t59 ^ 2;
t50 = t56 ^ 2;
t48 = t53 ^ 2;
t47 = t50 * t62;
t39 = t48 * t60 ^ 2;
t33 = Ifges(5,1) * t55 + t82;
t32 = Ifges(6,1) * t55 + t80;
t31 = Ifges(5,2) * t58 + t83;
t30 = Ifges(6,2) * t58 + t81;
t29 = t71 * t58;
t28 = -t59 * mrSges(4,1) + t56 * mrSges(4,2);
t25 = t71 * t55;
t23 = (pkin(4) * t55 + pkin(7)) * t56;
t22 = -t59 * mrSges(5,1) - mrSges(5,3) * t74;
t21 = -t59 * mrSges(6,1) - mrSges(6,3) * t74;
t20 = t59 * mrSges(5,2) - mrSges(5,3) * t75;
t19 = t59 * mrSges(6,2) - mrSges(6,3) * t75;
t18 = t58 * t24;
t16 = t63 * t56;
t11 = -Ifges(5,5) * t59 + (Ifges(5,1) * t58 - t83) * t56;
t10 = -Ifges(6,5) * t59 + (Ifges(6,1) * t58 - t81) * t56;
t9 = -Ifges(5,6) * t59 + (-Ifges(5,2) * t55 + t82) * t56;
t8 = -Ifges(6,6) * t59 + (-Ifges(6,2) * t55 + t80) * t56;
t6 = -t55 * t84 + t18;
t5 = -t55 * t65 + t7;
t4 = t14 * t58 - t55 * t76;
t3 = -t14 * t55 - t58 * t76;
t2 = -t58 * t65 + t18 + (-pkin(7) * t55 - pkin(4)) * t59;
t1 = [m(2) + m(4) * (t14 ^ 2 + t39 + t87) + m(3) * (t48 * t57 ^ 2 + t54 ^ 2 + t39) + 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * (t3 ^ 2 + t4 ^ 2 + t87); mrSges(4,3) * t78 + (t19 + t20) * t4 + (t21 + t22) * t3 + (-t57 * mrSges(3,2) + (mrSges(3,1) - t28) * t60) * t53 + (t56 * mrSges(4,3) + t15 + t16) * t12 + m(5) * (pkin(7) * t79 + t6 * t3 + t7 * t4) + m(6) * (t23 * t12 + t2 * t3 + t5 * t4) + m(4) * (pkin(2) * t76 + (t78 + t79) * pkin(7)); -0.2e1 * pkin(2) * t28 + 0.2e1 * t23 * t15 + 0.2e1 * t5 * t19 + 0.2e1 * t2 * t21 + 0.2e1 * t7 * t20 + 0.2e1 * t6 * t22 + Ifges(3,3) + (t50 + t52) * mrSges(4,3) * t92 + m(6) * (t2 ^ 2 + t23 ^ 2 + t5 ^ 2) + m(5) * (t6 ^ 2 + t7 ^ 2 + t47) + m(4) * (pkin(2) ^ 2 + t52 * t62 + t47) + ((Ifges(4,2) + t72) * t59 - t69) * t59 + (Ifges(4,1) * t56 + 0.2e1 * Ifges(4,4) * t59 + t16 * t92 + (t10 + t11) * t58 + (t73 * t59 - t8 - t9) * t55) * t56; -t14 * mrSges(4,2) + m(6) * (t25 * t3 - t29 * t4) + (t88 + t89) * t12 + (m(5) * pkin(8) + mrSges(5,3) + mrSges(6,3)) * (-t3 * t55 + t4 * t58); t25 * t21 + t23 * t26 - t29 * t19 + t41 * t15 - pkin(3) * t16 + m(6) * (t25 * t2 + t41 * t23 - t29 * t5) + (t8 / 0.2e1 + t9 / 0.2e1 + t7 * mrSges(5,3) + t5 * mrSges(6,3) + (m(5) * t7 + t20) * pkin(8)) * t58 + (t10 / 0.2e1 + t11 / 0.2e1 - t6 * mrSges(5,3) - t2 * mrSges(6,3) + (-m(5) * t6 - t22) * pkin(8)) * t55 + (Ifges(4,5) + t89 * pkin(7) + (t32 / 0.2e1 + t33 / 0.2e1) * t58 + (-t30 / 0.2e1 - t31 / 0.2e1) * t55) * t56 + (Ifges(4,6) - pkin(7) * mrSges(4,2) - t90 / 0.2e1) * t59; -0.2e1 * pkin(3) * t27 + 0.2e1 * t41 * t26 + Ifges(4,3) + 0.2e1 * t66 * pkin(8) * mrSges(5,3) + m(6) * (t25 ^ 2 + t29 ^ 2 + t41 ^ 2) + m(5) * (t66 * pkin(8) ^ 2 + pkin(3) ^ 2) + (t29 * t91 + t30 + t31) * t58 + (t25 * t91 + t32 + t33) * t55; (-mrSges(5,2) - mrSges(6,2)) * t4 + (mrSges(5,1) + mrSges(6,1) + t86) * t3; t6 * mrSges(5,1) + t2 * mrSges(6,1) - t7 * mrSges(5,2) - t5 * mrSges(6,2) - t72 * t59 - t73 * t75 + (m(6) * t2 + t21) * pkin(4) + t69; t25 * mrSges(6,1) + t29 * mrSges(6,2) - t63 * pkin(8) + (m(6) * t25 - t55 * mrSges(6,3)) * pkin(4) + t90; (0.2e1 * mrSges(6,1) + t86) * pkin(4) + t72; m(6) * t12; m(6) * t23 + t15; t88; 0; m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
