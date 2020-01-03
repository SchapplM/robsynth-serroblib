% Calculate joint inertia matrix for
% S5RPRRP13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-31 19:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP13_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP13_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP13_inertiaJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP13_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP13_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP13_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:58:32
% EndTime: 2019-12-31 18:58:33
% DurationCPUTime: 0.47s
% Computational Cost: add. (374->160), mult. (734->211), div. (0->0), fcn. (500->4), ass. (0->66)
t83 = Ifges(6,2) + Ifges(5,3);
t47 = sin(qJ(4));
t49 = cos(qJ(4));
t65 = t47 ^ 2 + t49 ^ 2;
t51 = (-pkin(1) - pkin(6));
t82 = -2 * t51;
t81 = (mrSges(6,2) + mrSges(5,3)) * t65;
t80 = 2 * qJ(2);
t78 = Ifges(5,4) * t47;
t77 = Ifges(5,4) * t49;
t76 = Ifges(6,5) * t47;
t75 = Ifges(6,5) * t49;
t48 = sin(qJ(3));
t74 = Ifges(5,6) * t48;
t73 = Ifges(6,6) * t48;
t50 = cos(qJ(3));
t72 = t47 * t50;
t71 = t48 * t51;
t18 = t48 * pkin(3) - t50 * pkin(7) + qJ(2);
t70 = t49 * t18;
t69 = t49 * t50;
t21 = -t49 * mrSges(5,1) + t47 * mrSges(5,2);
t68 = mrSges(4,1) - t21;
t4 = t47 * t18 + t49 * t71;
t67 = t65 * pkin(7) * t48;
t66 = t65 * pkin(7) ^ 2;
t43 = t48 ^ 2;
t45 = t50 ^ 2;
t64 = t45 + t43;
t63 = t64 * mrSges(4,3);
t16 = -t48 * mrSges(6,1) + mrSges(6,2) * t69;
t61 = Ifges(6,6) * t72 + (Ifges(6,4) + Ifges(5,5)) * t69 + t83 * t48;
t1 = t48 * qJ(5) + t4;
t2 = -t70 + (t47 * t51 - pkin(4)) * t48;
t60 = t1 * t49 + t2 * t47;
t3 = -t47 * t71 + t70;
t59 = -t3 * t47 + t4 * t49;
t58 = t47 * mrSges(5,1) + t49 * mrSges(5,2);
t57 = t47 * mrSges(6,1) - t49 * mrSges(6,3);
t56 = -pkin(4) * t47 + qJ(5) * t49;
t14 = -t48 * mrSges(5,2) - mrSges(5,3) * t72;
t15 = t48 * mrSges(5,1) - mrSges(5,3) * t69;
t17 = -mrSges(6,2) * t72 + t48 * mrSges(6,3);
t55 = (t14 + t17) * t49 + (-t15 + t16) * t47;
t54 = m(6) * t56 - t57 - t58;
t52 = qJ(2) ^ 2;
t46 = t51 ^ 2;
t39 = Ifges(6,4) * t47;
t38 = Ifges(5,5) * t47;
t36 = Ifges(5,6) * t49;
t34 = t45 * t51;
t33 = t45 * t46;
t25 = Ifges(5,1) * t47 + t77;
t24 = Ifges(6,1) * t47 - t75;
t23 = Ifges(5,2) * t49 + t78;
t22 = -Ifges(6,3) * t49 + t76;
t20 = -t49 * mrSges(6,1) - t47 * mrSges(6,3);
t19 = -t49 * pkin(4) - t47 * qJ(5) - pkin(3);
t12 = t58 * t50;
t11 = t57 * t50;
t9 = Ifges(5,5) * t48 + (Ifges(5,1) * t49 - t78) * t50;
t8 = Ifges(6,4) * t48 + (Ifges(6,1) * t49 + t76) * t50;
t7 = t74 + (-Ifges(5,2) * t47 + t77) * t50;
t6 = t73 + (Ifges(6,3) * t47 + t75) * t50;
t5 = (-t51 - t56) * t50;
t10 = [-(2 * pkin(1) * mrSges(3,2)) + (mrSges(3,3) * t80) + 0.2e1 * t1 * t17 + 0.2e1 * t5 * t11 + 0.2e1 * t4 * t14 + 0.2e1 * t3 * t15 + 0.2e1 * t2 * t16 + Ifges(3,1) + Ifges(2,3) + t63 * t82 + (mrSges(4,1) * t80 + Ifges(4,2) * t48 + t61) * t48 + m(6) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) + m(5) * (t3 ^ 2 + t4 ^ 2 + t33) + (m(3) * (pkin(1) ^ 2 + t52)) + m(4) * (t43 * t46 + t33 + t52) + ((mrSges(4,2) * t80) + Ifges(4,1) * t50 - 0.2e1 * Ifges(4,4) * t48 + t12 * t82 + (t8 + t9) * t49 + (t6 - t7 - t74) * t47) * t50; -(m(3) * pkin(1)) + mrSges(3,2) + (-t11 - t12) * t50 - t63 + t55 * t48 + m(6) * (t60 * t48 - t50 * t5) + m(5) * (t59 * t48 + t34) + m(4) * (t43 * t51 + t34); m(3) + m(4) * t64 + 0.2e1 * (m(6) / 0.2e1 + m(5) / 0.2e1) * (t65 * t43 + t45); -pkin(3) * t12 + t5 * t20 + (m(6) * t5 + t11) * t19 + (t39 / 0.2e1 + t38 / 0.2e1 + t36 / 0.2e1 - Ifges(4,6) - (t51 * mrSges(4,2))) * t48 + (-t73 / 0.2e1 - t6 / 0.2e1 + t7 / 0.2e1 + t1 * mrSges(6,2) + t4 * mrSges(5,3)) * t49 + (t8 / 0.2e1 + t9 / 0.2e1 + t2 * mrSges(6,2) - t3 * mrSges(5,3)) * t47 + (m(5) * t59 + m(6) * t60 + t55) * pkin(7) + (Ifges(4,5) + (t24 / 0.2e1 + t25 / 0.2e1) * t49 + (t22 / 0.2e1 - t23 / 0.2e1) * t47 + (m(5) * pkin(3) + t68) * t51) * t50; (-t20 + t68) * t50 + m(6) * (-t19 * t50 + t67) + m(5) * (pkin(3) * t50 + t67) + (-mrSges(4,2) + t81) * t48; -0.2e1 * pkin(3) * t21 + 0.2e1 * t19 * t20 + Ifges(4,3) + (-t22 + t23) * t49 + (t24 + t25) * t47 + m(6) * (t19 ^ 2 + t66) + m(5) * (pkin(3) ^ 2 + t66) + 0.2e1 * pkin(7) * t81; -Ifges(5,6) * t72 + m(6) * (-pkin(4) * t2 + qJ(5) * t1) + t1 * mrSges(6,3) + qJ(5) * t17 - t4 * mrSges(5,2) + t3 * mrSges(5,1) - t2 * mrSges(6,1) - pkin(4) * t16 + t61; t54 * t48; t56 * mrSges(6,2) - Ifges(6,6) * t49 + t54 * pkin(7) + t36 + t38 + t39; 0.2e1 * pkin(4) * mrSges(6,1) + 0.2e1 * qJ(5) * mrSges(6,3) + m(6) * (pkin(4) ^ 2 + qJ(5) ^ 2) + t83; m(6) * t2 + t16; m(6) * t47 * t48; (m(6) * pkin(7) + mrSges(6,2)) * t47; -m(6) * pkin(4) - mrSges(6,1); m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t10(1), t10(2), t10(4), t10(7), t10(11); t10(2), t10(3), t10(5), t10(8), t10(12); t10(4), t10(5), t10(6), t10(9), t10(13); t10(7), t10(8), t10(9), t10(10), t10(14); t10(11), t10(12), t10(13), t10(14), t10(15);];
Mq = res;
