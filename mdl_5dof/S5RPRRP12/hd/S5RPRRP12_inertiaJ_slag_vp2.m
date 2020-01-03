% Calculate joint inertia matrix for
% S5RPRRP12
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
% Datum: 2019-12-31 18:57
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP12_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP12_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP12_inertiaJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP12_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP12_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP12_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:56:16
% EndTime: 2019-12-31 18:56:18
% DurationCPUTime: 0.47s
% Computational Cost: add. (370->170), mult. (717->224), div. (0->0), fcn. (507->4), ass. (0->64)
t76 = Ifges(5,3) + Ifges(6,3);
t75 = -2 * mrSges(6,3);
t52 = (-pkin(1) - pkin(6));
t74 = -2 * t52;
t73 = 2 * qJ(2);
t72 = m(6) * pkin(4);
t48 = sin(qJ(4));
t71 = Ifges(5,4) * t48;
t50 = cos(qJ(4));
t70 = Ifges(5,4) * t50;
t69 = Ifges(6,4) * t48;
t68 = Ifges(6,4) * t50;
t51 = cos(qJ(3));
t67 = t48 * t51;
t49 = sin(qJ(3));
t66 = t49 * t52;
t65 = t50 * t51;
t22 = -t50 * mrSges(5,1) + t48 * mrSges(5,2);
t64 = mrSges(4,1) - t22;
t63 = -Ifges(5,6) - Ifges(6,6);
t62 = -qJ(5) - pkin(7);
t19 = t49 * pkin(3) - t51 * pkin(7) + qJ(2);
t4 = t48 * t19 + t50 * t66;
t10 = mrSges(6,1) * t67 + mrSges(6,2) * t65;
t61 = t48 ^ 2 + t50 ^ 2;
t44 = t49 ^ 2;
t46 = t51 ^ 2;
t60 = t46 + t44;
t59 = qJ(5) * t51;
t58 = t61 * mrSges(5,3);
t57 = t60 * mrSges(4,3);
t21 = -t50 * mrSges(6,1) + t48 * mrSges(6,2);
t56 = (Ifges(5,5) + Ifges(6,5)) * t65 + t76 * t49;
t55 = mrSges(5,1) * t48 + mrSges(5,2) * t50;
t53 = qJ(2) ^ 2;
t47 = t52 ^ 2;
t42 = Ifges(5,5) * t48;
t41 = Ifges(6,5) * t48;
t40 = Ifges(5,6) * t50;
t39 = Ifges(6,6) * t50;
t35 = t46 * t52;
t34 = t46 * t47;
t33 = -t50 * pkin(4) - pkin(3);
t27 = Ifges(5,1) * t48 + t70;
t26 = Ifges(6,1) * t48 + t68;
t25 = Ifges(5,2) * t50 + t71;
t24 = Ifges(6,2) * t50 + t69;
t23 = t62 * t50;
t20 = t62 * t48;
t18 = t49 * mrSges(5,1) - mrSges(5,3) * t65;
t17 = t49 * mrSges(6,1) - mrSges(6,3) * t65;
t16 = -t49 * mrSges(5,2) - mrSges(5,3) * t67;
t15 = -t49 * mrSges(6,2) - mrSges(6,3) * t67;
t14 = (pkin(4) * t48 - t52) * t51;
t13 = t50 * t19;
t11 = t55 * t51;
t8 = Ifges(5,5) * t49 + (Ifges(5,1) * t50 - t71) * t51;
t7 = Ifges(6,5) * t49 + (Ifges(6,1) * t50 - t69) * t51;
t6 = Ifges(5,6) * t49 + (-Ifges(5,2) * t48 + t70) * t51;
t5 = Ifges(6,6) * t49 + (-Ifges(6,2) * t48 + t68) * t51;
t3 = -t48 * t66 + t13;
t2 = -t48 * t59 + t4;
t1 = -t50 * t59 + t13 + (-t48 * t52 + pkin(4)) * t49;
t9 = [-(2 * pkin(1) * mrSges(3,2)) + (mrSges(3,3) * t73) + 0.2e1 * t1 * t17 + 0.2e1 * t14 * t10 + 0.2e1 * t2 * t15 + 0.2e1 * t4 * t16 + 0.2e1 * t3 * t18 + Ifges(3,1) + Ifges(2,3) + t57 * t74 + (mrSges(4,1) * t73 + Ifges(4,2) * t49 + t56) * t49 + m(6) * (t1 ^ 2 + t14 ^ 2 + t2 ^ 2) + m(5) * (t3 ^ 2 + t4 ^ 2 + t34) + (m(3) * (pkin(1) ^ 2 + t53)) + m(4) * (t44 * t47 + t34 + t53) + ((mrSges(4,2) * t73) + Ifges(4,1) * t51 - 0.2e1 * Ifges(4,4) * t49 + t11 * t74 + (t7 + t8) * t50 + (t63 * t49 - t5 - t6) * t48) * t51; -(m(3) * pkin(1)) + mrSges(3,2) + (-t10 - t11) * t51 - t57 + ((t15 + t16) * t50 + (-t17 - t18) * t48) * t49 + m(6) * (-t51 * t14 + (-t1 * t48 + t2 * t50) * t49) + m(5) * (t35 + (-t3 * t48 + t4 * t50) * t49) + m(4) * (t44 * t52 + t35); m(3) + m(4) * t60 + 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * (t61 * t44 + t46); t33 * t10 - pkin(3) * t11 + t20 * t17 + t14 * t21 - t23 * t15 + m(6) * (t20 * t1 + t33 * t14 - t23 * t2) + (t41 / 0.2e1 + t39 / 0.2e1 + t42 / 0.2e1 + t40 / 0.2e1 - Ifges(4,6) - (t52 * mrSges(4,2))) * t49 + (Ifges(4,5) + (m(5) * pkin(3) + t64) * t52) * t51 + (t5 / 0.2e1 + t6 / 0.2e1 + t2 * mrSges(6,3) + t4 * mrSges(5,3) + (t26 / 0.2e1 + t27 / 0.2e1) * t51 + (m(5) * t4 + t16) * pkin(7)) * t50 + (t7 / 0.2e1 + t8 / 0.2e1 - t1 * mrSges(6,3) - t3 * mrSges(5,3) + (-t24 / 0.2e1 - t25 / 0.2e1) * t51 + (-m(5) * t3 - t18) * pkin(7)) * t48; (-t21 + t64) * t51 + (t61 * mrSges(6,3) - mrSges(4,2) + t58) * t49 + m(5) * (t61 * t49 * pkin(7) + pkin(3) * t51) + m(6) * (-t33 * t51 + (-t20 * t48 - t23 * t50) * t49); -0.2e1 * pkin(3) * t22 + 0.2e1 * t33 * t21 + Ifges(4,3) + 0.2e1 * pkin(7) * t58 + m(6) * (t20 ^ 2 + t23 ^ 2 + t33 ^ 2) + m(5) * (t61 * pkin(7) ^ 2 + pkin(3) ^ 2) + (t23 * t75 + t24 + t25) * t50 + (t20 * t75 + t26 + t27) * t48; t3 * mrSges(5,1) + t1 * mrSges(6,1) - t4 * mrSges(5,2) - t2 * mrSges(6,2) + t63 * t67 + (m(6) * t1 + t17) * pkin(4) + t56; ((-mrSges(5,2) - mrSges(6,2)) * t50 + (-mrSges(5,1) - mrSges(6,1) - t72) * t48) * t49; t20 * mrSges(6,1) + t23 * mrSges(6,2) + t39 + t40 + t41 + t42 - t55 * pkin(7) + (m(6) * t20 - t48 * mrSges(6,3)) * pkin(4); (0.2e1 * mrSges(6,1) + t72) * pkin(4) + t76; m(6) * t14 + t10; -m(6) * t51; m(6) * t33 + t21; 0; m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t9(1), t9(2), t9(4), t9(7), t9(11); t9(2), t9(3), t9(5), t9(8), t9(12); t9(4), t9(5), t9(6), t9(9), t9(13); t9(7), t9(8), t9(9), t9(10), t9(14); t9(11), t9(12), t9(13), t9(14), t9(15);];
Mq = res;
