% Calculate joint inertia matrix for
% S5RRPRP10
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
% Datum: 2019-12-31 20:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP10_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP10_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP10_inertiaJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP10_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP10_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP10_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:09:26
% EndTime: 2019-12-31 20:09:27
% DurationCPUTime: 0.53s
% Computational Cost: add. (410->171), mult. (748->216), div. (0->0), fcn. (532->4), ass. (0->63)
t81 = pkin(3) + pkin(6);
t80 = Ifges(5,3) + Ifges(6,3);
t50 = sin(qJ(2));
t52 = cos(qJ(2));
t79 = t50 ^ 2 + t52 ^ 2;
t78 = m(4) * pkin(6) + mrSges(4,1);
t77 = -2 * mrSges(6,3);
t76 = m(6) * pkin(4);
t53 = -pkin(2) - pkin(7);
t75 = (m(5) * t53);
t58 = -t50 * qJ(3) - pkin(1);
t13 = t52 * t53 + t58;
t29 = t81 * t50;
t49 = sin(qJ(4));
t51 = cos(qJ(4));
t4 = t51 * t13 + t49 * t29;
t74 = Ifges(5,4) * t49;
t73 = Ifges(5,4) * t51;
t72 = Ifges(6,4) * t49;
t71 = Ifges(6,4) * t51;
t70 = t49 * t52;
t69 = t51 * mrSges(5,1);
t68 = t51 * t52;
t67 = -Ifges(5,6) - Ifges(6,6);
t23 = t49 * mrSges(6,1) + t51 * mrSges(6,2);
t66 = t80 * t50;
t65 = t79 * pkin(6) ^ 2;
t30 = t81 * t52;
t64 = t49 ^ 2 + t51 ^ 2;
t63 = qJ(5) * t52;
t62 = -qJ(5) + t53;
t61 = -m(4) * pkin(2) + mrSges(4,2);
t60 = -mrSges(5,3) + t75;
t59 = t64 * mrSges(5,3);
t11 = mrSges(6,1) * t68 - mrSges(6,2) * t70;
t56 = t67 * t51 + (-Ifges(5,5) - Ifges(6,5)) * t49;
t54 = qJ(3) ^ 2;
t39 = Ifges(5,5) * t51;
t38 = Ifges(6,5) * t51;
t33 = t49 * pkin(4) + qJ(3);
t28 = Ifges(5,1) * t51 - t74;
t27 = Ifges(6,1) * t51 - t72;
t26 = -Ifges(5,2) * t49 + t73;
t25 = -Ifges(6,2) * t49 + t71;
t24 = t49 * mrSges(5,1) + t51 * mrSges(5,2);
t22 = -t52 * pkin(2) + t58;
t21 = t62 * t51;
t20 = t62 * t49;
t19 = -t50 * mrSges(5,2) - mrSges(5,3) * t68;
t18 = -t50 * mrSges(6,2) - mrSges(6,3) * t68;
t17 = t50 * mrSges(5,1) + mrSges(5,3) * t70;
t16 = t50 * mrSges(6,1) + mrSges(6,3) * t70;
t15 = t51 * t29;
t12 = (-t49 * mrSges(5,2) + t69) * t52;
t10 = pkin(4) * t68 + t30;
t8 = Ifges(5,5) * t50 + (-Ifges(5,1) * t49 - t73) * t52;
t7 = Ifges(6,5) * t50 + (-Ifges(6,1) * t49 - t71) * t52;
t6 = Ifges(5,6) * t50 + (-Ifges(5,2) * t51 - t74) * t52;
t5 = Ifges(6,6) * t50 + (-Ifges(6,2) * t51 - t72) * t52;
t3 = -t49 * t13 + t15;
t2 = -t51 * t63 + t4;
t1 = t50 * pkin(4) + t15 + (-t13 + t63) * t49;
t9 = [0.2e1 * t1 * t16 + 0.2e1 * t10 * t11 + 0.2e1 * t30 * t12 + 0.2e1 * t3 * t17 + 0.2e1 * t2 * t18 + 0.2e1 * t4 * t19 + Ifges(2,3) + m(6) * (t1 ^ 2 + t10 ^ 2 + t2 ^ 2) + m(5) * (t3 ^ 2 + t30 ^ 2 + t4 ^ 2) + m(3) * (pkin(1) ^ 2 + t65) + m(4) * (t22 ^ 2 + t65) + (-0.2e1 * pkin(1) * mrSges(3,2) - 0.2e1 * t22 * mrSges(4,3) + (Ifges(4,2) + Ifges(3,1)) * t50 + t66) * t50 + (0.2e1 * t22 * mrSges(4,2) + 0.2e1 * pkin(1) * mrSges(3,1) + (Ifges(4,3) + Ifges(3,2)) * t52 + (-t5 - t6) * t51 + (-t7 - t8) * t49 + ((2 * Ifges(3,4)) + (2 * Ifges(4,6)) + t56) * t50) * t52 + 0.2e1 * (mrSges(4,1) + mrSges(3,3)) * pkin(6) * t79; t10 * t23 + t33 * t11 + t21 * t16 + t20 * t18 + t30 * t24 + m(6) * (t21 * t1 + t33 * t10 + t20 * t2) + (t7 / 0.2e1 + t8 / 0.2e1 - t1 * mrSges(6,3) + t53 * t17 + t60 * t3) * t51 + (-t5 / 0.2e1 - t6 / 0.2e1 - t2 * mrSges(6,3) + t53 * t19 + t60 * t4) * t49 + (t38 / 0.2e1 + t39 / 0.2e1 - Ifges(4,4) + Ifges(3,5) - pkin(2) * mrSges(4,1) + (-Ifges(6,6) / 0.2e1 - Ifges(5,6) / 0.2e1) * t49 + (-mrSges(3,1) + t61) * pkin(6)) * t50 + (-Ifges(4,5) + Ifges(3,6) + (-t25 / 0.2e1 - t26 / 0.2e1) * t51 + (-t27 / 0.2e1 - t28 / 0.2e1) * t49 + (-mrSges(3,2) + mrSges(4,3)) * pkin(6)) * t52 + (m(5) * t30 + t78 * t52 + t12) * qJ(3); -0.2e1 * pkin(2) * mrSges(4,2) + 0.2e1 * t33 * t23 + Ifges(4,1) + Ifges(3,3) + (t21 * t77 + t27 + t28) * t51 + (t20 * t77 - t25 - t26) * t49 + m(5) * (t53 ^ 2 * t64 + t54) + m(6) * (t20 ^ 2 + t21 ^ 2 + t33 ^ 2) + m(4) * (pkin(2) ^ 2 + t54) + 0.2e1 * (t24 + mrSges(4,3)) * qJ(3) - 0.2e1 * t53 * t59; (t16 + t17) * t51 + t78 * t50 + (t18 + t19) * t49 + m(6) * (t51 * t1 + t49 * t2) + m(5) * (t51 * t3 + t49 * t4); -t59 + m(6) * (t49 * t20 + t51 * t21) + t61 + (-mrSges(6,3) + t75) * t64; m(4) + 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * t64; t3 * mrSges(5,1) + t1 * mrSges(6,1) - t4 * mrSges(5,2) - t2 * mrSges(6,2) + (m(6) * t1 + t16) * pkin(4) + t56 * t52 + t66; t53 * t69 + t21 * mrSges(6,1) - t20 * mrSges(6,2) + t38 + t39 + (m(6) * t21 - t51 * mrSges(6,3)) * pkin(4) + (-mrSges(5,2) * t53 + t67) * t49; (-mrSges(5,2) - mrSges(6,2)) * t49 + (mrSges(5,1) + mrSges(6,1) + t76) * t51; (0.2e1 * mrSges(6,1) + t76) * pkin(4) + t80; m(6) * t10 + t11; m(6) * t33 + t23; 0; 0; m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t9(1), t9(2), t9(4), t9(7), t9(11); t9(2), t9(3), t9(5), t9(8), t9(12); t9(4), t9(5), t9(6), t9(9), t9(13); t9(7), t9(8), t9(9), t9(10), t9(14); t9(11), t9(12), t9(13), t9(14), t9(15);];
Mq = res;
