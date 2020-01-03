% Calculate joint inertia matrix for
% S5RRPPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5]';
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
% Datum: 2019-12-31 19:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR9_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR9_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPPR9_inertiaJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR9_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR9_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR9_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:40:42
% EndTime: 2019-12-31 19:40:43
% DurationCPUTime: 0.35s
% Computational Cost: add. (311->132), mult. (531->163), div. (0->0), fcn. (344->4), ass. (0->46)
t35 = sin(qJ(2));
t37 = cos(qJ(2));
t61 = t35 ^ 2 + t37 ^ 2;
t60 = pkin(6) - qJ(4);
t14 = t60 * t37;
t59 = t14 ^ 2;
t58 = -2 * mrSges(5,3);
t11 = -t37 * pkin(2) - t35 * qJ(3) - pkin(1);
t57 = -0.2e1 * t11;
t34 = sin(qJ(5));
t36 = cos(qJ(5));
t46 = t34 ^ 2 + t36 ^ 2;
t8 = m(6) * t46;
t56 = m(5) + t8;
t38 = -pkin(2) - pkin(3);
t55 = Ifges(6,4) * t34;
t54 = Ifges(6,4) * t36;
t53 = Ifges(6,5) * t34;
t52 = Ifges(6,6) * t36;
t51 = t34 * t37;
t50 = t36 * t37;
t49 = -mrSges(5,3) + mrSges(4,2);
t48 = Ifges(6,6) * t51 + Ifges(6,3) * t35;
t47 = t61 * pkin(6) ^ 2;
t45 = -m(4) * pkin(2) - mrSges(4,1);
t44 = t46 * mrSges(6,3);
t7 = t37 * pkin(3) - t11;
t13 = t60 * t35;
t3 = t35 * pkin(4) + t37 * pkin(7) + t7;
t1 = -t34 * t13 + t36 * t3;
t2 = t36 * t13 + t34 * t3;
t42 = -t34 * t1 + t36 * t2;
t41 = -t34 * mrSges(6,1) - t36 * mrSges(6,2);
t39 = qJ(3) ^ 2;
t33 = qJ(3) + pkin(4);
t28 = -pkin(7) + t38;
t21 = t35 * mrSges(5,1);
t17 = -Ifges(6,1) * t34 - t54;
t16 = -Ifges(6,2) * t36 - t55;
t12 = t36 * mrSges(6,1) - t34 * mrSges(6,2);
t10 = t35 * mrSges(6,1) + mrSges(6,3) * t50;
t9 = -t35 * mrSges(6,2) + mrSges(6,3) * t51;
t6 = t41 * t37;
t5 = Ifges(6,5) * t35 + (-Ifges(6,1) * t36 + t55) * t37;
t4 = Ifges(6,6) * t35 + (Ifges(6,2) * t34 - t54) * t37;
t15 = [0.2e1 * t1 * t10 + 0.2e1 * t14 * t6 + 0.2e1 * t2 * t9 + 0.2e1 * t7 * t21 + Ifges(2,3) + m(6) * (t1 ^ 2 + t2 ^ 2 + t59) + m(5) * (t13 ^ 2 + t7 ^ 2 + t59) + m(3) * (pkin(1) ^ 2 + t47) + m(4) * (t11 ^ 2 + t47) + (0.2e1 * pkin(1) * mrSges(3,1) + mrSges(4,1) * t57 - 0.2e1 * t7 * mrSges(5,2) + t14 * t58 + t34 * t4 - t36 * t5 + (Ifges(5,1) + Ifges(4,3) + Ifges(3,2)) * t37) * t37 + (-0.2e1 * pkin(1) * mrSges(3,2) + mrSges(4,3) * t57 + t13 * t58 + (Ifges(5,2) + Ifges(4,1) + Ifges(3,1)) * t35 + (-Ifges(6,5) * t36 + (2 * Ifges(3,4)) + (2 * Ifges(5,4)) - (2 * Ifges(4,5))) * t37 + t48) * t35 + 0.2e1 * (mrSges(4,2) + mrSges(3,3)) * pkin(6) * t61; t13 * mrSges(5,2) + t33 * t6 + (t12 + mrSges(5,1)) * t14 + (-t2 * mrSges(6,3) + t28 * t9 - t4 / 0.2e1) * t36 + (t1 * mrSges(6,3) - t28 * t10 - t5 / 0.2e1) * t34 + m(6) * (t33 * t14 + t42 * t28) + m(5) * (qJ(3) * t14 + t38 * t13) + (-pkin(2) * mrSges(4,2) - t38 * mrSges(5,3) + Ifges(5,6) + Ifges(4,4) + Ifges(3,5) - t53 / 0.2e1 - t52 / 0.2e1 + (-mrSges(3,1) + t45) * pkin(6)) * t35 + (t34 * t16 / 0.2e1 - t36 * t17 / 0.2e1 + Ifges(5,5) - Ifges(4,6) + Ifges(3,6) + t49 * qJ(3) + (m(4) * qJ(3) - mrSges(3,2) + mrSges(4,3)) * pkin(6)) * t37; 0.2e1 * pkin(2) * mrSges(4,1) + 0.2e1 * t38 * mrSges(5,2) + 0.2e1 * t33 * t12 - t36 * t16 - t34 * t17 + Ifges(4,2) + Ifges(3,3) + Ifges(5,3) + m(5) * (t38 ^ 2 + t39) + m(6) * (t46 * t28 ^ 2 + t33 ^ 2) + m(4) * (pkin(2) ^ 2 + t39) + 0.2e1 * (mrSges(5,1) + mrSges(4,3)) * qJ(3) - 0.2e1 * t28 * t44; -t34 * t10 + t36 * t9 + m(6) * t42 + m(5) * t13 + (m(4) * pkin(6) + t49) * t35; m(5) * t38 + t28 * t8 + mrSges(5,2) - t44 + t45; m(4) + t56; -t37 * mrSges(5,2) + t36 * t10 + t34 * t9 + t21 + m(6) * (t36 * t1 + t34 * t2) + m(5) * t7; 0; 0; t56; t1 * mrSges(6,1) - t2 * mrSges(6,2) - Ifges(6,5) * t50 + t48; t41 * t28 - t52 - t53; t41; t12; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t15(1), t15(2), t15(4), t15(7), t15(11); t15(2), t15(3), t15(5), t15(8), t15(12); t15(4), t15(5), t15(6), t15(9), t15(13); t15(7), t15(8), t15(9), t15(10), t15(14); t15(11), t15(12), t15(13), t15(14), t15(15);];
Mq = res;
