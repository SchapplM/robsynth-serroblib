% Calculate joint inertia matrix for
% S5RPPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-31 18:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR12_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR12_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR12_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR12_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR12_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR12_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:06:57
% EndTime: 2019-12-31 18:06:58
% DurationCPUTime: 0.32s
% Computational Cost: add. (454->114), mult. (808->162), div. (0->0), fcn. (767->6), ass. (0->52)
t40 = sin(qJ(5));
t42 = cos(qJ(5));
t37 = sin(pkin(8));
t38 = cos(pkin(8));
t41 = sin(qJ(4));
t56 = t41 * t38;
t61 = cos(qJ(4));
t19 = t37 * t61 + t56;
t49 = t61 * t38;
t18 = t41 * t37 - t49;
t58 = t18 * t40;
t7 = -mrSges(6,2) * t19 + mrSges(6,3) * t58;
t57 = t18 * t42;
t8 = mrSges(6,1) * t19 + mrSges(6,3) * t57;
t71 = -t40 * t8 + t42 * t7;
t22 = -mrSges(6,1) * t42 + mrSges(6,2) * t40;
t70 = -m(6) * pkin(4) - mrSges(5,1) + t22;
t39 = -pkin(1) - qJ(3);
t62 = -pkin(6) + t39;
t21 = t62 * t37;
t9 = t41 * t21 - t49 * t62;
t69 = t9 ^ 2;
t17 = t18 ^ 2;
t33 = t38 ^ 2;
t68 = -2 * mrSges(5,3);
t66 = t40 / 0.2e1;
t65 = t18 * t9;
t60 = Ifges(6,4) * t40;
t59 = Ifges(6,4) * t42;
t55 = -Ifges(6,5) * t57 + Ifges(6,3) * t19;
t54 = t37 * mrSges(4,1) + t38 * mrSges(4,2);
t53 = Ifges(6,5) * t40 + Ifges(6,6) * t42;
t52 = t37 ^ 2 + t33;
t51 = t40 ^ 2 + t42 ^ 2;
t25 = t37 * pkin(3) + qJ(2);
t50 = m(4) * t52;
t48 = t52 * mrSges(4,3);
t11 = t21 * t61 + t56 * t62;
t5 = pkin(4) * t19 + pkin(7) * t18 + t25;
t1 = -t11 * t40 + t42 * t5;
t2 = t11 * t42 + t40 * t5;
t46 = -t1 * t40 + t2 * t42;
t45 = -mrSges(6,1) * t40 - mrSges(6,2) * t42;
t43 = qJ(2) ^ 2;
t24 = Ifges(6,1) * t40 + t59;
t23 = Ifges(6,2) * t42 + t60;
t16 = t19 ^ 2;
t14 = t18 * mrSges(5,2);
t6 = t45 * t18;
t4 = Ifges(6,5) * t19 + (-Ifges(6,1) * t42 + t60) * t18;
t3 = Ifges(6,6) * t19 + (Ifges(6,2) * t40 - t59) * t18;
t10 = [Ifges(3,1) + Ifges(2,3) + 0.2e1 * t9 * t6 + 0.2e1 * t2 * t7 + 0.2e1 * t1 * t8 - 0.2e1 * t25 * t14 + Ifges(4,1) * t33 - (2 * pkin(1) * mrSges(3,2)) + (-0.2e1 * Ifges(4,4) * t38 + Ifges(4,2) * t37) * t37 + (0.2e1 * t25 * mrSges(5,1) + Ifges(5,2) * t19 + t11 * t68 + t55) * t19 - 0.2e1 * t39 * t48 + (t9 * t68 + Ifges(5,1) * t18 + t40 * t3 - t42 * t4 + (Ifges(6,6) * t40 + (2 * Ifges(5,4))) * t19) * t18 + m(6) * (t1 ^ 2 + t2 ^ 2 + t69) + m(5) * (t11 ^ 2 + t25 ^ 2 + t69) + m(4) * (t39 ^ 2 * t52 + t43) + m(3) * ((pkin(1) ^ 2) + t43) + 0.2e1 * (mrSges(3,3) + t54) * qJ(2); -m(3) * pkin(1) - t17 * mrSges(5,3) + t18 * t6 + mrSges(3,2) - t48 + (-mrSges(5,3) * t19 + t71) * t19 + m(6) * (t19 * t46 + t65) + m(5) * (t11 * t19 + t65) + t39 * t50; m(3) + t50 + m(5) * (t16 + t17) + m(6) * (t16 * t51 + t17); m(4) * qJ(2) + t19 * mrSges(5,1) + t40 * t7 + t42 * t8 - t14 + m(6) * (t1 * t42 + t2 * t40) + m(5) * t25 + t54; 0; m(6) * t51 + m(4) + m(5); -pkin(4) * t6 + t4 * t66 + t42 * t3 / 0.2e1 - t11 * mrSges(5,2) + t46 * mrSges(6,3) + (t23 * t66 - t42 * t24 / 0.2e1 - Ifges(5,5)) * t18 + (t53 / 0.2e1 - Ifges(5,6)) * t19 + t70 * t9 + (m(6) * t46 + t71) * pkin(7); t70 * t18 + (-mrSges(5,2) + (m(6) * pkin(7) + mrSges(6,3)) * t51) * t19; 0; Ifges(5,3) - 0.2e1 * pkin(4) * t22 + t42 * t23 + t40 * t24 + m(6) * (pkin(7) ^ 2 * t51 + pkin(4) ^ 2) + 0.2e1 * t51 * pkin(7) * mrSges(6,3); mrSges(6,1) * t1 - mrSges(6,2) * t2 + Ifges(6,6) * t58 + t55; t45 * t19; -t22; pkin(7) * t45 + t53; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t10(1), t10(2), t10(4), t10(7), t10(11); t10(2), t10(3), t10(5), t10(8), t10(12); t10(4), t10(5), t10(6), t10(9), t10(13); t10(7), t10(8), t10(9), t10(10), t10(14); t10(11), t10(12), t10(13), t10(14), t10(15);];
Mq = res;
