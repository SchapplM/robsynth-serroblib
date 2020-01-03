% Calculate joint inertia matrix for
% S5RRPRP5
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
% Datum: 2019-12-31 19:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP5_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP5_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP5_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP5_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP5_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:54:02
% EndTime: 2019-12-31 19:54:03
% DurationCPUTime: 0.39s
% Computational Cost: add. (638->128), mult. (1198->173), div. (0->0), fcn. (1217->6), ass. (0->42)
t60 = 2 * mrSges(6,3);
t37 = sin(pkin(8));
t38 = cos(pkin(8));
t40 = sin(qJ(2));
t41 = cos(qJ(2));
t27 = -t37 * t40 + t38 * t41;
t34 = -pkin(2) * t41 - pkin(1);
t20 = -pkin(3) * t27 + t34;
t59 = 0.2e1 * t20;
t58 = 0.2e1 * t27;
t57 = pkin(2) * t37;
t56 = cos(qJ(4));
t33 = pkin(2) * t38 + pkin(3);
t39 = sin(qJ(4));
t23 = t56 * t33 - t39 * t57;
t55 = t23 * mrSges(5,1);
t24 = t39 * t33 + t56 * t57;
t54 = t24 * mrSges(5,2);
t53 = mrSges(6,2) + mrSges(5,3);
t52 = Ifges(6,2) + Ifges(5,3);
t51 = -qJ(3) - pkin(6);
t47 = t51 * t40;
t48 = t51 * t41;
t19 = t37 * t47 - t38 * t48;
t50 = t40 ^ 2 + t41 ^ 2;
t18 = t37 * t48 + t38 * t47;
t28 = t37 * t41 + t38 * t40;
t44 = -t28 * pkin(7) + t18;
t9 = pkin(7) * t27 + t19;
t5 = t39 * t9 - t56 * t44;
t7 = t39 * t44 + t56 * t9;
t49 = t5 ^ 2 + t7 ^ 2;
t46 = -t27 * mrSges(4,1) + t28 * mrSges(4,2);
t16 = -t56 * t27 + t28 * t39;
t17 = t39 * t27 + t56 * t28;
t45 = (-mrSges(5,2) + mrSges(6,3)) * t7 + (-mrSges(5,1) - mrSges(6,1)) * t5 + (Ifges(6,4) + Ifges(5,5)) * t17 + (-Ifges(5,6) + Ifges(6,6)) * t16;
t22 = -pkin(4) - t23;
t21 = qJ(5) + t24;
t11 = t17 * mrSges(5,2);
t10 = t16 * mrSges(6,1);
t3 = pkin(4) * t16 - qJ(5) * t17 + t20;
t1 = [t40 * (Ifges(3,1) * t40 + Ifges(3,4) * t41) - 0.2e1 * pkin(1) * (-mrSges(3,1) * t41 + mrSges(3,2) * t40) + t41 * (Ifges(3,4) * t40 + Ifges(3,2) * t41) + 0.2e1 * t34 * t46 + Ifges(4,2) * t27 ^ 2 + 0.2e1 * t3 * t10 + t11 * t59 + Ifges(2,3) + t19 * mrSges(4,3) * t58 + 0.2e1 * t50 * pkin(6) * mrSges(3,3) + (-0.2e1 * t18 * mrSges(4,3) + Ifges(4,1) * t28 + Ifges(4,4) * t58) * t28 + (-0.2e1 * t3 * mrSges(6,3) + (Ifges(6,1) + Ifges(5,1)) * t17 + 0.2e1 * t53 * t5) * t17 + (mrSges(5,1) * t59 + (Ifges(6,3) + Ifges(5,2)) * t16 - 0.2e1 * t53 * t7 + 0.2e1 * (-Ifges(5,4) + Ifges(6,5)) * t17) * t16 + m(6) * (t3 ^ 2 + t49) + m(5) * (t20 ^ 2 + t49) + m(4) * (t18 ^ 2 + t19 ^ 2 + t34 ^ 2) + m(3) * (t50 * pkin(6) ^ 2 + pkin(1) ^ 2); t45 + (m(4) * (t18 * t38 + t19 * t37) + (t37 * t27 - t38 * t28) * mrSges(4,3)) * pkin(2) + m(6) * (t21 * t7 + t22 * t5) + m(5) * (-t23 * t5 + t24 * t7) + (-mrSges(3,1) * t40 - mrSges(3,2) * t41) * pkin(6) + (-t16 * t24 - t17 * t23) * mrSges(5,3) + (-t16 * t21 + t17 * t22) * mrSges(6,2) + Ifges(3,5) * t40 + Ifges(3,6) * t41 + Ifges(4,6) * t27 + Ifges(4,5) * t28 + t18 * mrSges(4,1) - t19 * mrSges(4,2); 0.2e1 * t55 - 0.2e1 * t22 * mrSges(6,1) - 0.2e1 * t54 + t21 * t60 + Ifges(3,3) + Ifges(4,3) + m(6) * (t21 ^ 2 + t22 ^ 2) + m(5) * (t23 ^ 2 + t24 ^ 2) + t52 + (0.2e1 * mrSges(4,1) * t38 - 0.2e1 * mrSges(4,2) * t37 + m(4) * (t37 ^ 2 + t38 ^ 2) * pkin(2)) * pkin(2); m(4) * t34 + m(5) * t20 + m(6) * t3 + t16 * mrSges(5,1) - t17 * mrSges(6,3) + t10 + t11 + t46; 0; m(4) + m(5) + m(6); m(6) * (-pkin(4) * t5 + qJ(5) * t7) + (-pkin(4) * t17 - qJ(5) * t16) * mrSges(6,2) + t45; m(6) * (-pkin(4) * t22 + qJ(5) * t21) - t54 + t55 + (t21 + qJ(5)) * mrSges(6,3) + (-t22 + pkin(4)) * mrSges(6,1) + t52; 0; 0.2e1 * pkin(4) * mrSges(6,1) + qJ(5) * t60 + m(6) * (pkin(4) ^ 2 + qJ(5) ^ 2) + t52; m(6) * t5 + t17 * mrSges(6,2); m(6) * t22 - mrSges(6,1); 0; -m(6) * pkin(4) - mrSges(6,1); m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
