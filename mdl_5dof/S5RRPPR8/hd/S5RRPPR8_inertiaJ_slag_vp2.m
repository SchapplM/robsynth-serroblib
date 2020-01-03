% Calculate joint inertia matrix for
% S5RRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
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
% Datum: 2019-12-31 19:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR8_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR8_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR8_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR8_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR8_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR8_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:38:08
% EndTime: 2019-12-31 19:38:09
% DurationCPUTime: 0.42s
% Computational Cost: add. (492->129), mult. (870->177), div. (0->0), fcn. (818->6), ass. (0->43)
t39 = sin(qJ(2));
t41 = cos(qJ(2));
t57 = t39 ^ 2 + t41 ^ 2;
t56 = -m(4) * pkin(2) - mrSges(4,1);
t36 = sin(pkin(8));
t37 = cos(pkin(8));
t16 = -t39 * t36 - t41 * t37;
t17 = -t41 * t36 + t39 * t37;
t38 = sin(qJ(5));
t40 = cos(qJ(5));
t6 = t40 * t16 - t38 * t17;
t55 = 0.2e1 * t6;
t54 = 0.2e1 * t16;
t25 = -t41 * pkin(2) - t39 * qJ(3) - pkin(1);
t53 = -0.2e1 * t25;
t42 = -pkin(2) - pkin(3);
t23 = -t36 * qJ(3) + t37 * t42;
t22 = -pkin(4) + t23;
t24 = t37 * qJ(3) + t36 * t42;
t8 = t40 * t22 - t38 * t24;
t52 = t8 * mrSges(6,1);
t9 = t38 * t22 + t40 * t24;
t51 = t9 * mrSges(6,2);
t50 = pkin(6) - qJ(4);
t26 = t50 * t39;
t27 = t50 * t41;
t12 = t36 * t26 + t37 * t27;
t49 = t57 * pkin(6) ^ 2;
t7 = t38 * t16 + t40 * t17;
t48 = -t6 * mrSges(6,1) + t7 * mrSges(6,2);
t47 = -t16 * mrSges(5,1) + t17 * mrSges(5,2);
t11 = t37 * t26 - t36 * t27;
t14 = t41 * pkin(3) - t25;
t15 = -t38 * t36 + t40 * t37;
t18 = t40 * t36 + t38 * t37;
t45 = t15 * mrSges(6,1) - t18 * mrSges(6,2);
t3 = -t17 * pkin(7) + t11;
t4 = t16 * pkin(7) + t12;
t1 = t40 * t3 - t38 * t4;
t2 = t38 * t3 + t40 * t4;
t44 = t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t7 + Ifges(6,6) * t6;
t10 = -t16 * pkin(4) + t14;
t5 = [0.2e1 * t14 * t47 + Ifges(5,2) * t16 ^ 2 + Ifges(6,2) * t6 ^ 2 + 0.2e1 * t10 * t48 + Ifges(2,3) + t2 * mrSges(6,3) * t55 + t12 * mrSges(5,3) * t54 + (0.2e1 * pkin(1) * mrSges(3,1) + mrSges(4,1) * t53 + (Ifges(3,2) + Ifges(4,3)) * t41) * t41 + (-0.2e1 * pkin(1) * mrSges(3,2) + mrSges(4,3) * t53 + (Ifges(3,1) + Ifges(4,1)) * t39 + 0.2e1 * (Ifges(3,4) - Ifges(4,5)) * t41) * t39 + (-0.2e1 * t1 * mrSges(6,3) + Ifges(6,1) * t7 + Ifges(6,4) * t55) * t7 + (-0.2e1 * t11 * mrSges(5,3) + Ifges(5,1) * t17 + Ifges(5,4) * t54) * t17 + m(6) * (t1 ^ 2 + t10 ^ 2 + t2 ^ 2) + m(5) * (t11 ^ 2 + t12 ^ 2 + t14 ^ 2) + m(3) * (pkin(1) ^ 2 + t49) + m(4) * (t25 ^ 2 + t49) + 0.2e1 * (mrSges(4,2) + mrSges(3,3)) * pkin(6) * t57; -t11 * mrSges(5,1) + t12 * mrSges(5,2) - Ifges(5,5) * t17 - Ifges(5,6) * t16 + (t9 * t6 - t8 * t7) * mrSges(6,3) + (t24 * t16 - t23 * t17) * mrSges(5,3) + m(6) * (t8 * t1 + t9 * t2) + m(5) * (t23 * t11 + t24 * t12) + (qJ(3) * mrSges(4,2) + Ifges(3,6) - Ifges(4,6)) * t41 + (-pkin(2) * mrSges(4,2) + Ifges(4,4) + Ifges(3,5)) * t39 + ((m(4) * qJ(3) - mrSges(3,2) + mrSges(4,3)) * t41 + (-mrSges(3,1) + t56) * t39) * pkin(6) - t44; 0.2e1 * pkin(2) * mrSges(4,1) - 0.2e1 * t23 * mrSges(5,1) - 0.2e1 * t52 + 0.2e1 * t24 * mrSges(5,2) + 0.2e1 * t51 + 0.2e1 * qJ(3) * mrSges(4,3) + Ifges(4,2) + Ifges(3,3) + Ifges(5,3) + Ifges(6,3) + m(6) * (t8 ^ 2 + t9 ^ 2) + m(4) * (pkin(2) ^ 2 + qJ(3) ^ 2) + m(5) * (t23 ^ 2 + t24 ^ 2); (m(4) * pkin(6) + mrSges(4,2)) * t39 + (-t15 * t7 + t18 * t6) * mrSges(6,3) + (t36 * t16 - t37 * t17) * mrSges(5,3) + m(6) * (t15 * t1 + t18 * t2) + m(5) * (t37 * t11 + t36 * t12); -t37 * mrSges(5,1) + t36 * mrSges(5,2) + m(6) * (t15 * t8 + t18 * t9) + m(5) * (t37 * t23 + t36 * t24) - t45 + t56; m(4) + m(5) * (t36 ^ 2 + t37 ^ 2) + m(6) * (t15 ^ 2 + t18 ^ 2); m(5) * t14 + m(6) * t10 + t47 + t48; 0; 0; m(5) + m(6); t44; -Ifges(6,3) - t51 + t52; t45; 0; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t5(1), t5(2), t5(4), t5(7), t5(11); t5(2), t5(3), t5(5), t5(8), t5(12); t5(4), t5(5), t5(6), t5(9), t5(13); t5(7), t5(8), t5(9), t5(10), t5(14); t5(11), t5(12), t5(13), t5(14), t5(15);];
Mq = res;
