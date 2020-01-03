% Calculate joint inertia matrix for
% S5RPRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-31 18:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR14_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR14_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR14_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR14_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR14_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR14_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:34:31
% EndTime: 2019-12-31 18:34:32
% DurationCPUTime: 0.39s
% Computational Cost: add. (500->131), mult. (888->191), div. (0->0), fcn. (826->6), ass. (0->53)
t43 = cos(qJ(3));
t71 = t43 ^ 2;
t40 = sin(qJ(5));
t42 = cos(qJ(5));
t38 = sin(pkin(8));
t39 = cos(pkin(8));
t41 = sin(qJ(3));
t19 = t38 * t43 + t39 * t41;
t18 = t38 * t41 - t39 * t43;
t60 = t18 * t40;
t7 = -t19 * mrSges(6,2) + mrSges(6,3) * t60;
t59 = t18 * t42;
t8 = t19 * mrSges(6,1) + mrSges(6,3) * t59;
t70 = -t40 * t8 + t42 * t7;
t44 = -pkin(1) - pkin(6);
t54 = -qJ(4) + t44;
t21 = t54 * t41;
t50 = t54 * t43;
t9 = t38 * t21 - t39 * t50;
t69 = t9 ^ 2;
t17 = t18 ^ 2;
t68 = -2 * mrSges(5,3);
t66 = t40 / 0.2e1;
t65 = t18 * t9;
t62 = Ifges(6,4) * t40;
t61 = Ifges(6,4) * t42;
t58 = -Ifges(6,5) * t59 + Ifges(6,3) * t19;
t57 = Ifges(6,5) * t40 + Ifges(6,6) * t42;
t56 = t40 ^ 2 + t42 ^ 2;
t55 = t41 ^ 2 + t71;
t29 = t41 * pkin(3) + qJ(2);
t53 = m(4) * t55;
t27 = t38 * pkin(3) + pkin(7);
t52 = t56 * t27;
t51 = t55 * mrSges(4,3);
t11 = t39 * t21 + t38 * t50;
t6 = t19 * pkin(4) + t18 * pkin(7) + t29;
t1 = -t40 * t11 + t42 * t6;
t2 = t42 * t11 + t40 * t6;
t49 = -t1 * t40 + t2 * t42;
t48 = -mrSges(6,1) * t40 - mrSges(6,2) * t42;
t47 = t18 * t39 - t19 * t38;
t45 = qJ(2) ^ 2;
t28 = -t39 * pkin(3) - pkin(4);
t24 = Ifges(6,1) * t40 + t61;
t23 = Ifges(6,2) * t42 + t62;
t22 = -t42 * mrSges(6,1) + t40 * mrSges(6,2);
t16 = t19 ^ 2;
t13 = t18 * mrSges(5,2);
t5 = t48 * t18;
t4 = Ifges(6,5) * t19 + (-Ifges(6,1) * t42 + t62) * t18;
t3 = Ifges(6,6) * t19 + (Ifges(6,2) * t40 - t61) * t18;
t10 = [Ifges(3,1) + Ifges(2,3) + 0.2e1 * t2 * t7 + 0.2e1 * t1 * t8 + 0.2e1 * t9 * t5 + Ifges(4,1) * t71 - 0.2e1 * t29 * t13 - (2 * pkin(1) * mrSges(3,2)) + (0.2e1 * t29 * mrSges(5,1) + Ifges(5,2) * t19 + t11 * t68 + t58) * t19 - 0.2e1 * t44 * t51 + (t9 * t68 + Ifges(5,1) * t18 + t40 * t3 - t42 * t4 + (Ifges(6,6) * t40 + (2 * Ifges(5,4))) * t19) * t18 + m(6) * (t1 ^ 2 + t2 ^ 2 + t69) + m(5) * (t11 ^ 2 + t29 ^ 2 + t69) + m(4) * (t55 * t44 ^ 2 + t45) + m(3) * ((pkin(1) ^ 2) + t45) + (-0.2e1 * Ifges(4,4) * t43 + Ifges(4,2) * t41) * t41 + 0.2e1 * (t41 * mrSges(4,1) + t43 * mrSges(4,2) + mrSges(3,3)) * qJ(2); -m(3) * pkin(1) - t17 * mrSges(5,3) + t18 * t5 + mrSges(3,2) - t51 + (-mrSges(5,3) * t19 + t70) * t19 + m(6) * (t49 * t19 + t65) + m(5) * (t19 * t11 + t65) + t44 * t53; m(3) + t53 + m(5) * (t16 + t17) + m(6) * (t56 * t16 + t17); t9 * t22 + t4 * t66 + t42 * t3 / 0.2e1 - t11 * mrSges(5,2) - t9 * mrSges(5,1) + (t44 * mrSges(4,1) + Ifges(4,5)) * t43 + (-t44 * mrSges(4,2) - Ifges(4,6)) * t41 + t49 * mrSges(6,3) + (-Ifges(5,5) - t42 * t24 / 0.2e1 + t23 * t66) * t18 + (m(5) * (t11 * t38 - t39 * t9) + t47 * mrSges(5,3)) * pkin(3) + (m(6) * t9 + t5) * t28 + (m(6) * t49 + t70) * t27 + (-Ifges(5,6) + t57 / 0.2e1) * t19; t43 * mrSges(4,1) - t41 * mrSges(4,2) + (-mrSges(5,1) + t22) * t18 + (t56 * mrSges(6,3) - mrSges(5,2)) * t19 + m(6) * (t28 * t18 + t19 * t52) - m(5) * t47 * pkin(3); 0.2e1 * t28 * t22 + t42 * t23 + t40 * t24 + Ifges(4,3) + Ifges(5,3) + m(6) * (t56 * t27 ^ 2 + t28 ^ 2) + m(5) * (t38 ^ 2 + t39 ^ 2) * pkin(3) ^ 2 + 0.2e1 * (t39 * mrSges(5,1) - t38 * mrSges(5,2)) * pkin(3) + 0.2e1 * mrSges(6,3) * t52; t19 * mrSges(5,1) + t40 * t7 + t42 * t8 - t13 + m(6) * (t42 * t1 + t40 * t2) + m(5) * t29; 0; 0; m(6) * t56 + m(5); t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,6) * t60 + t58; t48 * t19; t48 * t27 + t57; -t22; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t10(1), t10(2), t10(4), t10(7), t10(11); t10(2), t10(3), t10(5), t10(8), t10(12); t10(4), t10(5), t10(6), t10(9), t10(13); t10(7), t10(8), t10(9), t10(10), t10(14); t10(11), t10(12), t10(13), t10(14), t10(15);];
Mq = res;
