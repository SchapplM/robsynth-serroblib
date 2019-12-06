% Calculate joint inertia matrix for
% S5PRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRR3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR3_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR3_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR3_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR3_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR3_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:46:39
% EndTime: 2019-12-05 15:46:40
% DurationCPUTime: 0.25s
% Computational Cost: add. (276->87), mult. (597->139), div. (0->0), fcn. (559->8), ass. (0->36)
t31 = sin(qJ(4));
t34 = cos(qJ(4));
t21 = -t34 * mrSges(5,1) + t31 * mrSges(5,2);
t30 = sin(qJ(5));
t33 = cos(qJ(5));
t18 = -t30 * t31 + t33 * t34;
t19 = t30 * t34 + t33 * t31;
t7 = -t18 * mrSges(6,1) + t19 * mrSges(6,2);
t46 = t21 + t7;
t28 = sin(pkin(9));
t29 = cos(pkin(9));
t32 = sin(qJ(2));
t35 = cos(qJ(2));
t15 = t28 * t32 - t29 * t35;
t13 = t15 ^ 2;
t27 = t34 ^ 2;
t45 = m(6) * pkin(4);
t24 = t28 * pkin(2) + pkin(6);
t44 = pkin(7) + t24;
t43 = t31 ^ 2 + t27;
t17 = t28 * t35 + t29 * t32;
t2 = t19 * t17;
t3 = t18 * t17;
t42 = -t2 * mrSges(6,1) - t3 * mrSges(6,2);
t25 = -t29 * pkin(2) - pkin(3);
t41 = t43 * t24;
t40 = -t31 * mrSges(5,1) - t34 * mrSges(5,2);
t8 = t44 * t31;
t9 = t44 * t34;
t5 = -t30 * t9 - t33 * t8;
t6 = -t30 * t8 + t33 * t9;
t39 = t5 * mrSges(6,1) - t6 * mrSges(6,2) + Ifges(6,5) * t19 + Ifges(6,6) * t18;
t38 = (t33 * mrSges(6,1) - t30 * mrSges(6,2)) * pkin(4);
t20 = -t34 * pkin(4) + t25;
t14 = t17 ^ 2;
t1 = [m(2) + m(6) * (t2 ^ 2 + t3 ^ 2 + t13) + m(5) * (t43 * t14 + t13) + m(4) * (t14 + t13) + m(3) * (t32 ^ 2 + t35 ^ 2); t35 * mrSges(3,1) - t32 * mrSges(3,2) + (t3 * t18 + t2 * t19) * mrSges(6,3) + (t43 * mrSges(5,3) - mrSges(4,2)) * t17 + (-mrSges(4,1) + t46) * t15 + m(6) * (t20 * t15 - t5 * t2 + t6 * t3) + m(5) * (t25 * t15 + t17 * t41) + m(4) * (-t15 * t29 + t17 * t28) * pkin(2); Ifges(5,2) * t27 + 0.2e1 * t20 * t7 + 0.2e1 * t25 * t21 + Ifges(3,3) + Ifges(4,3) + (Ifges(5,1) * t31 + 0.2e1 * Ifges(5,4) * t34) * t31 + (-0.2e1 * t5 * mrSges(6,3) + Ifges(6,1) * t19) * t19 + (0.2e1 * t6 * mrSges(6,3) + 0.2e1 * Ifges(6,4) * t19 + Ifges(6,2) * t18) * t18 + m(6) * (t20 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(5) * (t43 * t24 ^ 2 + t25 ^ 2) + m(4) * (t28 ^ 2 + t29 ^ 2) * pkin(2) ^ 2 + 0.2e1 * (t29 * mrSges(4,1) - t28 * mrSges(4,2)) * pkin(2) + 0.2e1 * mrSges(5,3) * t41; m(6) * (-t18 * t2 + t19 * t3); m(6) * (t18 * t5 + t19 * t6); m(4) + m(5) * t43 + m(6) * (t18 ^ 2 + t19 ^ 2); t40 * t17 + (-t2 * t33 + t3 * t30) * t45 + t42; Ifges(5,5) * t31 + Ifges(5,6) * t34 + t40 * t24 + (m(6) * (t30 * t6 + t33 * t5) + (t30 * t18 - t33 * t19) * mrSges(6,3)) * pkin(4) + t39; (t18 * t33 + t19 * t30) * t45 - t46; Ifges(5,3) + Ifges(6,3) + m(6) * (t30 ^ 2 + t33 ^ 2) * pkin(4) ^ 2 + 0.2e1 * t38; t42; t39; -t7; Ifges(6,3) + t38; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
