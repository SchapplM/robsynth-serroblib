% Calculate joint inertia matrix for
% S5PPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRRR2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR2_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR2_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR2_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR2_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRR2_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:14:12
% EndTime: 2019-12-05 15:14:12
% DurationCPUTime: 0.24s
% Computational Cost: add. (224->77), mult. (525->122), div. (0->0), fcn. (501->8), ass. (0->34)
t28 = sin(qJ(4));
t31 = cos(qJ(4));
t18 = -t31 * mrSges(5,1) + t28 * mrSges(5,2);
t27 = sin(qJ(5));
t30 = cos(qJ(5));
t16 = -t27 * t28 + t30 * t31;
t17 = t27 * t31 + t30 * t28;
t4 = -t16 * mrSges(6,1) + t17 * mrSges(6,2);
t43 = t18 + t4;
t25 = sin(pkin(9));
t26 = cos(pkin(9));
t29 = sin(qJ(3));
t40 = cos(qJ(3));
t13 = t29 * t25 - t40 * t26;
t11 = t13 ^ 2;
t24 = t31 ^ 2;
t42 = m(6) * pkin(4);
t41 = -pkin(7) - pkin(6);
t39 = t28 ^ 2 + t24;
t15 = t40 * t25 + t29 * t26;
t2 = t17 * t15;
t3 = t16 * t15;
t38 = -t2 * mrSges(6,1) - t3 * mrSges(6,2);
t37 = t39 * mrSges(5,3);
t19 = t41 * t28;
t20 = t41 * t31;
t6 = t30 * t19 + t27 * t20;
t7 = t27 * t19 - t30 * t20;
t36 = t6 * mrSges(6,1) - t7 * mrSges(6,2) + Ifges(6,5) * t17 + Ifges(6,6) * t16;
t35 = -t28 * mrSges(5,1) - t31 * mrSges(5,2);
t34 = (t30 * mrSges(6,1) - t27 * mrSges(6,2)) * pkin(4);
t22 = -t31 * pkin(4) - pkin(3);
t12 = t15 ^ 2;
t1 = [m(2) + m(6) * (t2 ^ 2 + t3 ^ 2 + t11) + m(5) * (t39 * t12 + t11) + m(4) * (t12 + t11) + m(3) * (t25 ^ 2 + t26 ^ 2); m(6) * (-t16 * t2 + t17 * t3); m(3) + m(4) + m(5) * t39 + m(6) * (t16 ^ 2 + t17 ^ 2); (t3 * t16 + t2 * t17) * mrSges(6,3) + (-mrSges(4,2) + t37) * t15 + (-mrSges(4,1) + t43) * t13 + m(6) * (t22 * t13 - t6 * t2 + t7 * t3) + m(5) * (t39 * t15 * pkin(6) - pkin(3) * t13); m(6) * (t6 * t16 + t7 * t17); Ifges(5,2) * t24 - 0.2e1 * pkin(3) * t18 + 0.2e1 * t22 * t4 + Ifges(4,3) + (Ifges(5,1) * t28 + 0.2e1 * Ifges(5,4) * t31) * t28 + (-0.2e1 * t6 * mrSges(6,3) + Ifges(6,1) * t17) * t17 + 0.2e1 * pkin(6) * t37 + m(6) * (t22 ^ 2 + t6 ^ 2 + t7 ^ 2) + m(5) * (t39 * pkin(6) ^ 2 + pkin(3) ^ 2) + (0.2e1 * t7 * mrSges(6,3) + 0.2e1 * Ifges(6,4) * t17 + Ifges(6,2) * t16) * t16; t35 * t15 + (-t2 * t30 + t27 * t3) * t42 + t38; (t16 * t30 + t17 * t27) * t42 - t43; Ifges(5,5) * t28 + Ifges(5,6) * t31 + t35 * pkin(6) + (m(6) * (t27 * t7 + t30 * t6) + (t27 * t16 - t30 * t17) * mrSges(6,3)) * pkin(4) + t36; Ifges(5,3) + Ifges(6,3) + m(6) * (t27 ^ 2 + t30 ^ 2) * pkin(4) ^ 2 + 0.2e1 * t34; t38; -t4; t36; Ifges(6,3) + t34; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
