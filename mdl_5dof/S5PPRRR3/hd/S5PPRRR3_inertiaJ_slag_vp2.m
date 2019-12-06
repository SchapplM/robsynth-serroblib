% Calculate joint inertia matrix for
% S5PPRRR3
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
% Datum: 2019-12-05 15:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRRR3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR3_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR3_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR3_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR3_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRR3_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:16:23
% EndTime: 2019-12-05 15:16:24
% DurationCPUTime: 0.25s
% Computational Cost: add. (245->97), mult. (628->153), div. (0->0), fcn. (583->8), ass. (0->43)
t34 = cos(qJ(4));
t26 = t34 ^ 2;
t50 = m(6) * pkin(4);
t49 = -pkin(7) - pkin(6);
t28 = sin(pkin(9));
t32 = sin(qJ(3));
t48 = t28 * t32;
t35 = cos(qJ(3));
t47 = t28 * t35;
t31 = sin(qJ(4));
t46 = t31 ^ 2 + t26;
t17 = -t34 * mrSges(5,1) + t31 * mrSges(5,2);
t30 = sin(qJ(5));
t33 = cos(qJ(5));
t15 = -t30 * t31 + t33 * t34;
t16 = t30 * t34 + t33 * t31;
t4 = -t15 * mrSges(6,1) + t16 * mrSges(6,2);
t45 = mrSges(4,1) - t17 - t4;
t29 = cos(pkin(9));
t11 = -t29 * t34 - t31 * t47;
t12 = -t29 * t31 + t34 * t47;
t2 = t33 * t11 - t30 * t12;
t3 = t30 * t11 + t33 * t12;
t44 = t2 * mrSges(6,1) - t3 * mrSges(6,2);
t10 = t15 * t32;
t9 = t16 * t32;
t43 = -t9 * mrSges(6,1) - t10 * mrSges(6,2);
t42 = t46 * mrSges(5,3);
t41 = -t31 * mrSges(5,1) - t34 * mrSges(5,2);
t40 = -t11 * t31 + t12 * t34;
t18 = t49 * t31;
t19 = t49 * t34;
t6 = t33 * t18 + t30 * t19;
t7 = t30 * t18 - t33 * t19;
t39 = t6 * mrSges(6,1) - t7 * mrSges(6,2) + Ifges(6,5) * t16 + Ifges(6,6) * t15;
t38 = (t33 * mrSges(6,1) - t30 * mrSges(6,2)) * pkin(4);
t27 = t35 ^ 2;
t25 = t32 ^ 2;
t23 = t29 ^ 2;
t22 = t28 ^ 2;
t21 = -t34 * pkin(4) - pkin(3);
t20 = t25 * t22;
t1 = [m(2) + m(6) * (t2 ^ 2 + t3 ^ 2 + t20) + m(5) * (t11 ^ 2 + t12 ^ 2 + t20) + m(4) * (t27 * t22 + t20 + t23) + m(3) * (t22 + t23); m(5) * (t40 - t47) * t32 + (t10 * t3 - t9 * t2 - t47 * t32) * m(6); m(3) + m(4) * (t25 + t27) + m(5) * (t46 * t25 + t27) + m(6) * (t10 ^ 2 + t9 ^ 2 + t27); (t3 * t15 - t2 * t16) * mrSges(6,3) + t40 * mrSges(5,3) + (-t35 * mrSges(4,2) - t45 * t32) * t28 + m(6) * (t6 * t2 + t21 * t48 + t7 * t3) + m(5) * (-pkin(3) * t48 + t40 * pkin(6)); (t10 * t15 + t9 * t16) * mrSges(6,3) + t45 * t35 + (-mrSges(4,2) + t42) * t32 + m(5) * (t46 * t32 * pkin(6) + pkin(3) * t35) + m(6) * (t7 * t10 - t21 * t35 - t6 * t9); Ifges(5,2) * t26 - 0.2e1 * pkin(3) * t17 + 0.2e1 * t21 * t4 + Ifges(4,3) + (Ifges(5,1) * t31 + 0.2e1 * Ifges(5,4) * t34) * t31 + (-0.2e1 * t6 * mrSges(6,3) + Ifges(6,1) * t16) * t16 + 0.2e1 * pkin(6) * t42 + m(6) * (t21 ^ 2 + t6 ^ 2 + t7 ^ 2) + m(5) * (t46 * pkin(6) ^ 2 + pkin(3) ^ 2) + (0.2e1 * t7 * mrSges(6,3) + 0.2e1 * Ifges(6,4) * t16 + Ifges(6,2) * t15) * t15; t11 * mrSges(5,1) - t12 * mrSges(5,2) + (t2 * t33 + t3 * t30) * t50 + t44; t41 * t32 + (t10 * t30 - t33 * t9) * t50 + t43; Ifges(5,5) * t31 + Ifges(5,6) * t34 + t41 * pkin(6) + (m(6) * (t30 * t7 + t33 * t6) + (t30 * t15 - t33 * t16) * mrSges(6,3)) * pkin(4) + t39; Ifges(5,3) + Ifges(6,3) + m(6) * (t30 ^ 2 + t33 ^ 2) * pkin(4) ^ 2 + 0.2e1 * t38; t44; t43; t39; Ifges(6,3) + t38; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
