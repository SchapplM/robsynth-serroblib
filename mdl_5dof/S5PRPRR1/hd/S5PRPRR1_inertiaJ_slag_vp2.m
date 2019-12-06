% Calculate joint inertia matrix for
% S5PRPRR1
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
% Datum: 2019-12-05 15:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRR1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR1_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR1_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR1_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR1_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR1_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:42:40
% EndTime: 2019-12-05 15:42:40
% DurationCPUTime: 0.24s
% Computational Cost: add. (352->78), mult. (695->119), div. (0->0), fcn. (718->6), ass. (0->32)
t27 = sin(pkin(9));
t28 = cos(pkin(9));
t30 = sin(qJ(4));
t32 = cos(qJ(4));
t19 = -t30 * t27 + t32 * t28;
t43 = t19 ^ 2;
t26 = t28 ^ 2;
t42 = 0.2e1 * t19;
t41 = pkin(6) + qJ(3);
t21 = t41 * t27;
t22 = t41 * t28;
t13 = -t30 * t21 + t32 * t22;
t40 = t27 ^ 2 + t26;
t23 = -t28 * pkin(3) - pkin(2);
t20 = t32 * t27 + t30 * t28;
t29 = sin(qJ(5));
t31 = cos(qJ(5));
t10 = t31 * t19 - t29 * t20;
t11 = t29 * t19 + t31 * t20;
t4 = -t10 * mrSges(6,1) + t11 * mrSges(6,2);
t39 = -t28 * mrSges(4,1) + t27 * mrSges(4,2);
t38 = -t19 * mrSges(5,1) + t20 * mrSges(5,2);
t12 = -t32 * t21 - t30 * t22;
t5 = -t20 * pkin(7) + t12;
t6 = t19 * pkin(7) + t13;
t2 = -t29 * t6 + t31 * t5;
t3 = t29 * t5 + t31 * t6;
t37 = t2 * mrSges(6,1) - t3 * mrSges(6,2) + Ifges(6,5) * t11 + Ifges(6,6) * t10;
t36 = (t31 * mrSges(6,1) - t29 * mrSges(6,2)) * pkin(4);
t35 = -t38 - t4;
t14 = -t19 * pkin(4) + t23;
t1 = [m(2) + m(3) + m(4) * t40 + m(5) * (t20 ^ 2 + t43) + m(6) * (t10 ^ 2 + t11 ^ 2); m(5) * (t12 * t19 + t13 * t20) + m(6) * (t2 * t10 + t3 * t11); Ifges(3,3) + Ifges(5,2) * t43 + Ifges(4,2) * t26 + 0.2e1 * t14 * t4 + 0.2e1 * t23 * t38 - 0.2e1 * pkin(2) * t39 + t13 * mrSges(5,3) * t42 + (Ifges(4,1) * t27 + 0.2e1 * Ifges(4,4) * t28) * t27 + (-0.2e1 * t2 * mrSges(6,3) + Ifges(6,1) * t11) * t11 + 0.2e1 * t40 * qJ(3) * mrSges(4,3) + (-0.2e1 * t12 * mrSges(5,3) + Ifges(5,1) * t20 + Ifges(5,4) * t42) * t20 + (0.2e1 * t3 * mrSges(6,3) + 0.2e1 * Ifges(6,4) * t11 + Ifges(6,2) * t10) * t10 + m(6) * (t14 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(5) * (t12 ^ 2 + t13 ^ 2 + t23 ^ 2) + m(4) * (t40 * qJ(3) ^ 2 + pkin(2) ^ 2); 0; -m(4) * pkin(2) + m(5) * t23 + m(6) * t14 - t35 + t39; m(4) + m(5) + m(6); m(6) * (t10 * t31 + t11 * t29) * pkin(4) + t35; t12 * mrSges(5,1) - t13 * mrSges(5,2) + Ifges(5,5) * t20 + Ifges(5,6) * t19 + (m(6) * (t2 * t31 + t29 * t3) + (t29 * t10 - t31 * t11) * mrSges(6,3)) * pkin(4) + t37; 0; Ifges(5,3) + Ifges(6,3) + m(6) * (t29 ^ 2 + t31 ^ 2) * pkin(4) ^ 2 + 0.2e1 * t36; -t4; t37; 0; Ifges(6,3) + t36; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
