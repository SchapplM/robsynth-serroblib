% Calculate joint inertia matrix for
% S5PRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
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
% Datum: 2019-12-05 15:25
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPPR2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR2_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR2_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR2_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR2_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPPR2_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:23:59
% EndTime: 2019-12-05 15:24:00
% DurationCPUTime: 0.20s
% Computational Cost: add. (218->72), mult. (473->116), div. (0->0), fcn. (436->8), ass. (0->31)
t24 = sin(pkin(9));
t26 = cos(pkin(9));
t16 = -t26 * mrSges(5,1) + t24 * mrSges(5,2);
t28 = sin(qJ(5));
t30 = cos(qJ(5));
t12 = -t28 * t24 + t30 * t26;
t14 = t30 * t24 + t28 * t26;
t5 = -t12 * mrSges(6,1) + t14 * mrSges(6,2);
t36 = t16 + t5;
t25 = sin(pkin(8));
t27 = cos(pkin(8));
t29 = sin(qJ(2));
t31 = cos(qJ(2));
t10 = t25 * t29 - t27 * t31;
t8 = t10 ^ 2;
t23 = t26 ^ 2;
t18 = t25 * pkin(2) + qJ(4);
t35 = pkin(6) + t18;
t34 = t24 ^ 2 + t23;
t20 = -t27 * pkin(2) - pkin(3);
t33 = t34 * t18;
t15 = -t26 * pkin(4) + t20;
t13 = t25 * t31 + t27 * t29;
t9 = t13 ^ 2;
t7 = t35 * t26;
t6 = t35 * t24;
t4 = -t28 * t6 + t30 * t7;
t3 = -t28 * t7 - t30 * t6;
t2 = t12 * t13;
t1 = t14 * t13;
t11 = [m(2) + m(6) * (t1 ^ 2 + t2 ^ 2 + t8) + m(4) * (t9 + t8) + m(5) * (t34 * t9 + t8) + m(3) * (t29 ^ 2 + t31 ^ 2); t31 * mrSges(3,1) - t29 * mrSges(3,2) + (t1 * t14 + t2 * t12) * mrSges(6,3) + (t34 * mrSges(5,3) - mrSges(4,2)) * t13 + (-mrSges(4,1) + t36) * t10 + m(6) * (-t3 * t1 + t15 * t10 + t4 * t2) + m(5) * (t20 * t10 + t13 * t33) + m(4) * (-t10 * t27 + t13 * t25) * pkin(2); Ifges(5,2) * t23 + 0.2e1 * t15 * t5 + 0.2e1 * t20 * t16 + Ifges(3,3) + Ifges(4,3) + (Ifges(5,1) * t24 + 0.2e1 * Ifges(5,4) * t26) * t24 + (-0.2e1 * t3 * mrSges(6,3) + Ifges(6,1) * t14) * t14 + (0.2e1 * t4 * mrSges(6,3) + 0.2e1 * Ifges(6,4) * t14 + Ifges(6,2) * t12) * t12 + m(6) * (t15 ^ 2 + t3 ^ 2 + t4 ^ 2) + m(5) * (t34 * t18 ^ 2 + t20 ^ 2) + m(4) * (t25 ^ 2 + t27 ^ 2) * pkin(2) ^ 2 + 0.2e1 * (t27 * mrSges(4,1) - t25 * mrSges(4,2)) * pkin(2) + 0.2e1 * mrSges(5,3) * t33; m(6) * (-t12 * t1 + t14 * t2); m(6) * (t12 * t3 + t14 * t4); m(4) + m(5) * t34 + m(6) * (t12 ^ 2 + t14 ^ 2); 0.2e1 * (m(6) / 0.2e1 + m(5) / 0.2e1) * t10; m(5) * t20 + m(6) * t15 + t36; 0; m(5) + m(6); -t1 * mrSges(6,1) - t2 * mrSges(6,2); t3 * mrSges(6,1) - t4 * mrSges(6,2) + Ifges(6,5) * t14 + Ifges(6,6) * t12; -t5; 0; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t11(1), t11(2), t11(4), t11(7), t11(11); t11(2), t11(3), t11(5), t11(8), t11(12); t11(4), t11(5), t11(6), t11(9), t11(13); t11(7), t11(8), t11(9), t11(10), t11(14); t11(11), t11(12), t11(13), t11(14), t11(15);];
Mq = res;
