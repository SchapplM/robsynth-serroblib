% Calculate joint inertia matrix for
% S5PRPRR5
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
% Datum: 2019-12-05 15:55
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRR5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR5_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR5_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR5_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR5_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR5_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:53:27
% EndTime: 2019-12-05 15:53:28
% DurationCPUTime: 0.28s
% Computational Cost: add. (405->97), mult. (872->149), div. (0->0), fcn. (900->8), ass. (0->41)
t36 = sin(pkin(9));
t37 = cos(pkin(9));
t39 = sin(qJ(4));
t42 = cos(qJ(4));
t25 = -t39 * t36 + t42 * t37;
t26 = t42 * t36 + t39 * t37;
t15 = -t25 * mrSges(5,1) + t26 * mrSges(5,2);
t28 = -t37 * mrSges(4,1) + t36 * mrSges(4,2);
t38 = sin(qJ(5));
t41 = cos(qJ(5));
t13 = t41 * t25 - t38 * t26;
t14 = t38 * t25 + t41 * t26;
t4 = -t13 * mrSges(6,1) + t14 * mrSges(6,2);
t53 = -t15 - t28 - t4;
t33 = t37 ^ 2;
t52 = pkin(6) + qJ(3);
t27 = t52 * t36;
t29 = t52 * t37;
t17 = -t39 * t27 + t42 * t29;
t51 = t36 ^ 2 + t33;
t50 = m(4) + m(5) + m(6);
t40 = sin(qJ(2));
t19 = t26 * t40;
t20 = t25 * t40;
t6 = -t41 * t19 - t38 * t20;
t7 = -t38 * t19 + t41 * t20;
t49 = t6 * mrSges(6,1) - t7 * mrSges(6,2);
t30 = -t37 * pkin(3) - pkin(2);
t16 = -t42 * t27 - t39 * t29;
t48 = qJ(3) * t51;
t8 = -t26 * pkin(7) + t16;
t9 = t25 * pkin(7) + t17;
t2 = -t38 * t9 + t41 * t8;
t3 = t38 * t8 + t41 * t9;
t47 = t2 * mrSges(6,1) - t3 * mrSges(6,2) + Ifges(6,5) * t14 + Ifges(6,6) * t13;
t46 = (t41 * mrSges(6,1) - t38 * mrSges(6,2)) * pkin(4);
t43 = cos(qJ(2));
t35 = t43 ^ 2;
t34 = t40 ^ 2;
t18 = -t25 * pkin(4) + t30;
t1 = [m(2) + m(6) * (t6 ^ 2 + t7 ^ 2 + t35) + m(5) * (t19 ^ 2 + t20 ^ 2 + t35) + m(4) * (t51 * t34 + t35) + m(3) * (t34 + t35); (t7 * t13 - t6 * t14) * mrSges(6,3) + (t19 * t26 + t20 * t25) * mrSges(5,3) + (t51 * mrSges(4,3) - mrSges(3,2)) * t40 + (mrSges(3,1) + t53) * t43 + m(6) * (-t18 * t43 + t2 * t6 + t3 * t7) + m(5) * (-t16 * t19 + t17 * t20 - t30 * t43) + m(4) * (t43 * pkin(2) + t40 * t48); Ifges(4,2) * t33 - 0.2e1 * pkin(2) * t28 + 0.2e1 * t30 * t15 + 0.2e1 * t18 * t4 + Ifges(3,3) + (Ifges(4,1) * t36 + 0.2e1 * Ifges(4,4) * t37) * t36 + (-0.2e1 * t16 * mrSges(5,3) + Ifges(5,1) * t26) * t26 + (-0.2e1 * t2 * mrSges(6,3) + Ifges(6,1) * t14) * t14 + 0.2e1 * mrSges(4,3) * t48 + (0.2e1 * t17 * mrSges(5,3) + 0.2e1 * Ifges(5,4) * t26 + Ifges(5,2) * t25) * t25 + (0.2e1 * t3 * mrSges(6,3) + 0.2e1 * Ifges(6,4) * t14 + Ifges(6,2) * t13) * t13 + m(6) * (t18 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(5) * (t16 ^ 2 + t17 ^ 2 + t30 ^ 2) + m(4) * (t51 * qJ(3) ^ 2 + pkin(2) ^ 2); -t50 * t43; -m(4) * pkin(2) + m(5) * t30 + m(6) * t18 - t53; t50; -t19 * mrSges(5,1) - t20 * mrSges(5,2) + m(6) * (t38 * t7 + t41 * t6) * pkin(4) + t49; t16 * mrSges(5,1) - t17 * mrSges(5,2) + Ifges(5,5) * t26 + Ifges(5,6) * t25 + (m(6) * (t2 * t41 + t3 * t38) + (t38 * t13 - t41 * t14) * mrSges(6,3)) * pkin(4) + t47; 0; Ifges(5,3) + Ifges(6,3) + m(6) * (t38 ^ 2 + t41 ^ 2) * pkin(4) ^ 2 + 0.2e1 * t46; t49; t47; 0; Ifges(6,3) + t46; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
