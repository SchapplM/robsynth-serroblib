% Calculate joint inertia matrix for
% S5RPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% m [6x1]
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
% Datum: 2022-01-23 09:13
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRP1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP1_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP1_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP1_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:12:12
% EndTime: 2022-01-23 09:12:12
% DurationCPUTime: 0.34s
% Computational Cost: add. (261->99), mult. (515->131), div. (0->0), fcn. (391->6), ass. (0->47)
t58 = -Ifges(5,5) - Ifges(6,5);
t32 = sin(pkin(7));
t24 = t32 * pkin(1) + qJ(3);
t57 = 0.2e1 * t24;
t56 = m(5) + m(6);
t31 = sin(pkin(8));
t35 = sin(qJ(4));
t7 = (pkin(4) * t35 + t24) * t31;
t36 = cos(qJ(4));
t48 = t36 * t31;
t49 = t31 * t35;
t9 = mrSges(6,1) * t49 + mrSges(6,2) * t48;
t55 = m(6) * t7 + t9;
t30 = t36 ^ 2;
t54 = m(6) * pkin(4);
t34 = cos(pkin(7));
t25 = -t34 * pkin(1) - pkin(2);
t33 = cos(pkin(8));
t11 = -t33 * pkin(3) - t31 * pkin(6) + t25;
t51 = t24 * t33;
t4 = t35 * t11 + t36 * t51;
t52 = mrSges(5,2) * t36;
t50 = t24 * t35;
t47 = Ifges(5,6) + Ifges(6,6);
t46 = Ifges(5,3) + Ifges(6,3);
t12 = t33 * mrSges(6,2) - mrSges(6,3) * t49;
t13 = t33 * mrSges(5,2) - mrSges(5,3) * t49;
t45 = t12 + t13;
t14 = -t33 * mrSges(6,1) - mrSges(6,3) * t48;
t15 = -t33 * mrSges(5,1) - mrSges(5,3) * t48;
t44 = t14 + t15;
t43 = t58 * t48;
t27 = t31 ^ 2;
t28 = t33 ^ 2;
t42 = t27 + t28;
t41 = t35 ^ 2 + t30;
t40 = qJ(5) * t31;
t38 = -mrSges(5,1) - t54;
t26 = t31 * mrSges(4,2);
t23 = t24 ^ 2;
t17 = t27 * t23;
t10 = (mrSges(5,1) * t35 + t52) * t31;
t6 = t36 * t11;
t3 = -t33 * t50 + t6;
t2 = -t35 * t40 + t4;
t1 = -t36 * t40 + t6 + (-pkin(4) - t50) * t33;
t5 = [0.2e1 * t1 * t14 + 0.2e1 * t2 * t12 + 0.2e1 * t4 * t13 + 0.2e1 * t3 * t15 + 0.2e1 * t25 * t26 + 0.2e1 * t7 * t9 + Ifges(2,3) + Ifges(3,3) + 0.2e1 * (t34 * mrSges(3,1) - t32 * mrSges(3,2)) * pkin(1) + t42 * mrSges(4,3) * t57 + (-0.2e1 * t25 * mrSges(4,1) + (Ifges(4,2) + t46) * t33 + t43) * t33 + m(4) * (t28 * t23 + t25 ^ 2 + t17) + m(6) * (t1 ^ 2 + t2 ^ 2 + t7 ^ 2) + m(5) * (t3 ^ 2 + t4 ^ 2 + t17) + m(3) * (t32 ^ 2 + t34 ^ 2) * pkin(1) ^ 2 + (t10 * t57 + ((Ifges(6,2) + Ifges(5,2)) * t35 + 0.2e1 * (-Ifges(5,4) - Ifges(6,4)) * t36) * t49 + (0.2e1 * t47 * t35 + t58 * t36 + (2 * Ifges(4,4))) * t33 + (Ifges(4,1) + (Ifges(5,1) + Ifges(6,1)) * t30) * t31) * t31; (-t10 - t55) * t33 + (t45 * t36 - t44 * t35 + m(5) * (-t3 * t35 + t36 * t4 - t51) + m(6) * (-t1 * t35 + t2 * t36)) * t31; m(3) + m(4) * t42 + (t41 * t27 + t28) * t56; -t33 * mrSges(4,1) + t26 + t44 * t36 + t45 * t35 + m(6) * (t36 * t1 + t35 * t2) + m(5) * (t36 * t3 + t35 * t4) + m(4) * t25; 0; t41 * t56 + m(4); t3 * mrSges(5,1) + t1 * mrSges(6,1) - t4 * mrSges(5,2) - t2 * mrSges(6,2) - t46 * t33 - t47 * t49 + (m(6) * t1 + t14) * pkin(4) - t43; (t38 * t35 - t52) * t31 - t9; (-mrSges(5,2) - mrSges(6,2)) * t35 + (mrSges(6,1) - t38) * t36; (0.2e1 * mrSges(6,1) + t54) * pkin(4) + t46; t55; -m(6) * t33; 0; 0; m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t5(1), t5(2), t5(4), t5(7), t5(11); t5(2), t5(3), t5(5), t5(8), t5(12); t5(4), t5(5), t5(6), t5(9), t5(13); t5(7), t5(8), t5(9), t5(10), t5(14); t5(11), t5(12), t5(13), t5(14), t5(15);];
Mq = res;
