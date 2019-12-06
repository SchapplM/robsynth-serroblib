% Calculate joint inertia matrix for
% S5PRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:18
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPR2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR2_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR2_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR2_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR2_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR2_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:17:19
% EndTime: 2019-12-05 16:17:21
% DurationCPUTime: 0.33s
% Computational Cost: add. (285->98), mult. (587->132), div. (0->0), fcn. (430->6), ass. (0->48)
t33 = sin(pkin(9));
t37 = cos(qJ(5));
t52 = t33 * t37;
t60 = Ifges(6,5) * t52;
t29 = t33 ^ 2;
t34 = cos(pkin(9));
t30 = t34 ^ 2;
t50 = t29 + t30;
t59 = t50 * mrSges(5,3);
t36 = sin(qJ(3));
t25 = t36 * pkin(2) + qJ(4);
t58 = qJ(4) + t25;
t35 = sin(qJ(5));
t53 = t33 * t35;
t15 = t34 * mrSges(6,2) - mrSges(6,3) * t53;
t57 = 0.2e1 * t15;
t16 = -t34 * mrSges(6,1) - mrSges(6,3) * t52;
t56 = 0.2e1 * t16;
t38 = cos(qJ(3));
t55 = t38 * pkin(2);
t54 = t25 * t34;
t51 = t35 * t16;
t49 = t35 ^ 2 + t37 ^ 2;
t48 = qJ(4) * t25;
t47 = qJ(4) * t34;
t11 = -mrSges(6,1) * t53 - mrSges(6,2) * t52;
t46 = t34 * t11 + t15 * t52;
t19 = -t34 * mrSges(5,1) + t33 * mrSges(5,2);
t45 = t35 * t15 + t37 * t16 + t19;
t18 = -t34 * pkin(4) - t33 * pkin(7) - pkin(3);
t6 = -Ifges(6,6) * t53 - Ifges(6,3) * t34 + t60;
t44 = Ifges(4,3) + ((Ifges(6,1) * t37 - Ifges(6,4) * t35) * t52 + Ifges(5,1) * t33) * t33 + (0.2e1 * Ifges(5,4) * t33 + Ifges(5,2) * t34 - t6 - t60) * t34;
t43 = (t38 * mrSges(4,1) - t36 * mrSges(4,2)) * pkin(2);
t42 = -0.2e1 * t33 * t11 + 0.2e1 * t59;
t7 = -Ifges(6,6) * t34 + (Ifges(6,4) * t37 - Ifges(6,2) * t35) * t33;
t41 = -t7 * t53 + t44;
t39 = qJ(4) ^ 2;
t27 = t29 * t39;
t26 = -pkin(3) - t55;
t24 = t25 ^ 2;
t20 = t29 * t24;
t17 = t29 * t48;
t12 = t18 - t55;
t5 = t35 * t18 + t37 * t47;
t4 = t37 * t18 - t35 * t47;
t2 = t35 * t12 + t37 * t54;
t1 = t37 * t12 - t35 * t54;
t3 = [m(2) + m(3) + m(4) + m(5) * t50 + m(6) * (t49 * t29 + t30); (m(6) * (-t1 * t35 + t2 * t37 - t54) - t51) * t33 + t46; t1 * t56 + t2 * t57 + 0.2e1 * t26 * t19 + Ifges(3,3) + 0.2e1 * t43 + t42 * t25 + m(6) * (t1 ^ 2 + t2 ^ 2 + t20) + m(5) * (t30 * t24 + t26 ^ 2 + t20) + m(4) * (t36 ^ 2 + t38 ^ 2) * pkin(2) ^ 2 + t41; (m(6) * (-t35 * t4 + t37 * t5 - t47) - t51) * t33 + t46; (t26 - pkin(3)) * t19 + (t1 + t4) * t16 + (t2 + t5) * t15 + t43 + m(6) * (t4 * t1 + t5 * t2 + t17) + m(5) * (-pkin(3) * t26 + t30 * t48 + t17) + (-t58 * t11 - t35 * t7) * t33 + t44 + t58 * t59; -0.2e1 * pkin(3) * t19 + t5 * t57 + t4 * t56 + m(6) * (t4 ^ 2 + t5 ^ 2 + t27) + m(5) * (pkin(3) ^ 2 + t30 * t39 + t27) + t42 * qJ(4) + t41; 0; m(6) * (t37 * t1 + t35 * t2) + m(5) * t26 + t45; m(6) * (t35 * t5 + t37 * t4) - m(5) * pkin(3) + t45; m(6) * t49 + m(5); t11; t1 * mrSges(6,1) - t2 * mrSges(6,2) + t6; t4 * mrSges(6,1) - t5 * mrSges(6,2) + t6; t37 * mrSges(6,1) - t35 * mrSges(6,2); Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t3(1), t3(2), t3(4), t3(7), t3(11); t3(2), t3(3), t3(5), t3(8), t3(12); t3(4), t3(5), t3(6), t3(9), t3(13); t3(7), t3(8), t3(9), t3(10), t3(14); t3(11), t3(12), t3(13), t3(14), t3(15);];
Mq = res;
