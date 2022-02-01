% Calculate joint inertia matrix for
% S5RPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2022-01-23 09:30
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP3_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP3_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP3_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:29:33
% EndTime: 2022-01-23 09:29:33
% DurationCPUTime: 0.32s
% Computational Cost: add. (348->98), mult. (639->138), div. (0->0), fcn. (555->6), ass. (0->39)
t37 = cos(qJ(3));
t56 = t37 ^ 2;
t34 = sin(qJ(4));
t36 = cos(qJ(4));
t55 = (t36 * mrSges(5,1) + (-mrSges(5,2) - mrSges(6,2)) * t34) * pkin(3);
t54 = 2 * mrSges(6,1);
t53 = m(6) * pkin(4);
t35 = sin(qJ(3));
t22 = -t34 * t35 + t36 * t37;
t52 = t22 * pkin(4);
t51 = t36 * pkin(3);
t32 = sin(pkin(8));
t26 = t32 * pkin(1) + pkin(6);
t50 = pkin(7) + t26;
t13 = t50 * t35;
t14 = t50 * t37;
t6 = -t34 * t13 + t36 * t14;
t49 = t34 * t22;
t47 = Ifges(5,3) + Ifges(6,3);
t23 = t34 * t37 + t36 * t35;
t46 = -t22 * mrSges(6,1) + t23 * mrSges(6,2);
t45 = t35 ^ 2 + t56;
t33 = cos(pkin(8));
t27 = -t33 * pkin(1) - pkin(2);
t5 = -t36 * t13 - t34 * t14;
t2 = -t23 * qJ(5) + t5;
t43 = m(6) * t2 - t23 * mrSges(6,3);
t42 = -t37 * mrSges(4,1) + t35 * mrSges(4,2);
t16 = t22 * mrSges(5,1);
t41 = -t23 * mrSges(5,2) + t16 - t46;
t24 = -t37 * pkin(3) + t27;
t3 = t22 * qJ(5) + t6;
t40 = t5 * mrSges(5,1) + t2 * mrSges(6,1) - t6 * mrSges(5,2) - t3 * mrSges(6,2) + (Ifges(5,5) + Ifges(6,5)) * t23 + (Ifges(5,6) + Ifges(6,6)) * t22;
t38 = pkin(3) ^ 2;
t29 = t34 ^ 2 * t38;
t28 = pkin(4) + t51;
t12 = t34 * pkin(3) * t23;
t8 = t24 - t52;
t1 = [Ifges(2,3) + Ifges(3,3) - 0.2e1 * t24 * t16 + 0.2e1 * t8 * t46 + 0.2e1 * t27 * t42 + Ifges(4,2) * t56 + (0.2e1 * t24 * mrSges(5,2) - 0.2e1 * t5 * mrSges(5,3) - 0.2e1 * t2 * mrSges(6,3) + (Ifges(6,1) + Ifges(5,1)) * t23) * t23 + m(6) * (t2 ^ 2 + t3 ^ 2 + t8 ^ 2) + m(5) * (t24 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(4) * (t45 * t26 ^ 2 + t27 ^ 2) + m(3) * (t32 ^ 2 + t33 ^ 2) * pkin(1) ^ 2 + (Ifges(4,1) * t35 + 0.2e1 * Ifges(4,4) * t37) * t35 + 0.2e1 * (t33 * mrSges(3,1) - t32 * mrSges(3,2)) * pkin(1) + 0.2e1 * t45 * t26 * mrSges(4,3) + (0.2e1 * t6 * mrSges(5,3) + 0.2e1 * t3 * mrSges(6,3) + 0.2e1 * (Ifges(5,4) + Ifges(6,4)) * t23 + (Ifges(6,2) + Ifges(5,2)) * t22) * t22; m(5) * (t5 * t22 + t6 * t23) + m(6) * (t2 * t22 + t3 * t23); m(3) + m(4) * t45 + 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * (t22 ^ 2 + t23 ^ 2); Ifges(4,5) * t35 + Ifges(4,6) * t37 + t43 * t28 + (-t35 * mrSges(4,1) - t37 * mrSges(4,2)) * t26 + (mrSges(6,3) * t49 + (-t36 * t23 + t49) * mrSges(5,3) + m(6) * t3 * t34 + m(5) * (t34 * t6 + t36 * t5)) * pkin(3) + t40; m(5) * (t22 * t51 + t12) + m(6) * (t28 * t22 + t12) + t41 - t42; t28 * t54 + Ifges(4,3) + m(5) * (t36 ^ 2 * t38 + t29) + m(6) * (t28 ^ 2 + t29) + 0.2e1 * t55 + t47; t43 * pkin(4) + t40; m(6) * t52 + t41; t28 * t53 + (pkin(4) + t28) * mrSges(6,1) + t55 + t47; (t54 + t53) * pkin(4) + t47; m(6) * t8 + t46; 0; 0; 0; m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
