% Calculate joint inertia matrix for
% S5PRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRP2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP2_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP2_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP2_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP2_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP2_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:30:31
% EndTime: 2019-12-05 15:30:33
% DurationCPUTime: 0.32s
% Computational Cost: add. (201->93), mult. (440->120), div. (0->0), fcn. (331->4), ass. (0->43)
t53 = -Ifges(5,5) - Ifges(6,5);
t52 = m(5) + m(6);
t51 = 2 * qJ(3);
t28 = sin(pkin(8));
t30 = sin(qJ(4));
t14 = (pkin(4) * t30 + qJ(3)) * t28;
t45 = t30 * t28;
t31 = cos(qJ(4));
t46 = t28 * t31;
t6 = mrSges(6,1) * t45 + mrSges(6,2) * t46;
t50 = m(6) * t14 + t6;
t27 = t31 ^ 2;
t49 = m(6) * pkin(4);
t29 = cos(pkin(8));
t15 = -t29 * pkin(3) - t28 * pkin(6) - pkin(2);
t37 = qJ(3) * t29;
t4 = t30 * t15 + t31 * t37;
t47 = mrSges(5,2) * t31;
t44 = Ifges(5,6) + Ifges(6,6);
t43 = Ifges(5,3) + Ifges(6,3);
t10 = t29 * mrSges(6,2) - mrSges(6,3) * t45;
t11 = t29 * mrSges(5,2) - mrSges(5,3) * t45;
t42 = t10 + t11;
t12 = -t29 * mrSges(6,1) - mrSges(6,3) * t46;
t13 = -t29 * mrSges(5,1) - mrSges(5,3) * t46;
t41 = t12 + t13;
t40 = t53 * t46;
t24 = t28 ^ 2;
t25 = t29 ^ 2;
t39 = t24 + t25;
t38 = t30 ^ 2 + t27;
t36 = qJ(3) * t30;
t35 = qJ(5) * t28;
t33 = -mrSges(5,1) - t49;
t32 = qJ(3) ^ 2;
t23 = t28 * mrSges(4,2);
t22 = t24 * t32;
t9 = t31 * t15;
t7 = (mrSges(5,1) * t30 + t47) * t28;
t3 = -t29 * t36 + t9;
t2 = -t30 * t35 + t4;
t1 = -t31 * t35 + t9 + (-pkin(4) - t36) * t29;
t5 = [m(2) + m(3) + m(4) * t39 + (t38 * t24 + t25) * t52; (-t7 - t50) * t29 + (t42 * t31 - t41 * t30 + m(5) * (-t3 * t30 + t31 * t4 - t37) + m(6) * (-t1 * t30 + t2 * t31)) * t28; -0.2e1 * pkin(2) * t23 + 0.2e1 * t1 * t12 + 0.2e1 * t2 * t10 + 0.2e1 * t4 * t11 + 0.2e1 * t3 * t13 + 0.2e1 * t14 * t6 + Ifges(3,3) + t39 * mrSges(4,3) * t51 + m(6) * (t1 ^ 2 + t14 ^ 2 + t2 ^ 2) + m(5) * (t3 ^ 2 + t4 ^ 2 + t22) + m(4) * (pkin(2) ^ 2 + t25 * t32 + t22) + (0.2e1 * pkin(2) * mrSges(4,1) + (Ifges(4,2) + t43) * t29 + t40) * t29 + (t7 * t51 + ((Ifges(5,2) + Ifges(6,2)) * t30 + 0.2e1 * (-Ifges(5,4) - Ifges(6,4)) * t31) * t45 + (0.2e1 * t44 * t30 + t53 * t31 + (2 * Ifges(4,4))) * t29 + (Ifges(4,1) + (Ifges(6,1) + Ifges(5,1)) * t27) * t28) * t28; 0; -m(4) * pkin(2) - t29 * mrSges(4,1) + t23 + t41 * t31 + t42 * t30 + m(6) * (t31 * t1 + t30 * t2) + m(5) * (t31 * t3 + t30 * t4); t38 * t52 + m(4); (t33 * t30 - t47) * t28 - t6; t3 * mrSges(5,1) + t1 * mrSges(6,1) - t4 * mrSges(5,2) - t2 * mrSges(6,2) - t43 * t29 - t44 * t45 + (m(6) * t1 + t12) * pkin(4) - t40; (-mrSges(5,2) - mrSges(6,2)) * t30 + (mrSges(6,1) - t33) * t31; (0.2e1 * mrSges(6,1) + t49) * pkin(4) + t43; -m(6) * t29; t50; 0; 0; m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t5(1), t5(2), t5(4), t5(7), t5(11); t5(2), t5(3), t5(5), t5(8), t5(12); t5(4), t5(5), t5(6), t5(9), t5(13); t5(7), t5(8), t5(9), t5(10), t5(14); t5(11), t5(12), t5(13), t5(14), t5(15);];
Mq = res;
