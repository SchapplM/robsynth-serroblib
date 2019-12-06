% Calculate joint inertia matrix for
% S5RPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-05 18:00
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_inertiaJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP1_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP1_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP1_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:59:14
% EndTime: 2019-12-05 17:59:16
% DurationCPUTime: 0.34s
% Computational Cost: add. (359->105), mult. (602->137), div. (0->0), fcn. (516->4), ass. (0->38)
t38 = cos(qJ(3));
t62 = t38 ^ 2;
t35 = sin(qJ(4));
t36 = sin(qJ(3));
t37 = cos(qJ(4));
t21 = -t35 * t38 - t37 * t36;
t23 = -t35 * t36 + t37 * t38;
t61 = t21 ^ 2 + t23 ^ 2;
t51 = mrSges(5,2) + mrSges(6,2);
t60 = (t37 * mrSges(5,1) - t51 * t35) * pkin(3);
t58 = 2 * mrSges(6,1);
t56 = m(6) * pkin(4);
t39 = -pkin(1) - pkin(6);
t55 = -pkin(7) + t39;
t54 = t35 * t21;
t52 = t37 * t23;
t49 = Ifges(5,3) + Ifges(6,3);
t24 = t55 * t36;
t25 = t55 * t38;
t6 = t37 * t24 + t35 * t25;
t48 = t36 ^ 2 + t62;
t26 = t36 * pkin(3) + qJ(2);
t47 = m(4) * t48;
t45 = t48 * mrSges(4,3);
t5 = -t35 * t24 + t37 * t25;
t2 = -t23 * qJ(5) + t5;
t44 = m(6) * t2 - t23 * mrSges(6,3);
t43 = t51 * t21 + (mrSges(5,1) + mrSges(6,1)) * t23;
t3 = t21 * qJ(5) + t6;
t42 = t5 * mrSges(5,1) + t2 * mrSges(6,1) - t6 * mrSges(5,2) - t3 * mrSges(6,2) + (Ifges(5,5) + Ifges(6,5)) * t23 + (Ifges(5,6) + Ifges(6,6)) * t21;
t41 = pkin(3) ^ 2;
t40 = qJ(2) ^ 2;
t30 = t35 ^ 2 * t41;
t29 = t37 * pkin(3) + pkin(4);
t10 = t23 * mrSges(6,2);
t9 = pkin(3) * t54;
t8 = -t21 * pkin(4) + t26;
t1 = [Ifges(3,1) + Ifges(2,3) + 0.2e1 * t8 * t10 + Ifges(4,1) * t62 - (2 * pkin(1) * mrSges(3,2)) + (0.2e1 * t26 * mrSges(5,2) - 0.2e1 * t5 * mrSges(5,3) - 0.2e1 * t2 * mrSges(6,3) + (Ifges(6,1) + Ifges(5,1)) * t23) * t23 + m(6) * (t2 ^ 2 + t3 ^ 2 + t8 ^ 2) + m(5) * (t26 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(4) * (t48 * t39 ^ 2 + t40) + m(3) * ((pkin(1) ^ 2) + t40) - 0.2e1 * t39 * t45 + (-0.2e1 * Ifges(4,4) * t38 + Ifges(4,2) * t36) * t36 + 0.2e1 * (t36 * mrSges(4,1) + t38 * mrSges(4,2) + mrSges(3,3)) * qJ(2) + (-0.2e1 * t26 * mrSges(5,1) - 0.2e1 * t8 * mrSges(6,1) + 0.2e1 * t6 * mrSges(5,3) + 0.2e1 * t3 * mrSges(6,3) + 0.2e1 * (Ifges(5,4) + Ifges(6,4)) * t23 + (Ifges(6,2) + Ifges(5,2)) * t21) * t21; -m(3) * pkin(1) + mrSges(3,2) - t45 + m(6) * (t23 * t2 - t21 * t3) + m(5) * (-t21 * t6 + t23 * t5) + t39 * t47 + t61 * (-mrSges(5,3) - mrSges(6,3)); m(3) + t47 + 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * t61; (t39 * mrSges(4,1) + Ifges(4,5)) * t38 + (-t39 * mrSges(4,2) - Ifges(4,6)) * t36 + t44 * t29 + (mrSges(6,3) * t54 + (-t52 + t54) * mrSges(5,3) + m(6) * t3 * t35 + m(5) * (t35 * t6 + t37 * t5)) * pkin(3) + t42; t38 * mrSges(4,1) - t36 * mrSges(4,2) + m(5) * (pkin(3) * t52 - t9) + m(6) * (t29 * t23 - t9) + t43; t29 * t58 + Ifges(4,3) + m(5) * (t37 ^ 2 * t41 + t30) + m(6) * (t29 ^ 2 + t30) + 0.2e1 * t60 + t49; t44 * pkin(4) + t42; t23 * t56 + t43; t29 * t56 + (pkin(4) + t29) * mrSges(6,1) + t60 + t49; (t58 + t56) * pkin(4) + t49; m(6) * t8 - t21 * mrSges(6,1) + t10; 0; 0; 0; m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
