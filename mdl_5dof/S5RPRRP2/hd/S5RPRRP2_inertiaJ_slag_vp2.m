% Calculate joint inertia matrix for
% S5RPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-05 18:02
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP2_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP2_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP2_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:01:24
% EndTime: 2019-12-05 18:01:25
% DurationCPUTime: 0.27s
% Computational Cost: add. (266->95), mult. (490->118), div. (0->0), fcn. (323->6), ass. (0->39)
t35 = cos(qJ(4));
t55 = t35 ^ 2;
t33 = sin(qJ(4));
t45 = t33 ^ 2 + t55;
t54 = 0.2e1 * t45;
t24 = t33 * mrSges(6,2);
t14 = -t35 * mrSges(6,1) + t24;
t53 = 0.2e1 * t14;
t52 = m(6) * pkin(4);
t31 = sin(pkin(8));
t51 = pkin(1) * t31;
t50 = t35 * pkin(4);
t32 = cos(pkin(8));
t21 = t32 * pkin(1) + pkin(2);
t34 = sin(qJ(3));
t36 = cos(qJ(3));
t7 = t36 * t21 - t34 * t51;
t49 = t7 * mrSges(4,1);
t8 = t34 * t21 + t36 * t51;
t48 = t8 * mrSges(4,2);
t47 = t33 * mrSges(5,2);
t46 = t33 * mrSges(6,3);
t44 = 0.2e1 * mrSges(6,3);
t6 = pkin(7) + t8;
t43 = t45 * t6;
t42 = (Ifges(5,6) + Ifges(6,6)) * t35 + (Ifges(5,5) + Ifges(6,5)) * t33;
t5 = -pkin(3) - t7;
t41 = Ifges(4,3) + (Ifges(6,2) + Ifges(5,2)) * t55 + ((Ifges(5,1) + Ifges(6,1)) * t33 + 0.2e1 * (Ifges(5,4) + Ifges(6,4)) * t35) * t33;
t40 = -mrSges(5,1) * t33 - mrSges(5,2) * t35;
t39 = mrSges(5,3) * t54;
t23 = t35 * qJ(5);
t22 = -pkin(3) - t50;
t16 = t35 * pkin(7) + t23;
t15 = -t35 * mrSges(5,1) + t47;
t13 = (-qJ(5) - pkin(7)) * t33;
t3 = t5 - t50;
t2 = t35 * t6 + t23;
t1 = (-qJ(5) - t6) * t33;
t4 = [0.2e1 * t49 - 0.2e1 * t48 + t3 * t53 + 0.2e1 * t5 * t15 + Ifges(2,3) + Ifges(3,3) + (-t1 * t33 + t2 * t35) * t44 + t6 * t39 + m(6) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(5) * (t45 * t6 ^ 2 + t5 ^ 2) + m(4) * (t7 ^ 2 + t8 ^ 2) + t41 + (0.2e1 * t32 * mrSges(3,1) - 0.2e1 * t31 * mrSges(3,2) + m(3) * (t31 ^ 2 + t32 ^ 2) * pkin(1)) * pkin(1); m(6) * (t1 * t35 + t2 * t33); m(3) + m(4) + (m(5) / 0.2e1 + m(6) / 0.2e1) * t54; t49 - t48 + (t5 - pkin(3)) * t15 + (t3 + t22) * t14 + m(6) * (t13 * t1 + t16 * t2 + t22 * t3) + m(5) * (-pkin(3) * t5 + pkin(7) * t43) + ((t16 + t2) * t35 + (-t1 - t13) * t33) * mrSges(6,3) + (t45 * pkin(7) + t43) * mrSges(5,3) + t41; m(6) * (t13 * t35 + t16 * t33); -0.2e1 * pkin(3) * t15 + t22 * t53 + (-t13 * t33 + t16 * t35) * t44 + pkin(7) * t39 + m(6) * (t13 ^ 2 + t16 ^ 2 + t22 ^ 2) + m(5) * (t45 * pkin(7) ^ 2 + pkin(3) ^ 2) + t41; t1 * mrSges(6,1) - t2 * mrSges(6,2) + t40 * t6 + (m(6) * t1 - t46) * pkin(4) + t42; -t47 - t24 + (mrSges(5,1) + mrSges(6,1) + t52) * t35; t13 * mrSges(6,1) - t16 * mrSges(6,2) + t40 * pkin(7) + (m(6) * t13 - t46) * pkin(4) + t42; Ifges(5,3) + Ifges(6,3) + (0.2e1 * mrSges(6,1) + t52) * pkin(4); m(6) * t3 + t14; 0; m(6) * t22 + t14; 0; m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t4(1), t4(2), t4(4), t4(7), t4(11); t4(2), t4(3), t4(5), t4(8), t4(12); t4(4), t4(5), t4(6), t4(9), t4(13); t4(7), t4(8), t4(9), t4(10), t4(14); t4(11), t4(12), t4(13), t4(14), t4(15);];
Mq = res;
