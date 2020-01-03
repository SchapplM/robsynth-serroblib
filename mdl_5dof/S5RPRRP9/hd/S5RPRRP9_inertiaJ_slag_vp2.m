% Calculate joint inertia matrix for
% S5RPRRP9
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
% Datum: 2019-12-31 18:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP9_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP9_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP9_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP9_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP9_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP9_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:48:37
% EndTime: 2019-12-31 18:48:38
% DurationCPUTime: 0.34s
% Computational Cost: add. (575->113), mult. (1081->151), div. (0->0), fcn. (1118->6), ass. (0->39)
t34 = cos(pkin(8));
t32 = t34 ^ 2;
t55 = 2 * mrSges(6,3);
t33 = sin(pkin(8));
t36 = sin(qJ(3));
t38 = cos(qJ(3));
t25 = -t36 * t33 + t38 * t34;
t28 = -t34 * pkin(2) - pkin(1);
t20 = -t25 * pkin(3) + t28;
t54 = 0.2e1 * t20;
t53 = 0.2e1 * t25;
t52 = mrSges(6,2) + mrSges(5,3);
t51 = Ifges(6,2) + Ifges(5,3);
t50 = pkin(6) + qJ(2);
t46 = t50 * t34;
t47 = t50 * t33;
t19 = -t36 * t47 + t38 * t46;
t49 = t33 ^ 2 + t32;
t35 = sin(qJ(4));
t37 = cos(qJ(4));
t18 = -t36 * t46 - t38 * t47;
t26 = t38 * t33 + t36 * t34;
t41 = -t26 * pkin(7) + t18;
t9 = t25 * pkin(7) + t19;
t5 = t35 * t9 - t37 * t41;
t7 = t35 * t41 + t37 * t9;
t48 = t5 ^ 2 + t7 ^ 2;
t45 = -t34 * mrSges(3,1) + t33 * mrSges(3,2);
t44 = -t25 * mrSges(4,1) + t26 * mrSges(4,2);
t43 = (t37 * mrSges(5,1) - t35 * mrSges(5,2)) * pkin(3);
t16 = -t37 * t25 + t35 * t26;
t17 = t35 * t25 + t37 * t26;
t42 = (-mrSges(5,2) + mrSges(6,3)) * t7 + (-mrSges(5,1) - mrSges(6,1)) * t5 + (Ifges(6,4) + Ifges(5,5)) * t17 + (-Ifges(5,6) + Ifges(6,6)) * t16;
t29 = -t37 * pkin(3) - pkin(4);
t27 = t35 * pkin(3) + qJ(5);
t11 = t17 * mrSges(5,2);
t10 = t16 * mrSges(6,1);
t4 = t16 * pkin(4) - t17 * qJ(5) + t20;
t1 = [Ifges(3,2) * t32 - 0.2e1 * pkin(1) * t45 + Ifges(4,2) * t25 ^ 2 + 0.2e1 * t28 * t44 + 0.2e1 * t4 * t10 + t11 * t54 + t19 * mrSges(4,3) * t53 + Ifges(2,3) + (Ifges(3,1) * t33 + 0.2e1 * Ifges(3,4) * t34) * t33 + 0.2e1 * t49 * qJ(2) * mrSges(3,3) + (-0.2e1 * t18 * mrSges(4,3) + Ifges(4,1) * t26 + Ifges(4,4) * t53) * t26 + (-0.2e1 * t4 * mrSges(6,3) + (Ifges(6,1) + Ifges(5,1)) * t17 + 0.2e1 * t52 * t5) * t17 + (mrSges(5,1) * t54 + (Ifges(6,3) + Ifges(5,2)) * t16 - 0.2e1 * t52 * t7 + 0.2e1 * (-Ifges(5,4) + Ifges(6,5)) * t17) * t16 + m(6) * (t4 ^ 2 + t48) + m(5) * (t20 ^ 2 + t48) + m(4) * (t18 ^ 2 + t19 ^ 2 + t28 ^ 2) + m(3) * (t49 * qJ(2) ^ 2 + pkin(1) ^ 2); -m(3) * pkin(1) + m(4) * t28 + m(5) * t20 + m(6) * t4 + t16 * mrSges(5,1) - t17 * mrSges(6,3) + t10 + t11 + t44 + t45; m(3) + m(4) + m(5) + m(6); m(6) * (t27 * t7 + t29 * t5) - t19 * mrSges(4,2) + Ifges(4,5) * t26 + Ifges(4,6) * t25 + t18 * mrSges(4,1) + (-t27 * t16 + t29 * t17) * mrSges(6,2) + (m(5) * (t35 * t7 - t37 * t5) + (-t35 * t16 - t37 * t17) * mrSges(5,3)) * pkin(3) + t42; 0; -0.2e1 * t29 * mrSges(6,1) + t27 * t55 + Ifges(4,3) + 0.2e1 * t43 + m(6) * (t27 ^ 2 + t29 ^ 2) + m(5) * (t35 ^ 2 + t37 ^ 2) * pkin(3) ^ 2 + t51; m(6) * (-pkin(4) * t5 + qJ(5) * t7) + (-pkin(4) * t17 - qJ(5) * t16) * mrSges(6,2) + t42; 0; m(6) * (-pkin(4) * t29 + qJ(5) * t27) + t43 + (t27 + qJ(5)) * mrSges(6,3) + (-t29 + pkin(4)) * mrSges(6,1) + t51; 0.2e1 * pkin(4) * mrSges(6,1) + qJ(5) * t55 + m(6) * (pkin(4) ^ 2 + qJ(5) ^ 2) + t51; m(6) * t5 + t17 * mrSges(6,2); 0; m(6) * t29 - mrSges(6,1); -m(6) * pkin(4) - mrSges(6,1); m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
