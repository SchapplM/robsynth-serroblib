% Calculate joint inertia matrix for
% S4RRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPR7_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR7_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR7_inertiaJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR7_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR7_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR7_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:05:58
% EndTime: 2019-12-31 17:05:59
% DurationCPUTime: 0.32s
% Computational Cost: add. (327->102), mult. (650->161), div. (0->0), fcn. (607->6), ass. (0->43)
t38 = cos(qJ(2));
t48 = -qJ(3) - pkin(5);
t19 = t48 * t38;
t33 = sin(pkin(7));
t34 = cos(pkin(7));
t36 = sin(qJ(2));
t43 = t48 * t36;
t9 = -t33 * t19 - t34 * t43;
t56 = t9 ^ 2;
t55 = 0.2e1 * t9;
t26 = -t38 * pkin(2) - pkin(1);
t54 = 0.2e1 * t26;
t37 = cos(qJ(4));
t53 = t37 / 0.2e1;
t35 = sin(qJ(4));
t52 = Ifges(5,4) * t35;
t51 = Ifges(5,4) * t37;
t16 = t33 * t38 + t34 * t36;
t50 = t16 * t35;
t49 = t16 * t37;
t15 = t33 * t36 - t34 * t38;
t47 = Ifges(5,5) * t49 + Ifges(5,3) * t15;
t46 = Ifges(5,5) * t35 + Ifges(5,6) * t37;
t45 = t35 ^ 2 + t37 ^ 2;
t44 = t36 ^ 2 + t38 ^ 2;
t11 = -t34 * t19 + t33 * t43;
t6 = t15 * pkin(3) - t16 * pkin(6) + t26;
t1 = -t35 * t11 + t37 * t6;
t2 = t37 * t11 + t35 * t6;
t42 = -t1 * t35 + t2 * t37;
t41 = mrSges(5,1) * t35 + mrSges(5,2) * t37;
t25 = -t34 * pkin(2) - pkin(3);
t24 = t33 * pkin(2) + pkin(6);
t21 = Ifges(5,1) * t35 + t51;
t20 = Ifges(5,2) * t37 + t52;
t18 = -t37 * mrSges(5,1) + t35 * mrSges(5,2);
t13 = t16 * mrSges(4,2);
t8 = t15 * mrSges(5,1) - mrSges(5,3) * t49;
t7 = -t15 * mrSges(5,2) - mrSges(5,3) * t50;
t5 = t41 * t16;
t4 = Ifges(5,5) * t15 + (Ifges(5,1) * t37 - t52) * t16;
t3 = Ifges(5,6) * t15 + (-Ifges(5,2) * t35 + t51) * t16;
t10 = [Ifges(2,3) + 0.2e1 * t2 * t7 + t5 * t55 + 0.2e1 * t1 * t8 + t13 * t54 - 0.2e1 * pkin(1) * (-t38 * mrSges(3,1) + t36 * mrSges(3,2)) + t36 * (Ifges(3,1) * t36 + Ifges(3,4) * t38) + t38 * (Ifges(3,4) * t36 + Ifges(3,2) * t38) + (mrSges(4,1) * t54 - 0.2e1 * t11 * mrSges(4,3) + Ifges(4,2) * t15 + t47) * t15 + 0.2e1 * t44 * pkin(5) * mrSges(3,3) + m(5) * (t1 ^ 2 + t2 ^ 2 + t56) + m(4) * (t11 ^ 2 + t26 ^ 2 + t56) + m(3) * (t44 * pkin(5) ^ 2 + pkin(1) ^ 2) + (mrSges(4,3) * t55 + Ifges(4,1) * t16 - t35 * t3 + t37 * t4 + (-Ifges(5,6) * t35 - (2 * Ifges(4,4))) * t15) * t16; -Ifges(4,6) * t15 + Ifges(3,5) * t36 + Ifges(3,6) * t38 + t9 * t18 + t15 * t46 / 0.2e1 + t35 * t4 / 0.2e1 + t3 * t53 - t11 * mrSges(4,2) - t9 * mrSges(4,1) + (-t36 * mrSges(3,1) - t38 * mrSges(3,2)) * pkin(5) + t42 * mrSges(5,3) + (Ifges(4,5) + t21 * t53 - t35 * t20 / 0.2e1) * t16 + (m(4) * (t11 * t33 - t34 * t9) + (-t33 * t15 - t34 * t16) * mrSges(4,3)) * pkin(2) + (m(5) * t9 + t5) * t25 + (m(5) * t42 - t35 * t8 + t37 * t7) * t24; 0.2e1 * t25 * t18 + t37 * t20 + t35 * t21 + Ifges(3,3) + Ifges(4,3) + m(5) * (t45 * t24 ^ 2 + t25 ^ 2) + m(4) * (t33 ^ 2 + t34 ^ 2) * pkin(2) ^ 2 + 0.2e1 * (t34 * mrSges(4,1) - t33 * mrSges(4,2)) * pkin(2) + 0.2e1 * t45 * t24 * mrSges(5,3); t15 * mrSges(4,1) + t35 * t7 + t37 * t8 + t13 + m(5) * (t37 * t1 + t35 * t2) + m(4) * t26; 0; m(5) * t45 + m(4); t1 * mrSges(5,1) - t2 * mrSges(5,2) - Ifges(5,6) * t50 + t47; -t24 * t41 + t46; -t18; Ifges(5,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t10(1), t10(2), t10(4), t10(7); t10(2), t10(3), t10(5), t10(8); t10(4), t10(5), t10(6), t10(9); t10(7), t10(8), t10(9), t10(10);];
Mq = res;
