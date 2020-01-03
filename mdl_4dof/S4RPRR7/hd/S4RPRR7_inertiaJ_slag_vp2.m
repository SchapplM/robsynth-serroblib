% Calculate joint inertia matrix for
% S4RPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 16:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRR7_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR7_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR7_inertiaJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR7_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR7_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR7_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:53:43
% EndTime: 2019-12-31 16:53:43
% DurationCPUTime: 0.25s
% Computational Cost: add. (296->89), mult. (587->136), div. (0->0), fcn. (563->6), ass. (0->43)
t32 = cos(pkin(7));
t46 = pkin(5) + qJ(2);
t18 = t46 * t32;
t34 = sin(qJ(3));
t31 = sin(pkin(7));
t41 = t46 * t31;
t51 = cos(qJ(3));
t9 = t34 * t18 + t51 * t41;
t55 = t9 ^ 2;
t54 = 0.2e1 * t9;
t28 = t32 ^ 2;
t23 = -t32 * pkin(2) - pkin(1);
t53 = 0.2e1 * t23;
t35 = cos(qJ(4));
t52 = t35 / 0.2e1;
t33 = sin(qJ(4));
t50 = Ifges(5,4) * t33;
t49 = Ifges(5,4) * t35;
t17 = t51 * t31 + t34 * t32;
t48 = t17 * t33;
t47 = t17 * t35;
t16 = t34 * t31 - t51 * t32;
t45 = Ifges(5,5) * t47 + Ifges(5,3) * t16;
t44 = Ifges(5,5) * t33 + Ifges(5,6) * t35;
t43 = t31 ^ 2 + t28;
t42 = t33 ^ 2 + t35 ^ 2;
t40 = -t32 * mrSges(3,1) + t31 * mrSges(3,2);
t11 = t51 * t18 - t34 * t41;
t6 = t16 * pkin(3) - t17 * pkin(6) + t23;
t1 = -t33 * t11 + t35 * t6;
t2 = t35 * t11 + t33 * t6;
t39 = -t1 * t33 + t2 * t35;
t38 = mrSges(5,1) * t33 + mrSges(5,2) * t35;
t21 = Ifges(5,1) * t33 + t49;
t20 = Ifges(5,2) * t35 + t50;
t19 = -t35 * mrSges(5,1) + t33 * mrSges(5,2);
t13 = t17 * mrSges(4,2);
t8 = t16 * mrSges(5,1) - mrSges(5,3) * t47;
t7 = -t16 * mrSges(5,2) - mrSges(5,3) * t48;
t5 = t38 * t17;
t4 = Ifges(5,5) * t16 + (Ifges(5,1) * t35 - t50) * t17;
t3 = Ifges(5,6) * t16 + (-Ifges(5,2) * t33 + t49) * t17;
t10 = [Ifges(2,3) + 0.2e1 * t2 * t7 + 0.2e1 * t1 * t8 + t5 * t54 + t13 * t53 + Ifges(3,2) * t28 - 0.2e1 * pkin(1) * t40 + (Ifges(3,1) * t31 + 0.2e1 * Ifges(3,4) * t32) * t31 + (mrSges(4,1) * t53 - 0.2e1 * t11 * mrSges(4,3) + Ifges(4,2) * t16 + t45) * t16 + 0.2e1 * t43 * qJ(2) * mrSges(3,3) + m(5) * (t1 ^ 2 + t2 ^ 2 + t55) + m(4) * (t11 ^ 2 + t23 ^ 2 + t55) + m(3) * (t43 * qJ(2) ^ 2 + pkin(1) ^ 2) + (mrSges(4,3) * t54 + Ifges(4,1) * t17 - t33 * t3 + t35 * t4 + (-Ifges(5,6) * t33 - (2 * Ifges(4,4))) * t16) * t17; -m(3) * pkin(1) + t16 * mrSges(4,1) + t33 * t7 + t35 * t8 + t13 + m(5) * (t35 * t1 + t33 * t2) + m(4) * t23 + t40; m(5) * t42 + m(3) + m(4); t33 * t4 / 0.2e1 + t3 * t52 - pkin(3) * t5 - t11 * mrSges(4,2) + t39 * mrSges(5,3) + (t21 * t52 - t33 * t20 / 0.2e1 + Ifges(4,5)) * t17 + (t44 / 0.2e1 - Ifges(4,6)) * t16 + (-m(5) * pkin(3) - mrSges(4,1) + t19) * t9 + (m(5) * t39 - t33 * t8 + t35 * t7) * pkin(6); 0; Ifges(4,3) + t35 * t20 - 0.2e1 * pkin(3) * t19 + t33 * t21 + m(5) * (t42 * pkin(6) ^ 2 + pkin(3) ^ 2) + 0.2e1 * t42 * pkin(6) * mrSges(5,3); t1 * mrSges(5,1) - t2 * mrSges(5,2) - Ifges(5,6) * t48 + t45; -t19; -t38 * pkin(6) + t44; Ifges(5,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t10(1), t10(2), t10(4), t10(7); t10(2), t10(3), t10(5), t10(8); t10(4), t10(5), t10(6), t10(9); t10(7), t10(8), t10(9), t10(10);];
Mq = res;
