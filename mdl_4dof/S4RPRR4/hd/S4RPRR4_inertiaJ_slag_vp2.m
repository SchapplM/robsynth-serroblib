% Calculate joint inertia matrix for
% S4RPRR4
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
% Datum: 2019-12-31 16:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRR4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR4_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR4_inertiaJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR4_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR4_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR4_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:50:14
% EndTime: 2019-12-31 16:50:15
% DurationCPUTime: 0.22s
% Computational Cost: add. (175->83), mult. (372->126), div. (0->0), fcn. (259->6), ass. (0->41)
t27 = sin(qJ(4));
t29 = cos(qJ(4));
t9 = -t29 * mrSges(5,1) + t27 * mrSges(5,2);
t46 = mrSges(4,1) - t9;
t45 = t29 / 0.2e1;
t30 = cos(qJ(3));
t44 = t30 * pkin(3);
t43 = Ifges(5,4) * t27;
t42 = Ifges(5,4) * t29;
t41 = Ifges(5,6) * t30;
t25 = sin(pkin(7));
t17 = t25 * pkin(1) + pkin(5);
t40 = t17 * t30;
t28 = sin(qJ(3));
t39 = t27 * t28;
t38 = t28 * t29;
t37 = t27 ^ 2 + t29 ^ 2;
t22 = t28 ^ 2;
t24 = t30 ^ 2;
t36 = t22 + t24;
t26 = cos(pkin(7));
t18 = -t26 * pkin(1) - pkin(2);
t35 = t37 * t28;
t6 = -t28 * pkin(6) + t18 - t44;
t1 = -t27 * t40 + t29 * t6;
t2 = t27 * t6 + t29 * t40;
t34 = -t1 * t27 + t2 * t29;
t7 = t30 * mrSges(5,2) - mrSges(5,3) * t39;
t8 = -t30 * mrSges(5,1) - mrSges(5,3) * t38;
t33 = -t27 * t8 + t29 * t7;
t20 = Ifges(5,5) * t27;
t19 = Ifges(5,6) * t29;
t16 = t17 ^ 2;
t15 = Ifges(5,5) * t38;
t12 = t22 * t16;
t11 = Ifges(5,1) * t27 + t42;
t10 = Ifges(5,2) * t29 + t43;
t5 = -mrSges(5,1) * t39 - mrSges(5,2) * t38;
t4 = -Ifges(5,5) * t30 + (Ifges(5,1) * t29 - t43) * t28;
t3 = -t41 + (-Ifges(5,2) * t27 + t42) * t28;
t13 = [0.2e1 * t1 * t8 + 0.2e1 * t2 * t7 + Ifges(2,3) + Ifges(3,3) + (-0.2e1 * t18 * mrSges(4,1) - t15 + (Ifges(5,3) + Ifges(4,2)) * t30) * t30 + (0.2e1 * t18 * mrSges(4,2) + Ifges(4,1) * t28 + 0.2e1 * Ifges(4,4) * t30 - 0.2e1 * t17 * t5 + t29 * t4 + (-t3 + t41) * t27) * t28 + m(5) * (t1 ^ 2 + t2 ^ 2 + t12) + m(4) * (t24 * t16 + t18 ^ 2 + t12) + m(3) * (t25 ^ 2 + t26 ^ 2) * pkin(1) ^ 2 + 0.2e1 * (t26 * mrSges(3,1) - t25 * mrSges(3,2)) * pkin(1) + 0.2e1 * t36 * t17 * mrSges(4,3); t30 * t5 + (m(5) * (t34 - t40) + t33) * t28; m(3) + m(4) * t36 + m(5) * (t37 * t22 + t24); t27 * t4 / 0.2e1 + t3 * t45 + pkin(3) * t5 + (-t17 * mrSges(4,2) - t20 / 0.2e1 - t19 / 0.2e1 + Ifges(4,6)) * t30 + t34 * mrSges(5,3) + (m(5) * t34 + t33) * pkin(6) + (t11 * t45 - t27 * t10 / 0.2e1 + Ifges(4,5) + (-m(5) * pkin(3) - t46) * t17) * t28; -t28 * mrSges(4,2) + m(5) * (pkin(6) * t35 + t44) + mrSges(5,3) * t35 + t46 * t30; Ifges(4,3) - 0.2e1 * pkin(3) * t9 + m(5) * (t37 * pkin(6) ^ 2 + pkin(3) ^ 2) + t27 * t11 + t29 * t10 + 0.2e1 * t37 * pkin(6) * mrSges(5,3); t1 * mrSges(5,1) - t2 * mrSges(5,2) - Ifges(5,6) * t39 - Ifges(5,3) * t30 + t15; t5; t20 + t19 + (-mrSges(5,1) * t27 - mrSges(5,2) * t29) * pkin(6); Ifges(5,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t13(1), t13(2), t13(4), t13(7); t13(2), t13(3), t13(5), t13(8); t13(4), t13(5), t13(6), t13(9); t13(7), t13(8), t13(9), t13(10);];
Mq = res;
