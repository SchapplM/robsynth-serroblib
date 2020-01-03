% Calculate joint inertia matrix for
% S4RPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% Datum: 2019-12-31 16:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRR9_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR9_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR9_inertiaJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR9_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR9_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR9_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:56:09
% EndTime: 2019-12-31 16:56:09
% DurationCPUTime: 0.21s
% Computational Cost: add. (169->85), mult. (342->121), div. (0->0), fcn. (221->4), ass. (0->42)
t27 = (-pkin(1) - pkin(5));
t47 = -2 * t27;
t23 = sin(qJ(4));
t25 = cos(qJ(4));
t9 = -t25 * mrSges(5,1) + t23 * mrSges(5,2);
t46 = m(5) * pkin(3) + mrSges(4,1) - t9;
t45 = 2 * qJ(2);
t44 = t25 / 0.2e1;
t43 = Ifges(5,4) * t23;
t42 = Ifges(5,4) * t25;
t24 = sin(qJ(3));
t41 = Ifges(5,6) * t24;
t26 = cos(qJ(3));
t40 = t23 * t26;
t39 = t24 * t27;
t38 = t25 * t26;
t37 = Ifges(5,5) * t38 + Ifges(5,3) * t24;
t36 = t23 ^ 2 + t25 ^ 2;
t19 = t24 ^ 2;
t21 = t26 ^ 2;
t35 = t21 + t19;
t33 = t35 * mrSges(4,3);
t8 = t24 * pkin(3) - t26 * pkin(6) + qJ(2);
t1 = -t23 * t39 + t25 * t8;
t2 = t23 * t8 + t25 * t39;
t32 = -t1 * t23 + t2 * t25;
t6 = -t24 * mrSges(5,2) - mrSges(5,3) * t40;
t7 = t24 * mrSges(5,1) - mrSges(5,3) * t38;
t31 = -t23 * t7 + t25 * t6;
t30 = mrSges(5,1) * t23 + mrSges(5,2) * t25;
t28 = qJ(2) ^ 2;
t22 = t27 ^ 2;
t17 = Ifges(5,5) * t23;
t16 = Ifges(5,6) * t25;
t14 = t21 * t27;
t13 = t21 * t22;
t11 = Ifges(5,1) * t23 + t42;
t10 = Ifges(5,2) * t25 + t43;
t5 = t30 * t26;
t4 = Ifges(5,5) * t24 + (Ifges(5,1) * t25 - t43) * t26;
t3 = t41 + (-Ifges(5,2) * t23 + t42) * t26;
t12 = [-(2 * pkin(1) * mrSges(3,2)) + (mrSges(3,3) * t45) + 0.2e1 * t1 * t7 + 0.2e1 * t2 * t6 + Ifges(3,1) + Ifges(2,3) + t33 * t47 + (mrSges(4,1) * t45 + Ifges(4,2) * t24 + t37) * t24 + ((mrSges(4,2) * t45) + Ifges(4,1) * t26 - 0.2e1 * Ifges(4,4) * t24 + t25 * t4 + t5 * t47 + (-t3 - t41) * t23) * t26 + m(5) * (t1 ^ 2 + t2 ^ 2 + t13) + m(4) * (t19 * t22 + t13 + t28) + (m(3) * (pkin(1) ^ 2 + t28)); -(m(3) * pkin(1)) - t26 * t5 + mrSges(3,2) + t31 * t24 - t33 + m(5) * (t32 * t24 + t14) + m(4) * (t19 * t27 + t14); m(3) + m(4) * t35 + m(5) * (t36 * t19 + t21); -pkin(3) * t5 + t23 * t4 / 0.2e1 + t3 * t44 + (-(t27 * mrSges(4,2)) + t17 / 0.2e1 + t16 / 0.2e1 - Ifges(4,6)) * t24 + t32 * mrSges(5,3) + (m(5) * t32 + t31) * pkin(6) + (t11 * t44 - t23 * t10 / 0.2e1 + Ifges(4,5) + t46 * t27) * t26; t46 * t26 + (-mrSges(4,2) + (m(5) * pkin(6) + mrSges(5,3)) * t36) * t24; Ifges(4,3) - 0.2e1 * pkin(3) * t9 + m(5) * (t36 * pkin(6) ^ 2 + pkin(3) ^ 2) + t23 * t11 + t25 * t10 + 0.2e1 * t36 * pkin(6) * mrSges(5,3); t1 * mrSges(5,1) - t2 * mrSges(5,2) - Ifges(5,6) * t40 + t37; -t30 * t24; -t30 * pkin(6) + t16 + t17; Ifges(5,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t12(1), t12(2), t12(4), t12(7); t12(2), t12(3), t12(5), t12(8); t12(4), t12(5), t12(6), t12(9); t12(7), t12(8), t12(9), t12(10);];
Mq = res;
