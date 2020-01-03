% Calculate joint inertia matrix for
% S5RPPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
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
% Datum: 2019-12-31 18:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR11_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR11_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR11_inertiaJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR11_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR11_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR11_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:05:33
% EndTime: 2019-12-31 18:05:33
% DurationCPUTime: 0.25s
% Computational Cost: add. (219->104), mult. (395->137), div. (0->0), fcn. (247->4), ass. (0->44)
t23 = (qJ(2) - pkin(6));
t50 = -2 * t23;
t26 = sin(qJ(5));
t28 = cos(qJ(5));
t9 = -t28 * mrSges(6,1) + t26 * mrSges(6,2);
t49 = m(6) * pkin(4) + mrSges(5,1) - t9;
t24 = (pkin(1) + qJ(3));
t48 = t24 ^ 2;
t47 = 2 * t24;
t46 = t28 / 0.2e1;
t45 = Ifges(6,4) * t26;
t44 = Ifges(6,4) * t28;
t27 = sin(qJ(4));
t43 = Ifges(6,6) * t27;
t42 = t23 * t27;
t29 = cos(qJ(4));
t41 = t26 * t29;
t40 = t28 * t29;
t39 = Ifges(6,5) * t40 + Ifges(6,3) * t27;
t38 = t26 ^ 2 + t28 ^ 2;
t20 = t27 ^ 2;
t22 = t29 ^ 2;
t37 = t22 + t20;
t35 = t37 * mrSges(5,3);
t8 = t27 * pkin(4) - t29 * pkin(7) + t24;
t1 = -t26 * t42 + t28 * t8;
t2 = t26 * t8 + t28 * t42;
t34 = -t1 * t26 + t2 * t28;
t6 = -t27 * mrSges(6,2) - mrSges(6,3) * t41;
t7 = t27 * mrSges(6,1) - mrSges(6,3) * t40;
t33 = -t26 * t7 + t28 * t6;
t32 = mrSges(6,1) * t26 + mrSges(6,2) * t28;
t30 = qJ(2) ^ 2;
t18 = t23 ^ 2;
t17 = Ifges(6,5) * t26;
t16 = Ifges(6,6) * t28;
t14 = t22 * t23;
t13 = t22 * t18;
t11 = Ifges(6,1) * t26 + t44;
t10 = Ifges(6,2) * t28 + t45;
t5 = t32 * t29;
t4 = Ifges(6,5) * t27 + (Ifges(6,1) * t28 - t45) * t29;
t3 = t43 + (-Ifges(6,2) * t26 + t44) * t29;
t12 = [-(2 * pkin(1) * mrSges(3,2)) + (mrSges(4,3) * t47) + 0.2e1 * t1 * t7 + 0.2e1 * t2 * t6 + Ifges(3,1) + Ifges(4,1) + Ifges(2,3) + (mrSges(5,1) * t47 + Ifges(5,2) * t27 + t39) * t27 + ((mrSges(5,2) * t47) + Ifges(5,1) * t29 - 0.2e1 * Ifges(5,4) * t27 + t5 * t50 + t28 * t4 + (-t3 - t43) * t26) * t29 + m(6) * (t1 ^ 2 + t2 ^ 2 + t13) + m(5) * (t20 * t18 + t13 + t48) + m(4) * (t30 + t48) + m(3) * ((pkin(1) ^ 2) + t30) + (2 * (mrSges(3,3) + mrSges(4,2)) * qJ(2)) + t35 * t50; -m(3) * pkin(1) - t27 * mrSges(5,1) - t29 * mrSges(5,2) - t26 * t6 - t28 * t7 + mrSges(3,2) - mrSges(4,3) + m(6) * (-t28 * t1 - t26 * t2) + (-m(5) / 0.2e1 - m(4) / 0.2e1) * t47; m(6) * t38 + m(3) + m(4) + m(5); m(4) * qJ(2) - t29 * t5 + mrSges(4,2) + t33 * t27 - t35 + m(6) * (t34 * t27 + t14) + m(5) * (t20 * t23 + t14); 0; m(4) + m(5) * t37 + m(6) * (t38 * t20 + t22); t26 * t4 / 0.2e1 + t3 * t46 - pkin(4) * t5 + (-(t23 * mrSges(5,2)) + t17 / 0.2e1 + t16 / 0.2e1 - Ifges(5,6)) * t27 + t34 * mrSges(6,3) + (m(6) * t34 + t33) * pkin(7) + (t11 * t46 - t26 * t10 / 0.2e1 + Ifges(5,5) + t49 * t23) * t29; 0; t49 * t29 + (-mrSges(5,2) + (m(6) * pkin(7) + mrSges(6,3)) * t38) * t27; Ifges(5,3) + m(6) * (t38 * pkin(7) ^ 2 + pkin(4) ^ 2) - 0.2e1 * pkin(4) * t9 + t26 * t11 + t28 * t10 + 0.2e1 * t38 * pkin(7) * mrSges(6,3); t1 * mrSges(6,1) - t2 * mrSges(6,2) - Ifges(6,6) * t41 + t39; t9; -t32 * t27; -t32 * pkin(7) + t16 + t17; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t12(1), t12(2), t12(4), t12(7), t12(11); t12(2), t12(3), t12(5), t12(8), t12(12); t12(4), t12(5), t12(6), t12(9), t12(13); t12(7), t12(8), t12(9), t12(10), t12(14); t12(11), t12(12), t12(13), t12(14), t12(15);];
Mq = res;
