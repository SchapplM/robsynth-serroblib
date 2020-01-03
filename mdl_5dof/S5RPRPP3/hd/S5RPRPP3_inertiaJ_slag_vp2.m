% Calculate joint inertia matrix for
% S5RPRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
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
% Datum: 2019-12-31 18:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPP3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP3_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP3_inertiaJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP3_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP3_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPP3_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:12:21
% EndTime: 2019-12-31 18:12:21
% DurationCPUTime: 0.28s
% Computational Cost: add. (284->97), mult. (535->111), div. (0->0), fcn. (469->4), ass. (0->33)
t23 = cos(pkin(7));
t21 = t23 ^ 2;
t40 = -2 * mrSges(5,2);
t18 = -t23 * pkin(2) - pkin(1);
t39 = 0.2e1 * t18;
t38 = m(5) + m(6);
t37 = cos(qJ(3));
t36 = mrSges(4,1) - mrSges(5,2);
t35 = mrSges(5,1) + mrSges(6,1);
t34 = mrSges(5,1) + mrSges(4,3);
t33 = pkin(3) + qJ(5);
t32 = pkin(6) + qJ(2);
t22 = sin(pkin(7));
t31 = t22 ^ 2 + t21;
t15 = t32 * t22;
t16 = t32 * t23;
t25 = sin(qJ(3));
t5 = t37 * t15 + t25 * t16;
t7 = -t25 * t15 + t37 * t16;
t30 = t5 ^ 2 + t7 ^ 2;
t29 = -t23 * mrSges(3,1) + t22 * mrSges(3,2);
t14 = t37 * t22 + t25 * t23;
t28 = -t14 * qJ(4) + t18;
t26 = qJ(4) ^ 2;
t13 = t25 * t22 - t37 * t23;
t11 = t14 * mrSges(5,3);
t10 = t14 * mrSges(4,2);
t9 = t13 * mrSges(6,3);
t4 = t13 * pkin(3) + t28;
t3 = -t13 * pkin(4) + t7;
t2 = t14 * pkin(4) + t5;
t1 = t33 * t13 + t28;
t6 = [Ifges(2,3) - 0.2e1 * t4 * t11 + 0.2e1 * t1 * t9 + t10 * t39 - 0.2e1 * pkin(1) * t29 + Ifges(3,2) * t21 + (Ifges(3,1) * t22 + 0.2e1 * Ifges(3,4) * t23) * t22 + 0.2e1 * t31 * qJ(2) * mrSges(3,3) + m(6) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(5) * (t4 ^ 2 + t30) + m(4) * (t18 ^ 2 + t30) + m(3) * (t31 * qJ(2) ^ 2 + pkin(1) ^ 2) + (mrSges(4,1) * t39 - 0.2e1 * t3 * mrSges(6,1) + t4 * t40 + (Ifges(6,2) + Ifges(4,2) + Ifges(5,3)) * t13 - 0.2e1 * t34 * t7) * t13 + (0.2e1 * t2 * mrSges(6,1) - 0.2e1 * t1 * mrSges(6,2) + (Ifges(6,3) + Ifges(4,1) + Ifges(5,2)) * t14 + 0.2e1 * t34 * t5 + 0.2e1 * (-Ifges(4,4) - Ifges(5,6) + Ifges(6,6)) * t13) * t14; -m(3) * pkin(1) + m(4) * t18 + m(5) * t4 + m(6) * t1 - t14 * mrSges(6,2) + t36 * t13 + t10 - t11 + t29 + t9; m(3) + m(4) + t38; t3 * mrSges(6,2) - t2 * mrSges(6,3) + (-mrSges(4,2) + mrSges(5,3)) * t7 - t36 * t5 + m(6) * (qJ(4) * t3 - t2 * t33) + m(5) * (-pkin(3) * t5 + qJ(4) * t7) + (-pkin(3) * mrSges(5,1) - mrSges(6,1) * t33 - Ifges(5,4) + Ifges(4,5) + Ifges(6,5)) * t14 + (-t35 * qJ(4) + Ifges(6,4) + Ifges(5,5) - Ifges(4,6)) * t13; 0; pkin(3) * t40 + 0.2e1 * t33 * mrSges(6,3) + Ifges(5,1) + Ifges(6,1) + Ifges(4,3) + 0.2e1 * (mrSges(5,3) + mrSges(6,2)) * qJ(4) + m(5) * (pkin(3) ^ 2 + t26) + m(6) * (t33 ^ 2 + t26); m(5) * t5 + m(6) * t2 + t35 * t14; 0; -m(5) * pkin(3) - m(6) * t33 + mrSges(5,2) - mrSges(6,3); t38; m(6) * t3 - t13 * mrSges(6,1); 0; m(6) * qJ(4) + mrSges(6,2); 0; m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t6(1), t6(2), t6(4), t6(7), t6(11); t6(2), t6(3), t6(5), t6(8), t6(12); t6(4), t6(5), t6(6), t6(9), t6(13); t6(7), t6(8), t6(9), t6(10), t6(14); t6(11), t6(12), t6(13), t6(14), t6(15);];
Mq = res;
