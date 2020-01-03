% Calculate joint inertia matrix for
% S5PRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1]';
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
% Datum: 2019-12-31 17:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPP4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP4_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRRPP4_inertiaJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP4_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPP4_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPP4_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:40:47
% EndTime: 2019-12-31 17:40:47
% DurationCPUTime: 0.20s
% Computational Cost: add. (119->79), mult. (242->83), div. (0->0), fcn. (110->2), ass. (0->22)
t14 = sin(qJ(3));
t15 = cos(qJ(3));
t32 = t14 ^ 2 + t15 ^ 2;
t31 = 0.2e1 * t32;
t30 = -mrSges(4,2) + mrSges(5,3);
t29 = -m(5) * pkin(3) - mrSges(5,1);
t6 = t14 * qJ(4);
t23 = t15 * pkin(3) + t6;
t2 = -pkin(2) - t23;
t28 = -0.2e1 * t2;
t27 = -2 * mrSges(6,3);
t26 = m(5) + m(6);
t25 = t15 * mrSges(6,1) + t14 * mrSges(6,2);
t24 = t32 * pkin(6) ^ 2;
t22 = mrSges(5,2) - mrSges(6,3);
t21 = pkin(6) - qJ(5);
t17 = qJ(4) ^ 2;
t16 = -pkin(3) - pkin(4);
t4 = t21 * t15;
t3 = t21 * t14;
t1 = t15 * pkin(4) - t2;
t5 = [m(2) + m(3) + (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * t31; m(6) * (t4 * t14 - t3 * t15); Ifges(3,3) + 0.2e1 * t1 * t25 + m(5) * (t2 ^ 2 + t24) + m(6) * (t1 ^ 2 + t3 ^ 2 + t4 ^ 2) + m(4) * (pkin(2) ^ 2 + t24) + (0.2e1 * pkin(2) * mrSges(4,1) + mrSges(5,1) * t28 + t4 * t27 + (Ifges(5,3) + Ifges(4,2) + Ifges(6,2)) * t15) * t15 + (-0.2e1 * pkin(2) * mrSges(4,2) + mrSges(5,3) * t28 + t3 * t27 + (Ifges(4,1) + Ifges(6,1) + Ifges(5,1)) * t14 + 0.2e1 * (Ifges(4,4) - Ifges(6,4) - Ifges(5,5)) * t15) * t14 + (mrSges(5,2) + mrSges(4,3)) * pkin(6) * t31; (mrSges(4,1) + mrSges(5,1)) * t15 + t30 * t14 + m(5) * t23 + m(6) * (-t16 * t15 + t6) + t25; m(6) * (qJ(4) * t4 + t16 * t3) - t3 * mrSges(6,1) + t4 * mrSges(6,2) + (t22 * qJ(4) + Ifges(4,6) - Ifges(5,6) + Ifges(6,6)) * t15 + (-pkin(3) * mrSges(5,2) - t16 * mrSges(6,3) + Ifges(5,4) + Ifges(4,5) - Ifges(6,5)) * t14 + ((m(5) * qJ(4) + t30) * t15 + (-mrSges(4,1) + t29) * t14) * pkin(6); 0.2e1 * pkin(3) * mrSges(5,1) - 0.2e1 * t16 * mrSges(6,1) + Ifges(5,2) + Ifges(4,3) + Ifges(6,3) + 0.2e1 * (mrSges(5,3) + mrSges(6,2)) * qJ(4) + m(5) * (pkin(3) ^ 2 + t17) + m(6) * (t16 ^ 2 + t17); -t26 * t15; m(6) * t3 + (m(5) * pkin(6) + t22) * t14; m(6) * t16 - mrSges(6,1) + t29; t26; 0; m(6) * t1 + t25; 0; 0; m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t5(1), t5(2), t5(4), t5(7), t5(11); t5(2), t5(3), t5(5), t5(8), t5(12); t5(4), t5(5), t5(6), t5(9), t5(13); t5(7), t5(8), t5(9), t5(10), t5(14); t5(11), t5(12), t5(13), t5(14), t5(15);];
Mq = res;
