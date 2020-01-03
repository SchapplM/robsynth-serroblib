% Calculate joint inertia matrix for
% S5PRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
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
% Datum: 2019-12-31 17:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPPR4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR4_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR4_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR4_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR4_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPPR4_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:36:46
% EndTime: 2019-12-31 17:36:47
% DurationCPUTime: 0.18s
% Computational Cost: add. (134->65), mult. (283->90), div. (0->0), fcn. (216->4), ass. (0->21)
t17 = sin(pkin(8));
t18 = cos(pkin(8));
t29 = t17 ^ 2 + t18 ^ 2;
t28 = 0.2e1 * t29;
t23 = t17 * qJ(4) + pkin(2);
t8 = -t18 * pkin(3) - t23;
t27 = -0.2e1 * t8;
t26 = -pkin(6) + qJ(3);
t25 = t29 * qJ(3) ^ 2;
t19 = sin(qJ(5));
t20 = cos(qJ(5));
t6 = -t17 * t19 - t18 * t20;
t7 = t17 * t20 - t18 * t19;
t1 = -t6 * mrSges(6,1) + t7 * mrSges(6,2);
t14 = t17 * mrSges(4,2);
t10 = t26 * t18;
t9 = t26 * t17;
t4 = (pkin(3) + pkin(4)) * t18 + t23;
t3 = t20 * t10 + t19 * t9;
t2 = -t19 * t10 + t20 * t9;
t5 = [m(2) + m(3) + m(6) * (t6 ^ 2 + t7 ^ 2) + (m(4) / 0.2e1 + m(5) / 0.2e1) * t28; m(6) * (t2 * t6 + t3 * t7); -0.2e1 * pkin(2) * t14 + 0.2e1 * t4 * t1 + Ifges(3,3) + (-0.2e1 * t2 * mrSges(6,3) + Ifges(6,1) * t7) * t7 + (0.2e1 * t3 * mrSges(6,3) + 0.2e1 * Ifges(6,4) * t7 + Ifges(6,2) * t6) * t6 + (0.2e1 * pkin(2) * mrSges(4,1) + mrSges(5,1) * t27 + (Ifges(5,3) + Ifges(4,2)) * t18) * t18 + (mrSges(5,3) * t27 + (Ifges(4,1) + Ifges(5,1)) * t17 + 0.2e1 * (Ifges(4,4) - Ifges(5,5)) * t18) * t17 + m(6) * (t2 ^ 2 + t3 ^ 2 + t4 ^ 2) + m(4) * (pkin(2) ^ 2 + t25) + m(5) * (t8 ^ 2 + t25) + (mrSges(5,2) + mrSges(4,3)) * qJ(3) * t28; 0; -m(4) * pkin(2) - t17 * mrSges(5,3) + t14 + (-mrSges(5,1) - mrSges(4,1)) * t18 + m(5) * t8 - m(6) * t4 - t1; m(4) + m(5) + m(6); -m(5) * t18 + m(6) * (t19 * t7 + t20 * t6); m(6) * (t19 * t3 + t20 * t2) + (m(5) * qJ(3) + mrSges(5,2)) * t17 + (t19 * t6 - t20 * t7) * mrSges(6,3); 0; m(5) + m(6) * (t19 ^ 2 + t20 ^ 2); -t1; t2 * mrSges(6,1) - t3 * mrSges(6,2) + Ifges(6,5) * t7 + Ifges(6,6) * t6; 0; t20 * mrSges(6,1) - t19 * mrSges(6,2); Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t5(1), t5(2), t5(4), t5(7), t5(11); t5(2), t5(3), t5(5), t5(8), t5(12); t5(4), t5(5), t5(6), t5(9), t5(13); t5(7), t5(8), t5(9), t5(10), t5(14); t5(11), t5(12), t5(13), t5(14), t5(15);];
Mq = res;
