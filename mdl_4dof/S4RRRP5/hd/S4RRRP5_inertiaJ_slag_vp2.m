% Calculate joint inertia matrix for
% S4RRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% Datum: 2019-12-31 17:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRP5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP5_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP5_inertiaJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP5_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP5_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP5_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:16:38
% EndTime: 2019-12-31 17:16:39
% DurationCPUTime: 0.23s
% Computational Cost: add. (219->85), mult. (431->115), div. (0->0), fcn. (337->4), ass. (0->25)
t38 = 2 * mrSges(5,1);
t37 = 2 * mrSges(5,3);
t25 = cos(qJ(2));
t19 = -t25 * pkin(2) - pkin(1);
t36 = 0.2e1 * t19;
t35 = -pkin(6) - pkin(5);
t34 = mrSges(5,2) + mrSges(4,3);
t33 = Ifges(5,2) + Ifges(4,3);
t23 = sin(qJ(2));
t32 = t23 ^ 2 + t25 ^ 2;
t15 = t35 * t25;
t22 = sin(qJ(3));
t24 = cos(qJ(3));
t30 = t35 * t23;
t5 = -t22 * t15 - t24 * t30;
t7 = -t24 * t15 + t22 * t30;
t31 = t5 ^ 2 + t7 ^ 2;
t29 = (t24 * mrSges(4,1) - t22 * mrSges(4,2)) * pkin(2);
t12 = t22 * t23 - t24 * t25;
t13 = t22 * t25 + t24 * t23;
t28 = (-mrSges(4,2) + mrSges(5,3)) * t7 + (-mrSges(4,1) - mrSges(5,1)) * t5 + (Ifges(5,4) + Ifges(4,5)) * t13 + (-Ifges(4,6) + Ifges(5,6)) * t12;
t18 = -t24 * pkin(2) - pkin(3);
t16 = t22 * pkin(2) + qJ(4);
t1 = t12 * pkin(3) - t13 * qJ(4) + t19;
t2 = [Ifges(2,3) - 0.2e1 * pkin(1) * (-t25 * mrSges(3,1) + t23 * mrSges(3,2)) + t23 * (Ifges(3,1) * t23 + Ifges(3,4) * t25) + t25 * (Ifges(3,4) * t23 + Ifges(3,2) * t25) + 0.2e1 * t32 * pkin(5) * mrSges(3,3) + m(5) * (t1 ^ 2 + t31) + m(4) * (t19 ^ 2 + t31) + m(3) * (t32 * pkin(5) ^ 2 + pkin(1) ^ 2) + (mrSges(4,2) * t36 - 0.2e1 * t1 * mrSges(5,3) + (Ifges(5,1) + Ifges(4,1)) * t13 + 0.2e1 * t34 * t5) * t13 + (mrSges(4,1) * t36 + t1 * t38 + (Ifges(4,2) + Ifges(5,3)) * t12 - 0.2e1 * t34 * t7 + 0.2e1 * (-Ifges(4,4) + Ifges(5,5)) * t13) * t12; m(5) * (t16 * t7 + t18 * t5) + Ifges(3,5) * t23 + Ifges(3,6) * t25 + (-t23 * mrSges(3,1) - t25 * mrSges(3,2)) * pkin(5) + (-t16 * t12 + t18 * t13) * mrSges(5,2) + (m(4) * (t22 * t7 - t24 * t5) + (-t22 * t12 - t24 * t13) * mrSges(4,3)) * pkin(2) + t28; -0.2e1 * t18 * mrSges(5,1) + t16 * t37 + Ifges(3,3) + 0.2e1 * t29 + m(5) * (t16 ^ 2 + t18 ^ 2) + m(4) * (t22 ^ 2 + t24 ^ 2) * pkin(2) ^ 2 + t33; m(5) * (-pkin(3) * t5 + qJ(4) * t7) + (-pkin(3) * t13 - qJ(4) * t12) * mrSges(5,2) + t28; m(5) * (-pkin(3) * t18 + qJ(4) * t16) + t29 + (t16 + qJ(4)) * mrSges(5,3) + (-t18 + pkin(3)) * mrSges(5,1) + t33; pkin(3) * t38 + qJ(4) * t37 + m(5) * (pkin(3) ^ 2 + qJ(4) ^ 2) + t33; m(5) * t5 + t13 * mrSges(5,2); m(5) * t18 - mrSges(5,1); -m(5) * pkin(3) - mrSges(5,1); m(5);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t2(1), t2(2), t2(4), t2(7); t2(2), t2(3), t2(5), t2(8); t2(4), t2(5), t2(6), t2(9); t2(7), t2(8), t2(9), t2(10);];
Mq = res;
