% Calculate joint inertia matrix for
% S4RPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
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
% Datum: 2019-12-31 16:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRP6_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP6_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP6_inertiaJ_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP6_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP6_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP6_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:45:59
% EndTime: 2019-12-31 16:46:00
% DurationCPUTime: 0.14s
% Computational Cost: add. (87->51), mult. (150->60), div. (0->0), fcn. (69->2), ass. (0->16)
t21 = -2 * mrSges(5,3);
t20 = 2 * qJ(2);
t19 = m(5) * pkin(3);
t10 = sin(qJ(3));
t11 = cos(qJ(3));
t18 = t10 * mrSges(5,1) + t11 * mrSges(5,2);
t17 = t10 ^ 2 + t11 ^ 2;
t12 = -pkin(1) - pkin(5);
t16 = -qJ(4) + t12;
t15 = t17 * mrSges(4,3);
t14 = mrSges(5,1) + t19;
t13 = qJ(2) ^ 2;
t4 = t10 * pkin(3) + qJ(2);
t2 = t16 * t11;
t1 = t16 * t10;
t3 = [Ifges(3,1) + Ifges(2,3) + 0.2e1 * t4 * t18 + (mrSges(3,3) * t20) - (2 * pkin(1) * mrSges(3,2)) + m(4) * (t17 * t12 ^ 2 + t13) + m(5) * (t1 ^ 2 + t2 ^ 2 + t4 ^ 2) + (m(3) * (pkin(1) ^ 2 + t13)) + ((mrSges(4,2) * t20) + t2 * t21 + (Ifges(5,1) + Ifges(4,1)) * t11) * t11 - 0.2e1 * t12 * t15 + ((mrSges(4,1) * t20) + t1 * t21 + 0.2e1 * (-Ifges(4,4) - Ifges(5,4)) * t11 + (Ifges(5,2) + Ifges(4,2)) * t10) * t10; -(m(3) * pkin(1)) + mrSges(3,2) - t15 + m(5) * (t10 * t1 + t11 * t2) + (m(4) * t12 - mrSges(5,3)) * t17; m(3) + 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1) * t17; -t1 * mrSges(5,2) + t14 * t2 + (-mrSges(4,2) * t12 - Ifges(4,6) - Ifges(5,6)) * t10 + (mrSges(4,1) * t12 - mrSges(5,3) * pkin(3) + Ifges(4,5) + Ifges(5,5)) * t11; (-mrSges(4,2) - mrSges(5,2)) * t10 + (mrSges(4,1) + t14) * t11; Ifges(4,3) + Ifges(5,3) + (0.2e1 * mrSges(5,1) + t19) * pkin(3); m(5) * t4 + t18; 0; 0; m(5);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t3(1), t3(2), t3(4), t3(7); t3(2), t3(3), t3(5), t3(8); t3(4), t3(5), t3(6), t3(9); t3(7), t3(8), t3(9), t3(10);];
Mq = res;
