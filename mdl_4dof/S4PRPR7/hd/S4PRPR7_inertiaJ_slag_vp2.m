% Calculate joint inertia matrix for
% S4PRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
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
% Datum: 2019-12-31 16:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRPR7_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR7_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR7_inertiaJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR7_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR7_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR7_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:25:37
% EndTime: 2019-12-31 16:25:37
% DurationCPUTime: 0.12s
% Computational Cost: add. (59->37), mult. (135->49), div. (0->0), fcn. (60->4), ass. (0->18)
t11 = cos(qJ(4));
t9 = sin(qJ(4));
t4 = t9 ^ 2;
t19 = t11 ^ 2 + t4;
t17 = m(5) * t19;
t21 = m(4) + t17;
t20 = t19 * mrSges(5,3) - mrSges(4,2);
t18 = t9 * mrSges(5,1) + t11 * mrSges(5,2) + mrSges(4,3);
t13 = -pkin(2) - pkin(5);
t16 = t19 * t13;
t15 = t11 * mrSges(5,1) - t9 * mrSges(5,2);
t14 = qJ(3) ^ 2;
t12 = cos(qJ(2));
t10 = sin(qJ(2));
t7 = t12 ^ 2;
t5 = t10 ^ 2;
t3 = t10 * qJ(3);
t1 = [m(2) + m(5) * (t19 * t7 + t5) + 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1) * (t5 + t7); (mrSges(3,1) + t20) * t12 + m(4) * (t12 * pkin(2) + t3) + m(5) * (-t12 * t16 + t3) + (-mrSges(3,2) + t18) * t10; t4 * Ifges(5,2) - 0.2e1 * pkin(2) * mrSges(4,2) + Ifges(4,1) + Ifges(3,3) + (Ifges(5,1) * t11 - 0.2e1 * Ifges(5,4) * t9) * t11 + m(5) * (t19 * t13 ^ 2 + t14) + m(4) * (pkin(2) ^ 2 + t14) - 0.2e1 * mrSges(5,3) * t16 + 0.2e1 * t18 * qJ(3); -t21 * t12; -m(4) * pkin(2) + t13 * t17 - t20; t21; -t15 * t12; (-mrSges(5,2) * t13 - Ifges(5,6)) * t9 + (mrSges(5,1) * t13 + Ifges(5,5)) * t11; t15; Ifges(5,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
