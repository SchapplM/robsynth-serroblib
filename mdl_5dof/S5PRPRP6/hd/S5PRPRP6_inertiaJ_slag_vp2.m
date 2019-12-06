% Calculate joint inertia matrix for
% S5PRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1]';
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
% Datum: 2019-12-05 15:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRP6_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP6_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRPRP6_inertiaJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP6_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP6_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP6_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:40:21
% EndTime: 2019-12-05 15:40:22
% DurationCPUTime: 0.20s
% Computational Cost: add. (131->66), mult. (285->78), div. (0->0), fcn. (130->4), ass. (0->26)
t17 = sin(qJ(4));
t12 = t17 ^ 2;
t19 = cos(qJ(4));
t14 = t19 ^ 2;
t31 = -t14 - t12;
t35 = (-mrSges(5,3) - mrSges(6,2)) * t31;
t38 = t35 - mrSges(4,2);
t30 = m(5) / 0.2e1 + m(6) / 0.2e1;
t37 = 0.2e1 * t30;
t23 = m(6) * (pkin(4) * t19 + qJ(5) * t17) + (mrSges(5,1) + mrSges(6,1)) * t19 + (-mrSges(5,2) + mrSges(6,3)) * t17;
t34 = m(6) * t19;
t21 = -pkin(2) - pkin(6);
t33 = t31 * t21 ^ 2;
t32 = t17 * mrSges(5,1) + t19 * mrSges(5,2) + mrSges(4,3);
t20 = cos(qJ(2));
t29 = t31 * t20;
t27 = t21 * t29;
t25 = -0.2e1 * t30 * t31;
t22 = qJ(3) ^ 2;
t18 = sin(qJ(2));
t15 = t20 ^ 2;
t13 = t18 ^ 2;
t11 = t18 * qJ(3);
t5 = t17 * mrSges(6,1) - t19 * mrSges(6,3);
t4 = t17 * pkin(4) - t19 * qJ(5) + qJ(3);
t1 = [m(2) + 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1) * (t13 + t15) + (-t31 * t15 + t13) * t37; (-mrSges(3,2) + t5 + t32) * t18 + (mrSges(3,1) + t38) * t20 + m(4) * (t20 * pkin(2) + t11) + m(5) * (t11 + t27) + m(6) * (t4 * t18 + t27); -0.2e1 * pkin(2) * mrSges(4,2) + 0.2e1 * t4 * t5 + Ifges(4,1) + Ifges(3,3) + (Ifges(5,1) + Ifges(6,1)) * t14 + m(6) * (t4 ^ 2 - t33) + m(5) * (t22 - t33) + m(4) * (pkin(2) ^ 2 + t22) + (Ifges(6,3) + Ifges(5,2)) * t12 + 0.2e1 * t32 * qJ(3) + 0.2e1 * (-Ifges(5,4) + Ifges(6,5)) * t19 * t17 - 0.2e1 * t21 * t35; -m(4) * t20 + t29 * t37; -m(4) * pkin(2) + t21 * t25 - t38; m(4) + t25; -t23 * t20; (-pkin(4) * mrSges(6,2) + Ifges(6,4) + Ifges(5,5)) * t19 + (-qJ(5) * mrSges(6,2) - Ifges(5,6) + Ifges(6,6)) * t17 + t23 * t21; t23; Ifges(6,2) + Ifges(5,3) + 0.2e1 * pkin(4) * mrSges(6,1) + 0.2e1 * qJ(5) * mrSges(6,3) + m(6) * (pkin(4) ^ 2 + qJ(5) ^ 2); t20 * t34; (-m(6) * t21 + mrSges(6,2)) * t19; -t34; -m(6) * pkin(4) - mrSges(6,1); m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
