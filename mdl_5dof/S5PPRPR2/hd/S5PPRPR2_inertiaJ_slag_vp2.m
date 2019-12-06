% Calculate joint inertia matrix for
% S5PPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRPR2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR2_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR2_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR2_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR2_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRPR2_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:03:00
% EndTime: 2019-12-05 15:03:01
% DurationCPUTime: 0.12s
% Computational Cost: add. (90->43), mult. (200->56), div. (0->0), fcn. (138->6), ass. (0->24)
t13 = sin(pkin(8));
t14 = cos(pkin(8));
t16 = sin(qJ(3));
t26 = cos(qJ(3));
t3 = t16 * t13 - t26 * t14;
t1 = t3 ^ 2;
t5 = t26 * t13 + t16 * t14;
t2 = t5 ^ 2;
t17 = cos(qJ(5));
t11 = t17 ^ 2;
t15 = sin(qJ(5));
t24 = t15 ^ 2 + t11;
t7 = m(6) * t24;
t28 = m(5) + t7;
t8 = t15 * mrSges(6,1) + t17 * mrSges(6,2);
t27 = t8 + mrSges(5,3);
t25 = t5 * qJ(4);
t23 = t24 * mrSges(6,3);
t18 = -pkin(3) - pkin(6);
t22 = t24 * t18;
t21 = mrSges(5,2) - t23;
t20 = t17 * mrSges(6,1) - t15 * mrSges(6,2);
t19 = qJ(4) ^ 2;
t4 = [m(2) + m(3) * (t13 ^ 2 + t14 ^ 2) + m(6) * (t24 * t1 + t2) + (m(4) + m(5)) * (t1 + t2); 0; m(3) + m(4) + t28; (-mrSges(4,2) + t27) * t5 + (-mrSges(4,1) + t21) * t3 + m(5) * (-pkin(3) * t3 + t25) + m(6) * (t3 * t22 + t25); 0; Ifges(6,1) * t11 - 0.2e1 * pkin(3) * mrSges(5,2) + Ifges(5,1) + Ifges(4,3) + (-0.2e1 * Ifges(6,4) * t17 + Ifges(6,2) * t15) * t15 + m(6) * (t24 * t18 ^ 2 + t19) + m(5) * (pkin(3) ^ 2 + t19) + 0.2e1 * t27 * qJ(4) - 0.2e1 * t18 * t23; 0.2e1 * (m(5) / 0.2e1 + t7 / 0.2e1) * t3; 0; -m(5) * pkin(3) + m(6) * t22 + t21; t28; t20 * t3; -t8; (mrSges(6,1) * t18 + Ifges(6,5)) * t17 + (-mrSges(6,2) * t18 - Ifges(6,6)) * t15; t20; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t4(1), t4(2), t4(4), t4(7), t4(11); t4(2), t4(3), t4(5), t4(8), t4(12); t4(4), t4(5), t4(6), t4(9), t4(13); t4(7), t4(8), t4(9), t4(10), t4(14); t4(11), t4(12), t4(13), t4(14), t4(15);];
Mq = res;
