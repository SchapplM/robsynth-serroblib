% Calculate joint inertia matrix for
% S5PRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRP1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP1_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP1_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP1_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP1_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP1_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:39:49
% EndTime: 2019-12-05 16:39:50
% DurationCPUTime: 0.26s
% Computational Cost: add. (186->87), mult. (361->107), div. (0->0), fcn. (205->4), ass. (0->32)
t27 = cos(qJ(4));
t45 = t27 ^ 2;
t25 = sin(qJ(4));
t38 = t25 ^ 2 + t45;
t44 = 0.2e1 * t38;
t18 = t25 * mrSges(6,2);
t9 = -t27 * mrSges(6,1) + t18;
t43 = 0.2e1 * t9;
t42 = m(6) * pkin(4);
t28 = cos(qJ(3));
t41 = t28 * pkin(2);
t40 = t25 * mrSges(5,2);
t39 = t25 * mrSges(6,3);
t37 = 0.2e1 * mrSges(6,3);
t16 = -t27 * pkin(4) - pkin(3);
t26 = sin(qJ(3));
t14 = t26 * pkin(2) + pkin(7);
t36 = t38 * t14;
t35 = (Ifges(5,6) + Ifges(6,6)) * t27 + (Ifges(5,5) + Ifges(6,5)) * t25;
t34 = Ifges(4,3) + (Ifges(6,2) + Ifges(5,2)) * t45 + ((Ifges(6,1) + Ifges(5,1)) * t25 + 0.2e1 * (Ifges(5,4) + Ifges(6,4)) * t27) * t25;
t33 = -mrSges(5,1) * t25 - mrSges(5,2) * t27;
t32 = mrSges(5,3) * t44;
t31 = (t28 * mrSges(4,1) - t26 * mrSges(4,2)) * pkin(2);
t17 = t27 * qJ(5);
t15 = -pkin(3) - t41;
t11 = t27 * pkin(7) + t17;
t10 = -t27 * mrSges(5,1) + t40;
t8 = (-qJ(5) - pkin(7)) * t25;
t7 = t16 - t41;
t2 = t27 * t14 + t17;
t1 = (-qJ(5) - t14) * t25;
t3 = [m(2) + m(3) + m(4) + (m(5) / 0.2e1 + m(6) / 0.2e1) * t44; m(6) * (t1 * t27 + t2 * t25); 0.2e1 * t15 * t10 + t7 * t43 + Ifges(3,3) + 0.2e1 * t31 + (-t1 * t25 + t2 * t27) * t37 + t14 * t32 + m(6) * (t1 ^ 2 + t2 ^ 2 + t7 ^ 2) + m(5) * (t38 * t14 ^ 2 + t15 ^ 2) + m(4) * (t26 ^ 2 + t28 ^ 2) * pkin(2) ^ 2 + t34; m(6) * (t11 * t25 + t8 * t27); (t7 + t16) * t9 + (t15 - pkin(3)) * t10 + t31 + m(6) * (t8 * t1 + t11 * t2 + t16 * t7) + m(5) * (-pkin(3) * t15 + pkin(7) * t36) + ((t11 + t2) * t27 + (-t1 - t8) * t25) * mrSges(6,3) + (t38 * pkin(7) + t36) * mrSges(5,3) + t34; -0.2e1 * pkin(3) * t10 + t16 * t43 + (t11 * t27 - t8 * t25) * t37 + pkin(7) * t32 + m(6) * (t11 ^ 2 + t16 ^ 2 + t8 ^ 2) + m(5) * (t38 * pkin(7) ^ 2 + pkin(3) ^ 2) + t34; -t40 - t18 + (mrSges(5,1) + mrSges(6,1) + t42) * t27; t1 * mrSges(6,1) - t2 * mrSges(6,2) + t33 * t14 + (m(6) * t1 - t39) * pkin(4) + t35; t8 * mrSges(6,1) - t11 * mrSges(6,2) + t33 * pkin(7) + (m(6) * t8 - t39) * pkin(4) + t35; Ifges(5,3) + Ifges(6,3) + (0.2e1 * mrSges(6,1) + t42) * pkin(4); 0; m(6) * t7 + t9; m(6) * t16 + t9; 0; m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t3(1), t3(2), t3(4), t3(7), t3(11); t3(2), t3(3), t3(5), t3(8), t3(12); t3(4), t3(5), t3(6), t3(9), t3(13); t3(7), t3(8), t3(9), t3(10), t3(14); t3(11), t3(12), t3(13), t3(14), t3(15);];
Mq = res;
