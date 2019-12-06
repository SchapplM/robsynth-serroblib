% Calculate joint inertia matrix for
% S5RPPRR1
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
% Datum: 2019-12-05 17:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR1_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR1_inertiaJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR1_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR1_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR1_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:37:59
% EndTime: 2019-12-05 17:38:00
% DurationCPUTime: 0.20s
% Computational Cost: add. (195->74), mult. (326->96), div. (0->0), fcn. (237->4), ass. (0->32)
t22 = sin(qJ(5));
t23 = sin(qJ(4));
t24 = cos(qJ(5));
t25 = cos(qJ(4));
t8 = -t22 * t25 - t24 * t23;
t45 = t8 ^ 2;
t44 = t25 ^ 2;
t43 = t23 * mrSges(5,1) + t25 * mrSges(5,2) + mrSges(4,3);
t20 = (pkin(1) + qJ(3));
t42 = t20 ^ 2;
t41 = 2 * t20;
t19 = qJ(2) - pkin(6);
t40 = -pkin(7) + t19;
t39 = t23 ^ 2 + t44;
t10 = -t22 * t23 + t24 * t25;
t38 = t10 ^ 2 + t45;
t37 = t10 * mrSges(6,1) + t8 * mrSges(6,2);
t36 = m(5) * t39;
t35 = t39 * mrSges(5,3);
t11 = t40 * t23;
t12 = t40 * t25;
t2 = -t22 * t11 + t24 * t12;
t3 = t24 * t11 + t22 * t12;
t34 = t10 * t2 - t8 * t3;
t33 = t2 * mrSges(6,1) - t3 * mrSges(6,2) + Ifges(6,5) * t10 + Ifges(6,6) * t8;
t32 = -t8 * mrSges(6,1) + t10 * mrSges(6,2);
t31 = t10 * t24 - t22 * t8;
t30 = t25 * mrSges(5,1) - t23 * mrSges(5,2);
t28 = (t24 * mrSges(6,1) - t22 * mrSges(6,2)) * pkin(4);
t26 = qJ(2) ^ 2;
t13 = t23 * pkin(4) + t20;
t1 = [Ifges(3,1) + Ifges(4,1) + Ifges(2,3) + 0.2e1 * t13 * t32 + Ifges(6,2) * t45 + Ifges(5,1) * t44 - (2 * pkin(1) * mrSges(3,2)) + m(6) * (t13 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(5) * (t39 * t19 ^ 2 + t42) + m(4) * (t26 + t42) + m(3) * ((pkin(1) ^ 2) + t26) + (-0.2e1 * Ifges(5,4) * t25 + Ifges(5,2) * t23) * t23 + (Ifges(6,1) * t10 + 0.2e1 * Ifges(6,4) * t8) * t10 + t43 * t41 + (2 * (mrSges(3,3) + mrSges(4,2)) * qJ(2)) - 0.2e1 * t34 * mrSges(6,3) - 0.2e1 * t19 * t35; -m(3) * pkin(1) + mrSges(3,2) - m(6) * t13 + (-m(5) / 0.2e1 - m(4) / 0.2e1) * t41 - t32 - t43; m(3) + m(4) + m(5) + m(6); m(4) * qJ(2) + m(6) * t34 - t38 * mrSges(6,3) + t19 * t36 + mrSges(4,2) - t35; 0; m(6) * t38 + m(4) + t36; Ifges(5,5) * t25 - Ifges(5,6) * t23 + t30 * t19 + (m(6) * (t2 * t24 + t22 * t3) - t31 * mrSges(6,3)) * pkin(4) + t33; 0; m(6) * t31 * pkin(4) + t30 + t37; Ifges(5,3) + Ifges(6,3) + m(6) * (t22 ^ 2 + t24 ^ 2) * pkin(4) ^ 2 + 0.2e1 * t28; t33; 0; t37; Ifges(6,3) + t28; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
