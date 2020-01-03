% Calculate joint inertia matrix for
% S5RPRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-31 18:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR10_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR10_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR10_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR10_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR10_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR10_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:25:55
% EndTime: 2019-12-31 18:25:56
% DurationCPUTime: 0.28s
% Computational Cost: add. (354->88), mult. (536->119), div. (0->0), fcn. (395->6), ass. (0->36)
t26 = sin(pkin(8));
t27 = cos(pkin(8));
t29 = sin(qJ(3));
t31 = cos(qJ(3));
t11 = t26 * t31 + t27 * t29;
t28 = sin(qJ(5));
t30 = cos(qJ(5));
t15 = -t30 * mrSges(6,1) + t28 * mrSges(6,2);
t51 = t30 ^ 2;
t41 = t28 ^ 2 + t51;
t39 = t41 * mrSges(6,3);
t9 = t26 * t29 - t27 * t31;
t52 = -(mrSges(5,2) - t39) * t11 + t31 * mrSges(4,1) - t29 * mrSges(4,2) - (mrSges(5,1) - t15) * t9;
t50 = m(5) * pkin(3);
t20 = t26 * pkin(3) + pkin(7);
t38 = t41 * t20;
t47 = t9 ^ 2;
t32 = -pkin(1) - pkin(2);
t13 = -t29 * qJ(2) + t31 * t32;
t12 = -pkin(3) + t13;
t14 = t31 * qJ(2) + t29 * t32;
t5 = t26 * t12 + t27 * t14;
t4 = t27 * t12 - t26 * t14;
t46 = t4 * mrSges(5,1);
t45 = t5 * mrSges(5,2);
t44 = t13 * mrSges(4,1);
t43 = t14 * mrSges(4,2);
t3 = -pkin(7) + t5;
t40 = t41 * t3;
t36 = t27 * mrSges(5,1) - t26 * mrSges(5,2);
t35 = -mrSges(6,1) * t28 - mrSges(6,2) * t30;
t34 = Ifges(6,2) * t51 + Ifges(4,3) + Ifges(5,3) + (Ifges(6,1) * t28 + 0.2e1 * Ifges(6,4) * t30) * t28;
t21 = -t27 * pkin(3) - pkin(4);
t8 = t11 ^ 2;
t2 = pkin(4) - t4;
t1 = [(2 * pkin(1) * mrSges(3,1)) - 0.2e1 * t44 - 0.2e1 * t46 + 0.2e1 * t43 + 0.2e1 * t45 + 0.2e1 * qJ(2) * mrSges(3,3) - 0.2e1 * t2 * t15 + Ifges(3,2) + Ifges(2,3) - 0.2e1 * t3 * t39 + m(6) * (t3 ^ 2 * t41 + t2 ^ 2) + m(5) * (t4 ^ 2 + t5 ^ 2) + m(4) * (t13 ^ 2 + t14 ^ 2) + m(3) * ((pkin(1) ^ 2) + qJ(2) ^ 2) + t34; -m(3) * pkin(1) - mrSges(3,1) + m(6) * (t11 * t40 + t9 * t2) + m(5) * (t11 * t5 - t9 * t4) + m(4) * (t31 * t13 + t29 * t14) - t52; m(3) + m(4) * (t29 ^ 2 + t31 ^ 2) + m(5) * (t8 + t47) + m(6) * (t41 * t8 + t47); m(6) * (t21 * t2 + t3 * t38) - t43 + t44 - t45 + t46 + (-t21 + t2) * t15 + (m(5) * (t26 * t5 + t27 * t4) - t36) * pkin(3) + (t40 - t38) * mrSges(6,3) - t34; m(6) * (t11 * t38 + t21 * t9) + (t11 * t26 - t27 * t9) * t50 + t52; 0.2e1 * t21 * t15 + m(6) * (t20 ^ 2 * t41 + t21 ^ 2) + t34 + 0.2e1 * mrSges(6,3) * t38 + (0.2e1 * t36 + (t26 ^ 2 + t27 ^ 2) * t50) * pkin(3); 0; 0; 0; m(6) * t41 + m(5); (-mrSges(6,2) * t3 - Ifges(6,6)) * t30 + (-mrSges(6,1) * t3 - Ifges(6,5)) * t28; t35 * t11; Ifges(6,5) * t28 + Ifges(6,6) * t30 + t20 * t35; -t15; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
