% Calculate joint inertia matrix for
% S5PRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:16
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPR1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR1_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR1_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR1_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR1_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR1_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:15:35
% EndTime: 2019-12-05 16:15:36
% DurationCPUTime: 0.26s
% Computational Cost: add. (281->83), mult. (544->113), div. (0->0), fcn. (448->6), ass. (0->36)
t29 = sin(pkin(9));
t30 = cos(pkin(9));
t31 = sin(qJ(5));
t33 = cos(qJ(5));
t16 = t33 * t29 + t31 * t30;
t48 = t16 ^ 2;
t47 = t30 ^ 2;
t15 = -t31 * t29 + t33 * t30;
t5 = -t15 * mrSges(6,1) + t16 * mrSges(6,2);
t46 = 0.2e1 * t5;
t34 = cos(qJ(3));
t45 = t34 * pkin(2);
t44 = Ifges(6,5) * t16 + Ifges(6,6) * t15;
t43 = t29 ^ 2 + t47;
t42 = 2 * mrSges(6,3);
t23 = -t30 * pkin(4) - pkin(3);
t19 = -t30 * mrSges(5,1) + t29 * mrSges(5,2);
t41 = t43 * qJ(4);
t40 = Ifges(6,1) * t48 + Ifges(5,2) * t47 + Ifges(4,3) + (Ifges(5,1) * t29 + 0.2e1 * Ifges(5,4) * t30) * t29 + (0.2e1 * Ifges(6,4) * t16 + Ifges(6,2) * t15) * t15;
t39 = 0.2e1 * t43 * mrSges(5,3);
t32 = sin(qJ(3));
t38 = (t34 * mrSges(4,1) - t32 * mrSges(4,2)) * pkin(2);
t37 = t19 + t5;
t26 = t30 * pkin(7);
t24 = -pkin(3) - t45;
t22 = t32 * pkin(2) + qJ(4);
t20 = t30 * qJ(4) + t26;
t18 = (-pkin(7) - qJ(4)) * t29;
t17 = t23 - t45;
t9 = t30 * t22 + t26;
t8 = (-pkin(7) - t22) * t29;
t7 = t31 * t18 + t33 * t20;
t6 = t33 * t18 - t31 * t20;
t4 = t31 * t8 + t33 * t9;
t3 = -t31 * t9 + t33 * t8;
t1 = [m(2) + m(3) + m(4) + m(5) * t43 + m(6) * (t15 ^ 2 + t48); m(6) * (t3 * t15 + t4 * t16); t17 * t46 + 0.2e1 * t24 * t19 + Ifges(3,3) + 0.2e1 * t38 + (t4 * t15 - t3 * t16) * t42 + t22 * t39 + m(6) * (t17 ^ 2 + t3 ^ 2 + t4 ^ 2) + m(5) * (t43 * t22 ^ 2 + t24 ^ 2) + m(4) * (t32 ^ 2 + t34 ^ 2) * pkin(2) ^ 2 + t40; m(6) * (t6 * t15 + t7 * t16); (t17 + t23) * t5 + (t24 - pkin(3)) * t19 + t38 + m(6) * (t23 * t17 + t6 * t3 + t7 * t4) + m(5) * (-pkin(3) * t24 + t22 * t41) + ((-t3 - t6) * t16 + (t4 + t7) * t15) * mrSges(6,3) + (t43 * t22 + t41) * mrSges(5,3) + t40; -0.2e1 * pkin(3) * t19 + t23 * t46 + (t7 * t15 - t6 * t16) * t42 + qJ(4) * t39 + m(6) * (t23 ^ 2 + t6 ^ 2 + t7 ^ 2) + m(5) * (t43 * qJ(4) ^ 2 + pkin(3) ^ 2) + t40; 0; m(5) * t24 + m(6) * t17 + t37; -m(5) * pkin(3) + m(6) * t23 + t37; m(5) + m(6); -t5; t3 * mrSges(6,1) - t4 * mrSges(6,2) + t44; t6 * mrSges(6,1) - t7 * mrSges(6,2) + t44; 0; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
