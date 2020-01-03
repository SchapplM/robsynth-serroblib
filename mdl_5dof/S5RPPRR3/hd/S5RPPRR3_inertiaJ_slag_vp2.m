% Calculate joint inertia matrix for
% S5RPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2020-01-03 11:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR3_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR3_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR3_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:27:29
% EndTime: 2020-01-03 11:27:31
% DurationCPUTime: 0.27s
% Computational Cost: add. (433->84), mult. (793->130), div. (0->0), fcn. (799->8), ass. (0->36)
t30 = sin(pkin(9));
t32 = cos(pkin(9));
t35 = sin(qJ(4));
t37 = cos(qJ(4));
t21 = -t35 * t30 + t37 * t32;
t48 = t21 ^ 2;
t29 = t32 ^ 2;
t47 = 0.2e1 * t21;
t31 = sin(pkin(8));
t25 = t31 * pkin(1) + qJ(3);
t46 = pkin(6) + t25;
t18 = t46 * t30;
t19 = t46 * t32;
t8 = -t35 * t18 + t37 * t19;
t45 = t30 ^ 2 + t29;
t33 = cos(pkin(8));
t26 = -t33 * pkin(1) - pkin(2);
t22 = t37 * t30 + t35 * t32;
t34 = sin(qJ(5));
t36 = cos(qJ(5));
t12 = t36 * t21 - t34 * t22;
t13 = t34 * t21 + t36 * t22;
t4 = -t12 * mrSges(6,1) + t13 * mrSges(6,2);
t44 = -t32 * mrSges(4,1) + t30 * mrSges(4,2);
t43 = -t21 * mrSges(5,1) + t22 * mrSges(5,2);
t7 = -t37 * t18 - t35 * t19;
t5 = -t22 * pkin(7) + t7;
t6 = t21 * pkin(7) + t8;
t2 = -t34 * t6 + t36 * t5;
t3 = t34 * t5 + t36 * t6;
t42 = t2 * mrSges(6,1) - t3 * mrSges(6,2) + Ifges(6,5) * t13 + Ifges(6,6) * t12;
t23 = -t32 * pkin(3) + t26;
t41 = (t36 * mrSges(6,1) - t34 * mrSges(6,2)) * pkin(4);
t40 = -t4 - t43;
t14 = -t21 * pkin(4) + t23;
t1 = [Ifges(2,3) + Ifges(3,3) + t8 * mrSges(5,3) * t47 + 0.2e1 * t14 * t4 + 0.2e1 * t23 * t43 + Ifges(5,2) * t48 + 0.2e1 * t26 * t44 + Ifges(4,2) * t29 + (Ifges(4,1) * t30 + 0.2e1 * Ifges(4,4) * t32) * t30 + (-0.2e1 * t2 * mrSges(6,3) + Ifges(6,1) * t13) * t13 + (-0.2e1 * t7 * mrSges(5,3) + Ifges(5,1) * t22 + Ifges(5,4) * t47) * t22 + (0.2e1 * t3 * mrSges(6,3) + 0.2e1 * Ifges(6,4) * t13 + Ifges(6,2) * t12) * t12 + m(6) * (t14 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(5) * (t23 ^ 2 + t7 ^ 2 + t8 ^ 2) + m(4) * (t45 * t25 ^ 2 + t26 ^ 2) + m(3) * (t31 ^ 2 + t33 ^ 2) * pkin(1) ^ 2 + 0.2e1 * (t33 * mrSges(3,1) - t31 * mrSges(3,2)) * pkin(1) + 0.2e1 * t45 * t25 * mrSges(4,3); m(5) * (t7 * t21 + t8 * t22) + m(6) * (t2 * t12 + t3 * t13); m(3) + m(4) * t45 + m(5) * (t22 ^ 2 + t48) + m(6) * (t12 ^ 2 + t13 ^ 2); m(4) * t26 + m(5) * t23 + m(6) * t14 - t40 + t44; 0; m(4) + m(5) + m(6); t7 * mrSges(5,1) - t8 * mrSges(5,2) + Ifges(5,5) * t22 + Ifges(5,6) * t21 + (m(6) * (t2 * t36 + t3 * t34) + (t34 * t12 - t36 * t13) * mrSges(6,3)) * pkin(4) + t42; m(6) * (t12 * t36 + t13 * t34) * pkin(4) + t40; 0; Ifges(5,3) + Ifges(6,3) + m(6) * (t34 ^ 2 + t36 ^ 2) * pkin(4) ^ 2 + 0.2e1 * t41; t42; -t4; 0; Ifges(6,3) + t41; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
