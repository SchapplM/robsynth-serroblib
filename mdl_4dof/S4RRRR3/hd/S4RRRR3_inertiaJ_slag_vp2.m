% Calculate joint inertia matrix for
% S4RRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 17:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRR3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR3_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR3_inertiaJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR3_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR3_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR3_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:24:17
% EndTime: 2019-12-31 17:24:18
% DurationCPUTime: 0.27s
% Computational Cost: add. (390->94), mult. (766->146), div. (0->0), fcn. (740->6), ass. (0->37)
t50 = -pkin(6) - pkin(5);
t32 = sin(qJ(3));
t49 = pkin(2) * t32;
t35 = cos(qJ(3));
t26 = t35 * pkin(2) + pkin(3);
t31 = sin(qJ(4));
t34 = cos(qJ(4));
t16 = t31 * t26 + t34 * t49;
t48 = t16 * mrSges(5,2);
t47 = t31 * mrSges(5,2);
t46 = Ifges(4,3) + Ifges(5,3);
t33 = sin(qJ(2));
t24 = t50 * t33;
t36 = cos(qJ(2));
t25 = t50 * t36;
t12 = t32 * t24 - t35 * t25;
t45 = t33 ^ 2 + t36 ^ 2;
t44 = pkin(3) * t47;
t27 = -t36 * pkin(2) - pkin(1);
t11 = t35 * t24 + t32 * t25;
t15 = t34 * t26 - t31 * t49;
t14 = t15 * mrSges(5,1);
t43 = Ifges(5,3) + t14 - t48;
t20 = t32 * t36 + t35 * t33;
t4 = -t20 * pkin(7) + t11;
t19 = -t32 * t33 + t35 * t36;
t5 = t19 * pkin(7) + t12;
t2 = -t31 * t5 + t34 * t4;
t3 = t31 * t4 + t34 * t5;
t8 = t34 * t19 - t31 * t20;
t9 = t31 * t19 + t34 * t20;
t42 = t2 * mrSges(5,1) - t3 * mrSges(5,2) + Ifges(5,5) * t9 + Ifges(5,6) * t8;
t41 = (t35 * mrSges(4,1) - t32 * mrSges(4,2)) * pkin(2);
t40 = t11 * mrSges(4,1) - t12 * mrSges(4,2) + Ifges(4,5) * t20 + Ifges(4,6) * t19 + t42;
t28 = t34 * pkin(3) * mrSges(5,1);
t13 = -t19 * pkin(3) + t27;
t1 = [Ifges(2,3) + 0.2e1 * t13 * (-t8 * mrSges(5,1) + t9 * mrSges(5,2)) + t9 * (Ifges(5,1) * t9 + Ifges(5,4) * t8) + t8 * (Ifges(5,4) * t9 + Ifges(5,2) * t8) + 0.2e1 * t27 * (-t19 * mrSges(4,1) + t20 * mrSges(4,2)) + t19 * (Ifges(4,4) * t20 + Ifges(4,2) * t19) + t20 * (Ifges(4,1) * t20 + Ifges(4,4) * t19) - 0.2e1 * pkin(1) * (-t36 * mrSges(3,1) + t33 * mrSges(3,2)) + t33 * (Ifges(3,1) * t33 + Ifges(3,4) * t36) + t36 * (Ifges(3,4) * t33 + Ifges(3,2) * t36) + m(5) * (t13 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(4) * (t11 ^ 2 + t12 ^ 2 + t27 ^ 2) + m(3) * (t45 * pkin(5) ^ 2 + pkin(1) ^ 2) + 0.2e1 * (-t2 * t9 + t3 * t8) * mrSges(5,3) + 0.2e1 * (-t11 * t20 + t12 * t19) * mrSges(4,3) + 0.2e1 * t45 * pkin(5) * mrSges(3,3); m(5) * (t15 * t2 + t16 * t3) + Ifges(3,5) * t33 + Ifges(3,6) * t36 + (-t33 * mrSges(3,1) - t36 * mrSges(3,2)) * pkin(5) + (-t15 * t9 + t16 * t8) * mrSges(5,3) + (m(4) * (t11 * t35 + t12 * t32) + (t32 * t19 - t35 * t20) * mrSges(4,3)) * pkin(2) + t40; -0.2e1 * t48 + Ifges(3,3) + 0.2e1 * t14 + 0.2e1 * t41 + m(5) * (t15 ^ 2 + t16 ^ 2) + m(4) * (t32 ^ 2 + t35 ^ 2) * pkin(2) ^ 2 + t46; (m(5) * (t2 * t34 + t3 * t31) + (t31 * t8 - t34 * t9) * mrSges(5,3)) * pkin(3) + t40; Ifges(4,3) + t28 + t41 + (m(5) * (t15 * t34 + t16 * t31) - t47) * pkin(3) + t43; -0.2e1 * t44 + 0.2e1 * t28 + m(5) * (t31 ^ 2 + t34 ^ 2) * pkin(3) ^ 2 + t46; t42; t43; Ifges(5,3) + t28 - t44; Ifges(5,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
