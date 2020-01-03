% Calculate joint inertia matrix for
% S5RRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR5_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR5_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR5_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR5_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR5_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:28:54
% EndTime: 2019-12-31 19:28:55
% DurationCPUTime: 0.34s
% Computational Cost: add. (486->125), mult. (893->168), div. (0->0), fcn. (860->6), ass. (0->40)
t35 = cos(qJ(2));
t27 = -t35 * pkin(2) - pkin(1);
t49 = 0.2e1 * t27;
t31 = cos(pkin(8));
t26 = -t31 * pkin(2) - pkin(3);
t23 = -pkin(4) + t26;
t30 = sin(pkin(8));
t24 = t30 * pkin(2) + qJ(4);
t32 = sin(qJ(5));
t34 = cos(qJ(5));
t13 = t34 * t23 - t32 * t24;
t48 = t13 * mrSges(6,1);
t14 = t32 * t23 + t34 * t24;
t47 = t14 * mrSges(6,2);
t46 = mrSges(5,2) + mrSges(4,3);
t45 = -qJ(3) - pkin(6);
t22 = t45 * t35;
t33 = sin(qJ(2));
t42 = t45 * t33;
t12 = -t31 * t22 + t30 * t42;
t44 = t33 ^ 2 + t35 ^ 2;
t10 = -t30 * t22 - t31 * t42;
t43 = t10 ^ 2 + t12 ^ 2;
t17 = t30 * t33 - t31 * t35;
t18 = t30 * t35 + t31 * t33;
t7 = t34 * t17 - t32 * t18;
t8 = t32 * t17 + t34 * t18;
t41 = -t7 * mrSges(6,1) + t8 * mrSges(6,2);
t40 = t34 * mrSges(6,1) - t32 * mrSges(6,2);
t39 = t18 * qJ(4) - t27;
t4 = -t18 * pkin(7) + t10;
t5 = t17 * pkin(7) + t12;
t1 = -t32 * t5 + t34 * t4;
t2 = t32 * t4 + t34 * t5;
t38 = t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t8 + Ifges(6,6) * t7;
t16 = t18 * mrSges(4,2);
t15 = t17 * mrSges(5,1);
t6 = t17 * pkin(3) - t39;
t3 = (-pkin(3) - pkin(4)) * t17 + t39;
t9 = [t33 * (Ifges(3,1) * t33 + Ifges(3,4) * t35) + t35 * (Ifges(3,4) * t33 + Ifges(3,2) * t35) - 0.2e1 * pkin(1) * (-t35 * mrSges(3,1) + t33 * mrSges(3,2)) + 0.2e1 * t6 * t15 + t16 * t49 + 0.2e1 * t3 * t41 + t8 * (Ifges(6,1) * t8 + Ifges(6,4) * t7) + t7 * (Ifges(6,4) * t8 + Ifges(6,2) * t7) + Ifges(2,3) + 0.2e1 * (-t1 * t8 + t2 * t7) * mrSges(6,3) + 0.2e1 * t44 * pkin(6) * mrSges(3,3) + (-0.2e1 * t6 * mrSges(5,3) + (Ifges(5,1) + Ifges(4,1)) * t18 + 0.2e1 * t46 * t10) * t18 + (mrSges(4,1) * t49 + (Ifges(5,3) + Ifges(4,2)) * t17 + 0.2e1 * (-Ifges(4,4) + Ifges(5,5)) * t18 - 0.2e1 * t46 * t12) * t17 + m(6) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(5) * (t6 ^ 2 + t43) + m(4) * (t27 ^ 2 + t43) + m(3) * (t44 * pkin(6) ^ 2 + pkin(1) ^ 2); Ifges(3,5) * t33 + Ifges(3,6) * t35 + (-mrSges(4,2) + mrSges(5,3)) * t12 + (-mrSges(4,1) - mrSges(5,1)) * t10 + (-t33 * mrSges(3,1) - t35 * mrSges(3,2)) * pkin(6) + (-t13 * t8 + t14 * t7) * mrSges(6,3) + m(6) * (t13 * t1 + t14 * t2) + m(5) * (t26 * t10 + t24 * t12) + (t26 * mrSges(5,2) + Ifges(5,4) + Ifges(4,5)) * t18 + (-t24 * mrSges(5,2) - Ifges(4,6) + Ifges(5,6)) * t17 + (m(4) * (-t10 * t31 + t12 * t30) + (-t30 * t17 - t31 * t18) * mrSges(4,3)) * pkin(2) - t38; -0.2e1 * t26 * mrSges(5,1) - 0.2e1 * t48 + 0.2e1 * t47 + 0.2e1 * t24 * mrSges(5,3) + Ifges(5,2) + Ifges(3,3) + Ifges(4,3) + Ifges(6,3) + m(6) * (t13 ^ 2 + t14 ^ 2) + m(5) * (t24 ^ 2 + t26 ^ 2) + (0.2e1 * t31 * mrSges(4,1) - 0.2e1 * t30 * mrSges(4,2) + m(4) * (t30 ^ 2 + t31 ^ 2) * pkin(2)) * pkin(2); m(4) * t27 + m(5) * t6 - m(6) * t3 + t17 * mrSges(4,1) - t18 * mrSges(5,3) + t15 + t16 - t41; 0; m(4) + m(5) + m(6); t18 * mrSges(5,2) + (t32 * t7 - t34 * t8) * mrSges(6,3) + m(6) * (t34 * t1 + t32 * t2) + m(5) * t10; -mrSges(5,1) + m(6) * (t34 * t13 + t32 * t14) + m(5) * t26 - t40; 0; m(5) + m(6) * (t32 ^ 2 + t34 ^ 2); t38; -Ifges(6,3) - t47 + t48; 0; t40; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t9(1), t9(2), t9(4), t9(7), t9(11); t9(2), t9(3), t9(5), t9(8), t9(12); t9(4), t9(5), t9(6), t9(9), t9(13); t9(7), t9(8), t9(9), t9(10), t9(14); t9(11), t9(12), t9(13), t9(14), t9(15);];
Mq = res;
