% Calculate joint inertia matrix for
% S5RRPPR3
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
% Datum: 2019-12-31 19:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR3_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR3_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR3_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR3_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR3_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:26:28
% EndTime: 2019-12-31 19:26:29
% DurationCPUTime: 0.23s
% Computational Cost: add. (241->71), mult. (404->87), div. (0->0), fcn. (235->6), ass. (0->39)
t28 = cos(qJ(5));
t53 = t28 ^ 2;
t26 = sin(qJ(5));
t42 = t26 ^ 2 + t53;
t41 = t42 * mrSges(6,3);
t10 = -t26 * mrSges(6,1) - t28 * mrSges(6,2);
t52 = mrSges(5,3) - t10;
t29 = cos(qJ(2));
t18 = t29 * pkin(1) + pkin(2);
t24 = sin(pkin(8));
t25 = cos(pkin(8));
t27 = sin(qJ(2));
t47 = pkin(1) * t27;
t7 = t24 * t18 + t25 * t47;
t3 = qJ(4) + t7;
t51 = t3 ^ 2;
t15 = t24 * pkin(2) + qJ(4);
t50 = t15 ^ 2;
t49 = 2 * mrSges(5,2);
t9 = m(6) * t42;
t48 = m(5) + t9;
t46 = t15 * t3;
t6 = t25 * t18 - t24 * t47;
t45 = t6 * mrSges(4,1);
t44 = t7 * mrSges(4,2);
t17 = -t25 * pkin(2) - pkin(3);
t14 = -pkin(7) + t17;
t40 = t42 * t14;
t39 = Ifges(6,5) * t28 - Ifges(6,6) * t26;
t38 = 0.2e1 * t52;
t5 = -pkin(3) - t6;
t37 = mrSges(5,2) - t41;
t36 = t25 * mrSges(4,1) - t24 * mrSges(4,2);
t35 = t28 * mrSges(6,1) - t26 * mrSges(6,2);
t34 = -0.2e1 * t41;
t33 = Ifges(6,1) * t53 + Ifges(5,1) + Ifges(3,3) + Ifges(4,3) + (-0.2e1 * Ifges(6,4) * t28 + Ifges(6,2) * t26) * t26;
t32 = (t29 * mrSges(3,1) - t27 * mrSges(3,2)) * pkin(1);
t2 = -pkin(7) + t5;
t1 = [0.2e1 * t45 - 0.2e1 * t44 + t5 * t49 + Ifges(2,3) + t3 * t38 + 0.2e1 * t32 + t2 * t34 + m(6) * (t42 * t2 ^ 2 + t51) + m(4) * (t6 ^ 2 + t7 ^ 2) + m(5) * (t5 ^ 2 + t51) + m(3) * (t27 ^ 2 + t29 ^ 2) * pkin(1) ^ 2 + t33; t45 - t44 + t32 + (t5 + t17) * mrSges(5,2) + m(6) * (t2 * t40 + t46) + m(5) * (t17 * t5 + t46) + (m(4) * (t24 * t7 + t25 * t6) + t36) * pkin(2) + t33 + t52 * (t3 + t15) + (-t14 - t2) * t41; t17 * t49 + t15 * t38 + t14 * t34 + m(5) * (t17 ^ 2 + t50) + m(6) * (t42 * t14 ^ 2 + t50) + t33 + (0.2e1 * t36 + m(4) * (t24 ^ 2 + t25 ^ 2) * pkin(2)) * pkin(2); 0; 0; m(4) + t48; m(5) * t5 + t2 * t9 + t37; m(5) * t17 + m(6) * t40 + t37; 0; t48; t35 * t2 + t39; t35 * t14 + t39; t10; t35; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
