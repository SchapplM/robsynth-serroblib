% Calculate joint inertia matrix for
% S5PRRPR3
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
% Datum: 2019-12-05 16:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPR3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR3_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR3_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR3_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR3_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR3_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:19:03
% EndTime: 2019-12-05 16:19:05
% DurationCPUTime: 0.28s
% Computational Cost: add. (399->95), mult. (792->144), div. (0->0), fcn. (796->6), ass. (0->38)
t34 = cos(qJ(3));
t49 = t34 ^ 2;
t48 = m(5) * pkin(3);
t29 = sin(pkin(9));
t30 = cos(pkin(9));
t32 = sin(qJ(3));
t18 = -t29 * t32 + t30 * t34;
t47 = t18 ^ 2;
t46 = 0.2e1 * t18;
t45 = pkin(3) * t29;
t25 = t30 * pkin(3) + pkin(4);
t31 = sin(qJ(5));
t33 = cos(qJ(5));
t15 = t33 * t25 - t31 * t45;
t44 = t15 * mrSges(6,1);
t16 = t31 * t25 + t33 * t45;
t43 = t16 * mrSges(6,2);
t42 = -qJ(4) - pkin(6);
t23 = t42 * t32;
t24 = t42 * t34;
t13 = t29 * t23 - t30 * t24;
t41 = t32 ^ 2 + t49;
t26 = -t34 * pkin(3) - pkin(2);
t19 = t29 * t34 + t30 * t32;
t10 = t33 * t18 - t31 * t19;
t11 = t31 * t18 + t33 * t19;
t4 = -t10 * mrSges(6,1) + t11 * mrSges(6,2);
t40 = -t18 * mrSges(5,1) + t19 * mrSges(5,2);
t12 = t30 * t23 + t29 * t24;
t5 = -t19 * pkin(7) + t12;
t6 = t18 * pkin(7) + t13;
t2 = -t31 * t6 + t33 * t5;
t3 = t31 * t5 + t33 * t6;
t39 = t2 * mrSges(6,1) - t3 * mrSges(6,2) + Ifges(6,5) * t11 + Ifges(6,6) * t10;
t38 = -t34 * mrSges(4,1) + t32 * mrSges(4,2);
t37 = -t4 - t40;
t14 = -t18 * pkin(4) + t26;
t1 = [m(2) + m(3) + m(4) * t41 + m(5) * (t19 ^ 2 + t47) + m(6) * (t10 ^ 2 + t11 ^ 2); m(5) * (t12 * t18 + t13 * t19) + m(6) * (t2 * t10 + t3 * t11); Ifges(3,3) + Ifges(5,2) * t47 + Ifges(4,2) * t49 + 0.2e1 * t14 * t4 - 0.2e1 * pkin(2) * t38 + 0.2e1 * t26 * t40 + t13 * mrSges(5,3) * t46 + (-0.2e1 * t2 * mrSges(6,3) + Ifges(6,1) * t11) * t11 + 0.2e1 * t41 * pkin(6) * mrSges(4,3) + (-0.2e1 * t12 * mrSges(5,3) + Ifges(5,1) * t19 + Ifges(5,4) * t46) * t19 + (0.2e1 * t3 * mrSges(6,3) + 0.2e1 * Ifges(6,4) * t11 + Ifges(6,2) * t10) * t10 + m(6) * (t14 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(5) * (t12 ^ 2 + t13 ^ 2 + t26 ^ 2) + m(4) * (t41 * pkin(6) ^ 2 + pkin(2) ^ 2) + (Ifges(4,1) * t32 + 0.2e1 * Ifges(4,4) * t34) * t32; m(6) * (t15 * t10 + t16 * t11) + (t18 * t30 + t19 * t29) * t48 + t37 - t38; Ifges(5,5) * t19 + Ifges(5,6) * t18 + Ifges(4,5) * t32 + Ifges(4,6) * t34 + m(6) * (t15 * t2 + t16 * t3) + t12 * mrSges(5,1) - t13 * mrSges(5,2) + (-t32 * mrSges(4,1) - t34 * mrSges(4,2)) * pkin(6) + (t16 * t10 - t15 * t11) * mrSges(6,3) + (m(5) * (t12 * t30 + t13 * t29) + (t29 * t18 - t30 * t19) * mrSges(5,3)) * pkin(3) + t39; 0.2e1 * t44 - 0.2e1 * t43 + Ifges(4,3) + Ifges(5,3) + Ifges(6,3) + m(6) * (t15 ^ 2 + t16 ^ 2) + (0.2e1 * t30 * mrSges(5,1) - 0.2e1 * t29 * mrSges(5,2) + (t29 ^ 2 + t30 ^ 2) * t48) * pkin(3); 0; m(5) * t26 + m(6) * t14 - t37; 0; m(5) + m(6); -t4; t39; Ifges(6,3) - t43 + t44; 0; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
