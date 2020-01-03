% Calculate joint inertia matrix for
% S5RRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2020-01-03 11:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR1_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR1_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR1_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:55:44
% EndTime: 2020-01-03 11:55:45
% DurationCPUTime: 0.32s
% Computational Cost: add. (438->100), mult. (790->137), div. (0->0), fcn. (656->8), ass. (0->48)
t38 = sin(pkin(9));
t40 = cos(pkin(9));
t42 = sin(qJ(5));
t44 = cos(qJ(5));
t24 = t38 * t44 + t40 * t42;
t64 = t24 ^ 2;
t63 = t40 ^ 2;
t23 = -t38 * t42 + t40 * t44;
t7 = -t23 * mrSges(6,1) + t24 * mrSges(6,2);
t62 = 0.2e1 * t7;
t26 = -t40 * mrSges(5,1) + t38 * mrSges(5,2);
t61 = 0.2e1 * t26;
t43 = sin(qJ(2));
t60 = pkin(1) * t43;
t59 = t40 * pkin(4);
t45 = cos(qJ(2));
t33 = pkin(1) * t45 + pkin(2);
t39 = sin(pkin(8));
t41 = cos(pkin(8));
t14 = t33 * t41 - t39 * t60;
t58 = t14 * mrSges(4,1);
t15 = t39 * t33 + t41 * t60;
t57 = t15 * mrSges(4,2);
t56 = Ifges(6,5) * t24 + Ifges(6,6) * t23;
t55 = t38 ^ 2 + t63;
t54 = 2 * mrSges(6,3);
t32 = -pkin(2) * t41 - pkin(3);
t31 = pkin(2) * t39 + qJ(4);
t53 = t55 * t31;
t13 = -pkin(3) - t14;
t52 = t41 * mrSges(4,1) - t39 * mrSges(4,2);
t51 = 0.2e1 * t55 * mrSges(5,3);
t50 = Ifges(6,1) * t64 + Ifges(5,2) * t63 + Ifges(3,3) + Ifges(4,3) + (Ifges(5,1) * t38 + 0.2e1 * Ifges(5,4) * t40) * t38 + (0.2e1 * Ifges(6,4) * t24 + Ifges(6,2) * t23) * t23;
t49 = (mrSges(3,1) * t45 - mrSges(3,2) * t43) * pkin(1);
t48 = t26 + t7;
t35 = t40 * pkin(7);
t25 = t32 - t59;
t17 = t31 * t40 + t35;
t16 = (-pkin(7) - t31) * t38;
t12 = qJ(4) + t15;
t10 = t13 - t59;
t9 = t12 * t40 + t35;
t8 = (-pkin(7) - t12) * t38;
t6 = t16 * t42 + t17 * t44;
t5 = t16 * t44 - t17 * t42;
t2 = t42 * t8 + t44 * t9;
t1 = -t42 * t9 + t44 * t8;
t3 = [t50 + m(5) * (t12 ^ 2 * t55 + t13 ^ 2) + (-t1 * t24 + t2 * t23) * t54 + m(3) * (t43 ^ 2 + t45 ^ 2) * pkin(1) ^ 2 + t12 * t51 + m(6) * (t1 ^ 2 + t10 ^ 2 + t2 ^ 2) + m(4) * (t14 ^ 2 + t15 ^ 2) + Ifges(2,3) + 0.2e1 * t49 + t10 * t62 + 0.2e1 * t58 - 0.2e1 * t57 + t13 * t61; t58 - t57 + (t10 + t25) * t7 + (t13 + t32) * t26 + t49 + m(6) * (t1 * t5 + t10 * t25 + t2 * t6) + m(5) * (t12 * t53 + t32 * t13) + (m(4) * (t14 * t41 + t15 * t39) + t52) * pkin(2) + ((-t1 - t5) * t24 + (t2 + t6) * t23) * mrSges(6,3) + (t12 * t55 + t53) * mrSges(5,3) + t50; t25 * t62 + t32 * t61 + (t6 * t23 - t5 * t24) * t54 + t31 * t51 + m(6) * (t25 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(5) * (t31 ^ 2 * t55 + t32 ^ 2) + t50 + (0.2e1 * t52 + m(4) * (t39 ^ 2 + t41 ^ 2) * pkin(2)) * pkin(2); m(6) * (t1 * t23 + t2 * t24); m(6) * (t23 * t5 + t24 * t6); m(4) + m(5) * t55 + m(6) * (t23 ^ 2 + t64); m(5) * t13 + m(6) * t10 + t48; m(5) * t32 + m(6) * t25 + t48; 0; m(5) + m(6); mrSges(6,1) * t1 - mrSges(6,2) * t2 + t56; mrSges(6,1) * t5 - mrSges(6,2) * t6 + t56; -t7; 0; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t3(1), t3(2), t3(4), t3(7), t3(11); t3(2), t3(3), t3(5), t3(8), t3(12); t3(4), t3(5), t3(6), t3(9), t3(13); t3(7), t3(8), t3(9), t3(10), t3(14); t3(11), t3(12), t3(13), t3(14), t3(15);];
Mq = res;
