% Calculate joint inertia matrix for
% S5RPPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-31 18:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR9_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR9_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR9_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR9_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR9_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR9_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:02:19
% EndTime: 2019-12-31 18:02:20
% DurationCPUTime: 0.29s
% Computational Cost: add. (350->121), mult. (611->175), div. (0->0), fcn. (460->6), ass. (0->57)
t40 = sin(qJ(5));
t42 = cos(qJ(5));
t18 = -mrSges(6,1) * t42 + mrSges(6,2) * t40;
t66 = -m(6) * pkin(4) - mrSges(5,1) + t18;
t65 = t40 / 0.2e1;
t64 = Ifges(6,4) * t40;
t63 = Ifges(6,4) * t42;
t43 = cos(qJ(4));
t62 = Ifges(6,5) * t43;
t38 = sin(pkin(8));
t39 = cos(pkin(8));
t44 = -pkin(1) - pkin(2);
t17 = t39 * qJ(2) + t38 * t44;
t13 = -pkin(6) + t17;
t61 = t13 * t38;
t59 = t38 * t43;
t41 = sin(qJ(4));
t58 = t40 * t41;
t57 = t40 * t43;
t56 = t41 * t42;
t55 = t42 * t43;
t54 = Ifges(6,6) * t58 + Ifges(6,3) * t43;
t53 = t40 ^ 2 + t42 ^ 2;
t35 = t41 ^ 2;
t37 = t43 ^ 2;
t52 = t35 + t37;
t51 = t53 * t41;
t50 = t52 * mrSges(5,3);
t19 = t43 * mrSges(5,1) - t41 * mrSges(5,2);
t16 = -t38 * qJ(2) + t39 * t44;
t12 = pkin(3) - t16;
t31 = t43 * pkin(4);
t3 = pkin(7) * t41 + t12 + t31;
t1 = -t13 * t57 + t3 * t42;
t2 = t13 * t55 + t3 * t40;
t49 = -t1 * t40 + t2 * t42;
t7 = -t38 * t57 - t39 * t42;
t8 = t38 * t55 - t39 * t40;
t48 = -t40 * t7 + t42 * t8;
t47 = -mrSges(6,1) * t40 - mrSges(6,2) * t42;
t14 = -mrSges(6,2) * t43 + mrSges(6,3) * t58;
t15 = mrSges(6,1) * t43 + mrSges(6,3) * t56;
t46 = t42 * t14 - t40 * t15;
t33 = t39 ^ 2;
t32 = t38 ^ 2;
t30 = Ifges(6,5) * t40;
t29 = Ifges(6,6) * t42;
t23 = t35 * t32;
t21 = Ifges(6,1) * t40 + t63;
t20 = Ifges(6,2) * t42 + t64;
t11 = t13 ^ 2;
t10 = t47 * t41;
t9 = t35 * t11;
t6 = t35 * t61;
t5 = t62 + (-Ifges(6,1) * t42 + t64) * t41;
t4 = Ifges(6,6) * t43 + (Ifges(6,2) * t40 - t63) * t41;
t22 = [(2 * pkin(1) * mrSges(3,1)) - 0.2e1 * t16 * mrSges(4,1) + 0.2e1 * t17 * mrSges(4,2) + 0.2e1 * qJ(2) * mrSges(3,3) + 0.2e1 * t1 * t15 + 0.2e1 * t12 * t19 + 0.2e1 * t2 * t14 + Ifges(3,2) + Ifges(2,3) + Ifges(4,3) + (Ifges(5,2) * t43 + t54) * t43 - 0.2e1 * t13 * t50 + (Ifges(5,1) * t41 + 0.2e1 * Ifges(5,4) * t43 + 0.2e1 * t13 * t10 + t40 * t4 + (-t5 - t62) * t42) * t41 + m(6) * (t1 ^ 2 + t2 ^ 2 + t9) + m(5) * (t11 * t37 + t12 ^ 2 + t9) + m(4) * (t16 ^ 2 + t17 ^ 2) + m(3) * ((pkin(1) ^ 2) + qJ(2) ^ 2); -m(3) * pkin(1) + t8 * t14 + t7 * t15 - mrSges(3,1) + (-mrSges(4,1) - t19) * t39 + (t41 * t10 + mrSges(4,2) - t50) * t38 + m(6) * (t1 * t7 + t2 * t8 + t6) + m(5) * (-t12 * t39 + t37 * t61 + t6) + m(4) * (t16 * t39 + t17 * t38); m(3) + m(4) * (t32 + t33) + m(5) * (t32 * t37 + t23 + t33) + m(6) * (t7 ^ 2 + t8 ^ 2 + t23); -t43 * t10 + (m(6) * (-t13 * t43 + t49) + t46) * t41; m(6) * (t48 - t59) * t41; m(4) + m(5) * t52 + m(6) * (t35 * t53 + t37); -pkin(4) * t10 + t5 * t65 + t42 * t4 / 0.2e1 + (-t13 * mrSges(5,2) + t30 / 0.2e1 + t29 / 0.2e1 - Ifges(5,6)) * t43 + t49 * mrSges(6,3) + (m(6) * t49 + t46) * pkin(7) + (-t42 * t21 / 0.2e1 + t20 * t65 - Ifges(5,5) + t66 * t13) * t41; -mrSges(5,2) * t59 + (m(6) * pkin(7) + mrSges(6,3)) * t48 + t66 * t38 * t41; -t43 * t18 + m(6) * (pkin(7) * t51 + t31) + mrSges(6,3) * t51 + t19; Ifges(5,3) + m(6) * (pkin(7) ^ 2 * t53 + pkin(4) ^ 2) + t40 * t21 + t42 * t20 - 0.2e1 * pkin(4) * t18 + 0.2e1 * t53 * pkin(7) * mrSges(6,3); mrSges(6,1) * t1 - mrSges(6,2) * t2 - Ifges(6,5) * t56 + t54; mrSges(6,1) * t7 - mrSges(6,2) * t8; t10; pkin(7) * t47 + t29 + t30; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t22(1), t22(2), t22(4), t22(7), t22(11); t22(2), t22(3), t22(5), t22(8), t22(12); t22(4), t22(5), t22(6), t22(9), t22(13); t22(7), t22(8), t22(9), t22(10), t22(14); t22(11), t22(12), t22(13), t22(14), t22(15);];
Mq = res;
