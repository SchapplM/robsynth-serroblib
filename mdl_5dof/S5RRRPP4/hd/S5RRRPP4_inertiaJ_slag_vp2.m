% Calculate joint inertia matrix for
% S5RRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
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
% Datum: 2019-12-31 20:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPP4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP4_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP4_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP4_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP4_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP4_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:54:53
% EndTime: 2019-12-31 20:54:54
% DurationCPUTime: 0.42s
% Computational Cost: add. (704->139), mult. (1329->192), div. (0->0), fcn. (1327->6), ass. (0->45)
t66 = -2 * mrSges(6,1);
t65 = 2 * mrSges(6,3);
t43 = sin(qJ(3));
t44 = sin(qJ(2));
t45 = cos(qJ(3));
t46 = cos(qJ(2));
t29 = -t43 * t44 + t45 * t46;
t38 = -pkin(2) * t46 - pkin(1);
t21 = -pkin(3) * t29 + t38;
t64 = 0.2e1 * t21;
t63 = -pkin(7) - pkin(6);
t62 = pkin(2) * t43;
t37 = pkin(2) * t45 + pkin(3);
t41 = sin(pkin(8));
t42 = cos(pkin(8));
t24 = t37 * t42 - t41 * t62;
t61 = t24 * mrSges(5,1);
t25 = t41 * t37 + t42 * t62;
t60 = t25 * mrSges(5,2);
t59 = mrSges(6,2) + mrSges(5,3);
t54 = t63 * t44;
t55 = t63 * t46;
t20 = t43 * t54 - t45 * t55;
t58 = t44 ^ 2 + t46 ^ 2;
t19 = t43 * t55 + t45 * t54;
t30 = t43 * t46 + t44 * t45;
t51 = -t30 * qJ(4) + t19;
t9 = qJ(4) * t29 + t20;
t5 = t41 * t9 - t42 * t51;
t7 = t41 * t51 + t42 * t9;
t57 = t5 ^ 2 + t7 ^ 2;
t56 = Ifges(6,2) + Ifges(4,3) + Ifges(5,3);
t53 = t42 * mrSges(5,1) - t41 * mrSges(5,2);
t52 = (mrSges(4,1) * t45 - mrSges(4,2) * t43) * pkin(2);
t16 = -t42 * t29 + t30 * t41;
t17 = t29 * t41 + t30 * t42;
t50 = t19 * mrSges(4,1) - t20 * mrSges(4,2) + Ifges(4,5) * t30 + Ifges(4,6) * t29 + (-mrSges(5,2) + mrSges(6,3)) * t7 + (-mrSges(5,1) - mrSges(6,1)) * t5 + (Ifges(6,4) + Ifges(5,5)) * t17 + (-Ifges(5,6) + Ifges(6,6)) * t16;
t36 = -pkin(3) * t42 - pkin(4);
t35 = pkin(3) * t41 + qJ(5);
t23 = -pkin(4) - t24;
t22 = qJ(5) + t25;
t11 = t17 * mrSges(5,2);
t10 = t16 * mrSges(6,1);
t1 = pkin(4) * t16 - qJ(5) * t17 + t21;
t2 = [t46 * (Ifges(3,4) * t44 + Ifges(3,2) * t46) - 0.2e1 * pkin(1) * (-mrSges(3,1) * t46 + mrSges(3,2) * t44) + t44 * (Ifges(3,1) * t44 + Ifges(3,4) * t46) + 0.2e1 * t38 * (-mrSges(4,1) * t29 + mrSges(4,2) * t30) + t29 * (Ifges(4,4) * t30 + Ifges(4,2) * t29) + t30 * (Ifges(4,1) * t30 + Ifges(4,4) * t29) + t11 * t64 + 0.2e1 * t1 * t10 + Ifges(2,3) + 0.2e1 * (-t19 * t30 + t20 * t29) * mrSges(4,3) + 0.2e1 * t58 * pkin(6) * mrSges(3,3) + (-0.2e1 * t1 * mrSges(6,3) + (Ifges(5,1) + Ifges(6,1)) * t17 + 0.2e1 * t59 * t5) * t17 + (mrSges(5,1) * t64 + (Ifges(6,3) + Ifges(5,2)) * t16 - 0.2e1 * t59 * t7 + 0.2e1 * (-Ifges(5,4) + Ifges(6,5)) * t17) * t16 + m(6) * (t1 ^ 2 + t57) + m(5) * (t21 ^ 2 + t57) + m(4) * (t19 ^ 2 + t20 ^ 2 + t38 ^ 2) + m(3) * (t58 * pkin(6) ^ 2 + pkin(1) ^ 2); t50 + (m(4) * (t19 * t45 + t20 * t43) + (t43 * t29 - t45 * t30) * mrSges(4,3)) * pkin(2) + m(6) * (t22 * t7 + t23 * t5) + m(5) * (-t24 * t5 + t25 * t7) + (-mrSges(3,1) * t44 - mrSges(3,2) * t46) * pkin(6) + (-t16 * t25 - t17 * t24) * mrSges(5,3) + (-t16 * t22 + t17 * t23) * mrSges(6,2) + Ifges(3,6) * t46 + Ifges(3,5) * t44; 0.2e1 * t61 + t23 * t66 - 0.2e1 * t60 + t22 * t65 + Ifges(3,3) + 0.2e1 * t52 + m(5) * (t24 ^ 2 + t25 ^ 2) + m(6) * (t22 ^ 2 + t23 ^ 2) + m(4) * (t43 ^ 2 + t45 ^ 2) * pkin(2) ^ 2 + t56; t50 + (m(5) * (t41 * t7 - t42 * t5) + (-t16 * t41 - t17 * t42) * mrSges(5,3)) * pkin(3) + m(6) * (t35 * t7 + t36 * t5) + (-t16 * t35 + t17 * t36) * mrSges(6,2); m(6) * (t22 * t35 + t23 * t36) - t60 + t61 + t52 + (t35 + t22) * mrSges(6,3) + (-t36 - t23) * mrSges(6,1) + (m(5) * (t24 * t42 + t25 * t41) + t53) * pkin(3) + t56; t36 * t66 + t35 * t65 + m(6) * (t35 ^ 2 + t36 ^ 2) + t56 + (0.2e1 * t53 + m(5) * (t41 ^ 2 + t42 ^ 2) * pkin(3)) * pkin(3); m(5) * t21 + m(6) * t1 + t16 * mrSges(5,1) - t17 * mrSges(6,3) + t10 + t11; 0; 0; m(5) + m(6); m(6) * t5 + t17 * mrSges(6,2); m(6) * t23 - mrSges(6,1); m(6) * t36 - mrSges(6,1); 0; m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
