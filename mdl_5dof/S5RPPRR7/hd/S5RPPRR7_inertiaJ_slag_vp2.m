% Calculate joint inertia matrix for
% S5RPPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
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
% Datum: 2019-12-31 18:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR7_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR7_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR7_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR7_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR7_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR7_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:59:31
% EndTime: 2019-12-31 17:59:32
% DurationCPUTime: 0.28s
% Computational Cost: add. (257->103), mult. (503->149), div. (0->0), fcn. (338->6), ass. (0->50)
t32 = sin(qJ(5));
t34 = cos(qJ(5));
t12 = -t34 * mrSges(6,1) + t32 * mrSges(6,2);
t60 = -mrSges(5,1) + t12;
t59 = m(6) * pkin(4) - t60;
t30 = sin(pkin(8));
t20 = t30 * pkin(1) + qJ(3);
t58 = t20 ^ 2;
t57 = 0.2e1 * t20;
t56 = t34 / 0.2e1;
t33 = sin(qJ(4));
t27 = t33 ^ 2;
t35 = cos(qJ(4));
t29 = t35 ^ 2;
t44 = t29 + t27;
t55 = m(5) * t44 + m(4);
t54 = t33 * pkin(4);
t53 = Ifges(6,4) * t32;
t52 = Ifges(6,4) * t34;
t51 = Ifges(6,6) * t33;
t31 = cos(pkin(8));
t22 = -t31 * pkin(1) - pkin(2);
t16 = -pkin(6) + t22;
t50 = t16 * t33;
t49 = t32 * t35;
t48 = t34 * t35;
t47 = t35 * mrSges(6,3);
t46 = Ifges(6,5) * t48 + Ifges(6,3) * t33;
t45 = t32 ^ 2 + t34 ^ 2;
t42 = t45 * t35;
t41 = t44 * mrSges(5,3);
t5 = -t35 * pkin(7) + t20 + t54;
t1 = -t32 * t50 + t34 * t5;
t2 = t32 * t5 + t34 * t50;
t40 = -t1 * t32 + t2 * t34;
t8 = -t33 * mrSges(6,2) - t32 * t47;
t9 = t33 * mrSges(6,1) - t34 * t47;
t39 = -t32 * t9 + t34 * t8;
t38 = -mrSges(6,1) * t32 - mrSges(6,2) * t34;
t25 = Ifges(6,5) * t32;
t24 = Ifges(6,6) * t34;
t15 = t16 ^ 2;
t14 = Ifges(6,1) * t32 + t52;
t13 = Ifges(6,2) * t34 + t53;
t11 = t29 * t16;
t10 = t29 * t15;
t6 = -mrSges(6,1) * t49 - mrSges(6,2) * t48;
t4 = Ifges(6,5) * t33 + (Ifges(6,1) * t34 - t53) * t35;
t3 = t51 + (-Ifges(6,2) * t32 + t52) * t35;
t7 = [0.2e1 * t22 * mrSges(4,2) + mrSges(4,3) * t57 + 0.2e1 * t1 * t9 + 0.2e1 * t2 * t8 + Ifges(4,1) + Ifges(2,3) + Ifges(3,3) + (mrSges(5,1) * t57 + Ifges(5,2) * t33 + t46) * t33 + (mrSges(5,2) * t57 + Ifges(5,1) * t35 - 0.2e1 * Ifges(5,4) * t33 + 0.2e1 * t16 * t6 + t34 * t4 + (-t3 - t51) * t32) * t35 + m(6) * (t1 ^ 2 + t2 ^ 2 + t10) + m(5) * (t27 * t15 + t10 + t58) + m(4) * (t22 ^ 2 + t58) + m(3) * (t30 ^ 2 + t31 ^ 2) * pkin(1) ^ 2 + 0.2e1 * (t31 * mrSges(3,1) - t30 * mrSges(3,2)) * pkin(1) - 0.2e1 * t16 * t41; -t33 * t6 + (m(6) * (t40 - t50) + t39) * t35; m(3) + m(6) * (t45 * t29 + t27) + t55; t35 * t6 + mrSges(4,2) + t39 * t33 - t41 + m(6) * (t40 * t33 + t11) + m(5) * (t27 * t16 + t11) + m(4) * t22; m(6) * (-0.1e1 + t45) * t35 * t33; m(6) * (t45 * t27 + t29) + t55; t3 * t56 + t32 * t4 / 0.2e1 + pkin(4) * t6 + (-t16 * mrSges(5,2) + t25 / 0.2e1 + t24 / 0.2e1 - Ifges(5,6)) * t33 + t40 * mrSges(6,3) + (m(6) * t40 + t39) * pkin(7) + (t14 * t56 - t32 * t13 / 0.2e1 + Ifges(5,5) + t59 * t16) * t35; -t35 * mrSges(5,2) + m(6) * (pkin(7) * t42 - t54) + mrSges(6,3) * t42 + t60 * t33; t59 * t35 + (-mrSges(5,2) + (m(6) * pkin(7) + mrSges(6,3)) * t45) * t33; Ifges(5,3) - 0.2e1 * pkin(4) * t12 + t32 * t14 + t34 * t13 + m(6) * (t45 * pkin(7) ^ 2 + pkin(4) ^ 2) + 0.2e1 * t45 * pkin(7) * mrSges(6,3); t1 * mrSges(6,1) - t2 * mrSges(6,2) - Ifges(6,6) * t49 + t46; t6; t38 * t33; t38 * pkin(7) + t24 + t25; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t7(1), t7(2), t7(4), t7(7), t7(11); t7(2), t7(3), t7(5), t7(8), t7(12); t7(4), t7(5), t7(6), t7(9), t7(13); t7(7), t7(8), t7(9), t7(10), t7(14); t7(11), t7(12), t7(13), t7(14), t7(15);];
Mq = res;
