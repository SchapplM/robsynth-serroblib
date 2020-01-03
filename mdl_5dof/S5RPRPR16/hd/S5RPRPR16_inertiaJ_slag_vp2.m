% Calculate joint inertia matrix for
% S5RPRPR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5]';
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
% Datum: 2019-12-31 18:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR16_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR16_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPR16_inertiaJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR16_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR16_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR16_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:38:45
% EndTime: 2019-12-31 18:38:46
% DurationCPUTime: 0.38s
% Computational Cost: add. (290->120), mult. (519->153), div. (0->0), fcn. (324->4), ass. (0->53)
t31 = sin(qJ(5));
t33 = cos(qJ(5));
t48 = t31 ^ 2 + t33 ^ 2;
t45 = m(6) * t48;
t67 = m(5) + t45;
t32 = sin(qJ(3));
t26 = t32 ^ 2;
t34 = cos(qJ(3));
t28 = t34 ^ 2;
t47 = t28 + t26;
t66 = (mrSges(5,1) + mrSges(4,3)) * t47;
t21 = qJ(4) * t32;
t64 = m(5) * (pkin(3) * t34 + t21);
t63 = t48 * mrSges(6,3) - mrSges(5,2);
t10 = t32 * pkin(3) - t34 * qJ(4) + qJ(2);
t62 = -0.2e1 * t10;
t61 = 0.2e1 * qJ(2);
t60 = -t31 / 0.2e1;
t59 = t33 / 0.2e1;
t52 = t32 * t33;
t9 = -t34 * mrSges(6,2) + mrSges(6,3) * t52;
t57 = t31 * t9;
t36 = -pkin(1) - pkin(6);
t56 = pkin(4) - t36;
t55 = Ifges(6,4) * t31;
t54 = Ifges(6,4) * t33;
t53 = t31 * t32;
t35 = -pkin(3) - pkin(7);
t51 = t33 * t35;
t13 = t31 * mrSges(6,1) + t33 * mrSges(6,2);
t50 = t13 + mrSges(5,3);
t49 = t47 * t36 ^ 2;
t46 = Ifges(6,5) * t53 + Ifges(6,6) * t52 + Ifges(6,3) * t34;
t44 = t48 * t35;
t12 = t56 * t34;
t5 = t32 * pkin(7) + t10;
t1 = t33 * t12 - t31 * t5;
t2 = t31 * t12 + t33 * t5;
t42 = t33 * t1 + t31 * t2;
t8 = t34 * mrSges(6,1) - mrSges(6,3) * t53;
t41 = -t33 * t8 - t57;
t40 = t33 * mrSges(6,1) - t31 * mrSges(6,2);
t39 = 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1) * t47;
t38 = qJ(2) ^ 2;
t37 = qJ(4) ^ 2;
t23 = Ifges(6,5) * t33;
t15 = Ifges(6,1) * t33 - t55;
t14 = -Ifges(6,2) * t31 + t54;
t11 = t56 * t32;
t6 = t40 * t32;
t4 = Ifges(6,5) * t34 + (Ifges(6,1) * t31 + t54) * t32;
t3 = Ifges(6,6) * t34 + (Ifges(6,2) * t33 + t55) * t32;
t7 = [-(2 * pkin(1) * mrSges(3,2)) + mrSges(3,3) * t61 + 0.2e1 * t1 * t8 + 0.2e1 * t11 * t6 + 0.2e1 * t2 * t9 + Ifges(3,1) + Ifges(2,3) + (mrSges(4,2) * t61 + mrSges(5,3) * t62 + (Ifges(4,1) + Ifges(5,2)) * t34 + t46) * t34 + (mrSges(4,1) * t61 + mrSges(5,2) * t62 + t33 * t3 + t31 * t4 + (Ifges(4,2) + Ifges(5,3)) * t32 + 0.2e1 * (-Ifges(4,4) - Ifges(5,6)) * t34) * t32 + m(6) * (t1 ^ 2 + t11 ^ 2 + t2 ^ 2) + m(5) * (t10 ^ 2 + t49) + m(4) * (t38 + t49) + m(3) * ((pkin(1) ^ 2) + t38) - 0.2e1 * t36 * t66; -m(3) * pkin(1) - t32 * t6 + mrSges(3,2) + t41 * t34 + m(6) * (-t32 * t11 - t42 * t34) + t36 * t39 - t66; m(3) + m(6) * (t48 * t28 + t26) + t39; t35 * t57 + t8 * t51 + t3 * t60 - t11 * t13 + t4 * t59 - qJ(4) * t6 + m(6) * (-qJ(4) * t11 + t42 * t35) - t42 * mrSges(6,3) + (-Ifges(5,4) + Ifges(4,5) - pkin(3) * mrSges(5,1) + Ifges(6,6) * t60 + t23 / 0.2e1) * t34 + (Ifges(5,5) - Ifges(4,6) - qJ(4) * mrSges(5,1) + t14 * t59 + t31 * t15 / 0.2e1) * t32 + (t64 + (mrSges(4,1) - mrSges(5,2)) * t34 + (-mrSges(4,2) + mrSges(5,3)) * t32) * t36; (mrSges(4,1) + t63) * t34 + t64 + m(6) * (-t34 * t44 + t21) + (-mrSges(4,2) + t50) * t32; -0.2e1 * pkin(3) * mrSges(5,2) - t31 * t14 + t33 * t15 + Ifges(5,1) + Ifges(4,3) + m(6) * (t48 * t35 ^ 2 + t37) + m(5) * (pkin(3) ^ 2 + t37) - 0.2e1 * mrSges(6,3) * t44 + 0.2e1 * t50 * qJ(4); m(6) * t42 + (-m(5) * t36 + mrSges(5,1)) * t34 - t41; -t67 * t34; -m(5) * pkin(3) + t35 * t45 - t63; t67; t1 * mrSges(6,1) - t2 * mrSges(6,2) + t46; -t40 * t34; mrSges(6,1) * t51 + t23 + (-mrSges(6,2) * t35 - Ifges(6,6)) * t31; t40; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t7(1), t7(2), t7(4), t7(7), t7(11); t7(2), t7(3), t7(5), t7(8), t7(12); t7(4), t7(5), t7(6), t7(9), t7(13); t7(7), t7(8), t7(9), t7(10), t7(14); t7(11), t7(12), t7(13), t7(14), t7(15);];
Mq = res;
