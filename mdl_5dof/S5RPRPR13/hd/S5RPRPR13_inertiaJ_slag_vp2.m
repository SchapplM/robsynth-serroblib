% Calculate joint inertia matrix for
% S5RPRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% Datum: 2019-12-31 18:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR13_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR13_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR13_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR13_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR13_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR13_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:31:58
% EndTime: 2019-12-31 18:31:59
% DurationCPUTime: 0.35s
% Computational Cost: add. (496->126), mult. (926->169), div. (0->0), fcn. (904->6), ass. (0->55)
t40 = cos(pkin(8));
t35 = t40 ^ 2;
t68 = -2 * mrSges(5,2);
t31 = -t40 * pkin(2) - pkin(1);
t67 = 0.2e1 * t31;
t66 = pkin(3) + pkin(7);
t41 = sin(qJ(5));
t43 = cos(qJ(5));
t54 = t41 ^ 2 + t43 ^ 2;
t24 = m(6) * t54;
t65 = m(5) + t24;
t64 = cos(qJ(3));
t63 = Ifges(6,4) * t41;
t62 = Ifges(6,4) * t43;
t39 = sin(pkin(8));
t42 = sin(qJ(3));
t22 = t42 * t39 - t64 * t40;
t61 = t22 * t41;
t60 = t22 * t43;
t59 = t43 * mrSges(6,1);
t58 = mrSges(5,1) + mrSges(4,3);
t57 = -mrSges(5,2) + mrSges(4,1);
t56 = pkin(6) + qJ(2);
t55 = t39 ^ 2 + t35;
t23 = t64 * t39 + t42 * t40;
t53 = Ifges(6,5) * t61 + Ifges(6,6) * t60 + Ifges(6,3) * t23;
t25 = t56 * t39;
t26 = t56 * t40;
t12 = t64 * t25 + t42 * t26;
t14 = -t42 * t25 + t64 * t26;
t52 = t12 ^ 2 + t14 ^ 2;
t51 = t54 * mrSges(6,3);
t50 = -t40 * mrSges(3,1) + t39 * mrSges(3,2);
t47 = -t23 * qJ(4) + t31;
t3 = t66 * t22 + t47;
t6 = t23 * pkin(4) + t12;
t1 = -t41 * t3 + t43 * t6;
t2 = t43 * t3 + t41 * t6;
t49 = t43 * t1 + t41 * t2;
t48 = -t41 * mrSges(6,2) + t59;
t45 = qJ(4) ^ 2;
t33 = Ifges(6,5) * t43;
t29 = Ifges(6,1) * t43 - t63;
t28 = -Ifges(6,2) * t41 + t62;
t27 = t41 * mrSges(6,1) + t43 * mrSges(6,2);
t19 = t23 * mrSges(5,3);
t18 = t23 * mrSges(4,2);
t11 = -t23 * mrSges(6,2) + mrSges(6,3) * t60;
t10 = t23 * mrSges(6,1) - mrSges(6,3) * t61;
t9 = t48 * t22;
t8 = t22 * pkin(3) + t47;
t7 = -t22 * pkin(4) + t14;
t5 = Ifges(6,5) * t23 + (Ifges(6,1) * t41 + t62) * t22;
t4 = Ifges(6,6) * t23 + (Ifges(6,2) * t43 + t63) * t22;
t13 = [Ifges(3,2) * t35 - 0.2e1 * pkin(1) * t50 + t18 * t67 - 0.2e1 * t8 * t19 - 0.2e1 * t7 * t9 + 0.2e1 * t1 * t10 + 0.2e1 * t2 * t11 + Ifges(2,3) + (Ifges(3,1) * t39 + 0.2e1 * Ifges(3,4) * t40) * t39 + 0.2e1 * t55 * qJ(2) * mrSges(3,3) + (mrSges(4,1) * t67 + t8 * t68 + t43 * t4 + t41 * t5 + (Ifges(4,2) + Ifges(5,3)) * t22 - 0.2e1 * t58 * t14) * t22 + m(6) * (t1 ^ 2 + t2 ^ 2 + t7 ^ 2) + m(5) * (t8 ^ 2 + t52) + m(4) * (t31 ^ 2 + t52) + m(3) * (t55 * qJ(2) ^ 2 + pkin(1) ^ 2) + ((Ifges(4,1) + Ifges(5,2)) * t23 + t53 + 0.2e1 * (-Ifges(4,4) - Ifges(5,6)) * t22 + 0.2e1 * t58 * t12) * t23; -m(3) * pkin(1) - t41 * t10 + t43 * t11 + t18 - t19 + t57 * t22 + m(6) * (-t41 * t1 + t43 * t2) + m(5) * t8 + m(4) * t31 + t50; m(3) + m(4) + t65; -qJ(4) * t9 + t7 * t27 + (-mrSges(4,2) + mrSges(5,3)) * t14 - t57 * t12 + (-t66 * t10 - t1 * mrSges(6,3) + t5 / 0.2e1) * t43 + (-t66 * t11 - t2 * mrSges(6,3) - t4 / 0.2e1) * t41 + m(6) * (qJ(4) * t7 - t49 * t66) + m(5) * (-pkin(3) * t12 + qJ(4) * t14) + (-Ifges(5,4) + Ifges(4,5) - pkin(3) * mrSges(5,1) - Ifges(6,6) * t41 / 0.2e1 + t33 / 0.2e1) * t23 + (Ifges(5,5) - Ifges(4,6) + t41 * t29 / 0.2e1 + t43 * t28 / 0.2e1 - qJ(4) * mrSges(5,1)) * t22; 0; pkin(3) * t68 - t41 * t28 + t43 * t29 + Ifges(5,1) + Ifges(4,3) + m(6) * (t54 * t66 ^ 2 + t45) + m(5) * (pkin(3) ^ 2 + t45) + 0.2e1 * (t27 + mrSges(5,3)) * qJ(4) + 0.2e1 * t66 * t51; m(5) * t12 + m(6) * t49 + t23 * mrSges(5,1) + t43 * t10 + t41 * t11; 0; -m(5) * pkin(3) - t24 * t66 + mrSges(5,2) - t51; t65; t1 * mrSges(6,1) - t2 * mrSges(6,2) + t53; -t27; -t66 * t59 + t33 + (mrSges(6,2) * t66 - Ifges(6,6)) * t41; t48; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t13(1), t13(2), t13(4), t13(7), t13(11); t13(2), t13(3), t13(5), t13(8), t13(12); t13(4), t13(5), t13(6), t13(9), t13(13); t13(7), t13(8), t13(9), t13(10), t13(14); t13(11), t13(12), t13(13), t13(14), t13(15);];
Mq = res;
