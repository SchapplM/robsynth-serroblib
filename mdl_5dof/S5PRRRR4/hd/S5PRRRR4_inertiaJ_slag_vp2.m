% Calculate joint inertia matrix for
% S5PRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:08
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRR4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR4_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR4_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR4_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR4_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR4_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:07:41
% EndTime: 2019-12-05 17:07:42
% DurationCPUTime: 0.29s
% Computational Cost: add. (342->99), mult. (672->140), div. (0->0), fcn. (556->6), ass. (0->41)
t33 = sin(qJ(5));
t34 = sin(qJ(4));
t36 = cos(qJ(5));
t37 = cos(qJ(4));
t19 = t33 * t37 + t36 * t34;
t58 = t19 ^ 2;
t57 = t37 ^ 2;
t18 = -t33 * t34 + t36 * t37;
t56 = t33 * pkin(4) * t18 * mrSges(6,3) + Ifges(5,5) * t34 + Ifges(5,6) * t37;
t6 = -t18 * mrSges(6,1) + t19 * mrSges(6,2);
t55 = 0.2e1 * t6;
t38 = cos(qJ(3));
t54 = t38 * pkin(2);
t53 = Ifges(6,5) * t19 + Ifges(6,6) * t18;
t52 = t34 ^ 2 + t57;
t51 = 0.2e1 * mrSges(6,3);
t50 = t36 * t19 * mrSges(6,3);
t27 = -t37 * pkin(4) - pkin(3);
t35 = sin(qJ(3));
t25 = t35 * pkin(2) + pkin(7);
t49 = t52 * t25;
t48 = Ifges(6,1) * t58 + Ifges(5,2) * t57 + Ifges(4,3) + (Ifges(5,1) * t34 + 0.2e1 * Ifges(5,4) * t37) * t34 + (0.2e1 * Ifges(6,4) * t19 + Ifges(6,2) * t18) * t18;
t21 = -t37 * mrSges(5,1) + t34 * mrSges(5,2);
t47 = -t34 * mrSges(5,1) - t37 * mrSges(5,2);
t14 = (-pkin(8) - t25) * t34;
t30 = t37 * pkin(8);
t15 = t37 * t25 + t30;
t4 = t36 * t14 - t33 * t15;
t5 = t33 * t14 + t36 * t15;
t46 = t4 * mrSges(6,1) - t5 * mrSges(6,2) + t53;
t22 = (-pkin(8) - pkin(7)) * t34;
t23 = t37 * pkin(7) + t30;
t8 = t36 * t22 - t33 * t23;
t9 = t33 * t22 + t36 * t23;
t45 = t8 * mrSges(6,1) - t9 * mrSges(6,2) + t53;
t44 = 0.2e1 * t52 * mrSges(5,3);
t43 = (t38 * mrSges(4,1) - t35 * mrSges(4,2)) * pkin(2);
t42 = (t36 * mrSges(6,1) - t33 * mrSges(6,2)) * pkin(4);
t26 = -pkin(3) - t54;
t20 = t27 - t54;
t1 = [m(2) + m(3) + m(4) + m(5) * t52 + m(6) * (t18 ^ 2 + t58); m(6) * (t4 * t18 + t5 * t19); t20 * t55 + 0.2e1 * t26 * t21 + Ifges(3,3) + 0.2e1 * t43 + (t5 * t18 - t4 * t19) * t51 + t25 * t44 + m(6) * (t20 ^ 2 + t4 ^ 2 + t5 ^ 2) + m(5) * (t52 * t25 ^ 2 + t26 ^ 2) + m(4) * (t35 ^ 2 + t38 ^ 2) * pkin(2) ^ 2 + t48; m(6) * (t8 * t18 + t9 * t19); (t20 + t27) * t6 + (t26 - pkin(3)) * t21 + t43 + m(6) * (t27 * t20 + t8 * t4 + t9 * t5) + m(5) * (-pkin(3) * t26 + pkin(7) * t49) + ((-t4 - t8) * t19 + (t5 + t9) * t18) * mrSges(6,3) + (t52 * pkin(7) + t49) * mrSges(5,3) + t48; -0.2e1 * pkin(3) * t21 + t27 * t55 + (t9 * t18 - t8 * t19) * t51 + pkin(7) * t44 + m(6) * (t27 ^ 2 + t8 ^ 2 + t9 ^ 2) + m(5) * (t52 * pkin(7) ^ 2 + pkin(3) ^ 2) + t48; m(6) * (t18 * t36 + t19 * t33) * pkin(4) - t21 - t6; t47 * t25 + (-t50 + m(6) * (t33 * t5 + t36 * t4)) * pkin(4) + t46 + t56; t47 * pkin(7) + (-t50 + m(6) * (t33 * t9 + t36 * t8)) * pkin(4) + t45 + t56; Ifges(5,3) + Ifges(6,3) + m(6) * (t33 ^ 2 + t36 ^ 2) * pkin(4) ^ 2 + 0.2e1 * t42; -t6; t46; t45; Ifges(6,3) + t42; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
