% Calculate joint inertia matrix for
% S5RRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-05 18:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR4_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR4_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR4_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:31:50
% EndTime: 2019-12-05 18:31:51
% DurationCPUTime: 0.38s
% Computational Cost: add. (523->116), mult. (944->164), div. (0->0), fcn. (796->8), ass. (0->53)
t44 = sin(qJ(5));
t45 = sin(qJ(4));
t47 = cos(qJ(5));
t48 = cos(qJ(4));
t27 = t44 * t48 + t47 * t45;
t74 = t27 ^ 2;
t73 = t48 ^ 2;
t26 = -t44 * t45 + t47 * t48;
t72 = t44 * pkin(4) * t26 * mrSges(6,3) + Ifges(5,5) * t45 + Ifges(5,6) * t48;
t9 = -t26 * mrSges(6,1) + t27 * mrSges(6,2);
t71 = 0.2e1 * t9;
t29 = -t48 * mrSges(5,1) + t45 * mrSges(5,2);
t70 = 0.2e1 * t29;
t46 = sin(qJ(2));
t69 = pkin(1) * t46;
t68 = t48 * pkin(4);
t49 = cos(qJ(2));
t36 = t49 * pkin(1) + pkin(2);
t42 = sin(pkin(9));
t43 = cos(pkin(9));
t17 = t43 * t36 - t42 * t69;
t67 = t17 * mrSges(4,1);
t18 = t42 * t36 + t43 * t69;
t66 = t18 * mrSges(4,2);
t65 = Ifges(6,5) * t27 + Ifges(6,6) * t26;
t64 = t45 ^ 2 + t73;
t63 = 0.2e1 * mrSges(6,3);
t62 = t47 * t27 * mrSges(6,3);
t35 = -t43 * pkin(2) - pkin(3);
t34 = t42 * pkin(2) + pkin(7);
t61 = t64 * t34;
t14 = -pkin(3) - t17;
t60 = t43 * mrSges(4,1) - t42 * mrSges(4,2);
t59 = -t45 * mrSges(5,1) - t48 * mrSges(5,2);
t15 = pkin(7) + t18;
t10 = (-pkin(8) - t15) * t45;
t39 = t48 * pkin(8);
t11 = t48 * t15 + t39;
t2 = t47 * t10 - t44 * t11;
t3 = t44 * t10 + t47 * t11;
t58 = t2 * mrSges(6,1) - t3 * mrSges(6,2) + t65;
t19 = (-pkin(8) - t34) * t45;
t20 = t48 * t34 + t39;
t7 = t47 * t19 - t44 * t20;
t8 = t44 * t19 + t47 * t20;
t57 = t7 * mrSges(6,1) - t8 * mrSges(6,2) + t65;
t56 = 0.2e1 * t64 * mrSges(5,3);
t55 = Ifges(6,1) * t74 + Ifges(5,2) * t73 + Ifges(3,3) + Ifges(4,3) + (Ifges(5,1) * t45 + 0.2e1 * Ifges(5,4) * t48) * t45 + (0.2e1 * Ifges(6,4) * t27 + Ifges(6,2) * t26) * t26;
t54 = (t49 * mrSges(3,1) - t46 * mrSges(3,2)) * pkin(1);
t53 = (t47 * mrSges(6,1) - t44 * mrSges(6,2)) * pkin(4);
t28 = t35 - t68;
t12 = t14 - t68;
t1 = [t55 + (-t2 * t27 + t3 * t26) * t63 + m(5) * (t64 * t15 ^ 2 + t14 ^ 2) + t15 * t56 + m(3) * (t46 ^ 2 + t49 ^ 2) * pkin(1) ^ 2 + m(6) * (t12 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(4) * (t17 ^ 2 + t18 ^ 2) + t14 * t70 + t12 * t71 + 0.2e1 * t67 - 0.2e1 * t66 + Ifges(2,3) + 0.2e1 * t54; t67 - t66 + (t28 + t12) * t9 + (t14 + t35) * t29 + t54 + m(6) * (t28 * t12 + t7 * t2 + t8 * t3) + m(5) * (t35 * t14 + t15 * t61) + (m(4) * (t17 * t43 + t18 * t42) + t60) * pkin(2) + ((-t2 - t7) * t27 + (t3 + t8) * t26) * mrSges(6,3) + (t64 * t15 + t61) * mrSges(5,3) + t55; t28 * t71 + t35 * t70 + (t8 * t26 - t7 * t27) * t63 + t34 * t56 + m(6) * (t28 ^ 2 + t7 ^ 2 + t8 ^ 2) + m(5) * (t64 * t34 ^ 2 + t35 ^ 2) + t55 + (0.2e1 * t60 + m(4) * (t42 ^ 2 + t43 ^ 2) * pkin(2)) * pkin(2); m(6) * (t2 * t26 + t3 * t27); m(6) * (t7 * t26 + t8 * t27); m(4) + m(5) * t64 + m(6) * (t26 ^ 2 + t74); t59 * t15 + (-t62 + m(6) * (t2 * t47 + t3 * t44)) * pkin(4) + t58 + t72; t59 * t34 + (-t62 + m(6) * (t44 * t8 + t47 * t7)) * pkin(4) + t57 + t72; m(6) * (t26 * t47 + t27 * t44) * pkin(4) - t29 - t9; Ifges(5,3) + Ifges(6,3) + m(6) * (t44 ^ 2 + t47 ^ 2) * pkin(4) ^ 2 + 0.2e1 * t53; t58; t57; -t9; Ifges(6,3) + t53; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
