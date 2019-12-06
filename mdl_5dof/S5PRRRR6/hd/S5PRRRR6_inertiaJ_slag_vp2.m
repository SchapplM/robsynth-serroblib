% Calculate joint inertia matrix for
% S5PRRRR6
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
% Datum: 2019-12-05 17:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRR6_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR6_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR6_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR6_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR6_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR6_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:09:35
% EndTime: 2019-12-05 17:09:36
% DurationCPUTime: 0.39s
% Computational Cost: add. (458->123), mult. (954->177), div. (0->0), fcn. (876->8), ass. (0->54)
t52 = cos(qJ(4));
t46 = t52 ^ 2;
t48 = sin(qJ(4));
t71 = t48 ^ 2 + t46;
t79 = mrSges(5,3) * t71;
t47 = sin(qJ(5));
t51 = cos(qJ(5));
t28 = -t47 * t48 + t51 * t52;
t74 = t28 * mrSges(6,3);
t78 = t47 * pkin(4) * t74 + Ifges(5,5) * t48 + Ifges(5,6) * t52;
t49 = sin(qJ(3));
t50 = sin(qJ(2));
t53 = cos(qJ(3));
t54 = cos(qJ(2));
t29 = t49 * t50 - t53 * t54;
t24 = t29 ^ 2;
t31 = t47 * t52 + t51 * t48;
t11 = -t28 * mrSges(6,1) + t31 * mrSges(6,2);
t77 = 0.2e1 * t11;
t76 = t53 * pkin(2);
t73 = t31 * mrSges(6,3);
t72 = Ifges(6,5) * t31 + Ifges(6,6) * t28;
t70 = 0.2e1 * mrSges(6,3);
t69 = t51 * t73;
t32 = t49 * t54 + t53 * t50;
t6 = t31 * t32;
t7 = t28 * t32;
t68 = -t6 * mrSges(6,1) - t7 * mrSges(6,2);
t41 = -t52 * pkin(4) - pkin(3);
t67 = t71 * pkin(7);
t39 = t49 * pkin(2) + pkin(7);
t66 = t71 * t39;
t65 = Ifges(6,1) * t31 ^ 2 + Ifges(5,2) * t46 + Ifges(4,3) + (Ifges(5,1) * t48 + 0.2e1 * Ifges(5,4) * t52) * t48 + (0.2e1 * Ifges(6,4) * t31 + Ifges(6,2) * t28) * t28;
t64 = -t48 * mrSges(5,1) - t52 * mrSges(5,2);
t63 = 0.2e1 * t79;
t22 = (-pkin(8) - t39) * t48;
t44 = t52 * pkin(8);
t23 = t52 * t39 + t44;
t10 = t47 * t22 + t51 * t23;
t9 = t51 * t22 - t47 * t23;
t62 = t9 * mrSges(6,1) - t10 * mrSges(6,2) + t72;
t35 = (-pkin(8) - pkin(7)) * t48;
t36 = t52 * pkin(7) + t44;
t13 = t51 * t35 - t47 * t36;
t14 = t47 * t35 + t51 * t36;
t61 = t13 * mrSges(6,1) - t14 * mrSges(6,2) + t72;
t60 = (t53 * mrSges(4,1) - t49 * mrSges(4,2)) * pkin(2);
t59 = (t51 * mrSges(6,1) - t47 * mrSges(6,2)) * pkin(4);
t34 = -t52 * mrSges(5,1) + t48 * mrSges(5,2);
t58 = t6 * t73 + t7 * t74 + (-mrSges(4,1) + t11 + t34) * t29 + (-mrSges(4,2) + t79) * t32;
t40 = -pkin(3) - t76;
t33 = t41 - t76;
t25 = t32 ^ 2;
t1 = [m(2) + m(6) * (t6 ^ 2 + t7 ^ 2 + t24) + m(5) * (t71 * t25 + t24) + m(4) * (t25 + t24) + m(3) * (t50 ^ 2 + t54 ^ 2); t54 * mrSges(3,1) - t50 * mrSges(3,2) + m(6) * (t10 * t7 + t33 * t29 - t9 * t6) + m(5) * (t40 * t29 + t32 * t66) + m(4) * (-t29 * t53 + t32 * t49) * pkin(2) + t58; t33 * t77 + 0.2e1 * t40 * t34 + Ifges(3,3) + 0.2e1 * t60 + (t10 * t28 - t9 * t31) * t70 + t39 * t63 + m(6) * (t10 ^ 2 + t33 ^ 2 + t9 ^ 2) + m(5) * (t71 * t39 ^ 2 + t40 ^ 2) + m(4) * (t49 ^ 2 + t53 ^ 2) * pkin(2) ^ 2 + t65; m(6) * (-t13 * t6 + t14 * t7 + t41 * t29) + m(5) * (-pkin(3) * t29 + t32 * t67) + t58; (-pkin(3) + t40) * t34 + (t41 + t33) * t11 + t60 + m(5) * (-pkin(3) * t40 + pkin(7) * t66) + m(6) * (t14 * t10 + t13 * t9 + t41 * t33) + ((-t13 - t9) * t31 + (t10 + t14) * t28) * mrSges(6,3) + (t66 + t67) * mrSges(5,3) + t65; -0.2e1 * pkin(3) * t34 + t41 * t77 + (-t13 * t31 + t14 * t28) * t70 + pkin(7) * t63 + m(6) * (t13 ^ 2 + t14 ^ 2 + t41 ^ 2) + m(5) * (t71 * pkin(7) ^ 2 + pkin(3) ^ 2) + t65; t64 * t32 + m(6) * (t47 * t7 - t51 * t6) * pkin(4) + t68; t64 * t39 + (-t69 + m(6) * (t10 * t47 + t51 * t9)) * pkin(4) + t62 + t78; t64 * pkin(7) + (-t69 + m(6) * (t13 * t51 + t14 * t47)) * pkin(4) + t61 + t78; Ifges(5,3) + Ifges(6,3) + m(6) * (t47 ^ 2 + t51 ^ 2) * pkin(4) ^ 2 + 0.2e1 * t59; t68; t62; t61; Ifges(6,3) + t59; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
