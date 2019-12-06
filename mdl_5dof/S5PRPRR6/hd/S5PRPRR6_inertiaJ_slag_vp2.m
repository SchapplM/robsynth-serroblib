% Calculate joint inertia matrix for
% S5PRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRR6_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR6_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR6_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR6_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR6_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR6_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:56:29
% EndTime: 2019-12-05 15:56:30
% DurationCPUTime: 0.37s
% Computational Cost: add. (483->130), mult. (1064->198), div. (0->0), fcn. (1109->10), ass. (0->62)
t43 = sin(pkin(10));
t45 = cos(pkin(10));
t48 = sin(qJ(4));
t70 = cos(qJ(4));
t25 = t48 * t43 - t70 * t45;
t26 = t70 * t43 + t48 * t45;
t14 = t25 * mrSges(5,1) + t26 * mrSges(5,2);
t27 = -t45 * mrSges(4,1) + t43 * mrSges(4,2);
t78 = -t14 - t27;
t47 = sin(qJ(5));
t50 = cos(qJ(5));
t29 = -t50 * mrSges(6,1) + t47 * mrSges(6,2);
t77 = -m(6) * pkin(4) - mrSges(5,1) + t29;
t46 = cos(pkin(5));
t44 = sin(pkin(5));
t49 = sin(qJ(2));
t65 = t44 * t49;
t20 = -t43 * t65 + t46 * t45;
t21 = t46 * t43 + t45 * t65;
t7 = -t70 * t20 + t48 * t21;
t76 = t7 ^ 2;
t63 = pkin(7) + qJ(3);
t28 = t63 * t45;
t58 = t63 * t43;
t15 = t48 * t28 + t70 * t58;
t75 = t15 ^ 2;
t40 = t45 ^ 2;
t74 = 0.2e1 * t15;
t73 = t50 / 0.2e1;
t72 = -m(4) - m(5);
t71 = t15 * t7;
t69 = Ifges(6,4) * t47;
t68 = Ifges(6,4) * t50;
t67 = t26 * t47;
t66 = t26 * t50;
t51 = cos(qJ(2));
t64 = t44 * t51;
t62 = Ifges(6,5) * t66 + Ifges(6,3) * t25;
t61 = Ifges(6,5) * t47 + Ifges(6,6) * t50;
t60 = t43 ^ 2 + t40;
t59 = t47 ^ 2 + t50 ^ 2;
t34 = -t45 * pkin(3) - pkin(2);
t11 = t25 * pkin(4) - t26 * pkin(8) + t34;
t17 = t70 * t28 - t48 * t58;
t1 = t50 * t11 - t47 * t17;
t2 = t47 * t11 + t50 * t17;
t57 = -t1 * t47 + t2 * t50;
t55 = mrSges(6,1) * t47 + mrSges(6,2) * t50;
t54 = -t20 * t43 + t21 * t45;
t39 = t44 ^ 2;
t32 = t39 * t51 ^ 2;
t31 = Ifges(6,1) * t47 + t68;
t30 = Ifges(6,2) * t50 + t69;
t13 = t25 * mrSges(6,1) - mrSges(6,3) * t66;
t12 = -t25 * mrSges(6,2) - mrSges(6,3) * t67;
t10 = t55 * t26;
t9 = t48 * t20 + t70 * t21;
t6 = Ifges(6,5) * t25 + (Ifges(6,1) * t50 - t69) * t26;
t5 = Ifges(6,6) * t25 + (-Ifges(6,2) * t47 + t68) * t26;
t4 = -t47 * t64 + t50 * t9;
t3 = -t47 * t9 - t50 * t64;
t8 = [m(2) + m(6) * (t3 ^ 2 + t4 ^ 2 + t76) + m(5) * (t9 ^ 2 + t32 + t76) + m(4) * (t20 ^ 2 + t21 ^ 2 + t32) + m(3) * (t39 * t49 ^ 2 + t46 ^ 2 + t32); t7 * t10 + t4 * t12 + t3 * t13 + (-t9 * t25 + t7 * t26) * mrSges(5,3) + t54 * mrSges(4,3) + (-t49 * mrSges(3,2) + (mrSges(3,1) + t78) * t51) * t44 + m(6) * (t1 * t3 + t2 * t4 + t71) + m(5) * (t17 * t9 - t34 * t64 + t71) + m(4) * (pkin(2) * t64 + t54 * qJ(3)); Ifges(4,2) * t40 - 0.2e1 * pkin(2) * t27 + 0.2e1 * t1 * t13 + t10 * t74 + 0.2e1 * t2 * t12 + 0.2e1 * t34 * t14 + Ifges(3,3) + (Ifges(4,1) * t43 + 0.2e1 * Ifges(4,4) * t45) * t43 + 0.2e1 * t60 * qJ(3) * mrSges(4,3) + (-0.2e1 * t17 * mrSges(5,3) + Ifges(5,2) * t25 + t62) * t25 + m(6) * (t1 ^ 2 + t2 ^ 2 + t75) + m(5) * (t17 ^ 2 + t34 ^ 2 + t75) + m(4) * (t60 * qJ(3) ^ 2 + pkin(2) ^ 2) + (mrSges(5,3) * t74 + Ifges(5,1) * t26 - t47 * t5 + t50 * t6 + (-Ifges(6,6) * t47 - (2 * Ifges(5,4))) * t25) * t26; m(6) * (t50 * t3 + t47 * t4) + t72 * t64; -m(4) * pkin(2) + t47 * t12 + t50 * t13 + m(6) * (t50 * t1 + t47 * t2) + m(5) * t34 - t78; m(6) * t59 - t72; -t9 * mrSges(5,2) + (m(6) * pkin(8) + mrSges(6,3)) * (-t3 * t47 + t4 * t50) + t77 * t7; -pkin(4) * t10 + t47 * t6 / 0.2e1 + t5 * t73 - t17 * mrSges(5,2) + t57 * mrSges(6,3) + (t31 * t73 - t47 * t30 / 0.2e1 + Ifges(5,5)) * t26 + (t61 / 0.2e1 - Ifges(5,6)) * t25 + t77 * t15 + (m(6) * t57 + t50 * t12 - t47 * t13) * pkin(8); 0; Ifges(5,3) + t50 * t30 - 0.2e1 * pkin(4) * t29 + t47 * t31 + m(6) * (t59 * pkin(8) ^ 2 + pkin(4) ^ 2) + 0.2e1 * t59 * pkin(8) * mrSges(6,3); t3 * mrSges(6,1) - t4 * mrSges(6,2); t1 * mrSges(6,1) - t2 * mrSges(6,2) - Ifges(6,6) * t67 + t62; -t29; -t55 * pkin(8) + t61; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t8(1), t8(2), t8(4), t8(7), t8(11); t8(2), t8(3), t8(5), t8(8), t8(12); t8(4), t8(5), t8(6), t8(9), t8(13); t8(7), t8(8), t8(9), t8(10), t8(14); t8(11), t8(12), t8(13), t8(14), t8(15);];
Mq = res;
