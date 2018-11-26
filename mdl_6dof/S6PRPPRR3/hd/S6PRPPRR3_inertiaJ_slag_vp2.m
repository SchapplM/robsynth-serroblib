% Calculate joint inertia matrix for
% S6PRPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 14:54
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6PRPPRR3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR3_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPPRR3_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPPRR3_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 14:54:38
% EndTime: 2018-11-23 14:54:39
% DurationCPUTime: 0.46s
% Computational Cost: add. (562->173), mult. (1177->253), div. (0->0), fcn. (1080->10), ass. (0->79)
t91 = -m(4) * pkin(2) - mrSges(4,1);
t90 = m(7) * pkin(9) + mrSges(7,3);
t53 = sin(qJ(6));
t56 = cos(qJ(6));
t27 = -t56 * mrSges(7,1) + t53 * mrSges(7,2);
t89 = -m(7) * pkin(5) - mrSges(6,1) + t27;
t49 = sin(pkin(11));
t50 = sin(pkin(6));
t51 = cos(pkin(11));
t55 = sin(qJ(2));
t58 = cos(qJ(2));
t11 = (-t49 * t58 + t51 * t55) * t50;
t52 = cos(pkin(6));
t54 = sin(qJ(5));
t57 = cos(qJ(5));
t5 = t54 * t11 + t52 * t57;
t88 = t5 ^ 2;
t9 = (-t49 * t55 - t51 * t58) * t50;
t87 = t9 ^ 2;
t86 = t53 / 0.2e1;
t85 = t5 * t54;
t84 = t51 * t9;
t83 = t57 * t5;
t82 = Ifges(7,4) * t53;
t81 = Ifges(7,4) * t56;
t80 = Ifges(7,5) * t57;
t59 = -pkin(2) - pkin(3);
t26 = t51 * qJ(3) + t49 * t59;
t22 = -pkin(8) + t26;
t79 = t22 * t49;
t78 = t49 * t54;
t77 = t49 * t57;
t76 = t53 * t54;
t75 = t53 * t57;
t74 = t54 * t56;
t73 = t56 * t57;
t28 = t57 * mrSges(6,1) - t54 * mrSges(6,2);
t72 = mrSges(5,1) + t28;
t71 = Ifges(7,6) * t76 + Ifges(7,3) * t57;
t70 = t53 ^ 2 + t56 ^ 2;
t46 = t54 ^ 2;
t48 = t57 ^ 2;
t69 = t46 + t48;
t68 = t70 * t54;
t67 = t69 * mrSges(6,3);
t25 = -t49 * qJ(3) + t51 * t59;
t21 = pkin(4) - t25;
t7 = t57 * t11 - t52 * t54;
t1 = -t53 * t7 - t56 * t9;
t2 = -t53 * t9 + t56 * t7;
t66 = -t1 * t53 + t2 * t56;
t40 = t57 * pkin(5);
t8 = t54 * pkin(9) + t21 + t40;
t3 = -t22 * t75 + t56 * t8;
t4 = t22 * t73 + t53 * t8;
t65 = -t3 * t53 + t4 * t56;
t64 = t57 * t7 + t85;
t63 = -mrSges(7,1) * t53 - mrSges(7,2) * t56;
t16 = -t49 * t75 - t56 * t51;
t17 = t49 * t73 - t53 * t51;
t62 = -t16 * t53 + t17 * t56;
t23 = -t57 * mrSges(7,2) + mrSges(7,3) * t76;
t24 = t57 * mrSges(7,1) + mrSges(7,3) * t74;
t61 = t56 * t23 - t53 * t24;
t44 = t52 ^ 2;
t43 = t51 ^ 2;
t41 = t49 ^ 2;
t39 = Ifges(7,5) * t53;
t38 = Ifges(7,6) * t56;
t32 = t46 * t41;
t30 = Ifges(7,1) * t53 + t81;
t29 = Ifges(7,2) * t56 + t82;
t20 = t22 ^ 2;
t19 = t63 * t54;
t18 = t46 * t20;
t14 = t46 * t79;
t13 = t80 + (-Ifges(7,1) * t56 + t82) * t54;
t12 = Ifges(7,6) * t57 + (Ifges(7,2) * t53 - t81) * t54;
t6 = [m(2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t88) + m(6) * (t7 ^ 2 + t87 + t88) + m(5) * (t11 ^ 2 + t44 + t87) + 0.2e1 * (m(4) / 0.2e1 + m(3) / 0.2e1) * (t44 + (t55 ^ 2 + t58 ^ 2) * t50 ^ 2); t11 * mrSges(5,2) + t1 * t24 + t5 * t19 + t2 * t23 - t72 * t9 - t64 * mrSges(6,3) + m(7) * (t3 * t1 + t4 * t2 + t22 * t85) + m(6) * (-t21 * t9 + t22 * t64) + m(5) * (t26 * t11 + t25 * t9) + ((mrSges(3,1) - t91) * t58 + (m(4) * qJ(3) - mrSges(3,2) + mrSges(4,3)) * t55) * t50; (2 * pkin(2) * mrSges(4,1)) - 0.2e1 * t25 * mrSges(5,1) + 0.2e1 * t26 * mrSges(5,2) + 0.2e1 * qJ(3) * mrSges(4,3) + 0.2e1 * t21 * t28 + 0.2e1 * t4 * t23 + 0.2e1 * t3 * t24 + Ifges(4,2) + Ifges(3,3) + Ifges(5,3) + (Ifges(6,2) * t57 + t71) * t57 - 0.2e1 * t22 * t67 + (Ifges(6,1) * t54 + 0.2e1 * Ifges(6,4) * t57 + t53 * t12 + 0.2e1 * t22 * t19 + (-t13 - t80) * t56) * t54 + m(7) * (t3 ^ 2 + t4 ^ 2 + t18) + m(6) * (t48 * t20 + t21 ^ 2 + t18) + m(5) * (t25 ^ 2 + t26 ^ 2) + m(4) * ((pkin(2) ^ 2) + qJ(3) ^ 2); -m(4) * t50 * t58 + m(7) * (t16 * t1 + t17 * t2 + t5 * t78) + m(6) * (t49 * t64 + t84) + m(5) * (t49 * t11 + t84); t16 * t24 + t17 * t23 - t72 * t51 + (t54 * t19 + mrSges(5,2) - t67) * t49 + m(7) * (t16 * t3 + t17 * t4 + t14) + m(6) * (-t51 * t21 + t48 * t79 + t14) + m(5) * (t51 * t25 + t49 * t26) + t91; m(4) + m(5) * (t41 + t43) + m(6) * (t48 * t41 + t32 + t43) + m(7) * (t16 ^ 2 + t17 ^ 2 + t32); -m(5) * t52 + m(7) * (t54 * t66 - t83) + m(6) * (t54 * t7 - t83); -t57 * t19 + (m(7) * (-t22 * t57 + t65) + t61) * t54; m(7) * (t62 - t77) * t54; m(5) + m(6) * t69 + m(7) * (t46 * t70 + t48); -t7 * mrSges(6,2) + t89 * t5 + t90 * t66; t56 * t12 / 0.2e1 + t13 * t86 - pkin(5) * t19 + (-t22 * mrSges(6,2) + t39 / 0.2e1 + t38 / 0.2e1 - Ifges(6,6)) * t57 + t65 * mrSges(7,3) + (m(7) * t65 + t61) * pkin(9) + (-t56 * t30 / 0.2e1 + t29 * t86 - Ifges(6,5) + t89 * t22) * t54; -mrSges(6,2) * t77 + t90 * t62 + t89 * t78; -t57 * t27 + m(7) * (pkin(9) * t68 + t40) + mrSges(7,3) * t68 + t28; Ifges(6,3) - 0.2e1 * pkin(5) * t27 + m(7) * (pkin(9) ^ 2 * t70 + pkin(5) ^ 2) + t53 * t30 + t56 * t29 + 0.2e1 * t70 * pkin(9) * mrSges(7,3); t1 * mrSges(7,1) - t2 * mrSges(7,2); t3 * mrSges(7,1) - t4 * mrSges(7,2) - Ifges(7,5) * t74 + t71; t16 * mrSges(7,1) - t17 * mrSges(7,2); t19; pkin(9) * t63 + t38 + t39; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t6(1) t6(2) t6(4) t6(7) t6(11) t6(16); t6(2) t6(3) t6(5) t6(8) t6(12) t6(17); t6(4) t6(5) t6(6) t6(9) t6(13) t6(18); t6(7) t6(8) t6(9) t6(10) t6(14) t6(19); t6(11) t6(12) t6(13) t6(14) t6(15) t6(20); t6(16) t6(17) t6(18) t6(19) t6(20) t6(21);];
Mq  = res;
