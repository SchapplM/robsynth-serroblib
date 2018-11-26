% Calculate joint inertia matrix for
% S6PRPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
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
% Datum: 2018-11-23 14:55
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6PRPRPR1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_inertiaJ_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR1_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR1_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR1_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 14:55:05
% EndTime: 2018-11-23 14:55:06
% DurationCPUTime: 0.58s
% Computational Cost: add. (868->188), mult. (1912->284), div. (0->0), fcn. (2042->12), ass. (0->76)
t58 = cos(pkin(12));
t97 = t58 * pkin(4);
t61 = sin(qJ(6));
t64 = cos(qJ(6));
t35 = -t64 * mrSges(7,1) + t61 * mrSges(7,2);
t80 = -mrSges(6,1) + t35;
t55 = sin(pkin(12));
t62 = sin(qJ(4));
t65 = cos(qJ(4));
t31 = t55 * t62 - t58 * t65;
t33 = t55 * t65 + t58 * t62;
t84 = t33 * t61;
t15 = -t31 * mrSges(7,2) - mrSges(7,3) * t84;
t83 = t33 * t64;
t16 = t31 * mrSges(7,1) - mrSges(7,3) * t83;
t96 = t64 * t15 - t61 * t16;
t56 = sin(pkin(11));
t57 = sin(pkin(6));
t59 = cos(pkin(11));
t63 = sin(qJ(2));
t66 = cos(qJ(2));
t24 = (t56 * t66 + t59 * t63) * t57;
t60 = cos(pkin(6));
t19 = -t62 * t24 + t60 * t65;
t20 = t65 * t24 + t60 * t62;
t5 = -t58 * t19 + t55 * t20;
t95 = t5 ^ 2;
t44 = t56 * pkin(2) + pkin(8);
t75 = qJ(5) + t44;
t29 = t75 * t65;
t73 = t75 * t62;
t11 = t55 * t29 + t58 * t73;
t94 = t11 ^ 2;
t22 = (t56 * t63 - t59 * t66) * t57;
t21 = t22 ^ 2;
t93 = t31 ^ 2;
t54 = t65 ^ 2;
t92 = 0.2e1 * t11;
t91 = m(6) * pkin(4);
t90 = t64 / 0.2e1;
t89 = t11 * t5;
t88 = t31 * t5;
t87 = Ifges(7,4) * t61;
t86 = Ifges(7,4) * t64;
t85 = t31 * t11;
t79 = Ifges(7,5) * t83 + Ifges(7,3) * t31;
t78 = Ifges(7,5) * t61 + Ifges(7,6) * t64;
t77 = t61 ^ 2 + t64 ^ 2;
t76 = t62 ^ 2 + t54;
t46 = -t59 * pkin(2) - pkin(3);
t43 = t55 * pkin(4) + pkin(9);
t74 = t77 * t43;
t27 = t33 * mrSges(6,2);
t18 = t31 * mrSges(6,1) + t27;
t7 = t55 * t19 + t58 * t20;
t1 = t64 * t22 - t61 * t7;
t2 = t61 * t22 + t64 * t7;
t72 = -t1 * t61 + t2 * t64;
t34 = -t65 * pkin(4) + t46;
t10 = t31 * pkin(5) - t33 * pkin(9) + t34;
t13 = t58 * t29 - t55 * t73;
t3 = t64 * t10 - t61 * t13;
t4 = t61 * t10 + t64 * t13;
t71 = -t3 * t61 + t4 * t64;
t36 = -t65 * mrSges(5,1) + t62 * mrSges(5,2);
t70 = mrSges(7,1) * t61 + mrSges(7,2) * t64;
t69 = -t19 * t62 + t20 * t65;
t50 = t60 ^ 2;
t45 = -pkin(5) - t97;
t38 = Ifges(7,1) * t61 + t86;
t37 = Ifges(7,2) * t64 + t87;
t30 = t33 ^ 2;
t14 = t70 * t33;
t9 = Ifges(7,5) * t31 + (Ifges(7,1) * t64 - t87) * t33;
t8 = Ifges(7,6) * t31 + (-Ifges(7,2) * t61 + t86) * t33;
t6 = [m(2) + m(6) * (t7 ^ 2 + t21 + t95) + m(7) * (t1 ^ 2 + t2 ^ 2 + t95) + m(5) * (t19 ^ 2 + t20 ^ 2 + t21) + m(4) * (t24 ^ 2 + t21 + t50) + m(3) * (t50 + (t63 ^ 2 + t66 ^ 2) * t57 ^ 2); -t24 * mrSges(4,2) + t1 * t16 + t5 * t14 + t2 * t15 + (t66 * mrSges(3,1) - t63 * mrSges(3,2)) * t57 + (-t7 * t31 + t5 * t33) * mrSges(6,3) + t69 * mrSges(5,3) + (-mrSges(4,1) + t18 + t36) * t22 + m(6) * (t13 * t7 + t34 * t22 + t89) + m(7) * (t3 * t1 + t4 * t2 + t89) + m(5) * (t46 * t22 + t69 * t44) + m(4) * (-t22 * t59 + t24 * t56) * pkin(2); Ifges(5,2) * t54 + t14 * t92 + 0.2e1 * t4 * t15 + 0.2e1 * t3 * t16 + 0.2e1 * t34 * t18 + 0.2e1 * t46 * t36 + Ifges(3,3) + Ifges(4,3) + (Ifges(5,1) * t62 + 0.2e1 * Ifges(5,4) * t65) * t62 + (-0.2e1 * t13 * mrSges(6,3) + Ifges(6,2) * t31 + t79) * t31 + (mrSges(6,3) * t92 + Ifges(6,1) * t33 - t61 * t8 + t64 * t9 + (-Ifges(7,6) * t61 - (2 * Ifges(6,4))) * t31) * t33 + m(7) * (t3 ^ 2 + t4 ^ 2 + t94) + m(6) * (t13 ^ 2 + t34 ^ 2 + t94) + m(5) * (t76 * t44 ^ 2 + t46 ^ 2) + m(4) * (t56 ^ 2 + t59 ^ 2) * pkin(2) ^ 2 + 0.2e1 * (t59 * mrSges(4,1) - t56 * mrSges(4,2)) * pkin(2) + 0.2e1 * t76 * t44 * mrSges(5,3); m(4) * t60 + m(6) * (t33 * t7 + t88) + m(7) * (t72 * t33 + t88) + m(5) * (t65 * t19 + t62 * t20); t31 * t14 + t96 * t33 + m(7) * (t71 * t33 + t85) + m(6) * (t33 * t13 + t85); m(4) + m(5) * t76 + m(6) * (t30 + t93) + m(7) * (t77 * t30 + t93); t19 * mrSges(5,1) - t20 * mrSges(5,2) - t7 * mrSges(6,2) + t80 * t5 + t72 * mrSges(7,3) + m(7) * (t72 * t43 + t45 * t5) + (-t5 * t58 + t55 * t7) * t91; -Ifges(6,6) * t31 + Ifges(5,5) * t62 + Ifges(5,6) * t65 + t45 * t14 + t31 * t78 / 0.2e1 + t61 * t9 / 0.2e1 + t8 * t90 - t13 * mrSges(6,2) + (-t62 * mrSges(5,1) - t65 * mrSges(5,2)) * t44 + t71 * mrSges(7,3) + (Ifges(6,5) + t38 * t90 - t61 * t37 / 0.2e1) * t33 + (m(6) * t13 * t55 + (-t55 * t31 - t58 * t33) * mrSges(6,3)) * pkin(4) + (-m(6) * t97 + m(7) * t45 + t80) * t11 + (m(7) * t71 + t96) * t43; -t27 + t80 * t31 + t77 * t33 * mrSges(7,3) + m(7) * (t45 * t31 + t33 * t74) + (-t31 * t58 + t33 * t55) * t91 - t36; 0.2e1 * t45 * t35 + t64 * t37 + t61 * t38 + Ifges(5,3) + Ifges(6,3) + m(7) * (t77 * t43 ^ 2 + t45 ^ 2) + m(6) * (t55 ^ 2 + t58 ^ 2) * pkin(4) ^ 2 + 0.2e1 * (t58 * mrSges(6,1) - t55 * mrSges(6,2)) * pkin(4) + 0.2e1 * mrSges(7,3) * t74; m(6) * t22 + m(7) * (t64 * t1 + t61 * t2); t61 * t15 + t64 * t16 + m(7) * (t64 * t3 + t61 * t4) + m(6) * t34 + t18; 0; 0; m(7) * t77 + m(6); t1 * mrSges(7,1) - t2 * mrSges(7,2); t3 * mrSges(7,1) - t4 * mrSges(7,2) - Ifges(7,6) * t84 + t79; -t14; -t70 * t43 + t78; -t35; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t6(1) t6(2) t6(4) t6(7) t6(11) t6(16); t6(2) t6(3) t6(5) t6(8) t6(12) t6(17); t6(4) t6(5) t6(6) t6(9) t6(13) t6(18); t6(7) t6(8) t6(9) t6(10) t6(14) t6(19); t6(11) t6(12) t6(13) t6(14) t6(15) t6(20); t6(16) t6(17) t6(18) t6(19) t6(20) t6(21);];
Mq  = res;
