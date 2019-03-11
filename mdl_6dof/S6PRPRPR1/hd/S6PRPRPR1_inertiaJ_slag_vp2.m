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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRPR1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_inertiaJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR1_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR1_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR1_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:24:46
% EndTime: 2019-03-08 19:24:47
% DurationCPUTime: 0.58s
% Computational Cost: add. (868->188), mult. (1912->286), div. (0->0), fcn. (2042->12), ass. (0->74)
t60 = cos(pkin(12));
t97 = t60 * pkin(4);
t63 = sin(qJ(6));
t66 = cos(qJ(6));
t37 = -t66 * mrSges(7,1) + t63 * mrSges(7,2);
t81 = -mrSges(6,1) + t37;
t57 = sin(pkin(12));
t64 = sin(qJ(4));
t67 = cos(qJ(4));
t33 = t57 * t64 - t60 * t67;
t35 = t57 * t67 + t60 * t64;
t85 = t35 * t63;
t15 = -t33 * mrSges(7,2) - mrSges(7,3) * t85;
t84 = t35 * t66;
t16 = t33 * mrSges(7,1) - mrSges(7,3) * t84;
t96 = t66 * t15 - t63 * t16;
t58 = sin(pkin(11));
t59 = sin(pkin(6));
t61 = cos(pkin(11));
t65 = sin(qJ(2));
t68 = cos(qJ(2));
t24 = (t58 * t68 + t61 * t65) * t59;
t62 = cos(pkin(6));
t19 = -t24 * t64 + t62 * t67;
t20 = t24 * t67 + t62 * t64;
t5 = -t60 * t19 + t57 * t20;
t95 = t5 ^ 2;
t46 = t58 * pkin(2) + pkin(8);
t76 = qJ(5) + t46;
t31 = t76 * t67;
t74 = t76 * t64;
t11 = t57 * t31 + t60 * t74;
t94 = t11 ^ 2;
t22 = (t58 * t65 - t61 * t68) * t59;
t21 = t22 ^ 2;
t93 = t33 ^ 2;
t56 = t67 ^ 2;
t92 = m(6) * pkin(4);
t91 = t66 / 0.2e1;
t90 = t11 * t5;
t89 = t33 * t5;
t88 = Ifges(7,4) * t63;
t87 = Ifges(7,4) * t66;
t86 = t33 * t11;
t80 = Ifges(7,5) * t84 + Ifges(7,3) * t33;
t79 = Ifges(7,5) * t63 + Ifges(7,6) * t66;
t78 = t63 ^ 2 + t66 ^ 2;
t77 = t64 ^ 2 + t56;
t48 = -t61 * pkin(2) - pkin(3);
t45 = t57 * pkin(4) + pkin(9);
t75 = t78 * t45;
t29 = t35 * mrSges(6,2);
t18 = t33 * mrSges(6,1) + t29;
t7 = t57 * t19 + t60 * t20;
t1 = t22 * t66 - t63 * t7;
t2 = t22 * t63 + t66 * t7;
t73 = -t1 * t63 + t2 * t66;
t36 = -t67 * pkin(4) + t48;
t10 = t33 * pkin(5) - t35 * pkin(9) + t36;
t13 = t60 * t31 - t57 * t74;
t3 = t66 * t10 - t63 * t13;
t4 = t63 * t10 + t66 * t13;
t72 = -t3 * t63 + t4 * t66;
t38 = -t67 * mrSges(5,1) + t64 * mrSges(5,2);
t71 = -t19 * t64 + t20 * t67;
t52 = t62 ^ 2;
t47 = -pkin(5) - t97;
t40 = Ifges(7,1) * t63 + t87;
t39 = Ifges(7,2) * t66 + t88;
t32 = t35 ^ 2;
t14 = -mrSges(7,1) * t85 - mrSges(7,2) * t84;
t9 = Ifges(7,5) * t33 + (Ifges(7,1) * t66 - t88) * t35;
t8 = Ifges(7,6) * t33 + (-Ifges(7,2) * t63 + t87) * t35;
t6 = [m(2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t95) + m(5) * (t19 ^ 2 + t20 ^ 2 + t21) + m(6) * (t7 ^ 2 + t21 + t95) + m(4) * (t24 ^ 2 + t21 + t52) + m(3) * (t52 + (t65 ^ 2 + t68 ^ 2) * t59 ^ 2); -t24 * mrSges(4,2) + t1 * t16 - t5 * t14 + t2 * t15 + (t68 * mrSges(3,1) - t65 * mrSges(3,2)) * t59 + (-t7 * t33 + t5 * t35) * mrSges(6,3) + t71 * mrSges(5,3) + (-mrSges(4,1) + t18 + t38) * t22 + m(7) * (t3 * t1 + t4 * t2 + t90) + m(5) * (t48 * t22 + t71 * t46) + m(6) * (t13 * t7 + t36 * t22 + t90) + m(4) * (-t22 * t61 + t24 * t58) * pkin(2); Ifges(5,2) * t56 - 0.2e1 * t11 * t14 + 0.2e1 * t4 * t15 + 0.2e1 * t3 * t16 + 0.2e1 * t36 * t18 + 0.2e1 * t48 * t38 + Ifges(3,3) + Ifges(4,3) + (Ifges(5,1) * t64 + 0.2e1 * Ifges(5,4) * t67) * t64 + (-0.2e1 * t13 * mrSges(6,3) + Ifges(6,2) * t33 + t80) * t33 + (0.2e1 * t11 * mrSges(6,3) + Ifges(6,1) * t35 - t63 * t8 + t66 * t9 + (-Ifges(7,6) * t63 - (2 * Ifges(6,4))) * t33) * t35 + m(7) * (t3 ^ 2 + t4 ^ 2 + t94) + m(6) * (t13 ^ 2 + t36 ^ 2 + t94) + m(5) * (t77 * t46 ^ 2 + t48 ^ 2) + m(4) * (t58 ^ 2 + t61 ^ 2) * pkin(2) ^ 2 + 0.2e1 * (t61 * mrSges(4,1) - t58 * mrSges(4,2)) * pkin(2) + 0.2e1 * t77 * t46 * mrSges(5,3); m(4) * t62 + m(7) * (t73 * t35 + t89) + m(5) * (t67 * t19 + t64 * t20) + m(6) * (t35 * t7 + t89); -t33 * t14 + t96 * t35 + m(7) * (t72 * t35 + t86) + m(6) * (t35 * t13 + t86); m(4) + m(5) * t77 + m(6) * (t32 + t93) + m(7) * (t78 * t32 + t93); t19 * mrSges(5,1) - t20 * mrSges(5,2) - t7 * mrSges(6,2) + t81 * t5 + t73 * mrSges(7,3) + m(7) * (t73 * t45 + t47 * t5) + (-t5 * t60 + t57 * t7) * t92; Ifges(5,5) * t64 + Ifges(5,6) * t67 - Ifges(6,6) * t33 + t33 * t79 / 0.2e1 + t63 * t9 / 0.2e1 + t8 * t91 - t47 * t14 - t13 * mrSges(6,2) + (-t64 * mrSges(5,1) - t67 * mrSges(5,2)) * t46 + t72 * mrSges(7,3) + (Ifges(6,5) - t63 * t39 / 0.2e1 + t40 * t91) * t35 + (m(6) * t13 * t57 + (-t57 * t33 - t60 * t35) * mrSges(6,3)) * pkin(4) + (-m(6) * t97 + m(7) * t47 + t81) * t11 + (m(7) * t72 + t96) * t45; -t29 + t81 * t33 + t78 * t35 * mrSges(7,3) + m(7) * (t47 * t33 + t35 * t75) + (-t33 * t60 + t35 * t57) * t92 - t38; 0.2e1 * t47 * t37 + t66 * t39 + t63 * t40 + Ifges(5,3) + Ifges(6,3) + m(7) * (t78 * t45 ^ 2 + t47 ^ 2) + m(6) * (t57 ^ 2 + t60 ^ 2) * pkin(4) ^ 2 + 0.2e1 * (t60 * mrSges(6,1) - t57 * mrSges(6,2)) * pkin(4) + 0.2e1 * mrSges(7,3) * t75; m(7) * (t66 * t1 + t63 * t2) + m(6) * t22; t63 * t15 + t66 * t16 + m(7) * (t66 * t3 + t63 * t4) + m(6) * t36 + t18; 0; 0; m(7) * t78 + m(6); t1 * mrSges(7,1) - t2 * mrSges(7,2); t3 * mrSges(7,1) - t4 * mrSges(7,2) - Ifges(7,6) * t85 + t80; t14; (-mrSges(7,1) * t63 - mrSges(7,2) * t66) * t45 + t79; -t37; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t6(1) t6(2) t6(4) t6(7) t6(11) t6(16); t6(2) t6(3) t6(5) t6(8) t6(12) t6(17); t6(4) t6(5) t6(6) t6(9) t6(13) t6(18); t6(7) t6(8) t6(9) t6(10) t6(14) t6(19); t6(11) t6(12) t6(13) t6(14) t6(15) t6(20); t6(16) t6(17) t6(18) t6(19) t6(20) t6(21);];
Mq  = res;
