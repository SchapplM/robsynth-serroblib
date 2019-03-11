% Calculate joint inertia matrix for
% S6PRPRPR2
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
% Datum: 2019-03-08 19:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRPR2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR2_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR2_inertiaJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR2_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR2_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR2_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:29:37
% EndTime: 2019-03-08 19:29:39
% DurationCPUTime: 0.70s
% Computational Cost: add. (907->225), mult. (2080->328), div. (0->0), fcn. (2140->12), ass. (0->94)
t70 = sin(pkin(11));
t59 = t70 * pkin(2) + pkin(8);
t110 = 0.2e1 * t59;
t69 = sin(pkin(12));
t72 = cos(pkin(12));
t75 = sin(qJ(6));
t78 = cos(qJ(6));
t43 = t78 * t69 + t75 * t72;
t76 = sin(qJ(4));
t32 = t43 * t76;
t42 = -t75 * t69 + t78 * t72;
t33 = t42 * t76;
t10 = t32 * mrSges(7,1) + t33 * mrSges(7,2);
t95 = t72 * t76;
t96 = t69 * t76;
t37 = mrSges(6,1) * t96 + mrSges(6,2) * t95;
t109 = t10 + t37;
t36 = (pkin(5) * t69 + t59) * t76;
t108 = -m(7) * t36 - t109;
t71 = sin(pkin(6));
t73 = cos(pkin(11));
t77 = sin(qJ(2));
t80 = cos(qJ(2));
t31 = (t70 * t80 + t73 * t77) * t71;
t74 = cos(pkin(6));
t79 = cos(qJ(4));
t18 = t31 * t76 - t74 * t79;
t17 = t18 ^ 2;
t29 = (t70 * t77 - t73 * t80) * t71;
t107 = t29 ^ 2;
t106 = m(6) * pkin(4);
t105 = -t69 / 0.2e1;
t104 = t72 / 0.2e1;
t103 = m(6) + m(7);
t102 = t79 * pkin(4);
t101 = Ifges(6,4) * t69;
t100 = Ifges(6,4) * t72;
t20 = t31 * t79 + t74 * t76;
t99 = t20 * t79;
t98 = t59 * t76;
t97 = t59 * t79;
t94 = t76 * mrSges(6,3);
t93 = t79 * t18;
t92 = pkin(9) + qJ(5);
t91 = -Ifges(7,5) * t33 + Ifges(7,6) * t32;
t60 = -t73 * pkin(2) - pkin(3);
t41 = -t76 * qJ(5) - t102 + t60;
t13 = t69 * t41 + t72 * t97;
t48 = -t72 * mrSges(6,1) + t69 * mrSges(6,2);
t90 = t48 - mrSges(5,1);
t89 = t69 ^ 2 + t72 ^ 2;
t67 = t76 ^ 2;
t68 = t79 ^ 2;
t88 = t67 + t68;
t14 = -t42 * mrSges(7,1) + t43 * mrSges(7,2);
t87 = -t14 - t90;
t86 = qJ(5) * t89;
t5 = -t20 * t69 + t29 * t72;
t6 = t20 * t72 + t29 * t69;
t85 = -t5 * t69 + t6 * t72;
t35 = t72 * t41;
t12 = -t69 * t97 + t35;
t84 = -t12 * t69 + t13 * t72;
t44 = t79 * mrSges(6,2) - t69 * t94;
t45 = -t79 * mrSges(6,1) - t72 * t94;
t83 = t72 * t44 - t69 * t45;
t66 = t74 ^ 2;
t61 = -t72 * pkin(5) - pkin(4);
t57 = t59 ^ 2;
t53 = t67 * t57;
t52 = -t79 * mrSges(5,1) + t76 * mrSges(5,2);
t51 = Ifges(6,1) * t69 + t100;
t50 = Ifges(6,2) * t72 + t101;
t49 = t92 * t72;
t47 = t92 * t69;
t40 = Ifges(7,5) * t43;
t39 = Ifges(7,6) * t42;
t28 = -Ifges(6,5) * t79 + (Ifges(6,1) * t72 - t101) * t76;
t27 = -Ifges(6,6) * t79 + (-Ifges(6,2) * t69 + t100) * t76;
t24 = -t79 * mrSges(7,1) - t33 * mrSges(7,3);
t23 = t79 * mrSges(7,2) - t32 * mrSges(7,3);
t22 = -t75 * t47 + t78 * t49;
t21 = -t78 * t47 - t75 * t49;
t16 = Ifges(7,1) * t43 + Ifges(7,4) * t42;
t15 = Ifges(7,4) * t43 + Ifges(7,2) * t42;
t11 = -pkin(9) * t96 + t13;
t9 = -pkin(9) * t95 + t35 + (-t59 * t69 - pkin(5)) * t79;
t8 = Ifges(7,1) * t33 - Ifges(7,4) * t32 - Ifges(7,5) * t79;
t7 = Ifges(7,4) * t33 - Ifges(7,2) * t32 - Ifges(7,6) * t79;
t4 = t78 * t11 + t75 * t9;
t3 = -t75 * t11 + t78 * t9;
t2 = t75 * t5 + t78 * t6;
t1 = t78 * t5 - t75 * t6;
t19 = [m(2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t17) + m(5) * (t20 ^ 2 + t107 + t17) + m(6) * (t5 ^ 2 + t6 ^ 2 + t17) + m(4) * (t31 ^ 2 + t107 + t66) + m(3) * (t66 + (t77 ^ 2 + t80 ^ 2) * t71 ^ 2); mrSges(5,3) * t99 - t31 * mrSges(4,2) + t1 * t24 + t2 * t23 + t6 * t44 + t5 * t45 + (t80 * mrSges(3,1) - t77 * mrSges(3,2)) * t71 + (-mrSges(4,1) + t52) * t29 + (t76 * mrSges(5,3) + t109) * t18 + m(7) * (t3 * t1 + t36 * t18 + t4 * t2) + m(5) * (t60 * t29 + (t18 * t76 + t99) * t59) + m(6) * (t12 * t5 + t13 * t6 + t18 * t98) + m(4) * (-t29 * t73 + t31 * t70) * pkin(2); 0.2e1 * t36 * t10 + 0.2e1 * t12 * t45 + 0.2e1 * t13 * t44 + 0.2e1 * t4 * t23 + 0.2e1 * t3 * t24 - t32 * t7 + t33 * t8 + 0.2e1 * t60 * t52 + Ifges(3,3) + Ifges(4,3) + ((Ifges(6,3) + Ifges(5,2) + Ifges(7,3)) * t79 + t91) * t79 + m(5) * (t68 * t57 + t60 ^ 2 + t53) + m(7) * (t3 ^ 2 + t36 ^ 2 + t4 ^ 2) + m(6) * (t12 ^ 2 + t13 ^ 2 + t53) + m(4) * (t70 ^ 2 + t73 ^ 2) * pkin(2) ^ 2 + (Ifges(5,1) * t76 - t69 * t27 + t72 * t28 + t37 * t110 + (-Ifges(6,5) * t72 + Ifges(6,6) * t69 + (2 * Ifges(5,4))) * t79) * t76 + 0.2e1 * (t73 * mrSges(4,1) - t70 * mrSges(4,2)) * pkin(2) + t88 * mrSges(5,3) * t110; m(4) * t74 + m(7) * (-t32 * t1 + t33 * t2 - t93) + m(5) * (t76 * t20 - t93) + m(6) * (t85 * t76 - t93); m(7) * (-t32 * t3 + t33 * t4) - t32 * t24 + t33 * t23 + (m(6) * (t84 - t97) + t83) * t76 + t108 * t79; m(4) + m(5) * t88 + m(6) * (t89 * t67 + t68) + m(7) * (t32 ^ 2 + t33 ^ 2 + t68); -t20 * mrSges(5,2) + (-t1 * t43 + t2 * t42) * mrSges(7,3) + t85 * mrSges(6,3) - t87 * t18 + m(7) * (t21 * t1 + t61 * t18 + t22 * t2) + m(6) * (-pkin(4) * t18 + t85 * qJ(5)); t27 * t104 + t69 * t28 / 0.2e1 + t61 * t10 + t33 * t16 / 0.2e1 + t36 * t14 - pkin(4) * t37 + t42 * t7 / 0.2e1 + t43 * t8 / 0.2e1 + t21 * t24 - t32 * t15 / 0.2e1 + t22 * t23 + m(7) * (t21 * t3 + t22 * t4 + t61 * t36) + (-t59 * mrSges(5,2) + Ifges(6,5) * t105 - Ifges(6,6) * t72 / 0.2e1 - t40 / 0.2e1 - t39 / 0.2e1 + Ifges(5,6)) * t79 + (-t3 * t43 + t4 * t42) * mrSges(7,3) + t84 * mrSges(6,3) + (m(6) * t84 + t83) * qJ(5) + (t51 * t104 + t50 * t105 + Ifges(5,5) + (t90 - t106) * t59) * t76; (t32 * t43 + t33 * t42) * mrSges(7,3) + t87 * t79 + (t89 * mrSges(6,3) - mrSges(5,2)) * t76 + m(6) * (t76 * t86 + t102) + m(7) * (-t21 * t32 + t22 * t33 - t61 * t79); -0.2e1 * pkin(4) * t48 + 0.2e1 * t61 * t14 + t42 * t15 + t43 * t16 + t72 * t50 + t69 * t51 + Ifges(5,3) + m(7) * (t21 ^ 2 + t22 ^ 2 + t61 ^ 2) + m(6) * (t89 * qJ(5) ^ 2 + pkin(4) ^ 2) + 0.2e1 * (-t21 * t43 + t22 * t42) * mrSges(7,3) + 0.2e1 * mrSges(6,3) * t86; 0.2e1 * (m(7) / 0.2e1 + m(6) / 0.2e1) * t18; m(6) * t98 - t108; -t103 * t79; m(7) * t61 - t106 + t14 + t48; t103; t1 * mrSges(7,1) - t2 * mrSges(7,2); t3 * mrSges(7,1) - t4 * mrSges(7,2) - Ifges(7,3) * t79 - t91; -t10; t21 * mrSges(7,1) - t22 * mrSges(7,2) + t39 + t40; 0; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t19(1) t19(2) t19(4) t19(7) t19(11) t19(16); t19(2) t19(3) t19(5) t19(8) t19(12) t19(17); t19(4) t19(5) t19(6) t19(9) t19(13) t19(18); t19(7) t19(8) t19(9) t19(10) t19(14) t19(19); t19(11) t19(12) t19(13) t19(14) t19(15) t19(20); t19(16) t19(17) t19(18) t19(19) t19(20) t19(21);];
Mq  = res;
