% Calculate joint inertia matrix for
% S6PRPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
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
% Datum: 2018-11-23 15:06
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRR6_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR6_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR6_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR6_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR6_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR6_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:06:36
% EndTime: 2018-11-23 15:06:37
% DurationCPUTime: 0.77s
% Computational Cost: add. (883->243), mult. (1900->351), div. (0->0), fcn. (1836->10), ass. (0->101)
t79 = cos(qJ(5));
t80 = cos(qJ(4));
t105 = t79 * t80;
t76 = sin(qJ(4));
t118 = Ifges(6,5) * t105 + Ifges(6,3) * t76;
t72 = sin(pkin(6));
t81 = cos(qJ(2));
t109 = t72 * t81;
t73 = cos(pkin(6));
t35 = t80 * t109 + t73 * t76;
t34 = t35 ^ 2;
t117 = m(7) * pkin(5);
t116 = t79 / 0.2e1;
t115 = -pkin(10) - pkin(9);
t75 = sin(qJ(5));
t93 = t75 * mrSges(6,1) + t79 * mrSges(6,2);
t38 = t93 * t80;
t74 = sin(qJ(6));
t78 = cos(qJ(6));
t45 = t74 * t79 + t78 * t75;
t30 = t45 * t80;
t44 = -t74 * t75 + t78 * t79;
t32 = t44 * t80;
t9 = t30 * mrSges(7,1) + t32 * mrSges(7,2);
t114 = -t38 - t9;
t113 = Ifges(6,4) * t75;
t112 = Ifges(6,4) * t79;
t111 = Ifges(6,6) * t76;
t77 = sin(qJ(2));
t110 = t72 * t77;
t108 = t75 * t80;
t37 = -t76 * t109 + t73 * t80;
t107 = t76 * t37;
t82 = -pkin(2) - pkin(8);
t106 = t76 * t82;
t104 = t80 * t35;
t49 = -t79 * mrSges(6,1) + t75 * mrSges(6,2);
t103 = t49 - mrSges(5,1);
t102 = t76 * mrSges(5,1) + t80 * mrSges(5,2) + mrSges(4,3);
t48 = t76 * pkin(4) - t80 * pkin(9) + qJ(3);
t23 = t79 * t106 + t75 * t48;
t101 = t75 ^ 2 + t79 ^ 2;
t68 = t76 ^ 2;
t70 = t80 ^ 2;
t100 = -t70 - t68;
t12 = -t44 * mrSges(7,1) + t45 * mrSges(7,2);
t99 = -t12 - t103;
t98 = Ifges(7,5) * t32 - Ifges(7,6) * t30 + Ifges(7,3) * t76;
t15 = t79 * t110 - t75 * t37;
t16 = t75 * t110 + t79 * t37;
t5 = t78 * t15 - t74 * t16;
t6 = t74 * t15 + t78 * t16;
t97 = t5 * mrSges(7,1) - t6 * mrSges(7,2);
t96 = t101 * mrSges(6,3);
t95 = t100 * mrSges(5,3);
t29 = t45 * t76;
t31 = t44 * t76;
t94 = -t29 * mrSges(7,1) - t31 * mrSges(7,2);
t92 = -t15 * t75 + t16 * t79;
t42 = t79 * t48;
t22 = -t75 * t106 + t42;
t91 = -t22 * t75 + t23 * t79;
t90 = -t104 + t107;
t46 = -t76 * mrSges(6,2) - mrSges(6,3) * t108;
t47 = t76 * mrSges(6,1) - mrSges(6,3) * t105;
t89 = t79 * t46 - t75 * t47;
t53 = t115 * t75;
t54 = t115 * t79;
t18 = t78 * t53 + t74 * t54;
t19 = t74 * t53 - t78 * t54;
t39 = Ifges(7,6) * t44;
t40 = Ifges(7,5) * t45;
t88 = t18 * mrSges(7,1) - t19 * mrSges(7,2) + t39 + t40;
t10 = -pkin(10) * t105 + t42 + (-t75 * t82 + pkin(5)) * t76;
t11 = -pkin(10) * t108 + t23;
t2 = t78 * t10 - t74 * t11;
t3 = t74 * t10 + t78 * t11;
t87 = t2 * mrSges(7,1) - t3 * mrSges(7,2) + t98;
t86 = (t78 * mrSges(7,1) - t74 * mrSges(7,2)) * pkin(5);
t83 = qJ(3) ^ 2;
t71 = t82 ^ 2;
t66 = t72 ^ 2;
t65 = Ifges(6,5) * t75;
t64 = Ifges(6,6) * t79;
t61 = t70 * t82;
t60 = t70 * t71;
t59 = -t79 * pkin(5) - pkin(4);
t58 = t66 * t77 ^ 2;
t56 = qJ(3) * t110;
t52 = Ifges(6,1) * t75 + t112;
t51 = Ifges(6,2) * t79 + t113;
t43 = (pkin(5) * t75 - t82) * t80;
t28 = Ifges(6,5) * t76 + (Ifges(6,1) * t79 - t113) * t80;
t27 = t111 + (-Ifges(6,2) * t75 + t112) * t80;
t21 = t76 * mrSges(7,1) - t32 * mrSges(7,3);
t20 = -t76 * mrSges(7,2) - t30 * mrSges(7,3);
t14 = Ifges(7,1) * t45 + Ifges(7,4) * t44;
t13 = Ifges(7,4) * t45 + Ifges(7,2) * t44;
t8 = Ifges(7,1) * t32 - Ifges(7,4) * t30 + Ifges(7,5) * t76;
t7 = Ifges(7,4) * t32 - Ifges(7,2) * t30 + Ifges(7,6) * t76;
t1 = [m(2) + m(7) * (t5 ^ 2 + t6 ^ 2 + t34) + m(6) * (t15 ^ 2 + t16 ^ 2 + t34) + m(5) * (t37 ^ 2 + t34 + t58) + 0.2e1 * (m(4) / 0.2e1 + m(3) / 0.2e1) * (t66 * t81 ^ 2 + t73 ^ 2 + t58); -mrSges(5,3) * t107 + t15 * t47 + t16 * t46 + t6 * t20 + t5 * t21 + (t80 * mrSges(5,3) - t114) * t35 + ((mrSges(3,1) - mrSges(4,2)) * t81 + (-mrSges(3,2) + t102) * t77) * t72 + m(7) * (t2 * t5 + t3 * t6 + t43 * t35) + m(6) * (-t82 * t104 + t22 * t15 + t23 * t16) + m(5) * (t90 * t82 + t56) + m(4) * (pkin(2) * t109 + t56); -0.2e1 * pkin(2) * mrSges(4,2) + 0.2e1 * t2 * t21 + 0.2e1 * t3 * t20 + 0.2e1 * t22 * t47 + 0.2e1 * t23 * t46 - t30 * t7 + t32 * t8 + 0.2e1 * t43 * t9 + Ifges(4,1) + Ifges(3,3) + (Ifges(5,2) * t76 + t118 + t98) * t76 + (Ifges(5,1) * t80 - 0.2e1 * Ifges(5,4) * t76 + t79 * t28 - 0.2e1 * t82 * t38 + (-t27 - t111) * t75) * t80 + m(5) * (t68 * t71 + t60 + t83) + m(4) * (pkin(2) ^ 2 + t83) + m(6) * (t22 ^ 2 + t23 ^ 2 + t60) + m(7) * (t2 ^ 2 + t3 ^ 2 + t43 ^ 2) + 0.2e1 * t102 * qJ(3) + 0.2e1 * t82 * t95; -m(4) * t109 + m(7) * (-t29 * t5 + t31 * t6 - t104) + m(6) * (t92 * t76 - t104) + m(5) * t90; -m(4) * pkin(2) + t31 * t20 - t29 * t21 + mrSges(4,2) + t114 * t80 + t89 * t76 + t95 + m(7) * (-t29 * t2 + t31 * t3 - t80 * t43) + m(6) * (t91 * t76 + t61) + m(5) * (t68 * t82 + t61); m(4) - m(5) * t100 + m(6) * (t101 * t68 + t70) + m(7) * (t29 ^ 2 + t31 ^ 2 + t70); -t37 * mrSges(5,2) + (t6 * t44 - t5 * t45) * mrSges(7,3) + t92 * mrSges(6,3) - t99 * t35 + m(7) * (t18 * t5 + t19 * t6 + t59 * t35) + m(6) * (-pkin(4) * t35 + t92 * pkin(9)); m(7) * (t18 * t2 + t19 * t3 + t59 * t43) + t27 * t116 + t75 * t28 / 0.2e1 + t59 * t9 + t44 * t7 / 0.2e1 + t45 * t8 / 0.2e1 + t32 * t14 / 0.2e1 - pkin(4) * t38 + t43 * t12 - t30 * t13 / 0.2e1 + t19 * t20 + t18 * t21 + (-t82 * mrSges(5,2) + t65 / 0.2e1 + t64 / 0.2e1 + t40 / 0.2e1 + t39 / 0.2e1 - Ifges(5,6)) * t76 + (-t2 * t45 + t3 * t44) * mrSges(7,3) + t91 * mrSges(6,3) + (m(6) * t91 + t89) * pkin(9) + (t52 * t116 - t75 * t51 / 0.2e1 + Ifges(5,5) + (m(6) * pkin(4) - t103) * t82) * t80; (t29 * t45 + t31 * t44) * mrSges(7,3) + t99 * t80 + (-mrSges(5,2) + t96) * t76 + m(6) * (t101 * t76 * pkin(9) + pkin(4) * t80) + m(7) * (-t18 * t29 + t19 * t31 - t59 * t80); -0.2e1 * pkin(4) * t49 + 0.2e1 * t59 * t12 + t44 * t13 + t45 * t14 + t79 * t51 + t75 * t52 + Ifges(5,3) + m(7) * (t18 ^ 2 + t19 ^ 2 + t59 ^ 2) + m(6) * (t101 * pkin(9) ^ 2 + pkin(4) ^ 2) + 0.2e1 * (-t18 * t45 + t19 * t44) * mrSges(7,3) + 0.2e1 * pkin(9) * t96; t15 * mrSges(6,1) - t16 * mrSges(6,2) + (t5 * t78 + t6 * t74) * t117 + t97; -Ifges(6,6) * t108 + t22 * mrSges(6,1) - t23 * mrSges(6,2) + (m(7) * (t2 * t78 + t3 * t74) + t74 * t20 + t78 * t21) * pkin(5) + t87 + t118; -t93 * t76 + (-t29 * t78 + t31 * t74) * t117 + t94; t64 + t65 - t93 * pkin(9) + (m(7) * (t18 * t78 + t19 * t74) + (t74 * t44 - t78 * t45) * mrSges(7,3)) * pkin(5) + t88; Ifges(6,3) + Ifges(7,3) + m(7) * (t74 ^ 2 + t78 ^ 2) * pkin(5) ^ 2 + 0.2e1 * t86; t97; t87; t94; t88; Ifges(7,3) + t86; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
