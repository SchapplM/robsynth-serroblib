% Calculate joint inertia matrix for
% S6PRPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta5]';
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
% Datum: 2019-03-08 19:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRPR6_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR6_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR6_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR6_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR6_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR6_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:46:48
% EndTime: 2019-03-08 19:46:49
% DurationCPUTime: 0.66s
% Computational Cost: add. (759->223), mult. (1645->320), div. (0->0), fcn. (1581->10), ass. (0->95)
t79 = (-pkin(2) - pkin(8));
t110 = -2 * t79;
t70 = sin(pkin(6));
t78 = cos(qJ(2));
t100 = t70 * t78;
t72 = cos(pkin(6));
t74 = sin(qJ(4));
t77 = cos(qJ(4));
t32 = t77 * t100 + t72 * t74;
t31 = t32 ^ 2;
t109 = m(6) * pkin(4);
t69 = sin(pkin(11));
t108 = t69 / 0.2e1;
t71 = cos(pkin(11));
t107 = t71 / 0.2e1;
t106 = m(6) + m(7);
t102 = t69 * t77;
t99 = t71 * t77;
t35 = mrSges(6,1) * t102 + mrSges(6,2) * t99;
t73 = sin(qJ(6));
t76 = cos(qJ(6));
t43 = t76 * t69 + t73 * t71;
t27 = t43 * t77;
t42 = -t73 * t69 + t76 * t71;
t29 = t42 * t77;
t7 = t27 * mrSges(7,1) + t29 * mrSges(7,2);
t105 = -t35 - t7;
t104 = Ifges(6,4) * t69;
t103 = Ifges(6,4) * t71;
t75 = sin(qJ(2));
t101 = t70 * t75;
t34 = -t74 * t100 + t72 * t77;
t98 = t74 * t34;
t97 = t74 * t79;
t96 = t77 * t32;
t95 = t77 * t79;
t94 = pkin(9) + qJ(5);
t48 = -t71 * mrSges(6,1) + t69 * mrSges(6,2);
t93 = t48 - mrSges(5,1);
t92 = t74 * mrSges(5,1) + t77 * mrSges(5,2) + mrSges(4,3);
t46 = t74 * pkin(4) - t77 * qJ(5) + qJ(3);
t20 = t69 * t46 + t71 * t97;
t91 = t69 ^ 2 + t71 ^ 2;
t66 = t74 ^ 2;
t67 = t77 ^ 2;
t90 = t66 + t67;
t10 = -t42 * mrSges(7,1) + t43 * mrSges(7,2);
t89 = -t10 - t93;
t88 = Ifges(7,5) * t29 - Ifges(7,6) * t27 + Ifges(7,3) * t74;
t87 = t90 * mrSges(5,3);
t86 = qJ(5) * t91;
t13 = t71 * t101 - t34 * t69;
t14 = t69 * t101 + t34 * t71;
t85 = -t13 * t69 + t14 * t71;
t40 = t71 * t46;
t19 = -t69 * t97 + t40;
t84 = -t19 * t69 + t20 * t71;
t83 = -t96 + t98;
t44 = -t74 * mrSges(6,2) - mrSges(6,3) * t102;
t45 = t74 * mrSges(6,1) - mrSges(6,3) * t99;
t82 = t71 * t44 - t69 * t45;
t81 = qJ(3) ^ 2;
t68 = t79 ^ 2;
t64 = t70 ^ 2;
t60 = t67 * t79;
t59 = t67 * t68;
t58 = -t71 * pkin(5) - pkin(4);
t57 = t64 * t75 ^ 2;
t54 = qJ(3) * t101;
t51 = Ifges(6,1) * t69 + t103;
t50 = Ifges(6,2) * t71 + t104;
t49 = t94 * t71;
t47 = t94 * t69;
t41 = (pkin(5) * t69 - t79) * t77;
t38 = Ifges(7,5) * t43;
t37 = Ifges(7,6) * t42;
t28 = t42 * t74;
t26 = t43 * t74;
t25 = Ifges(6,5) * t74 + (Ifges(6,1) * t71 - t104) * t77;
t24 = Ifges(6,6) * t74 + (-Ifges(6,2) * t69 + t103) * t77;
t18 = t74 * mrSges(7,1) - t29 * mrSges(7,3);
t17 = -t74 * mrSges(7,2) - t27 * mrSges(7,3);
t16 = -t73 * t47 + t76 * t49;
t15 = -t76 * t47 - t73 * t49;
t12 = Ifges(7,1) * t43 + Ifges(7,4) * t42;
t11 = Ifges(7,4) * t43 + Ifges(7,2) * t42;
t9 = -pkin(9) * t102 + t20;
t8 = -pkin(9) * t99 + t40 + (-t69 * t79 + pkin(5)) * t74;
t6 = Ifges(7,1) * t29 - Ifges(7,4) * t27 + Ifges(7,5) * t74;
t5 = Ifges(7,4) * t29 - Ifges(7,2) * t27 + Ifges(7,6) * t74;
t4 = t73 * t13 + t76 * t14;
t3 = t76 * t13 - t73 * t14;
t2 = t73 * t8 + t76 * t9;
t1 = -t73 * t9 + t76 * t8;
t21 = [m(2) + m(6) * (t13 ^ 2 + t14 ^ 2 + t31) + m(7) * (t3 ^ 2 + t4 ^ 2 + t31) + m(5) * (t34 ^ 2 + t31 + t57) + 0.2e1 * (m(4) / 0.2e1 + m(3) / 0.2e1) * (t64 * t78 ^ 2 + t72 ^ 2 + t57); -mrSges(5,3) * t98 + t13 * t45 + t14 * t44 + t4 * t17 + t3 * t18 + (t77 * mrSges(5,3) - t105) * t32 + ((mrSges(3,1) - mrSges(4,2)) * t78 + (-mrSges(3,2) + t92) * t75) * t70 + m(6) * (t19 * t13 + t20 * t14 - t32 * t95) + m(7) * (t1 * t3 + t2 * t4 + t41 * t32) + m(5) * (t83 * t79 + t54) + m(4) * (pkin(2) * t100 + t54); -0.2e1 * pkin(2) * mrSges(4,2) + 0.2e1 * t1 * t18 + 0.2e1 * t2 * t17 + 0.2e1 * t19 * t45 + 0.2e1 * t20 * t44 - t27 * t5 + t29 * t6 + 0.2e1 * t41 * t7 + Ifges(4,1) + Ifges(3,3) + ((Ifges(6,3) + Ifges(5,2)) * t74 + t88) * t74 + m(5) * (t66 * t68 + t59 + t81) + m(4) * (pkin(2) ^ 2 + t81) + m(6) * (t19 ^ 2 + t20 ^ 2 + t59) + m(7) * (t1 ^ 2 + t2 ^ 2 + t41 ^ 2) + (Ifges(5,1) * t77 - t69 * t24 + t71 * t25 + t35 * t110 + (Ifges(6,5) * t71 - Ifges(6,6) * t69 - (2 * Ifges(5,4))) * t74) * t77 + 0.2e1 * t92 * qJ(3) + t87 * t110; -m(4) * t100 + m(6) * (t85 * t74 - t96) + m(7) * (-t26 * t3 + t28 * t4 - t96) + m(5) * t83; -m(4) * pkin(2) + t28 * t17 - t26 * t18 + mrSges(4,2) + t105 * t77 + t82 * t74 - t87 + m(7) * (-t26 * t1 + t28 * t2 - t77 * t41) + m(6) * (t84 * t74 + t60) + m(5) * (t66 * t79 + t60); m(4) + m(5) * t90 + m(6) * (t91 * t66 + t67) + m(7) * (t26 ^ 2 + t28 ^ 2 + t67); -t34 * mrSges(5,2) + (-t3 * t43 + t4 * t42) * mrSges(7,3) + t85 * mrSges(6,3) - t89 * t32 + m(6) * (-pkin(4) * t32 + t85 * qJ(5)) + m(7) * (t15 * t3 + t16 * t4 + t58 * t32); t25 * t108 + t24 * t107 + t58 * t7 + t41 * t10 + t42 * t5 / 0.2e1 + t43 * t6 / 0.2e1 - t27 * t11 / 0.2e1 + t29 * t12 / 0.2e1 - pkin(4) * t35 + t16 * t17 + t15 * t18 + m(7) * (t15 * t1 + t16 * t2 + t58 * t41) + (Ifges(6,5) * t108 + Ifges(6,6) * t107 + t38 / 0.2e1 + t37 / 0.2e1 - Ifges(5,6) - t79 * mrSges(5,2)) * t74 + (-t1 * t43 + t2 * t42) * mrSges(7,3) + t84 * mrSges(6,3) + (m(6) * t84 + t82) * qJ(5) + (Ifges(5,5) + t51 * t107 - t69 * t50 / 0.2e1 + (-t93 + t109) * t79) * t77; (t26 * t43 + t28 * t42) * mrSges(7,3) + t89 * t77 + (t91 * mrSges(6,3) - mrSges(5,2)) * t74 + m(6) * (pkin(4) * t77 + t74 * t86) + m(7) * (-t15 * t26 + t16 * t28 - t58 * t77); -0.2e1 * pkin(4) * t48 + 0.2e1 * t58 * t10 + t42 * t11 + t43 * t12 + t71 * t50 + t69 * t51 + Ifges(5,3) + m(7) * (t15 ^ 2 + t16 ^ 2 + t58 ^ 2) + m(6) * (t91 * qJ(5) ^ 2 + pkin(4) ^ 2) + 0.2e1 * (-t15 * t43 + t16 * t42) * mrSges(7,3) + 0.2e1 * mrSges(6,3) * t86; 0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * t32; -m(6) * t95 + m(7) * t41 - t105; -t106 * t77; m(7) * t58 + t10 - t109 + t48; t106; t3 * mrSges(7,1) - t4 * mrSges(7,2); t1 * mrSges(7,1) - t2 * mrSges(7,2) + t88; -t26 * mrSges(7,1) - t28 * mrSges(7,2); t15 * mrSges(7,1) - t16 * mrSges(7,2) + t37 + t38; 0; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t21(1) t21(2) t21(4) t21(7) t21(11) t21(16); t21(2) t21(3) t21(5) t21(8) t21(12) t21(17); t21(4) t21(5) t21(6) t21(9) t21(13) t21(18); t21(7) t21(8) t21(9) t21(10) t21(14) t21(19); t21(11) t21(12) t21(13) t21(14) t21(15) t21(20); t21(16) t21(17) t21(18) t21(19) t21(20) t21(21);];
Mq  = res;
