% Calculate joint inertia matrix for
% S6PPRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
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
% Datum: 2019-03-08 18:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PPRPRR1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_inertiaJ_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRPRR1_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRPRR1_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRPRR1_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:42:16
% EndTime: 2019-03-08 18:42:17
% DurationCPUTime: 0.55s
% Computational Cost: add. (844->187), mult. (2332->281), div. (0->0), fcn. (2645->14), ass. (0->89)
t52 = sin(pkin(13));
t41 = t52 * pkin(3) + pkin(9);
t102 = 0.2e1 * t41;
t101 = m(7) * pkin(10) + mrSges(7,3);
t60 = sin(qJ(6));
t63 = cos(qJ(6));
t32 = -t63 * mrSges(7,1) + t60 * mrSges(7,2);
t100 = -m(7) * pkin(5) - mrSges(6,1) + t32;
t54 = sin(pkin(7));
t55 = sin(pkin(6));
t57 = cos(pkin(12));
t58 = cos(pkin(7));
t59 = cos(pkin(6));
t27 = -t55 * t57 * t54 + t59 * t58;
t61 = sin(qJ(5));
t64 = cos(qJ(5));
t53 = sin(pkin(12));
t62 = sin(qJ(3));
t65 = cos(qJ(3));
t82 = t57 * t58;
t83 = t54 * t65;
t12 = t59 * t83 + (-t53 * t62 + t65 * t82) * t55;
t84 = t54 * t62;
t13 = t59 * t84 + (t53 * t65 + t62 * t82) * t55;
t56 = cos(pkin(13));
t8 = t52 * t12 + t56 * t13;
t3 = -t27 * t64 + t61 * t8;
t99 = t3 ^ 2;
t6 = -t56 * t12 + t52 * t13;
t98 = t6 ^ 2;
t23 = (t52 * t65 + t56 * t62) * t54;
t16 = t61 * t23 - t64 * t58;
t97 = t16 ^ 2;
t21 = t52 * t84 - t56 * t83;
t96 = t21 ^ 2;
t95 = m(5) * pkin(3);
t94 = t63 / 0.2e1;
t93 = t16 * t3;
t92 = t21 * t6;
t91 = t64 * pkin(5);
t90 = t64 * t3;
t89 = Ifges(7,4) * t60;
t88 = Ifges(7,4) * t63;
t87 = Ifges(7,6) * t64;
t86 = t41 * t61;
t85 = t41 * t64;
t81 = t60 * t61;
t80 = t61 * t63;
t79 = t64 * t16;
t33 = -t64 * mrSges(6,1) + t61 * mrSges(6,2);
t78 = -mrSges(5,1) + t33;
t77 = t60 ^ 2 + t63 ^ 2;
t49 = t61 ^ 2;
t51 = t64 ^ 2;
t76 = t49 + t51;
t42 = -t56 * pkin(3) - pkin(4);
t75 = t77 * t61;
t5 = t27 * t61 + t64 * t8;
t1 = -t60 * t5 + t63 * t6;
t2 = t63 * t5 + t60 * t6;
t74 = -t1 * t60 + t2 * t63;
t73 = t3 * t61 + t5 * t64;
t18 = t64 * t23 + t61 * t58;
t10 = t63 * t18 + t60 * t21;
t9 = -t60 * t18 + t63 * t21;
t72 = t10 * t63 - t60 * t9;
t71 = mrSges(7,1) * t60 + mrSges(7,2) * t63;
t29 = -t61 * pkin(10) + t42 - t91;
t14 = t63 * t29 - t60 * t85;
t15 = t60 * t29 + t63 * t85;
t70 = -t14 * t60 + t15 * t63;
t69 = t16 * t61 + t18 * t64;
t30 = t64 * mrSges(7,2) - mrSges(7,3) * t81;
t31 = -t64 * mrSges(7,1) - mrSges(7,3) * t80;
t68 = t63 * t30 - t60 * t31;
t47 = t58 ^ 2;
t44 = Ifges(7,5) * t60;
t43 = Ifges(7,6) * t63;
t39 = t41 ^ 2;
t38 = Ifges(7,5) * t80;
t36 = t49 * t39;
t35 = Ifges(7,1) * t60 + t88;
t34 = Ifges(7,2) * t63 + t89;
t28 = t71 * t61;
t26 = t27 ^ 2;
t25 = -Ifges(7,5) * t64 + (Ifges(7,1) * t63 - t89) * t61;
t24 = -t87 + (-Ifges(7,2) * t60 + t88) * t61;
t19 = t58 * t27;
t4 = [m(2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t99) + m(6) * (t5 ^ 2 + t98 + t99) + m(4) * (t12 ^ 2 + t13 ^ 2 + t26) + m(5) * (t8 ^ 2 + t26 + t98) + m(3) * (t59 ^ 2 + (t53 ^ 2 + t57 ^ 2) * t55 ^ 2); m(3) * t59 + m(7) * (t9 * t1 + t10 * t2 + t93) + m(6) * (t18 * t5 + t92 + t93) + m(4) * (t19 + (t12 * t65 + t13 * t62) * t54) + m(5) * (t23 * t8 + t19 + t92); m(3) + m(7) * (t10 ^ 2 + t9 ^ 2 + t97) + m(6) * (t18 ^ 2 + t96 + t97) + m(5) * (t23 ^ 2 + t47 + t96) + m(4) * (t47 + (t62 ^ 2 + t65 ^ 2) * t54 ^ 2); t12 * mrSges(4,1) - t13 * mrSges(4,2) - t8 * mrSges(5,2) + t1 * t31 + t2 * t30 + t3 * t28 + t78 * t6 + t73 * mrSges(6,3) + m(7) * (t14 * t1 + t15 * t2 + t3 * t86) + m(6) * (t73 * t41 + t42 * t6) + (t52 * t8 - t56 * t6) * t95; -t23 * mrSges(5,2) + t10 * t30 + t16 * t28 + t9 * t31 + (t65 * mrSges(4,1) - t62 * mrSges(4,2)) * t54 + t78 * t21 + t69 * mrSges(6,3) + m(7) * (t15 * t10 + t14 * t9 + t16 * t86) + m(6) * (t42 * t21 + t41 * t69) + (-t21 * t56 + t23 * t52) * t95; 0.2e1 * t14 * t31 + 0.2e1 * t15 * t30 + 0.2e1 * t42 * t33 + Ifges(4,3) + Ifges(5,3) + (-t38 + (Ifges(7,3) + Ifges(6,2)) * t64) * t64 + m(7) * (t14 ^ 2 + t15 ^ 2 + t36) + m(6) * (t51 * t39 + t42 ^ 2 + t36) + m(5) * (t52 ^ 2 + t56 ^ 2) * pkin(3) ^ 2 + (Ifges(6,1) * t61 + 0.2e1 * Ifges(6,4) * t64 + t63 * t25 + t28 * t102 + (-t24 + t87) * t60) * t61 + 0.2e1 * (t56 * mrSges(5,1) - t52 * mrSges(5,2)) * pkin(3) + t76 * mrSges(6,3) * t102; m(7) * (t74 * t61 - t90) + m(6) * (t61 * t5 - t90) + m(5) * t27; m(5) * t58 + m(7) * (t61 * t72 - t79) + m(6) * (t61 * t18 - t79); -t64 * t28 + (m(7) * (t70 - t85) + t68) * t61; m(5) + m(6) * t76 + m(7) * (t77 * t49 + t51); -t5 * mrSges(6,2) + t100 * t3 + t101 * t74; -t18 * mrSges(6,2) + t100 * t16 + t101 * t72; t60 * t25 / 0.2e1 + t24 * t94 - pkin(5) * t28 + (-t44 / 0.2e1 - t43 / 0.2e1 + Ifges(6,6) - t41 * mrSges(6,2)) * t64 + t70 * mrSges(7,3) + (m(7) * t70 + t68) * pkin(10) + (Ifges(6,5) - t60 * t34 / 0.2e1 + t35 * t94 + t100 * t41) * t61; -t64 * t32 + m(7) * (pkin(10) * t75 + t91) + mrSges(7,3) * t75 - t33; Ifges(6,3) + t63 * t34 + m(7) * (t77 * pkin(10) ^ 2 + pkin(5) ^ 2) + t60 * t35 - 0.2e1 * pkin(5) * t32 + 0.2e1 * t77 * pkin(10) * mrSges(7,3); t1 * mrSges(7,1) - t2 * mrSges(7,2); t9 * mrSges(7,1) - t10 * mrSges(7,2); t14 * mrSges(7,1) - t15 * mrSges(7,2) - Ifges(7,6) * t81 - Ifges(7,3) * t64 + t38; -t28; -pkin(10) * t71 + t43 + t44; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t4(1) t4(2) t4(4) t4(7) t4(11) t4(16); t4(2) t4(3) t4(5) t4(8) t4(12) t4(17); t4(4) t4(5) t4(6) t4(9) t4(13) t4(18); t4(7) t4(8) t4(9) t4(10) t4(14) t4(19); t4(11) t4(12) t4(13) t4(14) t4(15) t4(20); t4(16) t4(17) t4(18) t4(19) t4(20) t4(21);];
Mq  = res;
