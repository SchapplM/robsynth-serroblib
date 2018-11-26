% Calculate joint inertia matrix for
% S6PPRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
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
% Datum: 2018-11-23 14:50
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6PPRRPR2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_inertiaJ_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRPR2_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRPR2_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRPR2_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 14:50:18
% EndTime: 2018-11-23 14:50:19
% DurationCPUTime: 0.64s
% Computational Cost: add. (660->187), mult. (1701->265), div. (0->0), fcn. (1772->12), ass. (0->84)
t108 = m(6) + m(5);
t60 = sin(qJ(4));
t63 = cos(qJ(4));
t107 = t60 ^ 2 + t63 ^ 2;
t90 = mrSges(6,1) + mrSges(5,3);
t103 = -m(6) * pkin(4) + mrSges(6,2);
t53 = sin(pkin(12));
t55 = sin(pkin(6));
t58 = cos(pkin(6));
t61 = sin(qJ(3));
t64 = cos(qJ(3));
t56 = cos(pkin(12));
t57 = cos(pkin(7));
t95 = t56 * t57;
t54 = sin(pkin(7));
t97 = t54 * t61;
t12 = t58 * t97 + (t53 * t64 + t61 * t95) * t55;
t20 = -t54 * t55 * t56 + t58 * t57;
t7 = t63 * t12 + t60 * t20;
t4 = t7 ^ 2;
t96 = t54 * t64;
t10 = -t58 * t96 + (t53 * t61 - t64 * t95) * t55;
t9 = t10 ^ 2;
t24 = t60 * t57 + t63 * t97;
t21 = t24 ^ 2;
t59 = sin(qJ(6));
t102 = -t59 / 0.2e1;
t65 = -pkin(4) - pkin(10);
t101 = pkin(5) + pkin(9);
t3 = t24 * t7;
t100 = mrSges(7,3) * t63;
t99 = Ifges(7,4) * t59;
t62 = cos(qJ(6));
t98 = Ifges(7,4) * t62;
t30 = -t60 * mrSges(7,2) - t62 * t100;
t94 = t59 * t30;
t93 = t62 * t65;
t92 = t7 * qJ(5);
t91 = -mrSges(5,1) + mrSges(6,2);
t34 = t59 * mrSges(7,1) + t62 * mrSges(7,2);
t89 = t34 + mrSges(6,3);
t88 = t107 * pkin(9) ^ 2;
t87 = t59 ^ 2 + t62 ^ 2;
t86 = qJ(5) * t24;
t32 = t63 * mrSges(6,2) - t60 * mrSges(6,3);
t33 = -t63 * mrSges(5,1) + t60 * mrSges(5,2);
t83 = mrSges(4,1) - t32 - t33;
t82 = -mrSges(5,2) + t89;
t81 = m(7) * t87;
t80 = t90 * t60;
t79 = t87 * mrSges(7,3);
t78 = -t60 * qJ(5) - pkin(3);
t5 = t60 * t12 - t63 * t20;
t1 = -t59 * t10 + t62 * t5;
t2 = t62 * t10 + t59 * t5;
t76 = t62 * t1 + t59 * t2;
t74 = t62 * mrSges(7,1) - t59 * mrSges(7,2);
t26 = t74 * t63;
t75 = t90 * t63 + t26;
t73 = -Ifges(7,5) * t59 - Ifges(7,6) * t62;
t27 = t65 * t63 + t78;
t38 = t101 * t60;
t13 = -t59 * t27 + t62 * t38;
t14 = t62 * t27 + t59 * t38;
t72 = t62 * t13 + t59 * t14;
t22 = -t63 * t57 + t60 * t97;
t15 = t62 * t22 + t59 * t96;
t16 = t59 * t22 - t62 * t96;
t71 = t62 * t15 + t59 * t16;
t69 = (t5 * t60 + t7 * t63) * pkin(9);
t68 = (t22 * t60 + t24 * t63) * pkin(9);
t66 = qJ(5) ^ 2;
t46 = t54 ^ 2;
t43 = Ifges(7,5) * t62;
t42 = Ifges(7,3) * t60;
t40 = t46 * t64 ^ 2;
t39 = t101 * t63;
t36 = Ifges(7,1) * t62 - t99;
t35 = -Ifges(7,2) * t59 + t98;
t31 = -t63 * pkin(4) + t78;
t29 = t60 * mrSges(7,1) + t59 * t100;
t19 = Ifges(7,5) * t60 + (-Ifges(7,1) * t59 - t98) * t63;
t18 = Ifges(7,6) * t60 + (-Ifges(7,2) * t62 - t99) * t63;
t6 = [m(2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t4) + m(4) * (t12 ^ 2 + t20 ^ 2 + t9) + m(3) * (t58 ^ 2 + (t53 ^ 2 + t56 ^ 2) * t55 ^ 2) + t108 * (t5 ^ 2 + t4 + t9); m(3) * t58 + m(7) * (t15 * t1 + t16 * t2 + t3) + m(4) * (t57 * t20 + (-t10 * t64 + t12 * t61) * t54) + t108 * (-t10 * t96 + t22 * t5 + t3); m(3) + m(7) * (t15 ^ 2 + t16 ^ 2 + t21) + m(4) * (t46 * t61 ^ 2 + t57 ^ 2 + t40) + t108 * (t22 ^ 2 + t21 + t40); -t12 * mrSges(4,2) + t1 * t29 + t2 * t30 + t5 * t80 + t75 * t7 - t83 * t10 + m(7) * (t13 * t1 + t14 * t2 + t39 * t7) + m(6) * (t31 * t10 + t69) + m(5) * (-pkin(3) * t10 + t69); t15 * t29 + t16 * t30 + t22 * t80 + t75 * t24 + (-t61 * mrSges(4,2) + t83 * t64) * t54 + m(7) * (t13 * t15 + t14 * t16 + t39 * t24) + m(6) * (-t31 * t96 + t68) + m(5) * (pkin(3) * t96 + t68); -0.2e1 * pkin(3) * t33 + 0.2e1 * t13 * t29 + 0.2e1 * t14 * t30 + 0.2e1 * t39 * t26 + 0.2e1 * t31 * t32 + Ifges(4,3) + m(7) * (t13 ^ 2 + t14 ^ 2 + t39 ^ 2) + m(6) * (t31 ^ 2 + t88) + m(5) * (pkin(3) ^ 2 + t88) + (t42 + (Ifges(5,1) + Ifges(6,2)) * t60) * t60 + (-t59 * t19 - t62 * t18 + (Ifges(6,3) + Ifges(5,2)) * t63 + ((2 * Ifges(5,4)) + (2 * Ifges(6,6)) + t73) * t60) * t63 + 0.2e1 * t90 * pkin(9) * t107; t91 * t5 - t76 * mrSges(7,3) + t82 * t7 + m(7) * (t76 * t65 + t92) + m(6) * (-pkin(4) * t5 + t92); t91 * t22 - t71 * mrSges(7,3) + t82 * t24 + m(7) * (t71 * t65 + t86) + m(6) * (-pkin(4) * t22 + t86); t65 * t94 + t29 * t93 + t62 * t19 / 0.2e1 + t18 * t102 + t39 * t34 + qJ(5) * t26 + m(7) * (qJ(5) * t39 + t72 * t65) - t72 * mrSges(7,3) + (-pkin(4) * mrSges(6,1) - Ifges(6,4) + Ifges(5,5) + Ifges(7,6) * t102 + t43 / 0.2e1) * t60 + (t36 * t102 - t62 * t35 / 0.2e1 + qJ(5) * mrSges(6,1) - Ifges(6,5) + Ifges(5,6)) * t63 + ((m(6) * qJ(5) - mrSges(5,2) + mrSges(6,3)) * t63 + (-mrSges(5,1) + t103) * t60) * pkin(9); -0.2e1 * pkin(4) * mrSges(6,2) - t59 * t35 + t62 * t36 + Ifges(6,1) + Ifges(5,3) + m(6) * (pkin(4) ^ 2 + t66) + m(7) * (t87 * t65 ^ 2 + t66) + 0.2e1 * t89 * qJ(5) - 0.2e1 * t65 * t79; m(6) * t5 + m(7) * t76; m(6) * t22 + m(7) * t71; t94 + t62 * t29 + m(7) * t72 + (m(6) * pkin(9) + mrSges(6,1)) * t60; t65 * t81 + t103 - t79; m(6) + t81; t1 * mrSges(7,1) - t2 * mrSges(7,2); t15 * mrSges(7,1) - t16 * mrSges(7,2); t13 * mrSges(7,1) - t14 * mrSges(7,2) + t73 * t63 + t42; mrSges(7,1) * t93 + t43 + (-mrSges(7,2) * t65 - Ifges(7,6)) * t59; t74; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t6(1) t6(2) t6(4) t6(7) t6(11) t6(16); t6(2) t6(3) t6(5) t6(8) t6(12) t6(17); t6(4) t6(5) t6(6) t6(9) t6(13) t6(18); t6(7) t6(8) t6(9) t6(10) t6(14) t6(19); t6(11) t6(12) t6(13) t6(14) t6(15) t6(20); t6(16) t6(17) t6(18) t6(19) t6(20) t6(21);];
Mq  = res;
