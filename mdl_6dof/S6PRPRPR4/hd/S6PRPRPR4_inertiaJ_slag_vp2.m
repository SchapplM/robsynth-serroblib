% Calculate joint inertia matrix for
% S6PRPRPR4
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
% Datum: 2018-11-23 14:57
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6PRPRPR4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR4_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR4_inertiaJ_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR4_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR4_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR4_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 14:57:07
% EndTime: 2018-11-23 14:57:08
% DurationCPUTime: 0.68s
% Computational Cost: add. (1284->220), mult. (2732->329), div. (0->0), fcn. (3073->12), ass. (0->90)
t107 = cos(qJ(4));
t77 = sin(pkin(11));
t80 = cos(pkin(11));
t83 = sin(qJ(4));
t57 = t107 * t77 + t83 * t80;
t79 = cos(pkin(12));
t102 = t57 * t79;
t76 = sin(pkin(12));
t103 = t57 * t76;
t30 = mrSges(6,1) * t103 + mrSges(6,2) * t102;
t82 = sin(qJ(6));
t85 = cos(qJ(6));
t56 = t76 * t85 + t79 * t82;
t26 = t56 * t57;
t54 = -t76 * t82 + t79 * t85;
t27 = t54 * t57;
t9 = t26 * mrSges(7,1) + t27 * mrSges(7,2);
t114 = t30 + t9;
t55 = -t107 * t80 + t77 * t83;
t36 = t55 * mrSges(5,1) + t57 * mrSges(5,2);
t60 = -t80 * mrSges(4,1) + t77 * mrSges(4,2);
t113 = -t36 - t60;
t78 = sin(pkin(6));
t84 = sin(qJ(2));
t101 = t78 * t84;
t81 = cos(pkin(6));
t47 = -t77 * t101 + t80 * t81;
t48 = t80 * t101 + t77 * t81;
t23 = -t107 * t47 + t48 * t83;
t22 = t23 ^ 2;
t99 = pkin(8) + qJ(3);
t62 = t99 * t80;
t92 = t99 * t77;
t40 = t107 * t92 + t62 * t83;
t112 = t40 ^ 2;
t75 = t80 ^ 2;
t111 = 0.2e1 * t40;
t109 = t79 / 0.2e1;
t108 = -m(4) - m(5);
t106 = Ifges(6,4) * t76;
t105 = Ifges(6,4) * t79;
t104 = t23 * t40;
t86 = cos(qJ(2));
t100 = t78 * t86;
t98 = pkin(9) + qJ(5);
t68 = -pkin(3) * t80 - pkin(2);
t31 = pkin(4) * t55 - qJ(5) * t57 + t68;
t43 = t107 * t62 - t83 * t92;
t11 = t76 * t31 + t79 * t43;
t97 = Ifges(7,5) * t56 + Ifges(7,6) * t54;
t59 = -t79 * mrSges(6,1) + t76 * mrSges(6,2);
t96 = t59 - mrSges(5,1);
t95 = t76 ^ 2 + t79 ^ 2;
t94 = t77 ^ 2 + t75;
t93 = Ifges(7,5) * t27 - Ifges(7,6) * t26 + Ifges(7,3) * t55;
t10 = t79 * t31 - t43 * t76;
t35 = -t54 * mrSges(7,1) + t56 * mrSges(7,2);
t91 = -t10 * t76 + t11 * t79;
t25 = t107 * t48 + t83 * t47;
t14 = -t79 * t100 - t25 * t76;
t15 = -t76 * t100 + t25 * t79;
t90 = -t14 * t76 + t15 * t79;
t89 = -t47 * t77 + t48 * t80;
t73 = t78 ^ 2;
t67 = -pkin(5) * t79 - pkin(4);
t65 = t73 * t86 ^ 2;
t64 = Ifges(6,1) * t76 + t105;
t63 = Ifges(6,2) * t79 + t106;
t61 = t98 * t79;
t58 = t98 * t76;
t42 = -t58 * t82 + t61 * t85;
t39 = -t58 * t85 - t61 * t82;
t38 = Ifges(7,1) * t56 + Ifges(7,4) * t54;
t37 = Ifges(7,4) * t56 + Ifges(7,2) * t54;
t33 = mrSges(6,1) * t55 - mrSges(6,3) * t102;
t32 = -mrSges(6,2) * t55 - mrSges(6,3) * t103;
t18 = pkin(5) * t103 + t40;
t17 = Ifges(6,5) * t55 + (Ifges(6,1) * t79 - t106) * t57;
t16 = Ifges(6,6) * t55 + (-Ifges(6,2) * t76 + t105) * t57;
t13 = mrSges(7,1) * t55 - mrSges(7,3) * t27;
t12 = -mrSges(7,2) * t55 - mrSges(7,3) * t26;
t8 = -pkin(9) * t103 + t11;
t7 = Ifges(7,1) * t27 - Ifges(7,4) * t26 + Ifges(7,5) * t55;
t6 = Ifges(7,4) * t27 - Ifges(7,2) * t26 + Ifges(7,6) * t55;
t5 = pkin(5) * t55 - pkin(9) * t102 + t10;
t4 = t14 * t82 + t15 * t85;
t3 = t14 * t85 - t15 * t82;
t2 = t5 * t82 + t8 * t85;
t1 = t5 * t85 - t8 * t82;
t19 = [m(2) + m(6) * (t14 ^ 2 + t15 ^ 2 + t22) + m(7) * (t3 ^ 2 + t4 ^ 2 + t22) + m(5) * (t25 ^ 2 + t22 + t65) + m(4) * (t47 ^ 2 + t48 ^ 2 + t65) + m(3) * (t73 * t84 ^ 2 + t81 ^ 2 + t65); -t25 * t55 * mrSges(5,3) + t4 * t12 + t3 * t13 + t14 * t33 + t15 * t32 + t89 * mrSges(4,3) + (t57 * mrSges(5,3) + t114) * t23 + (-t84 * mrSges(3,2) + (mrSges(3,1) + t113) * t86) * t78 + m(6) * (t10 * t14 + t11 * t15 + t104) + m(7) * (t1 * t3 + t18 * t23 + t2 * t4) + m(5) * (-t68 * t100 + t25 * t43 + t104) + m(4) * (pkin(2) * t100 + t89 * qJ(3)); Ifges(4,2) * t75 - 0.2e1 * pkin(2) * t60 + 0.2e1 * t1 * t13 + 0.2e1 * t10 * t33 + 0.2e1 * t11 * t32 + 0.2e1 * t2 * t12 + 0.2e1 * t18 * t9 - t26 * t6 + t27 * t7 + t30 * t111 + 0.2e1 * t68 * t36 + Ifges(3,3) + (Ifges(4,1) * t77 + 0.2e1 * Ifges(4,4) * t80) * t77 + 0.2e1 * t94 * qJ(3) * mrSges(4,3) + (mrSges(5,3) * t111 + Ifges(5,1) * t57 - t76 * t16 + t79 * t17) * t57 + m(4) * (t94 * qJ(3) ^ 2 + pkin(2) ^ 2) + m(5) * (t43 ^ 2 + t68 ^ 2 + t112) + m(6) * (t10 ^ 2 + t11 ^ 2 + t112) + m(7) * (t1 ^ 2 + t18 ^ 2 + t2 ^ 2) + (-0.2e1 * t43 * mrSges(5,3) + (Ifges(6,3) + Ifges(5,2)) * t55 + (Ifges(6,5) * t79 - Ifges(6,6) * t76 - (2 * Ifges(5,4))) * t57 + t93) * t55; t108 * t100 + m(6) * (t14 * t79 + t15 * t76) + m(7) * (t3 * t54 + t4 * t56); -m(4) * pkin(2) + t56 * t12 + t54 * t13 + t76 * t32 + t79 * t33 + m(7) * (t1 * t54 + t2 * t56) + m(6) * (t10 * t79 + t11 * t76) + m(5) * t68 - t113; m(6) * t95 + m(7) * (t54 ^ 2 + t56 ^ 2) - t108; -t25 * mrSges(5,2) + (-t3 * t56 + t4 * t54) * mrSges(7,3) + t90 * mrSges(6,3) + (t35 + t96) * t23 + m(6) * (-pkin(4) * t23 + t90 * qJ(5)) + m(7) * (t23 * t67 + t3 * t39 + t4 * t42); t76 * t17 / 0.2e1 + t16 * t109 + t67 * t9 - Ifges(5,6) * t55 + t56 * t7 / 0.2e1 + t54 * t6 / 0.2e1 + t18 * t35 - t26 * t37 / 0.2e1 + t27 * t38 / 0.2e1 + t39 * t13 + t42 * t12 - t43 * mrSges(5,2) - pkin(4) * t30 + t96 * t40 + (t79 * t32 - t76 * t33) * qJ(5) + (-t1 * t56 + t2 * t54) * mrSges(7,3) + t91 * mrSges(6,3) + m(6) * (-pkin(4) * t40 + t91 * qJ(5)) + m(7) * (t1 * t39 + t18 * t67 + t2 * t42) + (Ifges(5,5) + t64 * t109 - t76 * t63 / 0.2e1) * t57 + (Ifges(6,5) * t76 + Ifges(6,6) * t79 + t97) * t55 / 0.2e1; m(7) * (t39 * t54 + t42 * t56); -0.2e1 * pkin(4) * t59 + 0.2e1 * t67 * t35 + t54 * t37 + t56 * t38 + t79 * t63 + t76 * t64 + Ifges(5,3) + m(7) * (t39 ^ 2 + t42 ^ 2 + t67 ^ 2) + m(6) * (t95 * qJ(5) ^ 2 + pkin(4) ^ 2) + 0.2e1 * (-t39 * t56 + t42 * t54) * mrSges(7,3) + 0.2e1 * t95 * qJ(5) * mrSges(6,3); 0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * t23; m(6) * t40 + m(7) * t18 + t114; 0; -m(6) * pkin(4) + m(7) * t67 + t35 + t59; m(6) + m(7); mrSges(7,1) * t3 - mrSges(7,2) * t4; mrSges(7,1) * t1 - mrSges(7,2) * t2 + t93; -t35; mrSges(7,1) * t39 - mrSges(7,2) * t42 + t97; 0; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t19(1) t19(2) t19(4) t19(7) t19(11) t19(16); t19(2) t19(3) t19(5) t19(8) t19(12) t19(17); t19(4) t19(5) t19(6) t19(9) t19(13) t19(18); t19(7) t19(8) t19(9) t19(10) t19(14) t19(19); t19(11) t19(12) t19(13) t19(14) t19(15) t19(20); t19(16) t19(17) t19(18) t19(19) t19(20) t19(21);];
Mq  = res;
