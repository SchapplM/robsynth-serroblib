% Calculate joint inertia matrix for
% S6PRPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
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
% Datum: 2018-11-23 15:02
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRP5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP5_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP5_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP5_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:01:58
% EndTime: 2018-11-23 15:01:58
% DurationCPUTime: 0.69s
% Computational Cost: add. (550->215), mult. (1193->286), div. (0->0), fcn. (995->8), ass. (0->90)
t108 = Ifges(6,3) + Ifges(7,3);
t107 = -2 * mrSges(7,3);
t106 = m(6) + m(7);
t105 = m(6) * pkin(4);
t65 = cos(qJ(5));
t44 = -t65 * pkin(5) - pkin(4);
t104 = m(7) * t44;
t61 = cos(pkin(6));
t63 = sin(qJ(4));
t66 = cos(qJ(4));
t60 = sin(pkin(6));
t67 = cos(qJ(2));
t96 = t60 * t67;
t15 = t61 * t63 + t66 * t96;
t103 = t15 ^ 2;
t102 = m(7) * pkin(5);
t62 = sin(qJ(5));
t101 = Ifges(6,4) * t62;
t100 = Ifges(6,4) * t65;
t99 = Ifges(7,4) * t62;
t98 = Ifges(7,4) * t65;
t64 = sin(qJ(2));
t97 = t60 * t64;
t95 = t62 * t66;
t17 = t61 * t66 - t63 * t96;
t94 = t63 * t17;
t68 = -pkin(2) - pkin(8);
t93 = t63 * t68;
t92 = t65 * t66;
t91 = t66 * t15;
t90 = -mrSges(6,2) - mrSges(7,2);
t89 = -Ifges(6,6) - Ifges(7,6);
t88 = -qJ(6) - pkin(9);
t18 = mrSges(7,1) * t95 + mrSges(7,2) * t92;
t72 = mrSges(6,1) * t62 + mrSges(6,2) * t65;
t19 = t72 * t66;
t87 = -t18 - t19;
t23 = -t63 * mrSges(7,2) - mrSges(7,3) * t95;
t24 = -t63 * mrSges(6,2) - mrSges(6,3) * t95;
t86 = t23 + t24;
t25 = t63 * mrSges(7,1) - mrSges(7,3) * t92;
t26 = t63 * mrSges(6,1) - mrSges(6,3) * t92;
t85 = t25 + t26;
t30 = -t65 * mrSges(6,1) + t62 * mrSges(6,2);
t84 = t30 - mrSges(5,1);
t83 = t63 * mrSges(5,1) + t66 * mrSges(5,2) + mrSges(4,3);
t27 = t63 * pkin(4) - t66 * pkin(9) + qJ(3);
t8 = t62 * t27 + t65 * t93;
t82 = t62 ^ 2 + t65 ^ 2;
t56 = t63 ^ 2;
t58 = t66 ^ 2;
t81 = -t58 - t56;
t80 = qJ(6) * t66;
t29 = -t65 * mrSges(7,1) + t62 * mrSges(7,2);
t78 = -t29 - t84;
t77 = t82 * mrSges(6,3);
t76 = t81 * mrSges(5,3);
t75 = (Ifges(6,5) + Ifges(7,5)) * t92 + t108 * t63;
t74 = mrSges(6,1) + mrSges(7,1) + t102;
t5 = -t62 * t17 + t65 * t97;
t6 = t65 * t17 + t62 * t97;
t73 = -t5 * t62 + t6 * t65;
t71 = -t91 + t94;
t69 = qJ(3) ^ 2;
t59 = t68 ^ 2;
t54 = t60 ^ 2;
t53 = Ifges(6,5) * t62;
t52 = Ifges(7,5) * t62;
t51 = Ifges(6,6) * t65;
t50 = Ifges(7,6) * t65;
t46 = t58 * t68;
t45 = t58 * t59;
t43 = t54 * t64 ^ 2;
t38 = qJ(3) * t97;
t36 = Ifges(6,1) * t62 + t100;
t35 = Ifges(7,1) * t62 + t98;
t34 = Ifges(6,2) * t65 + t101;
t33 = Ifges(7,2) * t65 + t99;
t31 = t88 * t65;
t28 = t88 * t62;
t22 = (pkin(5) * t62 - t68) * t66;
t21 = t65 * t27;
t12 = Ifges(6,5) * t63 + (Ifges(6,1) * t65 - t101) * t66;
t11 = Ifges(7,5) * t63 + (Ifges(7,1) * t65 - t99) * t66;
t10 = Ifges(6,6) * t63 + (-Ifges(6,2) * t62 + t100) * t66;
t9 = Ifges(7,6) * t63 + (-Ifges(7,2) * t62 + t98) * t66;
t7 = -t62 * t93 + t21;
t4 = -t62 * t80 + t8;
t3 = -t65 * t80 + t21 + (-t62 * t68 + pkin(5)) * t63;
t1 = [m(2) + m(5) * (t17 ^ 2 + t103 + t43) + (t5 ^ 2 + t6 ^ 2 + t103) * t106 + 0.2e1 * (m(4) / 0.2e1 + m(3) / 0.2e1) * (t54 * t67 ^ 2 + t61 ^ 2 + t43); -mrSges(5,3) * t94 + t86 * t6 + t85 * t5 + (t66 * mrSges(5,3) - t87) * t15 + ((mrSges(3,1) - mrSges(4,2)) * t67 + (-mrSges(3,2) + t83) * t64) * t60 + m(6) * (t7 * t5 + t8 * t6 - t68 * t91) + m(7) * (t22 * t15 + t3 * t5 + t4 * t6) + m(5) * (t71 * t68 + t38) + m(4) * (pkin(2) * t96 + t38); -0.2e1 * pkin(2) * mrSges(4,2) + 0.2e1 * t22 * t18 + 0.2e1 * t4 * t23 + 0.2e1 * t8 * t24 + 0.2e1 * t3 * t25 + 0.2e1 * t7 * t26 + Ifges(4,1) + Ifges(3,3) + (Ifges(5,2) * t63 + t75) * t63 + m(5) * (t56 * t59 + t45 + t69) + m(4) * (pkin(2) ^ 2 + t69) + m(6) * (t7 ^ 2 + t8 ^ 2 + t45) + m(7) * (t22 ^ 2 + t3 ^ 2 + t4 ^ 2) + (Ifges(5,1) * t66 - 0.2e1 * Ifges(5,4) * t63 - 0.2e1 * t68 * t19 + (t11 + t12) * t65 + (t89 * t63 - t10 - t9) * t62) * t66 + 0.2e1 * t83 * qJ(3) + 0.2e1 * t68 * t76; -m(4) * t96 + m(5) * t71 + (t73 * t63 - t91) * t106; -m(4) * pkin(2) + mrSges(4,2) + t87 * t66 + t76 + (-t85 * t62 + t86 * t65) * t63 + m(7) * (-t66 * t22 + (-t3 * t62 + t4 * t65) * t63) + m(6) * (t46 + (-t62 * t7 + t65 * t8) * t63) + m(5) * (t56 * t68 + t46); m(4) - m(5) * t81 + (t82 * t56 + t58) * t106; -t17 * mrSges(5,2) + m(7) * (t28 * t5 - t31 * t6) + (t104 - t78 - t105) * t15 + (m(6) * pkin(9) + mrSges(6,3) + mrSges(7,3)) * t73; m(7) * (t44 * t22 + t28 * t3 - t31 * t4) + t44 * t18 + t28 * t25 + t22 * t29 - t31 * t23 - pkin(4) * t19 + (-t68 * mrSges(5,2) + t52 / 0.2e1 + t50 / 0.2e1 + t53 / 0.2e1 + t51 / 0.2e1 - Ifges(5,6)) * t63 + (Ifges(5,5) + (-t84 + t105) * t68) * t66 + (t8 * mrSges(6,3) + t4 * mrSges(7,3) + t9 / 0.2e1 + t10 / 0.2e1 + (t35 / 0.2e1 + t36 / 0.2e1) * t66 + (m(6) * t8 + t24) * pkin(9)) * t65 + (-t7 * mrSges(6,3) - t3 * mrSges(7,3) + t11 / 0.2e1 + t12 / 0.2e1 + (-t33 / 0.2e1 - t34 / 0.2e1) * t66 + (-m(6) * t7 - t26) * pkin(9)) * t62; t78 * t66 + (t82 * mrSges(7,3) - mrSges(5,2) + t77) * t63 + m(7) * (-t44 * t66 + (-t28 * t62 - t31 * t65) * t63) + m(6) * (t82 * t63 * pkin(9) + pkin(4) * t66); -0.2e1 * pkin(4) * t30 + 0.2e1 * t44 * t29 + Ifges(5,3) + 0.2e1 * pkin(9) * t77 + m(7) * (t28 ^ 2 + t31 ^ 2 + t44 ^ 2) + m(6) * (t82 * pkin(9) ^ 2 + pkin(4) ^ 2) + (t31 * t107 + t33 + t34) * t65 + (t28 * t107 + t35 + t36) * t62; t74 * t5 + t90 * t6; t7 * mrSges(6,1) + t3 * mrSges(7,1) - t8 * mrSges(6,2) - t4 * mrSges(7,2) + t89 * t95 + (m(7) * t3 + t25) * pkin(5) + t75; (-t74 * t62 + t90 * t65) * t63; t28 * mrSges(7,1) + t31 * mrSges(7,2) + t50 + t51 + t52 + t53 - t72 * pkin(9) + (m(7) * t28 - t62 * mrSges(7,3)) * pkin(5); (0.2e1 * mrSges(7,1) + t102) * pkin(5) + t108; m(7) * t15; m(7) * t22 + t18; -m(7) * t66; t29 + t104; 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
