% Calculate joint inertia matrix for
% S6RRPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta5]';
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
% Datum: 2018-11-23 16:43
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRPPPR2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR2_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR2_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR2_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR2_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR2_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:42:48
% EndTime: 2018-11-23 16:42:49
% DurationCPUTime: 0.70s
% Computational Cost: add. (1364->218), mult. (2488->301), div. (0->0), fcn. (2665->8), ass. (0->85)
t77 = sin(pkin(9));
t107 = pkin(2) * t77;
t63 = qJ(4) + t107;
t111 = t63 ^ 2;
t79 = cos(pkin(9));
t81 = sin(qJ(2));
t83 = cos(qJ(2));
t48 = t77 * t81 - t79 * t83;
t50 = t77 * t83 + t79 * t81;
t69 = -pkin(2) * t83 - pkin(1);
t87 = -qJ(4) * t50 + t69;
t29 = pkin(3) * t48 + t87;
t110 = -0.2e1 * t29;
t109 = 0.2e1 * t69;
t78 = cos(pkin(10));
t108 = t78 / 0.2e1;
t106 = pkin(2) * t79;
t68 = -pkin(3) - t106;
t62 = -qJ(5) + t68;
t105 = -pkin(8) + t62;
t76 = sin(pkin(10));
t104 = Ifges(6,4) * t76;
t103 = Ifges(6,4) * t78;
t102 = t48 * t76;
t101 = t48 * t78;
t99 = mrSges(5,2) - mrSges(4,1);
t98 = -qJ(3) - pkin(7);
t16 = (pkin(3) + qJ(5)) * t48 + t87;
t59 = t98 * t81;
t60 = t98 * t83;
t35 = -t79 * t59 - t60 * t77;
t19 = pkin(4) * t50 + t35;
t6 = t78 * t16 + t76 * t19;
t56 = t76 * mrSges(6,1) + t78 * mrSges(6,2);
t97 = t76 ^ 2 + t78 ^ 2;
t96 = t81 ^ 2 + t83 ^ 2;
t80 = sin(qJ(6));
t82 = cos(qJ(6));
t52 = -t76 * t80 + t78 * t82;
t88 = t82 * t76 + t80 * t78;
t95 = t52 ^ 2 + t88 ^ 2;
t24 = t52 * t48;
t25 = t88 * t48;
t94 = Ifges(7,5) * t25 + Ifges(7,6) * t24 + Ifges(7,3) * t50;
t37 = t59 * t77 - t60 * t79;
t93 = t35 ^ 2 + t37 ^ 2;
t54 = m(6) * t97;
t92 = t97 * mrSges(6,3);
t9 = -t24 * mrSges(7,1) + t25 * mrSges(7,2);
t28 = -mrSges(6,1) * t101 + mrSges(6,2) * t102;
t18 = t78 * t19;
t3 = pkin(5) * t50 + t18 + (-pkin(8) * t48 - t16) * t76;
t4 = pkin(8) * t101 + t6;
t1 = t3 * t82 - t4 * t80;
t2 = t3 * t80 + t4 * t82;
t91 = t1 * t52 + t2 * t88;
t5 = -t16 * t76 + t18;
t90 = t5 * t78 + t6 * t76;
t32 = mrSges(7,1) * t88 + t52 * mrSges(7,2);
t40 = t105 * t76;
t41 = t105 * t78;
t26 = -t40 * t80 + t41 * t82;
t27 = t40 * t82 + t41 * t80;
t89 = t26 * t52 + t27 * t88;
t86 = m(7) * t95 + m(5) + t54;
t58 = Ifges(6,1) * t78 - t104;
t57 = -Ifges(6,2) * t76 + t103;
t55 = pkin(5) * t76 + t63;
t46 = Ifges(7,5) * t52;
t45 = Ifges(7,6) * t88;
t43 = t50 * mrSges(5,3);
t42 = t50 * mrSges(4,2);
t34 = Ifges(7,1) * t52 - Ifges(7,4) * t88;
t33 = Ifges(7,4) * t52 - Ifges(7,2) * t88;
t31 = -mrSges(6,2) * t50 + mrSges(6,3) * t101;
t30 = mrSges(6,1) * t50 - mrSges(6,3) * t102;
t20 = -pkin(4) * t48 + t37;
t15 = Ifges(6,5) * t50 + (Ifges(6,1) * t76 + t103) * t48;
t14 = Ifges(6,6) * t50 + (Ifges(6,2) * t78 + t104) * t48;
t12 = mrSges(7,1) * t50 - mrSges(7,3) * t25;
t11 = -mrSges(7,2) * t50 + mrSges(7,3) * t24;
t10 = (-pkin(5) * t78 - pkin(4)) * t48 + t37;
t8 = Ifges(7,1) * t25 + Ifges(7,4) * t24 + Ifges(7,5) * t50;
t7 = Ifges(7,4) * t25 + Ifges(7,2) * t24 + Ifges(7,6) * t50;
t13 = [-0.2e1 * pkin(1) * (-t83 * mrSges(3,1) + t81 * mrSges(3,2)) + t83 * (Ifges(3,4) * t81 + Ifges(3,2) * t83) + t81 * (Ifges(3,1) * t81 + Ifges(3,4) * t83) + t42 * t109 + t43 * t110 + 0.2e1 * t20 * t28 + 0.2e1 * t5 * t30 + 0.2e1 * t6 * t31 + t24 * t7 + t25 * t8 + 0.2e1 * t10 * t9 + 0.2e1 * t2 * t11 + 0.2e1 * t1 * t12 + Ifges(2,3) + (mrSges(4,1) * t109 + mrSges(5,2) * t110 + t78 * t14 + t76 * t15 + (Ifges(5,3) + Ifges(4,2)) * t48) * t48 + m(3) * (t96 * pkin(7) ^ 2 + pkin(1) ^ 2) + m(4) * (t69 ^ 2 + t93) + m(6) * (t20 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(5) * (t29 ^ 2 + t93) + m(7) * (t1 ^ 2 + t10 ^ 2 + t2 ^ 2) + ((Ifges(6,3) + Ifges(5,2) + Ifges(4,1)) * t50 + (Ifges(6,5) * t76 + Ifges(6,6) * t78 - (2 * Ifges(4,4)) - (2 * Ifges(5,6))) * t48 + t94) * t50 + 0.2e1 * t96 * pkin(7) * mrSges(3,3) + 0.2e1 * (t35 * t50 - t37 * t48) * (mrSges(5,1) + mrSges(4,3)); Ifges(3,6) * t83 + Ifges(3,5) * t81 + t63 * t28 + t52 * t8 / 0.2e1 + t55 * t9 + t20 * t56 - t88 * t7 / 0.2e1 + t10 * t32 + t24 * t33 / 0.2e1 + t25 * t34 / 0.2e1 + t26 * t12 + t27 * t11 + (-mrSges(4,2) + mrSges(5,3)) * t37 + t99 * t35 + (-t81 * mrSges(3,1) - t83 * mrSges(3,2)) * pkin(7) - t91 * mrSges(7,3) + (t15 / 0.2e1 + t62 * t30 - t5 * mrSges(6,3)) * t78 + (-t14 / 0.2e1 + t62 * t31 - t6 * mrSges(6,3)) * t76 + (Ifges(6,5) * t108 - Ifges(6,6) * t76 / 0.2e1 + t46 / 0.2e1 - t45 / 0.2e1 + Ifges(4,5) - Ifges(5,4) + t68 * mrSges(5,1) - mrSges(4,3) * t106) * t50 + m(6) * (t20 * t63 + t90 * t62) + m(5) * (t35 * t68 + t37 * t63) + m(7) * (t1 * t26 + t10 * t55 + t2 * t27) + m(4) * (-t35 * t79 + t37 * t77) * pkin(2) + (-Ifges(4,6) + Ifges(5,5) + t57 * t108 - t63 * mrSges(5,1) + t76 * t58 / 0.2e1 - mrSges(4,3) * t107) * t48; 0.2e1 * t68 * mrSges(5,2) + 0.2e1 * t55 * t32 - t88 * t33 + t52 * t34 - t76 * t57 + t78 * t58 + Ifges(5,1) + Ifges(3,3) + Ifges(4,3) + m(7) * (t26 ^ 2 + t27 ^ 2 + t55 ^ 2) + m(6) * (t97 * t62 ^ 2 + t111) + m(5) * (t68 ^ 2 + t111) + m(4) * (t77 ^ 2 + t79 ^ 2) * pkin(2) ^ 2 + 0.2e1 * (mrSges(5,3) + t56) * t63 + 0.2e1 * (mrSges(4,1) * t79 - mrSges(4,2) * t77) * pkin(2) - 0.2e1 * t89 * mrSges(7,3) - 0.2e1 * t62 * t92; t52 * t11 - t88 * t12 - t76 * t30 + t78 * t31 + t42 - t43 - t99 * t48 + m(7) * (-t1 * t88 + t2 * t52) + m(6) * (-t5 * t76 + t6 * t78) + m(5) * t29 + m(4) * t69; m(7) * (-t26 * t88 + t27 * t52); m(4) + t86; m(5) * t35 + m(6) * t90 + m(7) * t91 + t50 * mrSges(5,1) + t11 * t88 + t52 * t12 + t78 * t30 + t76 * t31; m(5) * t68 + m(7) * t89 - t95 * mrSges(7,3) + t62 * t54 + mrSges(5,2) - t92; 0; t86; m(6) * t20 + m(7) * t10 + t28 + t9; m(6) * t63 + m(7) * t55 + t32 + t56; 0; 0; m(6) + m(7); mrSges(7,1) * t1 - mrSges(7,2) * t2 + t94; mrSges(7,1) * t26 - t27 * mrSges(7,2) - t45 + t46; -t32; mrSges(7,1) * t52 - mrSges(7,2) * t88; 0; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t13(1) t13(2) t13(4) t13(7) t13(11) t13(16); t13(2) t13(3) t13(5) t13(8) t13(12) t13(17); t13(4) t13(5) t13(6) t13(9) t13(13) t13(18); t13(7) t13(8) t13(9) t13(10) t13(14) t13(19); t13(11) t13(12) t13(13) t13(14) t13(15) t13(20); t13(16) t13(17) t13(18) t13(19) t13(20) t13(21);];
Mq  = res;
