% Calculate joint inertia matrix for
% S6RPRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6]';
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
% Datum: 2018-11-23 16:10
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRR12_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR12_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRR12_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR12_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR12_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR12_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:10:09
% EndTime: 2018-11-23 16:10:10
% DurationCPUTime: 0.79s
% Computational Cost: add. (921->229), mult. (1608->309), div. (0->0), fcn. (1393->6), ass. (0->92)
t73 = sin(qJ(5));
t76 = cos(qJ(5));
t104 = t73 ^ 2 + t76 ^ 2;
t100 = m(6) * t104;
t122 = -m(5) - t100;
t74 = sin(qJ(3));
t67 = t74 ^ 2;
t77 = cos(qJ(3));
t69 = t77 ^ 2;
t103 = t69 + t67;
t121 = (-mrSges(4,3) - mrSges(5,1)) * t103;
t108 = t74 * t76;
t109 = t73 * t74;
t119 = Ifges(6,5) * t109 + Ifges(6,6) * t108 + Ifges(6,3) * t77;
t61 = qJ(4) * t74;
t118 = m(5) * (pkin(3) * t77 + t61);
t117 = t104 * mrSges(6,3) - mrSges(5,2);
t44 = t74 * pkin(3) - qJ(4) * t77 + qJ(2);
t116 = -0.2e1 * t44;
t115 = 0.2e1 * qJ(2);
t114 = m(7) * pkin(5);
t78 = -pkin(3) - pkin(8);
t112 = -pkin(9) + t78;
t111 = Ifges(6,4) * t73;
t110 = Ifges(6,4) * t76;
t107 = t76 * mrSges(6,1);
t30 = pkin(8) * t74 + t44;
t79 = -pkin(1) - pkin(7);
t48 = (pkin(4) - t79) * t77;
t10 = t76 * t30 + t73 * t48;
t49 = mrSges(6,1) * t73 + mrSges(6,2) * t76;
t106 = t49 + mrSges(5,3);
t105 = t103 * t79 ^ 2;
t72 = sin(qJ(6));
t75 = cos(qJ(6));
t87 = t72 * t73 - t75 * t76;
t88 = t72 * t76 + t75 * t73;
t102 = t87 ^ 2 + t88 ^ 2;
t25 = t87 * t74;
t27 = t88 * t74;
t101 = Ifges(7,5) * t27 - Ifges(7,6) * t25 + Ifges(7,3) * t77;
t99 = t104 * t78;
t26 = t87 * t77;
t28 = t88 * t77;
t98 = t26 * mrSges(7,1) + t28 * mrSges(7,2);
t97 = -mrSges(7,1) * t87 - mrSges(7,2) * t88;
t36 = t76 * t48;
t6 = t77 * pkin(5) + t36 + (-pkin(9) * t74 - t30) * t73;
t8 = pkin(9) * t108 + t10;
t2 = t6 * t75 - t72 * t8;
t3 = t6 * t72 + t75 * t8;
t95 = -t2 * t87 + t3 * t88;
t9 = -t30 * t73 + t36;
t94 = t10 * t73 + t76 * t9;
t93 = t73 * mrSges(6,2) - t107;
t45 = t112 * t73;
t47 = t112 * t76;
t15 = -t45 * t72 + t47 * t75;
t16 = t45 * t75 + t47 * t72;
t92 = -t15 * t87 + t16 * t88;
t91 = -t26 * t87 - t28 * t88;
t90 = -t72 * t88 + t75 * t87;
t42 = mrSges(6,1) * t77 - mrSges(6,3) * t109;
t43 = -mrSges(6,2) * t77 + mrSges(6,3) * t108;
t89 = -t76 * t42 - t73 * t43;
t33 = Ifges(7,6) * t88;
t34 = Ifges(7,5) * t87;
t86 = t15 * mrSges(7,1) - t16 * mrSges(7,2) - t33 - t34;
t85 = t2 * mrSges(7,1) - t3 * mrSges(7,2) + t101;
t84 = (mrSges(7,1) * t75 - mrSges(7,2) * t72) * pkin(5);
t83 = 0.2e1 * (m(5) / 0.2e1 + m(4) / 0.2e1) * t103;
t81 = qJ(2) ^ 2;
t80 = qJ(4) ^ 2;
t64 = Ifges(6,5) * t76;
t59 = t74 * t79;
t55 = pkin(5) * t73 + qJ(4);
t51 = Ifges(6,1) * t76 - t111;
t50 = -Ifges(6,2) * t73 + t110;
t46 = -pkin(4) * t74 + t59;
t31 = t93 * t74;
t29 = t59 + (-pkin(5) * t76 - pkin(4)) * t74;
t23 = Ifges(6,5) * t77 + (Ifges(6,1) * t73 + t110) * t74;
t22 = Ifges(6,6) * t77 + (Ifges(6,2) * t76 + t111) * t74;
t18 = mrSges(7,1) * t77 - mrSges(7,3) * t27;
t17 = -mrSges(7,2) * t77 - mrSges(7,3) * t25;
t14 = -Ifges(7,1) * t87 - Ifges(7,4) * t88;
t13 = -Ifges(7,4) * t87 - Ifges(7,2) * t88;
t12 = mrSges(7,1) * t88 - mrSges(7,2) * t87;
t7 = mrSges(7,1) * t25 + mrSges(7,2) * t27;
t5 = Ifges(7,1) * t27 - Ifges(7,4) * t25 + Ifges(7,5) * t77;
t4 = Ifges(7,4) * t27 - Ifges(7,2) * t25 + Ifges(7,6) * t77;
t1 = [-(2 * pkin(1) * mrSges(3,2)) + mrSges(3,3) * t115 + 0.2e1 * t10 * t43 + 0.2e1 * t3 * t17 + 0.2e1 * t2 * t18 - t25 * t4 + t27 * t5 + 0.2e1 * t29 * t7 + 0.2e1 * t46 * t31 + 0.2e1 * t9 * t42 + Ifges(3,1) + Ifges(2,3) + (mrSges(4,2) * t115 + mrSges(5,3) * t116 + (Ifges(5,2) + Ifges(4,1)) * t77 + t101 + t119) * t77 + (mrSges(4,1) * t115 + mrSges(5,2) * t116 + t76 * t22 + t73 * t23 + (Ifges(5,3) + Ifges(4,2)) * t74 + 0.2e1 * (-Ifges(4,4) - Ifges(5,6)) * t77) * t74 + m(4) * (t81 + t105) + m(3) * ((pkin(1) ^ 2) + t81) + m(6) * (t10 ^ 2 + t46 ^ 2 + t9 ^ 2) + m(5) * (t44 ^ 2 + t105) + m(7) * (t2 ^ 2 + t29 ^ 2 + t3 ^ 2) + 0.2e1 * t79 * t121; -m(3) * pkin(1) - t28 * t17 + t26 * t18 + mrSges(3,2) + t89 * t77 + (t31 + t7) * t74 + m(7) * (t2 * t26 - t28 * t3 + t29 * t74) + m(6) * (t74 * t46 - t94 * t77) + t79 * t83 + t121; m(3) + m(7) * (t26 ^ 2 + t28 ^ 2 + t67) + m(6) * (t104 * t69 + t67) + t83; t55 * t7 - t87 * t5 / 0.2e1 + t46 * t49 + qJ(4) * t31 - t88 * t4 / 0.2e1 - t25 * t13 / 0.2e1 + t27 * t14 / 0.2e1 + t29 * t12 + t16 * t17 + t15 * t18 - t95 * mrSges(7,3) + (t78 * t42 - t9 * mrSges(6,3) + t23 / 0.2e1) * t76 + (t78 * t43 - t10 * mrSges(6,3) - t22 / 0.2e1) * t73 + m(6) * (qJ(4) * t46 + t94 * t78) + m(7) * (t15 * t2 + t16 * t3 + t29 * t55) + (-pkin(3) * mrSges(5,1) - Ifges(6,6) * t73 / 0.2e1 + t64 / 0.2e1 - t34 / 0.2e1 - t33 / 0.2e1 - Ifges(5,4) + Ifges(4,5)) * t77 + (t76 * t50 / 0.2e1 + t73 * t51 / 0.2e1 - qJ(4) * mrSges(5,1) + Ifges(5,5) - Ifges(4,6)) * t74 + (t118 + (mrSges(4,1) - mrSges(5,2)) * t77 + (-mrSges(4,2) + mrSges(5,3)) * t74) * t79; -t91 * mrSges(7,3) + (mrSges(4,1) + t117) * t77 + (-mrSges(4,2) + t12 + t106) * t74 + m(7) * (t15 * t26 - t16 * t28 + t55 * t74) + m(6) * (-t77 * t99 + t61) + t118; -0.2e1 * pkin(3) * mrSges(5,2) + 0.2e1 * t55 * t12 - t88 * t13 - t87 * t14 - t73 * t50 + t76 * t51 + Ifges(5,1) + Ifges(4,3) + m(7) * (t15 ^ 2 + t16 ^ 2 + t55 ^ 2) + m(6) * (t104 * t78 ^ 2 + t80) + m(5) * (pkin(3) ^ 2 + t80) - 0.2e1 * mrSges(6,3) * t99 - 0.2e1 * t92 * mrSges(7,3) + 0.2e1 * t106 * qJ(4); t88 * t17 - t87 * t18 + (-m(5) * t79 + mrSges(5,1)) * t77 + m(7) * t95 + m(6) * t94 - t89; m(7) * t91 + t122 * t77; -m(5) * pkin(3) + m(7) * t92 - t102 * mrSges(7,3) + t78 * t100 - t117; m(7) * t102 - t122; t9 * mrSges(6,1) - t10 * mrSges(6,2) + (t72 * t17 + t75 * t18 + m(7) * (t2 * t75 + t3 * t72)) * pkin(5) + t85 + t119; t93 * t77 + (t26 * t75 - t28 * t72) * t114 + t98; t78 * t107 + t64 + (-mrSges(6,2) * t78 - Ifges(6,6)) * t73 + (m(7) * (t15 * t75 + t16 * t72) + t90 * mrSges(7,3)) * pkin(5) + t86; -t90 * t114 - t93 + t97; Ifges(6,3) + Ifges(7,3) + m(7) * (t72 ^ 2 + t75 ^ 2) * pkin(5) ^ 2 + 0.2e1 * t84; t85; t98; t86; t97; Ifges(7,3) + t84; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
