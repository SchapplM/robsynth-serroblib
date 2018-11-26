% Calculate joint inertia matrix for
% S6RRPPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta4]';
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
% Datum: 2018-11-23 16:44
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRPPPR4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR4_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR4_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR4_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR4_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR4_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:43:57
% EndTime: 2018-11-23 16:43:58
% DurationCPUTime: 0.83s
% Computational Cost: add. (855->237), mult. (1543->311), div. (0->0), fcn. (1309->6), ass. (0->87)
t80 = sin(pkin(9));
t81 = cos(pkin(9));
t103 = -t80 ^ 2 - t81 ^ 2;
t84 = sin(qJ(2));
t86 = cos(qJ(2));
t119 = t84 ^ 2 + t86 ^ 2;
t118 = (mrSges(6,2) + mrSges(5,3)) * t103;
t117 = pkin(3) + pkin(7);
t116 = -pkin(4) - pkin(5);
t82 = -pkin(2) - qJ(4);
t115 = pkin(8) + t82;
t114 = Ifges(5,4) * t80;
t113 = Ifges(5,4) * t81;
t112 = Ifges(6,5) * t80;
t111 = Ifges(6,5) * t81;
t110 = t80 * t86;
t109 = t81 * t86;
t83 = sin(qJ(6));
t85 = cos(qJ(6));
t90 = t83 * t80 + t85 * t81;
t28 = t90 * t86;
t42 = t80 * t85 - t81 * t83;
t29 = t42 * t86;
t108 = Ifges(7,5) * t29 - Ifges(7,6) * t28;
t99 = -qJ(3) * t84 - pkin(1);
t38 = t82 * t86 + t99;
t58 = t117 * t84;
t11 = t81 * t38 + t80 * t58;
t44 = mrSges(5,1) * t84 + mrSges(5,3) * t110;
t45 = -t84 * mrSges(6,1) - mrSges(6,2) * t110;
t107 = t44 - t45;
t46 = -mrSges(5,2) * t84 - mrSges(5,3) * t109;
t47 = -mrSges(6,2) * t109 + mrSges(6,3) * t84;
t106 = t46 + t47;
t34 = mrSges(6,1) * t109 + mrSges(6,3) * t110;
t105 = t103 * t82 ^ 2;
t53 = t80 * mrSges(5,1) + t81 * mrSges(5,2);
t104 = t119 * pkin(7) ^ 2;
t102 = -t81 * qJ(5) + qJ(3);
t101 = t42 ^ 2 + t90 ^ 2;
t8 = t84 * qJ(5) + t11;
t100 = -m(4) * pkin(2) + mrSges(4,2);
t52 = t80 * mrSges(6,1) - t81 * mrSges(6,3);
t10 = -t80 * t38 + t58 * t81;
t96 = -qJ(5) * t80 - t117;
t35 = mrSges(5,1) * t109 - mrSges(5,2) * t110;
t3 = pkin(8) * t110 + t116 * t84 - t10;
t6 = pkin(8) * t109 + t8;
t1 = t3 * t85 - t6 * t83;
t2 = t3 * t83 + t6 * t85;
t95 = -t1 * t90 + t2 * t42;
t9 = -pkin(4) * t84 - t10;
t94 = t8 * t80 - t9 * t81;
t7 = -t28 * mrSges(7,1) - t29 * mrSges(7,2);
t16 = -t42 * mrSges(7,1) + mrSges(7,2) * t90;
t93 = t10 * t81 + t11 * t80;
t48 = t115 * t80;
t49 = t115 * t81;
t14 = -t48 * t83 - t49 * t85;
t15 = t48 * t85 - t49 * t83;
t92 = -t14 * t90 + t15 * t42;
t91 = t42 * t83 - t85 * t90;
t89 = -0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * t103;
t87 = qJ(3) ^ 2;
t59 = t117 * t86;
t57 = Ifges(5,1) * t81 - t114;
t56 = Ifges(6,1) * t81 + t112;
t55 = -Ifges(5,2) * t80 + t113;
t54 = Ifges(6,3) * t80 + t111;
t51 = -pkin(2) * t86 + t99;
t50 = pkin(4) * t80 + t102;
t37 = Ifges(7,5) * t90;
t36 = Ifges(7,6) * t42;
t32 = -t116 * t80 + t102;
t27 = Ifges(5,5) * t84 + (-Ifges(5,1) * t80 - t113) * t86;
t26 = Ifges(6,4) * t84 + (-Ifges(6,1) * t80 + t111) * t86;
t25 = Ifges(5,6) * t84 + (-Ifges(5,2) * t81 - t114) * t86;
t24 = Ifges(6,6) * t84 + (Ifges(6,3) * t81 - t112) * t86;
t21 = (pkin(4) * t81 - t96) * t86;
t20 = -mrSges(7,1) * t84 + mrSges(7,3) * t29;
t19 = mrSges(7,2) * t84 + mrSges(7,3) * t28;
t18 = Ifges(7,1) * t90 + Ifges(7,4) * t42;
t17 = Ifges(7,4) * t90 + Ifges(7,2) * t42;
t12 = (t116 * t81 + t96) * t86;
t5 = -Ifges(7,1) * t29 + Ifges(7,4) * t28 - Ifges(7,5) * t84;
t4 = -Ifges(7,4) * t29 + Ifges(7,2) * t28 - Ifges(7,6) * t84;
t13 = [0.2e1 * t1 * t20 + 0.2e1 * t10 * t44 + 0.2e1 * t11 * t46 + 0.2e1 * t12 * t7 + 0.2e1 * t2 * t19 + 0.2e1 * t21 * t34 + t28 * t4 - t29 * t5 + 0.2e1 * t59 * t35 + 0.2e1 * t9 * t45 + 0.2e1 * t8 * t47 + Ifges(2,3) + m(4) * (t51 ^ 2 + t104) + m(3) * (pkin(1) ^ 2 + t104) + m(5) * (t10 ^ 2 + t11 ^ 2 + t59 ^ 2) + m(7) * (t1 ^ 2 + t12 ^ 2 + t2 ^ 2) + m(6) * (t21 ^ 2 + t8 ^ 2 + t9 ^ 2) + (0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * t51 * mrSges(4,2) + (Ifges(3,2) + Ifges(4,3)) * t86 + (t24 - t25) * t81 + (-t26 - t27) * t80) * t86 + (-0.2e1 * pkin(1) * mrSges(3,2) - 0.2e1 * t51 * mrSges(4,3) + (Ifges(6,2) + Ifges(5,3) + Ifges(3,1) + Ifges(4,2) + Ifges(7,3)) * t84 + ((2 * Ifges(3,4)) + (2 * Ifges(4,6)) + (-Ifges(5,6) + Ifges(6,6)) * t81 + (-Ifges(6,4) - Ifges(5,5)) * t80) * t86 + t108) * t84 + 0.2e1 * (mrSges(4,1) + mrSges(3,3)) * pkin(7) * t119; t50 * t34 + t21 * t52 + t59 * t53 + t90 * t5 / 0.2e1 + t42 * t4 / 0.2e1 + t28 * t17 / 0.2e1 - t29 * t18 / 0.2e1 - t32 * t7 + qJ(3) * t35 + t12 * t16 + t15 * t19 + t14 * t20 + t95 * mrSges(7,3) + (-t10 * mrSges(5,3) + t9 * mrSges(6,2) + t26 / 0.2e1 + t27 / 0.2e1 + t107 * t82) * t81 + (-t11 * mrSges(5,3) - t8 * mrSges(6,2) + t24 / 0.2e1 - t25 / 0.2e1 + t106 * t82) * t80 + m(6) * (t21 * t50 + t94 * t82) + m(5) * (qJ(3) * t59 + t93 * t82) + m(7) * (t1 * t14 - t12 * t32 + t15 * t2) + (-pkin(2) * mrSges(4,1) - t37 / 0.2e1 - t36 / 0.2e1 - Ifges(4,4) + Ifges(3,5) + (Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1) * t81 + (Ifges(6,6) / 0.2e1 - Ifges(5,6) / 0.2e1) * t80 + (-mrSges(3,1) + t100) * pkin(7)) * t84 + (qJ(3) * mrSges(4,1) - Ifges(4,5) + Ifges(3,6) + (t54 / 0.2e1 - t55 / 0.2e1) * t81 + (-t56 / 0.2e1 - t57 / 0.2e1) * t80 + (m(4) * qJ(3) - mrSges(3,2) + mrSges(4,3)) * pkin(7)) * t86; -0.2e1 * pkin(2) * mrSges(4,2) - 0.2e1 * t32 * t16 + t42 * t17 + t90 * t18 + 0.2e1 * t50 * t52 + Ifges(4,1) + Ifges(3,3) + (t57 + t56) * t81 + (-t55 + t54) * t80 + m(7) * (t14 ^ 2 + t15 ^ 2 + t32 ^ 2) + m(6) * (t50 ^ 2 - t105) + m(5) * (t87 - t105) + m(4) * (pkin(2) ^ 2 + t87) + 0.2e1 * t82 * t118 + 0.2e1 * (mrSges(4,3) + t53) * qJ(3) + 0.2e1 * t92 * mrSges(7,3); t42 * t19 - t90 * t20 + (m(4) * pkin(7) + mrSges(4,1)) * t84 + t107 * t81 + t106 * t80 + m(7) * t95 + m(6) * t94 + m(5) * t93; m(7) * t92 + t101 * mrSges(7,3) + t82 * t89 + t100 + t118; m(7) * t101 + m(4) + t89; m(5) * t59 + m(6) * t21 - m(7) * t12 + t34 + t35 - t7; m(5) * qJ(3) + m(6) * t50 + m(7) * t32 - t16 + t52 + t53; 0; m(5) + m(6) + m(7); t83 * t19 + t85 * t20 + m(7) * (t1 * t85 + t2 * t83) + m(6) * t9 + t45; m(7) * (t14 * t85 + t15 * t83) + (-m(6) * t82 + mrSges(6,2)) * t81 + t91 * mrSges(7,3); -m(6) * t81 + m(7) * t91; 0; m(6) + m(7) * (t83 ^ 2 + t85 ^ 2); mrSges(7,1) * t1 - mrSges(7,2) * t2 - Ifges(7,3) * t84 - t108; mrSges(7,1) * t14 - mrSges(7,2) * t15 + t36 + t37; -mrSges(7,1) * t90 - mrSges(7,2) * t42; 0; mrSges(7,1) * t85 - mrSges(7,2) * t83; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t13(1) t13(2) t13(4) t13(7) t13(11) t13(16); t13(2) t13(3) t13(5) t13(8) t13(12) t13(17); t13(4) t13(5) t13(6) t13(9) t13(13) t13(18); t13(7) t13(8) t13(9) t13(10) t13(14) t13(19); t13(11) t13(12) t13(13) t13(14) t13(15) t13(20); t13(16) t13(17) t13(18) t13(19) t13(20) t13(21);];
Mq  = res;
