% Calculate joint inertia matrix for
% S6RRPPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3]';
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

function Mq = S6RRPPPR5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR5_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR5_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR5_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR5_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR5_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:44:31
% EndTime: 2018-11-23 16:44:32
% DurationCPUTime: 1.02s
% Computational Cost: add. (855->280), mult. (1682->356), div. (0->0), fcn. (1432->6), ass. (0->102)
t126 = pkin(4) + qJ(3);
t93 = sin(pkin(9));
t94 = cos(pkin(9));
t125 = t93 ^ 2 + t94 ^ 2;
t124 = 2 * pkin(7);
t96 = sin(qJ(2));
t108 = t94 * t96;
t98 = cos(qJ(2));
t123 = pkin(4) * t108 + qJ(5) * t98;
t106 = pkin(3) + qJ(5);
t122 = t106 * t94;
t121 = -t106 * t93 - pkin(7);
t120 = -2 * mrSges(6,2);
t119 = m(4) * pkin(2);
t102 = qJ(4) * t93 + pkin(2);
t56 = -pkin(3) * t94 - t102;
t118 = m(5) * t56;
t117 = pkin(7) * t98;
t95 = sin(qJ(6));
t97 = cos(qJ(6));
t116 = m(6) + m(7) * (t95 ^ 2 + t97 ^ 2);
t115 = Ifges(4,4) * t93;
t114 = Ifges(4,4) * t94;
t113 = Ifges(6,5) * t93;
t112 = Ifges(6,5) * t94;
t111 = Ifges(5,6) * t93;
t110 = Ifges(5,6) * t94;
t109 = t93 * t96;
t107 = -mrSges(6,1) - mrSges(5,3);
t105 = -pkin(5) - qJ(4);
t48 = t93 * t95 - t94 * t97;
t35 = t48 * t96;
t47 = t93 * t97 + t94 * t95;
t36 = t47 * t96;
t104 = -Ifges(7,5) * t36 + Ifges(7,6) * t35;
t40 = mrSges(4,1) * t109 + mrSges(4,2) * t108;
t57 = -pkin(2) * t98 - qJ(3) * t96 - pkin(1);
t23 = t117 * t94 + t57 * t93;
t103 = t125 * qJ(3) ^ 2;
t58 = t126 * t93;
t61 = t126 * t94;
t50 = -t98 * mrSges(6,1) + mrSges(6,2) * t109;
t7 = t35 * mrSges(7,1) + mrSges(7,2) * t36;
t13 = -t47 * mrSges(7,1) + mrSges(7,2) * t48;
t54 = mrSges(5,1) * t108 - t98 * mrSges(5,2);
t74 = t93 * t117;
t22 = t57 * t94 - t74;
t52 = -mrSges(6,2) * t108 + mrSges(6,3) * t98;
t88 = t98 * pkin(3);
t19 = -t22 + t88;
t18 = qJ(4) * t98 - t23;
t100 = pkin(7) ^ 2;
t92 = t98 ^ 2;
t91 = t96 ^ 2;
t87 = t91 * t100;
t82 = t94 * mrSges(5,2);
t81 = t93 * mrSges(4,2);
t69 = qJ(4) * t108;
t68 = Ifges(4,1) * t93 + t114;
t67 = -Ifges(6,1) * t94 + t113;
t66 = Ifges(4,2) * t94 + t115;
t65 = Ifges(6,3) * t93 - t112;
t64 = -Ifges(5,2) * t93 - t110;
t63 = -Ifges(5,3) * t94 - t111;
t62 = t93 * mrSges(6,1) + t94 * mrSges(6,3);
t60 = -t94 * mrSges(4,1) + t81;
t59 = -t93 * mrSges(5,3) + t82;
t53 = mrSges(5,1) * t109 + mrSges(5,3) * t98;
t51 = -mrSges(4,1) * t98 - mrSges(4,3) * t108;
t49 = mrSges(4,2) * t98 - mrSges(4,3) * t109;
t45 = pkin(8) * t94 + t61;
t44 = pkin(8) * t93 + t58;
t43 = Ifges(7,5) * t48;
t42 = Ifges(7,6) * t47;
t39 = (-t93 * mrSges(5,2) - t94 * mrSges(5,3)) * t96;
t38 = (t94 * mrSges(6,1) - t93 * mrSges(6,3)) * t96;
t37 = t102 + t122;
t34 = -t69 + (pkin(3) * t93 + pkin(7)) * t96;
t33 = -Ifges(5,4) * t98 + (-Ifges(5,2) * t94 + t111) * t96;
t32 = -Ifges(5,5) * t98 + (Ifges(5,3) * t93 - t110) * t96;
t31 = -Ifges(4,5) * t98 + (Ifges(4,1) * t94 - t115) * t96;
t30 = Ifges(6,4) * t98 + (Ifges(6,1) * t93 + t112) * t96;
t29 = -Ifges(4,6) * t98 + (-Ifges(4,2) * t93 + t114) * t96;
t28 = Ifges(6,6) * t98 + (Ifges(6,3) * t94 + t113) * t96;
t24 = t105 * t93 - pkin(2) - t122;
t21 = -mrSges(7,1) * t98 - mrSges(7,3) * t36;
t20 = mrSges(7,2) * t98 - mrSges(7,3) * t35;
t16 = -t121 * t96 - t69;
t15 = Ifges(7,1) * t48 + Ifges(7,4) * t47;
t14 = Ifges(7,4) * t48 + Ifges(7,2) * t47;
t12 = -pkin(4) * t109 - t18;
t11 = -t69 + (-pkin(5) * t94 - t121) * t96;
t10 = t44 * t97 + t45 * t95;
t9 = -t44 * t95 + t45 * t97;
t8 = t19 + t123;
t6 = t105 * t98 + (-pkin(4) - pkin(8)) * t109 + t23;
t5 = t74 + t88 + (pkin(8) * t96 - t57) * t94 + t123;
t4 = Ifges(7,1) * t36 - Ifges(7,4) * t35 - Ifges(7,5) * t98;
t3 = Ifges(7,4) * t36 - Ifges(7,2) * t35 - Ifges(7,6) * t98;
t2 = t5 * t97 + t6 * t95;
t1 = -t5 * t95 + t6 * t97;
t17 = [0.2e1 * t1 * t21 + 0.2e1 * t11 * t7 + 0.2e1 * t12 * t50 - 0.2e1 * t16 * t38 + 0.2e1 * t18 * t53 + 0.2e1 * t19 * t54 + 0.2e1 * t2 * t20 + 0.2e1 * t22 * t51 + 0.2e1 * t23 * t49 - t35 * t3 + 0.2e1 * t34 * t39 + t36 * t4 + 0.2e1 * t8 * t52 + Ifges(2,3) + (t91 + t92) * mrSges(3,3) * t124 + m(3) * (pkin(1) ^ 2 + t100 * t92 + t87) + m(4) * (t22 ^ 2 + t23 ^ 2 + t87) + m(5) * (t18 ^ 2 + t19 ^ 2 + t34 ^ 2) + m(7) * (t1 ^ 2 + t11 ^ 2 + t2 ^ 2) + m(6) * (t12 ^ 2 + t16 ^ 2 + t8 ^ 2) + (0.2e1 * pkin(1) * mrSges(3,1) + (Ifges(3,2) + Ifges(5,1) + Ifges(4,3) + Ifges(6,2) + Ifges(7,3)) * t98 + t104) * t98 + (-0.2e1 * pkin(1) * mrSges(3,2) + Ifges(3,1) * t96 + t40 * t124 + (t28 + t31 - t33) * t94 + (-t29 + t30 + t32) * t93 + ((2 * Ifges(3,4)) + (Ifges(5,4) - Ifges(4,5) + Ifges(6,6)) * t94 + (Ifges(6,4) - Ifges(5,5) + Ifges(4,6)) * t93) * t98) * t96; m(6) * (t12 * t61 - t16 * t37 + t58 * t8) + m(7) * (t1 * t9 + t10 * t2 + t11 * t24) + (-t12 * mrSges(6,2) - t18 * mrSges(5,1) + t23 * mrSges(4,3) + t29 / 0.2e1 - t30 / 0.2e1 - t32 / 0.2e1 + (Ifges(5,5) / 0.2e1 - Ifges(4,6) / 0.2e1 - Ifges(6,4) / 0.2e1) * t98) * t94 + (t19 * mrSges(5,1) - t22 * mrSges(4,3) - t8 * mrSges(6,2) - t33 / 0.2e1 + t28 / 0.2e1 + t31 / 0.2e1 + (Ifges(5,4) / 0.2e1 - Ifges(4,5) / 0.2e1 + Ifges(6,6) / 0.2e1) * t98) * t93 + (-t1 * t48 + t2 * t47) * mrSges(7,3) + (Ifges(3,5) + (-t64 / 0.2e1 + t65 / 0.2e1 + t68 / 0.2e1) * t94 + (t63 / 0.2e1 - t66 / 0.2e1 + t67 / 0.2e1) * t93 + (-mrSges(3,1) + t60 - t119) * pkin(7)) * t96 + t47 * t3 / 0.2e1 + t48 * t4 / 0.2e1 + t56 * t39 + t58 * t52 + t61 * t50 - t16 * t62 - t35 * t14 / 0.2e1 + t36 * t15 / 0.2e1 + t37 * t38 - pkin(2) * t40 + t10 * t20 + t9 * t21 + t24 * t7 + t11 * t13 + (-pkin(7) * mrSges(3,2) - t43 / 0.2e1 - t42 / 0.2e1 + Ifges(3,6)) * t98 + (t59 + t118) * t34 + ((t49 - t53) * t94 + (-t51 + t54) * t93 + m(4) * (-t22 * t93 + t23 * t94) + m(5) * (-t18 * t94 + t19 * t93)) * qJ(3); -0.2e1 * pkin(2) * t60 + 0.2e1 * t24 * t13 + t47 * t14 + t48 * t15 + 0.2e1 * t37 * t62 + 0.2e1 * t56 * t59 + Ifges(3,3) + 0.2e1 * (t10 * t47 - t48 * t9) * mrSges(7,3) + (t120 * t61 - t63 + t66 - t67) * t94 + (t120 * t58 - t64 + t65 + t68) * t93 + m(5) * (t56 ^ 2 + t103) + m(4) * (pkin(2) ^ 2 + t103) + m(6) * (t37 ^ 2 + t58 ^ 2 + t61 ^ 2) + m(7) * (t10 ^ 2 + t24 ^ 2 + t9 ^ 2) + 0.2e1 * (mrSges(5,1) + mrSges(4,3)) * qJ(3) * t125; m(5) * t34 + m(6) * t16 + m(7) * t11 + (m(4) * pkin(7) + t107 * t94 + (-mrSges(5,2) + mrSges(6,3)) * t93) * t96 + t7 + t40; -t119 + t81 + t82 + (-mrSges(6,3) - mrSges(4,1)) * t94 + t107 * t93 + t118 - m(6) * t37 + m(7) * t24 + t13; m(4) + m(5) + m(6) + m(7); t97 * t20 - t95 * t21 + m(7) * (-t1 * t95 + t2 * t97) + m(6) * t8 + m(5) * t19 + t52 + t54; (t47 * t97 + t48 * t95) * mrSges(7,3) + m(7) * (t10 * t97 - t9 * t95) + m(6) * t58 + (m(5) * qJ(3) + mrSges(5,1) - mrSges(6,2)) * t93; 0; m(5) + t116; t95 * t20 + t97 * t21 + m(7) * (t1 * t97 + t2 * t95) + m(6) * t12 + t50; -t94 * mrSges(6,2) + (t47 * t95 - t48 * t97) * mrSges(7,3) + m(7) * (t10 * t95 + t9 * t97) + m(6) * t61; 0; 0; t116; mrSges(7,1) * t1 - mrSges(7,2) * t2 - Ifges(7,3) * t98 - t104; mrSges(7,1) * t9 - mrSges(7,2) * t10 + t42 + t43; 0; -mrSges(7,1) * t95 - mrSges(7,2) * t97; mrSges(7,1) * t97 - mrSges(7,2) * t95; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t17(1) t17(2) t17(4) t17(7) t17(11) t17(16); t17(2) t17(3) t17(5) t17(8) t17(12) t17(17); t17(4) t17(5) t17(6) t17(9) t17(13) t17(18); t17(7) t17(8) t17(9) t17(10) t17(14) t17(19); t17(11) t17(12) t17(13) t17(14) t17(15) t17(20); t17(16) t17(17) t17(18) t17(19) t17(20) t17(21);];
Mq  = res;
