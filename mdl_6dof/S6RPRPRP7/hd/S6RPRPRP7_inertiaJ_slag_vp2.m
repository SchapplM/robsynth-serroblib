% Calculate joint inertia matrix for
% S6RPRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
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
% Datum: 2018-11-23 16:01
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRP7_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP7_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP7_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP7_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP7_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP7_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:00:41
% EndTime: 2018-11-23 16:00:41
% DurationCPUTime: 0.71s
% Computational Cost: add. (973->217), mult. (1714->292), div. (0->0), fcn. (1673->6), ass. (0->86)
t110 = Ifges(6,3) + Ifges(7,3);
t69 = cos(qJ(3));
t109 = t69 ^ 2;
t108 = 2 * mrSges(7,3);
t107 = m(6) + m(7);
t67 = sin(qJ(3));
t70 = -pkin(1) - pkin(7);
t83 = -qJ(4) + t70;
t39 = t83 * t67;
t64 = sin(pkin(9));
t65 = cos(pkin(9));
t77 = t83 * t69;
t21 = t39 * t64 - t65 * t77;
t106 = t21 ^ 2;
t36 = t64 * t67 - t65 * t69;
t35 = t36 ^ 2;
t105 = -2 * mrSges(5,3);
t103 = m(5) * pkin(3);
t102 = m(7) * pkin(5);
t101 = pkin(3) * t64;
t100 = pkin(3) * t65;
t66 = sin(qJ(5));
t99 = Ifges(6,4) * t66;
t68 = cos(qJ(5));
t98 = Ifges(6,4) * t68;
t97 = Ifges(7,4) * t66;
t96 = Ifges(7,4) * t68;
t95 = t36 * t21;
t94 = t36 * t66;
t93 = t36 * t68;
t92 = t66 * mrSges(6,2);
t91 = t66 * mrSges(7,3);
t90 = Ifges(6,6) + Ifges(7,6);
t37 = t64 * t69 + t65 * t67;
t16 = -mrSges(7,2) * t37 + t36 * t91;
t17 = -mrSges(6,2) * t37 + mrSges(6,3) * t94;
t89 = t16 + t17;
t18 = mrSges(7,1) * t37 + mrSges(7,3) * t93;
t19 = mrSges(6,1) * t37 + mrSges(6,3) * t93;
t88 = t18 + t19;
t52 = t67 * pkin(3) + qJ(2);
t15 = pkin(4) * t37 + pkin(8) * t36 + t52;
t23 = t65 * t39 + t64 * t77;
t4 = t66 * t15 + t68 * t23;
t13 = -mrSges(7,1) * t94 - mrSges(7,2) * t93;
t42 = -mrSges(6,1) * t68 + t92;
t87 = t42 - mrSges(5,1);
t86 = t66 ^ 2 + t68 ^ 2;
t85 = t67 ^ 2 + t109;
t84 = qJ(6) * t36;
t50 = pkin(8) + t101;
t82 = qJ(6) + t50;
t51 = -pkin(4) - t100;
t80 = m(4) * t85;
t79 = t86 * t50;
t78 = t85 * mrSges(4,3);
t53 = t66 * mrSges(7,2);
t41 = -t68 * mrSges(7,1) + t53;
t3 = t68 * t15 - t23 * t66;
t76 = (-Ifges(6,5) - Ifges(7,5)) * t93 + t110 * t37;
t75 = mrSges(6,1) + mrSges(7,1) + t102;
t74 = -t3 * t66 + t4 * t68;
t73 = -mrSges(6,1) * t66 - mrSges(6,2) * t68;
t71 = qJ(2) ^ 2;
t57 = Ifges(6,5) * t66;
t56 = Ifges(7,5) * t66;
t55 = Ifges(6,6) * t68;
t54 = Ifges(7,6) * t68;
t46 = Ifges(6,1) * t66 + t98;
t45 = Ifges(7,1) * t66 + t96;
t44 = Ifges(6,2) * t68 + t99;
t43 = Ifges(7,2) * t68 + t97;
t40 = -pkin(5) * t68 + t51;
t34 = t37 ^ 2;
t32 = t82 * t68;
t31 = t82 * t66;
t28 = t36 * mrSges(5,2);
t14 = t73 * t36;
t10 = -pkin(5) * t94 + t21;
t8 = Ifges(6,5) * t37 + (-Ifges(6,1) * t68 + t99) * t36;
t7 = Ifges(7,5) * t37 + (-Ifges(7,1) * t68 + t97) * t36;
t6 = Ifges(6,6) * t37 + (Ifges(6,2) * t66 - t98) * t36;
t5 = Ifges(7,6) * t37 + (Ifges(7,2) * t66 - t96) * t36;
t2 = t66 * t84 + t4;
t1 = pkin(5) * t37 + t68 * t84 + t3;
t9 = [Ifges(4,1) * t109 - 0.2e1 * t52 * t28 + 0.2e1 * t10 * t13 + 0.2e1 * t2 * t16 + 0.2e1 * t4 * t17 + 0.2e1 * t1 * t18 + 0.2e1 * t3 * t19 + 0.2e1 * t21 * t14 - (2 * pkin(1) * mrSges(3,2)) + Ifges(3,1) + Ifges(2,3) - 0.2e1 * t70 * t78 + (0.2e1 * mrSges(5,1) * t52 + Ifges(5,2) * t37 + t23 * t105 + t76) * t37 + (t21 * t105 + Ifges(5,1) * t36 + 0.2e1 * Ifges(5,4) * t37 + (-t7 - t8) * t68 + (t90 * t37 + t5 + t6) * t66) * t36 + m(4) * (t85 * t70 ^ 2 + t71) + m(3) * ((pkin(1) ^ 2) + t71) + m(5) * (t23 ^ 2 + t52 ^ 2 + t106) + m(6) * (t3 ^ 2 + t4 ^ 2 + t106) + m(7) * (t1 ^ 2 + t10 ^ 2 + t2 ^ 2) + (-0.2e1 * Ifges(4,4) * t69 + Ifges(4,2) * t67) * t67 + 0.2e1 * (mrSges(4,1) * t67 + mrSges(4,2) * t69 + mrSges(3,3)) * qJ(2); -m(3) * pkin(1) - t35 * mrSges(5,3) + mrSges(3,2) + (t13 + t14) * t36 - t78 + (-mrSges(5,3) * t37 - t88 * t66 + t89 * t68) * t37 + m(7) * (t36 * t10 + (-t1 * t66 + t2 * t68) * t37) + m(6) * (t74 * t37 + t95) + m(5) * (t23 * t37 + t95) + t70 * t80; m(3) + m(5) * (t34 + t35) + t80 + (t86 * t34 + t35) * t107; -t23 * mrSges(5,2) + t10 * t41 + t40 * t13 + t51 * t14 + t32 * t16 - t31 * t18 + (mrSges(4,1) * t70 + Ifges(4,5)) * t69 + (-mrSges(4,2) * t70 - Ifges(4,6)) * t67 + (t57 / 0.2e1 + t55 / 0.2e1 + t56 / 0.2e1 + t54 / 0.2e1 - Ifges(5,6) - mrSges(5,3) * t101) * t37 + t87 * t21 + (t5 / 0.2e1 + t6 / 0.2e1 + t50 * t17 + t2 * mrSges(7,3) + t4 * mrSges(6,3)) * t68 + (t8 / 0.2e1 + t7 / 0.2e1 - t50 * t19 - t1 * mrSges(7,3) - t3 * mrSges(6,3)) * t66 + m(6) * (t51 * t21 + t74 * t50) + m(7) * (-t1 * t31 + t10 * t40 + t2 * t32) + (-t21 * t65 + t23 * t64) * t103 + (mrSges(5,3) * t100 - Ifges(5,5) + (-t45 / 0.2e1 - t46 / 0.2e1) * t68 + (t43 / 0.2e1 + t44 / 0.2e1) * t66) * t36; t69 * mrSges(4,1) - t67 * mrSges(4,2) + (t41 + t87) * t36 + (-mrSges(5,2) + (mrSges(6,3) + mrSges(7,3)) * t86) * t37 + m(7) * (t40 * t36 + (t31 * t66 + t32 * t68) * t37) + m(6) * (t51 * t36 + t37 * t79) + (-t36 * t65 + t37 * t64) * t103; 0.2e1 * t40 * t41 + 0.2e1 * t51 * t42 + Ifges(4,3) + Ifges(5,3) + (t32 * t108 + t43 + t44) * t68 + (t31 * t108 + t45 + t46) * t66 + m(6) * (t86 * t50 ^ 2 + t51 ^ 2) + m(7) * (t31 ^ 2 + t32 ^ 2 + t40 ^ 2) + m(5) * (t64 ^ 2 + t65 ^ 2) * pkin(3) ^ 2 + 0.2e1 * (mrSges(5,1) * t65 - mrSges(5,2) * t64) * pkin(3) + 0.2e1 * mrSges(6,3) * t79; t37 * mrSges(5,1) - t28 + t88 * t68 + t89 * t66 + m(7) * (t1 * t68 + t2 * t66) + m(6) * (t3 * t68 + t4 * t66) + m(5) * t52; 0; m(7) * (-t31 * t68 + t32 * t66); t86 * t107 + m(5); t3 * mrSges(6,1) + t1 * mrSges(7,1) - t4 * mrSges(6,2) - t2 * mrSges(7,2) + t90 * t94 + (m(7) * t1 + t18) * pkin(5) + t76; ((-mrSges(6,2) - mrSges(7,2)) * t68 - t75 * t66) * t37; -t31 * mrSges(7,1) - t32 * mrSges(7,2) + t54 + t55 + t56 + t57 + t73 * t50 + (-m(7) * t31 - t91) * pkin(5); t75 * t68 - t53 - t92; (0.2e1 * mrSges(7,1) + t102) * pkin(5) + t110; m(7) * t10 + t13; m(7) * t36; m(7) * t40 + t41; 0; 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t9(1) t9(2) t9(4) t9(7) t9(11) t9(16); t9(2) t9(3) t9(5) t9(8) t9(12) t9(17); t9(4) t9(5) t9(6) t9(9) t9(13) t9(18); t9(7) t9(8) t9(9) t9(10) t9(14) t9(19); t9(11) t9(12) t9(13) t9(14) t9(15) t9(20); t9(16) t9(17) t9(18) t9(19) t9(20) t9(21);];
Mq  = res;
