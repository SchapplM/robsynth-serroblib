% Calculate joint inertia matrix for
% S6RPPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
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
% Datum: 2019-03-09 02:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRP8_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP8_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP8_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP8_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP8_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP8_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:15:18
% EndTime: 2019-03-09 02:15:19
% DurationCPUTime: 0.63s
% Computational Cost: add. (920->195), mult. (1634->256), div. (0->0), fcn. (1601->6), ass. (0->80)
t109 = Ifges(7,2) + Ifges(6,3);
t65 = sin(qJ(5));
t67 = cos(qJ(5));
t83 = t65 ^ 2 + t67 ^ 2;
t108 = m(6) + m(7);
t107 = (mrSges(7,2) + mrSges(6,3)) * t83;
t64 = -pkin(1) - qJ(3);
t102 = -pkin(7) + t64;
t62 = sin(pkin(9));
t38 = t102 * t62;
t66 = sin(qJ(4));
t101 = cos(qJ(4));
t63 = cos(pkin(9));
t80 = t101 * t63;
t20 = -t102 * t80 + t66 * t38;
t106 = t20 ^ 2;
t35 = t66 * t62 - t80;
t34 = t35 ^ 2;
t58 = t63 ^ 2;
t105 = -2 * mrSges(5,3);
t100 = Ifges(6,4) * t65;
t99 = Ifges(6,4) * t67;
t98 = Ifges(7,5) * t65;
t97 = Ifges(7,5) * t67;
t91 = t66 * t63;
t36 = t101 * t62 + t91;
t96 = Ifges(6,6) * t36;
t95 = Ifges(7,6) * t36;
t94 = t35 * t20;
t93 = t35 * t65;
t92 = t35 * t67;
t47 = t62 * pkin(3) + qJ(2);
t12 = t36 * pkin(4) + t35 * pkin(8) + t47;
t22 = t101 * t38 + t102 * t91;
t4 = t65 * t12 + t67 * t22;
t16 = -t36 * mrSges(6,2) + mrSges(6,3) * t93;
t19 = mrSges(7,2) * t93 + t36 * mrSges(7,3);
t90 = t16 + t19;
t17 = t36 * mrSges(6,1) + mrSges(6,3) * t92;
t18 = -t36 * mrSges(7,1) - mrSges(7,2) * t92;
t89 = t17 - t18;
t88 = t83 * pkin(8) * t36;
t41 = -t67 * mrSges(6,1) + t65 * mrSges(6,2);
t87 = t41 - mrSges(5,1);
t86 = t62 * mrSges(4,1) + t63 * mrSges(4,2);
t85 = t83 * pkin(8) ^ 2;
t84 = t62 ^ 2 + t58;
t81 = m(4) * t84;
t79 = t84 * mrSges(4,3);
t77 = -Ifges(7,6) * t93 + (-Ifges(7,4) - Ifges(6,5)) * t92 + t109 * t36;
t1 = t36 * qJ(6) + t4;
t3 = t67 * t12 - t65 * t22;
t2 = -t36 * pkin(5) - t3;
t76 = t1 * t67 + t2 * t65;
t75 = -t3 * t65 + t4 * t67;
t74 = -t65 * mrSges(6,1) - t67 * mrSges(6,2);
t40 = -t67 * mrSges(7,1) - t65 * mrSges(7,3);
t73 = -t65 * mrSges(7,1) + t67 * mrSges(7,3);
t72 = t67 * pkin(5) + t65 * qJ(6);
t71 = -pkin(5) * t65 + qJ(6) * t67;
t70 = m(7) * t71 + t73 + t74;
t68 = qJ(2) ^ 2;
t54 = Ifges(7,4) * t65;
t53 = Ifges(6,5) * t65;
t52 = Ifges(6,6) * t67;
t45 = Ifges(6,1) * t65 + t99;
t44 = Ifges(7,1) * t65 - t97;
t43 = Ifges(6,2) * t67 + t100;
t42 = -Ifges(7,3) * t67 + t98;
t39 = -pkin(4) - t72;
t33 = t36 ^ 2;
t30 = t35 * mrSges(5,2);
t14 = t74 * t35;
t13 = t73 * t35;
t9 = Ifges(6,5) * t36 + (-Ifges(6,1) * t67 + t100) * t35;
t8 = Ifges(7,4) * t36 + (-Ifges(7,1) * t67 - t98) * t35;
t7 = t96 + (Ifges(6,2) * t65 - t99) * t35;
t6 = t95 + (-Ifges(7,3) * t65 - t97) * t35;
t5 = t71 * t35 + t20;
t10 = [Ifges(4,1) * t58 - 0.2e1 * t47 * t30 + 0.2e1 * t1 * t19 + 0.2e1 * t20 * t14 + 0.2e1 * t5 * t13 + 0.2e1 * t4 * t16 + 0.2e1 * t3 * t17 + 0.2e1 * t2 * t18 - (2 * pkin(1) * mrSges(3,2)) + Ifges(3,1) + Ifges(2,3) + (-0.2e1 * Ifges(4,4) * t63 + Ifges(4,2) * t62) * t62 - 0.2e1 * t64 * t79 + (0.2e1 * t47 * mrSges(5,1) + Ifges(5,2) * t36 + t22 * t105 + t77) * t36 + (t20 * t105 + Ifges(5,1) * t35 + 0.2e1 * Ifges(5,4) * t36 + (-t8 - t9) * t67 + (-t6 + t7 + t96) * t65) * t35 + m(4) * (t84 * t64 ^ 2 + t68) + m(3) * ((pkin(1) ^ 2) + t68) + m(5) * (t22 ^ 2 + t47 ^ 2 + t106) + m(6) * (t3 ^ 2 + t4 ^ 2 + t106) + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) + 0.2e1 * (t86 + mrSges(3,3)) * qJ(2); -m(3) * pkin(1) - t34 * mrSges(5,3) + mrSges(3,2) + (t13 + t14) * t35 - t79 + (-mrSges(5,3) * t36 - t89 * t65 + t90 * t67) * t36 + m(7) * (t35 * t5 + t76 * t36) + m(6) * (t75 * t36 + t94) + m(5) * (t36 * t22 + t94) + t64 * t81; m(3) + m(5) * (t33 + t34) + t81 + (t83 * t33 + t34) * t108; m(4) * qJ(2) + t36 * mrSges(5,1) - t30 + t89 * t67 + t90 * t65 + m(7) * (t65 * t1 - t67 * t2) + m(6) * (t67 * t3 + t65 * t4) + m(5) * t47 + t86; 0; t83 * t108 + m(4) + m(5); -t22 * mrSges(5,2) - pkin(4) * t14 + t39 * t13 + t5 * t40 + (t54 / 0.2e1 + t53 / 0.2e1 + t52 / 0.2e1 - Ifges(5,6)) * t36 + t87 * t20 + (t4 * mrSges(6,3) + t1 * mrSges(7,2) - t95 / 0.2e1 - t6 / 0.2e1 + t7 / 0.2e1 + t90 * pkin(8)) * t67 + (t2 * mrSges(7,2) - t3 * mrSges(6,3) + t8 / 0.2e1 + t9 / 0.2e1 - t89 * pkin(8)) * t65 + m(7) * (t76 * pkin(8) + t39 * t5) + m(6) * (-pkin(4) * t20 + t75 * pkin(8)) + (-Ifges(5,5) + (-t44 / 0.2e1 - t45 / 0.2e1) * t67 + (-t42 / 0.2e1 + t43 / 0.2e1) * t65) * t35; (t40 + t87) * t35 + m(7) * (t39 * t35 + t88) + m(6) * (-pkin(4) * t35 + t88) + (-mrSges(5,2) + t107) * t36; 0; -0.2e1 * pkin(4) * t41 + 0.2e1 * t39 * t40 + Ifges(5,3) + (t43 - t42) * t67 + (t44 + t45) * t65 + m(7) * (t39 ^ 2 + t85) + m(6) * (pkin(4) ^ 2 + t85) + 0.2e1 * pkin(8) * t107; Ifges(6,6) * t93 + m(7) * (-pkin(5) * t2 + qJ(6) * t1) + qJ(6) * t19 + t1 * mrSges(7,3) - pkin(5) * t18 + t3 * mrSges(6,1) - t2 * mrSges(7,1) - t4 * mrSges(6,2) + t77; t70 * t36; m(7) * t72 - t40 - t41; t71 * mrSges(7,2) - Ifges(7,6) * t67 + t70 * pkin(8) + t52 + t53 + t54; 0.2e1 * pkin(5) * mrSges(7,1) + 0.2e1 * qJ(6) * mrSges(7,3) + m(7) * (pkin(5) ^ 2 + qJ(6) ^ 2) + t109; m(7) * t2 + t18; m(7) * t65 * t36; -m(7) * t67; (m(7) * pkin(8) + mrSges(7,2)) * t65; -m(7) * pkin(5) - mrSges(7,1); m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t10(1) t10(2) t10(4) t10(7) t10(11) t10(16); t10(2) t10(3) t10(5) t10(8) t10(12) t10(17); t10(4) t10(5) t10(6) t10(9) t10(13) t10(18); t10(7) t10(8) t10(9) t10(10) t10(14) t10(19); t10(11) t10(12) t10(13) t10(14) t10(15) t10(20); t10(16) t10(17) t10(18) t10(19) t10(20) t10(21);];
Mq  = res;
