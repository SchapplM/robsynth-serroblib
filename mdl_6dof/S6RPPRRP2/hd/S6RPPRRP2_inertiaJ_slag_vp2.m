% Calculate joint inertia matrix for
% S6RPPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-03-09 02:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRP2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP2_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP2_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP2_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP2_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP2_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:00:13
% EndTime: 2019-03-09 02:00:14
% DurationCPUTime: 0.63s
% Computational Cost: add. (944->189), mult. (1760->256), div. (0->0), fcn. (1728->8), ass. (0->83)
t109 = Ifges(7,2) + Ifges(6,3);
t65 = sin(qJ(5));
t67 = cos(qJ(5));
t82 = t65 ^ 2 + t67 ^ 2;
t108 = m(6) + m(7);
t107 = (mrSges(7,2) + mrSges(6,3)) * t82;
t62 = sin(pkin(9));
t48 = t62 * pkin(1) + qJ(3);
t100 = pkin(7) + t48;
t63 = cos(pkin(10));
t30 = t100 * t63;
t66 = sin(qJ(4));
t61 = sin(pkin(10));
t99 = cos(qJ(4));
t80 = t99 * t61;
t14 = t100 * t80 + t66 * t30;
t106 = t14 ^ 2;
t89 = t66 * t61;
t35 = -t99 * t63 + t89;
t105 = t35 ^ 2;
t58 = t63 ^ 2;
t104 = 0.2e1 * t14;
t64 = cos(pkin(9));
t50 = -t64 * pkin(1) - pkin(2);
t38 = -t63 * pkin(3) + t50;
t103 = 0.2e1 * t38;
t102 = pkin(4) * t35;
t37 = t66 * t63 + t80;
t101 = t37 * pkin(8);
t98 = Ifges(6,4) * t65;
t97 = Ifges(6,4) * t67;
t96 = Ifges(7,5) * t65;
t95 = Ifges(7,5) * t67;
t94 = Ifges(6,6) * t35;
t93 = Ifges(7,6) * t35;
t92 = t14 * t35;
t91 = t37 * t65;
t90 = t37 * t67;
t13 = -t101 + t38 + t102;
t16 = -t100 * t89 + t99 * t30;
t4 = t65 * t13 + t67 * t16;
t19 = -t35 * mrSges(6,2) - mrSges(6,3) * t91;
t22 = -mrSges(7,2) * t91 + t35 * mrSges(7,3);
t88 = t19 + t22;
t20 = t35 * mrSges(6,1) - mrSges(6,3) * t90;
t21 = -t35 * mrSges(7,1) + mrSges(7,2) * t90;
t87 = t20 - t21;
t86 = t82 * t101;
t41 = -t67 * mrSges(6,1) + t65 * mrSges(6,2);
t85 = t41 - mrSges(5,1);
t84 = t82 * pkin(8) ^ 2;
t83 = t61 ^ 2 + t58;
t79 = -t63 * mrSges(4,1) + t61 * mrSges(4,2);
t77 = Ifges(7,6) * t91 + (Ifges(7,4) + Ifges(6,5)) * t90 + t109 * t35;
t1 = t35 * qJ(6) + t4;
t3 = t67 * t13 - t65 * t16;
t2 = -t35 * pkin(5) - t3;
t76 = t1 * t67 + t2 * t65;
t75 = -t3 * t65 + t4 * t67;
t74 = t65 * mrSges(6,1) + t67 * mrSges(6,2);
t40 = -t67 * mrSges(7,1) - t65 * mrSges(7,3);
t73 = t65 * mrSges(7,1) - t67 * mrSges(7,3);
t72 = t67 * pkin(5) + t65 * qJ(6);
t71 = pkin(5) * t65 - qJ(6) * t67;
t70 = -m(7) * t71 - t73 - t74;
t54 = Ifges(7,4) * t65;
t53 = Ifges(6,5) * t65;
t52 = Ifges(6,6) * t67;
t45 = Ifges(6,1) * t65 + t97;
t44 = Ifges(7,1) * t65 - t95;
t43 = Ifges(6,2) * t67 + t98;
t42 = -Ifges(7,3) * t67 + t96;
t39 = -pkin(4) - t72;
t34 = t37 ^ 2;
t31 = t37 * mrSges(5,2);
t18 = t74 * t37;
t17 = t73 * t37;
t9 = Ifges(6,5) * t35 + (Ifges(6,1) * t67 - t98) * t37;
t8 = Ifges(7,4) * t35 + (Ifges(7,1) * t67 + t96) * t37;
t7 = t94 + (-Ifges(6,2) * t65 + t97) * t37;
t6 = t93 + (Ifges(7,3) * t65 + t95) * t37;
t5 = t71 * t37 + t14;
t10 = [Ifges(4,2) * t58 + 0.2e1 * t50 * t79 + t31 * t103 + 0.2e1 * t5 * t17 + t18 * t104 + 0.2e1 * t4 * t19 + 0.2e1 * t3 * t20 + 0.2e1 * t2 * t21 + 0.2e1 * t1 * t22 + Ifges(2,3) + Ifges(3,3) + (Ifges(4,1) * t61 + 0.2e1 * Ifges(4,4) * t63) * t61 + (mrSges(5,1) * t103 - 0.2e1 * t16 * mrSges(5,3) + Ifges(5,2) * t35 + t77) * t35 + (mrSges(5,3) * t104 + Ifges(5,1) * t37 - 0.2e1 * Ifges(5,4) * t35 + (t8 + t9) * t67 + (t6 - t7 - t94) * t65) * t37 + m(4) * (t83 * t48 ^ 2 + t50 ^ 2) + m(5) * (t16 ^ 2 + t38 ^ 2 + t106) + m(6) * (t3 ^ 2 + t4 ^ 2 + t106) + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) + m(3) * (t62 ^ 2 + t64 ^ 2) * pkin(1) ^ 2 + 0.2e1 * (t64 * mrSges(3,1) - t62 * mrSges(3,2)) * pkin(1) + 0.2e1 * t83 * t48 * mrSges(4,3); (t17 + t18) * t35 + (-t87 * t65 + t88 * t67) * t37 + m(6) * (t75 * t37 + t92) + m(7) * (t5 * t35 + t76 * t37) + m(5) * (t16 * t37 + t92); m(3) + m(5) * (t34 + t105) + m(4) * t83 + (t82 * t34 + t105) * t108; t35 * mrSges(5,1) + t31 + t87 * t67 + t88 * t65 + m(7) * (t65 * t1 - t67 * t2) + m(6) * (t67 * t3 + t65 * t4) + m(5) * t38 + m(4) * t50 + t79; 0; t82 * t108 + m(4) + m(5); -t16 * mrSges(5,2) - pkin(4) * t18 + t39 * t17 + t5 * t40 + (t54 / 0.2e1 + t53 / 0.2e1 + t52 / 0.2e1 - Ifges(5,6)) * t35 + t85 * t14 + (-t93 / 0.2e1 - t6 / 0.2e1 + t7 / 0.2e1 + t4 * mrSges(6,3) + t1 * mrSges(7,2) + t88 * pkin(8)) * t67 + (t8 / 0.2e1 + t9 / 0.2e1 - t3 * mrSges(6,3) + t2 * mrSges(7,2) - t87 * pkin(8)) * t65 + m(6) * (-pkin(4) * t14 + t75 * pkin(8)) + m(7) * (t76 * pkin(8) + t39 * t5) + (Ifges(5,5) + (t44 / 0.2e1 + t45 / 0.2e1) * t67 + (t42 / 0.2e1 - t43 / 0.2e1) * t65) * t37; -t31 + (t40 + t85) * t35 + m(6) * (t86 - t102) + m(7) * (t39 * t35 + t86) + t37 * t107; 0; -0.2e1 * pkin(4) * t41 + 0.2e1 * t39 * t40 + Ifges(5,3) + (t43 - t42) * t67 + (t45 + t44) * t65 + m(7) * (t39 ^ 2 + t84) + m(6) * (pkin(4) ^ 2 + t84) + 0.2e1 * pkin(8) * t107; -Ifges(6,6) * t91 + qJ(6) * t22 + m(7) * (-pkin(5) * t2 + qJ(6) * t1) + t1 * mrSges(7,3) - t4 * mrSges(6,2) + t3 * mrSges(6,1) - t2 * mrSges(7,1) - pkin(5) * t21 + t77; t70 * t37; m(7) * t72 - t40 - t41; -t71 * mrSges(7,2) - Ifges(7,6) * t67 + t70 * pkin(8) + t52 + t53 + t54; 0.2e1 * pkin(5) * mrSges(7,1) + 0.2e1 * qJ(6) * mrSges(7,3) + m(7) * (pkin(5) ^ 2 + qJ(6) ^ 2) + t109; m(7) * t2 + t21; m(7) * t91; -m(7) * t67; (m(7) * pkin(8) + mrSges(7,2)) * t65; -m(7) * pkin(5) - mrSges(7,1); m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t10(1) t10(2) t10(4) t10(7) t10(11) t10(16); t10(2) t10(3) t10(5) t10(8) t10(12) t10(17); t10(4) t10(5) t10(6) t10(9) t10(13) t10(18); t10(7) t10(8) t10(9) t10(10) t10(14) t10(19); t10(11) t10(12) t10(13) t10(14) t10(15) t10(20); t10(16) t10(17) t10(18) t10(19) t10(20) t10(21);];
Mq  = res;
