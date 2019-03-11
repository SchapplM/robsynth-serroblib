% Calculate joint inertia matrix for
% S6RPPRRP7
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
% Datum: 2019-03-09 02:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRP7_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP7_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP7_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP7_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP7_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP7_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:12:42
% EndTime: 2019-03-09 02:12:44
% DurationCPUTime: 0.75s
% Computational Cost: add. (906->197), mult. (1612->262), div. (0->0), fcn. (1593->6), ass. (0->80)
t111 = Ifges(6,5) + Ifges(7,5);
t89 = Ifges(6,6) + Ifges(7,6);
t110 = Ifges(6,3) + Ifges(7,3);
t109 = -2 * mrSges(7,3);
t108 = m(6) + m(7);
t66 = sin(qJ(5));
t68 = cos(qJ(5));
t107 = t111 * t66 + t89 * t68;
t65 = -pkin(1) - qJ(3);
t101 = -pkin(7) + t65;
t63 = sin(pkin(9));
t37 = t101 * t63;
t67 = sin(qJ(4));
t100 = cos(qJ(4));
t64 = cos(pkin(9));
t76 = t100 * t64;
t21 = -t101 * t76 + t67 * t37;
t106 = t21 ^ 2;
t34 = t67 * t63 - t76;
t33 = t34 ^ 2;
t59 = t64 ^ 2;
t105 = -2 * mrSges(5,3);
t103 = m(7) * pkin(5);
t99 = Ifges(6,4) * t66;
t98 = Ifges(6,4) * t68;
t97 = Ifges(7,4) * t66;
t96 = Ifges(7,4) * t68;
t95 = t34 * t21;
t94 = t34 * t66;
t93 = t34 * t68;
t92 = t66 * mrSges(6,2);
t91 = t66 * mrSges(7,3);
t90 = t67 * t64;
t88 = -qJ(6) - pkin(8);
t35 = t100 * t63 + t90;
t47 = t63 * pkin(3) + qJ(2);
t13 = t35 * pkin(4) + t34 * pkin(8) + t47;
t23 = t100 * t37 + t101 * t90;
t4 = t66 * t13 + t68 * t23;
t17 = -t35 * mrSges(7,2) + t34 * t91;
t18 = -t35 * mrSges(6,2) + mrSges(6,3) * t94;
t87 = t17 + t18;
t19 = t35 * mrSges(7,1) + mrSges(7,3) * t93;
t20 = t35 * mrSges(6,1) + mrSges(6,3) * t93;
t86 = t19 + t20;
t14 = -mrSges(7,1) * t94 - mrSges(7,2) * t93;
t40 = -t68 * mrSges(6,1) + t92;
t85 = t40 - mrSges(5,1);
t84 = t63 * mrSges(4,1) + t64 * mrSges(4,2);
t81 = t63 ^ 2 + t59;
t80 = t66 ^ 2 + t68 ^ 2;
t79 = qJ(6) * t34;
t77 = m(4) * t81;
t75 = t81 * mrSges(4,3);
t74 = t80 * mrSges(6,3);
t53 = t66 * mrSges(7,2);
t39 = -t68 * mrSges(7,1) + t53;
t3 = t68 * t13 - t66 * t23;
t73 = t110 * t35 - t111 * t93;
t72 = mrSges(6,1) + mrSges(7,1) + t103;
t71 = -mrSges(6,1) * t66 - mrSges(6,2) * t68;
t69 = qJ(2) ^ 2;
t49 = -t68 * pkin(5) - pkin(4);
t45 = Ifges(6,1) * t66 + t98;
t44 = Ifges(7,1) * t66 + t96;
t43 = Ifges(6,2) * t68 + t99;
t42 = Ifges(7,2) * t68 + t97;
t41 = t88 * t68;
t38 = t88 * t66;
t32 = t35 ^ 2;
t29 = t34 * mrSges(5,2);
t15 = t71 * t34;
t9 = -pkin(5) * t94 + t21;
t8 = Ifges(6,5) * t35 + (-Ifges(6,1) * t68 + t99) * t34;
t7 = Ifges(7,5) * t35 + (-Ifges(7,1) * t68 + t97) * t34;
t6 = Ifges(6,6) * t35 + (Ifges(6,2) * t66 - t98) * t34;
t5 = Ifges(7,6) * t35 + (Ifges(7,2) * t66 - t96) * t34;
t2 = t66 * t79 + t4;
t1 = t35 * pkin(5) + t68 * t79 + t3;
t10 = [Ifges(4,1) * t59 - 0.2e1 * t47 * t29 + 0.2e1 * t9 * t14 + 0.2e1 * t2 * t17 + 0.2e1 * t4 * t18 + 0.2e1 * t1 * t19 + 0.2e1 * t3 * t20 + 0.2e1 * t21 * t15 - (2 * pkin(1) * mrSges(3,2)) + Ifges(3,1) + Ifges(2,3) + (-0.2e1 * Ifges(4,4) * t64 + Ifges(4,2) * t63) * t63 - 0.2e1 * t65 * t75 + (0.2e1 * t47 * mrSges(5,1) + Ifges(5,2) * t35 + t23 * t105 + t73) * t35 + (t21 * t105 + Ifges(5,1) * t34 + 0.2e1 * Ifges(5,4) * t35 + (-t7 - t8) * t68 + (t89 * t35 + t5 + t6) * t66) * t34 + m(4) * (t81 * t65 ^ 2 + t69) + m(3) * ((pkin(1) ^ 2) + t69) + m(5) * (t23 ^ 2 + t47 ^ 2 + t106) + m(6) * (t3 ^ 2 + t4 ^ 2 + t106) + m(7) * (t1 ^ 2 + t2 ^ 2 + t9 ^ 2) + 0.2e1 * (t84 + mrSges(3,3)) * qJ(2); -m(3) * pkin(1) - t33 * mrSges(5,3) + mrSges(3,2) + (t14 + t15) * t34 - t75 + (-mrSges(5,3) * t35 - t86 * t66 + t87 * t68) * t35 + m(6) * (t95 + (-t3 * t66 + t4 * t68) * t35) + m(7) * (t34 * t9 + (-t1 * t66 + t2 * t68) * t35) + m(5) * (t35 * t23 + t95) + t65 * t77; m(3) + m(5) * (t32 + t33) + t77 + (t80 * t32 + t33) * t108; m(4) * qJ(2) + t35 * mrSges(5,1) - t29 + t86 * t68 + t87 * t66 + m(6) * (t68 * t3 + t66 * t4) + m(7) * (t68 * t1 + t66 * t2) + m(5) * t47 + t84; 0; t80 * t108 + m(4) + m(5); m(7) * (t38 * t1 - t41 * t2 + t49 * t9) + t38 * t19 + t9 * t39 - t41 * t17 + t49 * t14 - t23 * mrSges(5,2) - pkin(4) * t15 + (-m(6) * pkin(4) + t85) * t21 + (t4 * mrSges(6,3) + t2 * mrSges(7,3) + t5 / 0.2e1 + t6 / 0.2e1 + (m(6) * t4 + t18) * pkin(8)) * t68 + (-t1 * mrSges(7,3) - t3 * mrSges(6,3) + t7 / 0.2e1 + t8 / 0.2e1 + (-m(6) * t3 - t20) * pkin(8)) * t66 + (-Ifges(5,5) + (-t44 / 0.2e1 - t45 / 0.2e1) * t68 + (t42 / 0.2e1 + t43 / 0.2e1) * t66) * t34 + (-Ifges(5,6) + t107 / 0.2e1) * t35; (t39 + t85) * t34 + (t80 * mrSges(7,3) - mrSges(5,2) + t74) * t35 + m(7) * (t49 * t34 + (-t38 * t66 - t41 * t68) * t35) + m(6) * (t80 * t35 * pkin(8) - pkin(4) * t34); m(7) * (t38 * t68 - t41 * t66); -0.2e1 * pkin(4) * t40 + 0.2e1 * t49 * t39 + Ifges(5,3) + 0.2e1 * pkin(8) * t74 + m(7) * (t38 ^ 2 + t41 ^ 2 + t49 ^ 2) + m(6) * (t80 * pkin(8) ^ 2 + pkin(4) ^ 2) + (t41 * t109 + t42 + t43) * t68 + (t38 * t109 + t44 + t45) * t66; t3 * mrSges(6,1) + t1 * mrSges(7,1) - t4 * mrSges(6,2) - t2 * mrSges(7,2) + t89 * t94 + (m(7) * t1 + t19) * pkin(5) + t73; ((-mrSges(6,2) - mrSges(7,2)) * t68 - t72 * t66) * t35; t72 * t68 - t53 - t92; t38 * mrSges(7,1) + t41 * mrSges(7,2) + t71 * pkin(8) + (m(7) * t38 - t91) * pkin(5) + t107; (0.2e1 * mrSges(7,1) + t103) * pkin(5) + t110; m(7) * t9 + t14; m(7) * t34; 0; m(7) * t49 + t39; 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t10(1) t10(2) t10(4) t10(7) t10(11) t10(16); t10(2) t10(3) t10(5) t10(8) t10(12) t10(17); t10(4) t10(5) t10(6) t10(9) t10(13) t10(18); t10(7) t10(8) t10(9) t10(10) t10(14) t10(19); t10(11) t10(12) t10(13) t10(14) t10(15) t10(20); t10(16) t10(17) t10(18) t10(19) t10(20) t10(21);];
Mq  = res;
