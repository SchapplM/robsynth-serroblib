% Calculate joint inertia matrix for
% S6RRPPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
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
% Datum: 2018-11-23 16:45
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRP2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP2_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP2_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP2_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP2_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRP2_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:45:36
% EndTime: 2018-11-23 16:45:37
% DurationCPUTime: 0.74s
% Computational Cost: add. (1019->221), mult. (1848->285), div. (0->0), fcn. (1845->6), ass. (0->85)
t108 = Ifges(6,6) + Ifges(7,6);
t107 = Ifges(6,3) + Ifges(7,3);
t68 = sin(pkin(9));
t101 = pkin(2) * t68;
t55 = qJ(4) + t101;
t106 = t55 ^ 2;
t105 = -2 * mrSges(7,3);
t69 = cos(pkin(9));
t71 = sin(qJ(2));
t73 = cos(qJ(2));
t38 = t68 * t71 - t69 * t73;
t39 = t68 * t73 + t69 * t71;
t59 = -pkin(2) * t73 - pkin(1);
t76 = -qJ(4) * t39 + t59;
t18 = pkin(3) * t38 + t76;
t104 = -0.2e1 * t18;
t103 = 0.2e1 * t59;
t102 = m(7) * pkin(5);
t70 = sin(qJ(5));
t72 = cos(qJ(5));
t53 = t70 ^ 2 + t72 ^ 2;
t42 = m(6) * t53;
t100 = pkin(2) * t69;
t11 = (pkin(3) + pkin(8)) * t38 + t76;
t88 = -qJ(3) - pkin(7);
t44 = t88 * t71;
t47 = t88 * t73;
t23 = -t69 * t44 - t47 * t68;
t14 = pkin(4) * t39 + t23;
t4 = t72 * t11 + t70 * t14;
t99 = mrSges(6,1) * t72;
t98 = mrSges(6,2) * t72;
t97 = Ifges(6,4) * t70;
t96 = Ifges(6,4) * t72;
t95 = Ifges(7,4) * t70;
t94 = Ifges(7,4) * t72;
t93 = t38 * t70;
t92 = t38 * t72;
t91 = t72 * mrSges(7,3);
t90 = mrSges(5,1) + mrSges(4,3);
t89 = mrSges(5,2) - mrSges(4,1);
t19 = mrSges(7,1) * t39 - mrSges(7,3) * t93;
t20 = mrSges(6,1) * t39 - mrSges(6,3) * t93;
t87 = t19 + t20;
t21 = -mrSges(7,2) * t39 + t38 * t91;
t22 = -mrSges(6,2) * t39 + mrSges(6,3) * t92;
t86 = t21 + t22;
t45 = t70 * mrSges(7,1) + t72 * mrSges(7,2);
t85 = t71 ^ 2 + t73 ^ 2;
t84 = qJ(6) * t38;
t58 = -pkin(3) - t100;
t54 = -pkin(8) + t58;
t83 = -qJ(6) + t54;
t82 = m(7) * t53 + m(5) + t42;
t25 = t44 * t68 - t47 * t69;
t81 = t23 ^ 2 + t25 ^ 2;
t80 = -mrSges(6,1) - t102;
t79 = t53 * mrSges(6,3);
t16 = -mrSges(7,1) * t92 + mrSges(7,2) * t93;
t13 = t72 * t14;
t3 = -t11 * t70 + t13;
t78 = t3 * t72 + t4 * t70;
t77 = (Ifges(6,5) + Ifges(7,5)) * t93 + t108 * t92 + t107 * t39;
t63 = Ifges(6,5) * t72;
t62 = Ifges(7,5) * t72;
t51 = Ifges(6,1) * t72 - t97;
t50 = Ifges(7,1) * t72 - t95;
t49 = -Ifges(6,2) * t70 + t96;
t48 = -Ifges(7,2) * t70 + t94;
t46 = mrSges(6,1) * t70 + t98;
t43 = pkin(5) * t70 + t55;
t35 = t39 * mrSges(5,3);
t34 = t39 * mrSges(4,2);
t33 = t83 * t72;
t32 = t83 * t70;
t17 = (mrSges(6,2) * t70 - t99) * t38;
t15 = -pkin(4) * t38 + t25;
t10 = Ifges(6,5) * t39 + (Ifges(6,1) * t70 + t96) * t38;
t9 = Ifges(7,5) * t39 + (Ifges(7,1) * t70 + t94) * t38;
t8 = Ifges(6,6) * t39 + (Ifges(6,2) * t72 + t97) * t38;
t7 = Ifges(7,6) * t39 + (Ifges(7,2) * t72 + t95) * t38;
t5 = (-pkin(5) * t72 - pkin(4)) * t38 + t25;
t2 = t72 * t84 + t4;
t1 = pkin(5) * t39 + t13 + (-t11 - t84) * t70;
t6 = [t34 * t103 + t35 * t104 - 0.2e1 * pkin(1) * (-t73 * mrSges(3,1) + t71 * mrSges(3,2)) + t71 * (Ifges(3,1) * t71 + Ifges(3,4) * t73) + t73 * (Ifges(3,4) * t71 + Ifges(3,2) * t73) + 0.2e1 * t5 * t16 + 0.2e1 * t15 * t17 + 0.2e1 * t1 * t19 + 0.2e1 * t3 * t20 + 0.2e1 * t2 * t21 + 0.2e1 * t4 * t22 + Ifges(2,3) + 0.2e1 * t85 * pkin(7) * mrSges(3,3) + m(3) * (t85 * pkin(7) ^ 2 + pkin(1) ^ 2) + m(4) * (t59 ^ 2 + t81) + m(6) * (t15 ^ 2 + t3 ^ 2 + t4 ^ 2) + m(5) * (t18 ^ 2 + t81) + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) + ((Ifges(5,2) + Ifges(4,1)) * t39 + 0.2e1 * t90 * t23 + t77) * t39 + (mrSges(4,1) * t103 + mrSges(5,2) * t104 + (Ifges(4,2) + Ifges(5,3)) * t38 + (t7 + t8) * t72 + (t10 + t9) * t70 + 0.2e1 * (-Ifges(4,4) - Ifges(5,6)) * t39 - 0.2e1 * t90 * t25) * t38; Ifges(3,5) * t71 + Ifges(3,6) * t73 + t15 * t46 + t43 * t16 + t55 * t17 + t33 * t19 + t32 * t21 + t5 * t45 + (-mrSges(4,2) + mrSges(5,3)) * t25 + t89 * t23 + (-mrSges(3,1) * t71 - mrSges(3,2) * t73) * pkin(7) + (t54 * t22 - t2 * mrSges(7,3) - t4 * mrSges(6,3) - t7 / 0.2e1 - t8 / 0.2e1) * t70 + (t58 * mrSges(5,1) + t62 / 0.2e1 + t63 / 0.2e1 - Ifges(5,4) + Ifges(4,5) - mrSges(4,3) * t100 + (-Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1) * t70) * t39 + (t54 * t20 - t1 * mrSges(7,3) - t3 * mrSges(6,3) + t10 / 0.2e1 + t9 / 0.2e1) * t72 + m(6) * (t15 * t55 + t78 * t54) + m(7) * (t1 * t33 + t2 * t32 + t43 * t5) + m(5) * (t23 * t58 + t25 * t55) + m(4) * (-t23 * t69 + t25 * t68) * pkin(2) + (-mrSges(4,3) * t101 - t55 * mrSges(5,1) + Ifges(5,5) - Ifges(4,6) + (t48 / 0.2e1 + t49 / 0.2e1) * t72 + (t50 / 0.2e1 + t51 / 0.2e1) * t70) * t38; 0.2e1 * t58 * mrSges(5,2) + 0.2e1 * t43 * t45 + Ifges(5,1) + Ifges(3,3) + Ifges(4,3) + (t33 * t105 + t50 + t51) * t72 + (t32 * t105 - t48 - t49) * t70 + m(7) * (t32 ^ 2 + t33 ^ 2 + t43 ^ 2) + m(6) * (t53 * t54 ^ 2 + t106) + m(5) * (t58 ^ 2 + t106) + m(4) * (t68 ^ 2 + t69 ^ 2) * pkin(2) ^ 2 + 0.2e1 * (t46 + mrSges(5,3)) * t55 + 0.2e1 * (mrSges(4,1) * t69 - mrSges(4,2) * t68) * pkin(2) - 0.2e1 * t54 * t79; t34 - t35 + t86 * t72 - t87 * t70 - t89 * t38 + m(6) * (-t3 * t70 + t4 * t72) + m(7) * (-t1 * t70 + t2 * t72) + m(5) * t18 + m(4) * t59; m(7) * (t32 * t72 - t33 * t70); m(4) + t82; t39 * mrSges(5,1) + t87 * t72 + t86 * t70 + m(6) * t78 + m(7) * (t1 * t72 + t2 * t70) + m(5) * t23; mrSges(5,2) - t53 * mrSges(7,3) - t79 + m(7) * (t32 * t70 + t33 * t72) + m(5) * t58 + t54 * t42; 0; t82; mrSges(6,1) * t3 + mrSges(7,1) * t1 - mrSges(6,2) * t4 - mrSges(7,2) * t2 + (m(7) * t1 + t19) * pkin(5) + t77; t54 * t99 + mrSges(7,1) * t33 - mrSges(7,2) * t32 + t62 + t63 + (m(7) * t33 - t91) * pkin(5) + (-mrSges(6,2) * t54 - t108) * t70; t80 * t70 - t45 - t98; (-mrSges(6,2) - mrSges(7,2)) * t70 + (mrSges(7,1) - t80) * t72; (0.2e1 * mrSges(7,1) + t102) * pkin(5) + t107; m(7) * t5 + t16; m(7) * t43 + t45; 0; 0; 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t6(1) t6(2) t6(4) t6(7) t6(11) t6(16); t6(2) t6(3) t6(5) t6(8) t6(12) t6(17); t6(4) t6(5) t6(6) t6(9) t6(13) t6(18); t6(7) t6(8) t6(9) t6(10) t6(14) t6(19); t6(11) t6(12) t6(13) t6(14) t6(15) t6(20); t6(16) t6(17) t6(18) t6(19) t6(20) t6(21);];
Mq  = res;
