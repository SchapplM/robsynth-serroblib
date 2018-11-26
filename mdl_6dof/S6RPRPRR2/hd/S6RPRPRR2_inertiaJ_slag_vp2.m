% Calculate joint inertia matrix for
% S6RPRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2018-11-23 16:03
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRR2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR2_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR2_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR2_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR2_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR2_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:03:19
% EndTime: 2018-11-23 16:03:20
% DurationCPUTime: 0.84s
% Computational Cost: add. (1553->232), mult. (2875->347), div. (0->0), fcn. (3054->10), ass. (0->85)
t73 = sin(pkin(11));
t75 = cos(pkin(11));
t79 = sin(qJ(3));
t82 = cos(qJ(3));
t52 = t73 * t82 + t75 * t79;
t77 = sin(qJ(6));
t78 = sin(qJ(5));
t80 = cos(qJ(6));
t81 = cos(qJ(5));
t54 = t77 * t81 + t78 * t80;
t23 = t54 * t52;
t53 = -t77 * t78 + t80 * t81;
t24 = t53 * t52;
t10 = t23 * mrSges(7,1) + t24 * mrSges(7,2);
t89 = mrSges(6,1) * t78 + mrSges(6,2) * t81;
t29 = t89 * t52;
t114 = t10 + t29;
t113 = t82 ^ 2;
t101 = t52 * t81;
t50 = t73 * t79 - t75 * t82;
t112 = Ifges(6,5) * t101 + Ifges(6,3) * t50;
t74 = sin(pkin(10));
t64 = pkin(1) * t74 + pkin(7);
t95 = qJ(4) + t64;
t42 = t95 * t82;
t92 = t95 * t79;
t26 = t42 * t73 + t75 * t92;
t111 = t26 ^ 2;
t48 = t50 ^ 2;
t110 = 0.2e1 * t26;
t76 = cos(pkin(10));
t66 = -pkin(1) * t76 - pkin(2);
t56 = -pkin(3) * t82 + t66;
t109 = 0.2e1 * t56;
t108 = m(7) * pkin(5);
t63 = pkin(3) * t73 + pkin(8);
t106 = pkin(9) + t63;
t105 = Ifges(6,4) * t78;
t104 = Ifges(6,4) * t81;
t103 = t26 * t50;
t102 = t52 * t78;
t22 = pkin(4) * t50 - pkin(8) * t52 + t56;
t28 = t75 * t42 - t73 * t92;
t9 = t78 * t22 + t81 * t28;
t100 = Ifges(7,5) * t54 + Ifges(7,6) * t53;
t57 = -t81 * mrSges(6,1) + t78 * mrSges(6,2);
t99 = t57 - mrSges(5,1);
t98 = Ifges(6,5) * t78 + Ifges(6,6) * t81;
t97 = t78 ^ 2 + t81 ^ 2;
t96 = t79 ^ 2 + t113;
t94 = Ifges(7,5) * t24 - Ifges(7,6) * t23 + Ifges(7,3) * t50;
t65 = -pkin(3) * t75 - pkin(4);
t93 = t97 * t63;
t34 = -t53 * mrSges(7,1) + t54 * mrSges(7,2);
t8 = t81 * t22 - t28 * t78;
t91 = -t8 * t78 + t9 * t81;
t90 = -t82 * mrSges(4,1) + t79 * mrSges(4,2);
t43 = t106 * t78;
t44 = t106 * t81;
t30 = -t43 * t80 - t44 * t77;
t31 = -t43 * t77 + t44 * t80;
t88 = t30 * mrSges(7,1) - t31 * mrSges(7,2) + t100;
t4 = pkin(5) * t50 - pkin(9) * t101 + t8;
t5 = -pkin(9) * t102 + t9;
t2 = t4 * t80 - t5 * t77;
t3 = t4 * t77 + t5 * t80;
t87 = t2 * mrSges(7,1) - t3 * mrSges(7,2) + t94;
t86 = (mrSges(7,1) * t80 - mrSges(7,2) * t77) * pkin(5);
t59 = Ifges(6,1) * t78 + t104;
t58 = Ifges(6,2) * t81 + t105;
t55 = -pkin(5) * t81 + t65;
t49 = t52 ^ 2;
t39 = t52 * mrSges(5,2);
t36 = Ifges(7,1) * t54 + Ifges(7,4) * t53;
t35 = Ifges(7,4) * t54 + Ifges(7,2) * t53;
t33 = mrSges(6,1) * t50 - mrSges(6,3) * t101;
t32 = -mrSges(6,2) * t50 - mrSges(6,3) * t102;
t15 = Ifges(6,5) * t50 + (Ifges(6,1) * t81 - t105) * t52;
t14 = Ifges(6,6) * t50 + (-Ifges(6,2) * t78 + t104) * t52;
t13 = pkin(5) * t102 + t26;
t12 = mrSges(7,1) * t50 - mrSges(7,3) * t24;
t11 = -mrSges(7,2) * t50 - mrSges(7,3) * t23;
t7 = Ifges(7,1) * t24 - Ifges(7,4) * t23 + Ifges(7,5) * t50;
t6 = Ifges(7,4) * t24 - Ifges(7,2) * t23 + Ifges(7,6) * t50;
t1 = [Ifges(4,2) * t113 + 0.2e1 * t66 * t90 + t39 * t109 + 0.2e1 * t9 * t32 + 0.2e1 * t8 * t33 - t23 * t6 + t24 * t7 + t29 * t110 + 0.2e1 * t3 * t11 + 0.2e1 * t2 * t12 + 0.2e1 * t13 * t10 + Ifges(2,3) + Ifges(3,3) + (mrSges(5,1) * t109 - 0.2e1 * mrSges(5,3) * t28 + Ifges(5,2) * t50 + t112 + t94) * t50 + (mrSges(5,3) * t110 + Ifges(5,1) * t52 - t78 * t14 + t81 * t15 + (-Ifges(6,6) * t78 - (2 * Ifges(5,4))) * t50) * t52 + m(4) * (t64 ^ 2 * t96 + t66 ^ 2) + m(5) * (t28 ^ 2 + t56 ^ 2 + t111) + m(6) * (t8 ^ 2 + t9 ^ 2 + t111) + m(7) * (t13 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(3) * (t74 ^ 2 + t76 ^ 2) * pkin(1) ^ 2 + (Ifges(4,1) * t79 + 0.2e1 * Ifges(4,4) * t82) * t79 + 0.2e1 * (mrSges(3,1) * t76 - mrSges(3,2) * t74) * pkin(1) + 0.2e1 * t96 * t64 * mrSges(4,3); t24 * t11 - t23 * t12 + (t81 * t32 - t78 * t33) * t52 + t114 * t50 + m(6) * (t52 * t91 + t103) + m(7) * (t13 * t50 - t2 * t23 + t24 * t3) + m(5) * (t28 * t52 + t103); m(3) + m(6) * (t49 * t97 + t48) + m(7) * (t23 ^ 2 + t24 ^ 2 + t48) + m(5) * (t49 + t48) + m(4) * t96; Ifges(4,6) * t82 + Ifges(4,5) * t79 + t65 * t29 + t53 * t6 / 0.2e1 + t54 * t7 / 0.2e1 + t55 * t10 - Ifges(5,6) * t50 + t31 * t11 + t13 * t34 - t23 * t35 / 0.2e1 + t24 * t36 / 0.2e1 - t28 * mrSges(5,2) + t30 * t12 + (-t79 * mrSges(4,1) - t82 * mrSges(4,2)) * t64 + t99 * t26 + (-t2 * t54 + t3 * t53) * mrSges(7,3) + (t63 * t32 + t9 * mrSges(6,3) + t14 / 0.2e1) * t81 + (-t63 * t33 - t8 * mrSges(6,3) + t15 / 0.2e1) * t78 + m(6) * (t26 * t65 + t63 * t91) + m(7) * (t13 * t55 + t2 * t30 + t3 * t31) + (t81 * t59 / 0.2e1 - t78 * t58 / 0.2e1 + Ifges(5,5)) * t52 + (m(5) * (-t26 * t75 + t28 * t73) + (-t73 * t50 - t75 * t52) * mrSges(5,3)) * pkin(3) + (t98 + t100) * t50 / 0.2e1; -t39 + (t23 * t54 + t24 * t53) * mrSges(7,3) + t97 * t52 * mrSges(6,3) + (t34 + t99) * t50 + m(6) * (t50 * t65 + t52 * t93) + m(7) * (-t23 * t30 + t24 * t31 + t50 * t55) + m(5) * (-t50 * t75 + t52 * t73) * pkin(3) - t90; 0.2e1 * t55 * t34 + t53 * t35 + t54 * t36 + 0.2e1 * t65 * t57 + t81 * t58 + t78 * t59 + Ifges(4,3) + Ifges(5,3) + m(7) * (t30 ^ 2 + t31 ^ 2 + t55 ^ 2) + m(6) * (t63 ^ 2 * t97 + t65 ^ 2) + m(5) * (t73 ^ 2 + t75 ^ 2) * pkin(3) ^ 2 + 0.2e1 * (mrSges(5,1) * t75 - mrSges(5,2) * t73) * pkin(3) + 0.2e1 * (-t30 * t54 + t31 * t53) * mrSges(7,3) + 0.2e1 * mrSges(6,3) * t93; t50 * mrSges(5,1) + t54 * t11 + t53 * t12 + t78 * t32 + t81 * t33 + t39 + m(7) * (t2 * t53 + t3 * t54) + m(6) * (t78 * t9 + t8 * t81) + m(5) * t56; m(7) * (-t23 * t53 + t24 * t54); m(7) * (t30 * t53 + t31 * t54); m(5) + m(6) * t97 + m(7) * (t53 ^ 2 + t54 ^ 2); -Ifges(6,6) * t102 + t8 * mrSges(6,1) - t9 * mrSges(6,2) + (m(7) * (t2 * t80 + t3 * t77) + t77 * t11 + t80 * t12) * pkin(5) + t87 + t112; (-t23 * t80 + t24 * t77) * t108 - t114; -t89 * t63 + (m(7) * (t30 * t80 + t31 * t77) + (t53 * t77 - t54 * t80) * mrSges(7,3)) * pkin(5) + t88 + t98; (t53 * t80 + t54 * t77) * t108 - t57 - t34; Ifges(6,3) + Ifges(7,3) + m(7) * (t77 ^ 2 + t80 ^ 2) * pkin(5) ^ 2 + 0.2e1 * t86; t87; -t10; t88; -t34; Ifges(7,3) + t86; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
