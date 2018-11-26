% Calculate joint inertia matrix for
% S6RRPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
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
% Datum: 2018-11-23 17:00
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPR2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR2_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR2_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR2_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR2_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR2_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:00:25
% EndTime: 2018-11-23 17:00:25
% DurationCPUTime: 0.74s
% Computational Cost: add. (1702->204), mult. (3121->274), div. (0->0), fcn. (3471->8), ass. (0->84)
t119 = mrSges(6,2) - mrSges(5,1);
t72 = sin(qJ(6));
t75 = cos(qJ(6));
t99 = t72 ^ 2 + t75 ^ 2;
t94 = t99 * mrSges(7,3);
t57 = mrSges(7,1) * t72 + mrSges(7,2) * t75;
t118 = mrSges(6,3) + t57;
t110 = cos(qJ(4));
t70 = sin(pkin(10));
t112 = pkin(2) * t70;
t71 = cos(pkin(10));
t62 = pkin(2) * t71 + pkin(3);
t73 = sin(qJ(4));
t46 = t110 * t112 + t73 * t62;
t42 = qJ(5) + t46;
t117 = t42 ^ 2;
t74 = sin(qJ(2));
t76 = cos(qJ(2));
t49 = -t70 * t74 + t71 * t76;
t50 = t70 * t76 + t71 * t74;
t35 = -t110 * t49 + t50 * t73;
t36 = t110 * t50 + t73 * t49;
t63 = -pkin(2) * t76 - pkin(1);
t39 = -pkin(3) * t49 + t63;
t82 = -qJ(5) * t36 + t39;
t13 = pkin(4) * t35 + t82;
t116 = -0.2e1 * t13;
t115 = 0.2e1 * t39;
t114 = 0.2e1 * t49;
t113 = pkin(4) + pkin(9);
t55 = m(7) * t99;
t111 = m(6) + t55;
t109 = mrSges(7,1) * t75;
t108 = Ifges(7,4) * t72;
t107 = Ifges(7,4) * t75;
t106 = t35 * t72;
t105 = t35 * t75;
t45 = t110 * t62 - t73 * t112;
t104 = t45 * mrSges(5,1);
t103 = t46 * mrSges(5,2);
t102 = mrSges(6,1) + mrSges(5,3);
t101 = -qJ(3) - pkin(7);
t56 = t101 * t74;
t58 = t101 * t76;
t38 = t70 * t56 - t71 * t58;
t100 = t74 ^ 2 + t76 ^ 2;
t98 = qJ(5) * t42;
t96 = Ifges(7,5) * t106 + Ifges(7,6) * t105 + Ifges(7,3) * t36;
t25 = pkin(8) * t49 + t38;
t37 = t71 * t56 + t58 * t70;
t83 = -pkin(8) * t50 + t37;
t14 = -t110 * t83 + t25 * t73;
t16 = t110 * t25 + t73 * t83;
t95 = t14 ^ 2 + t16 ^ 2;
t93 = t99 * t113;
t92 = -t49 * mrSges(4,1) + t50 * mrSges(4,2);
t64 = Ifges(7,5) * t75;
t91 = -Ifges(7,6) * t72 + t64;
t90 = 0.2e1 * t118;
t4 = t113 * t35 + t82;
t5 = pkin(5) * t36 + t14;
t1 = -t4 * t72 + t5 * t75;
t2 = t4 * t75 + t5 * t72;
t89 = t1 * t75 + t2 * t72;
t88 = mrSges(6,2) - t94;
t87 = -t72 * mrSges(7,2) + t109;
t19 = mrSges(7,1) * t36 - mrSges(7,3) * t106;
t20 = -mrSges(7,2) * t36 + mrSges(7,3) * t105;
t86 = t75 * t19 + t72 * t20;
t85 = -0.2e1 * t94;
t59 = -Ifges(7,2) * t72 + t107;
t60 = Ifges(7,1) * t75 - t108;
t84 = -t72 * t59 + t75 * t60 + Ifges(6,1) + Ifges(5,3);
t44 = -pkin(4) - t45;
t10 = Ifges(7,5) * t36 + (Ifges(7,1) * t72 + t107) * t35;
t6 = -t35 * pkin(5) + t16;
t9 = Ifges(7,6) * t36 + (Ifges(7,2) * t75 + t108) * t35;
t81 = t60 * t106 / 0.2e1 + t59 * t105 / 0.2e1 + t6 * t57 - t72 * t9 / 0.2e1 + t75 * t10 / 0.2e1 + (-mrSges(5,2) + mrSges(6,3)) * t16 - t89 * mrSges(7,3) + (Ifges(6,5) - Ifges(5,6)) * t35 + t119 * t14 + (t91 / 0.2e1 + Ifges(5,5) - Ifges(6,4)) * t36;
t78 = qJ(5) ^ 2;
t41 = -pkin(9) + t44;
t29 = t36 * mrSges(6,3);
t28 = t36 * mrSges(5,2);
t18 = t87 * t35;
t3 = [-0.2e1 * pkin(1) * (-t76 * mrSges(3,1) + t74 * mrSges(3,2)) + t76 * (Ifges(3,4) * t74 + Ifges(3,2) * t76) + t74 * (Ifges(3,1) * t74 + Ifges(3,4) * t76) + 0.2e1 * t63 * t92 + Ifges(4,2) * t49 ^ 2 + t28 * t115 + t29 * t116 - 0.2e1 * t6 * t18 + 0.2e1 * t1 * t19 + 0.2e1 * t2 * t20 + t38 * mrSges(4,3) * t114 + Ifges(2,3) + 0.2e1 * t100 * pkin(7) * mrSges(3,3) + (-0.2e1 * t37 * mrSges(4,3) + Ifges(4,1) * t50 + Ifges(4,4) * t114) * t50 + (mrSges(5,1) * t115 + mrSges(6,2) * t116 + t72 * t10 + t75 * t9 + (Ifges(5,2) + Ifges(6,3)) * t35 - 0.2e1 * t102 * t16) * t35 + m(3) * (pkin(7) ^ 2 * t100 + pkin(1) ^ 2) + m(4) * (t37 ^ 2 + t38 ^ 2 + t63 ^ 2) + m(5) * (t39 ^ 2 + t95) + m(7) * (t1 ^ 2 + t2 ^ 2 + t6 ^ 2) + m(6) * (t13 ^ 2 + t95) + ((Ifges(5,1) + Ifges(6,2)) * t36 + t96 + 0.2e1 * (-Ifges(5,4) - Ifges(6,6)) * t35 + 0.2e1 * t102 * t14) * t36; m(7) * (t41 * t89 + t42 * t6) + t81 + t86 * t41 + m(6) * (t14 * t44 + t16 * t42) + m(5) * (-t14 * t45 + t16 * t46) + Ifges(3,5) * t74 + Ifges(3,6) * t76 + Ifges(4,6) * t49 + Ifges(4,5) * t50 + t37 * mrSges(4,1) - t38 * mrSges(4,2) - t42 * t18 + (-t42 * t35 + t44 * t36) * mrSges(6,1) + (-t46 * t35 - t45 * t36) * mrSges(5,3) + (-t74 * mrSges(3,1) - t76 * mrSges(3,2)) * pkin(7) + ((t49 * t70 - t50 * t71) * mrSges(4,3) + m(4) * (t37 * t71 + t38 * t70)) * pkin(2); 0.2e1 * t104 - 0.2e1 * t103 + 0.2e1 * t44 * mrSges(6,2) + Ifges(3,3) + Ifges(4,3) + t42 * t90 + t41 * t85 + m(7) * (t41 ^ 2 * t99 + t117) + m(6) * (t44 ^ 2 + t117) + m(5) * (t45 ^ 2 + t46 ^ 2) + t84 + (0.2e1 * mrSges(4,1) * t71 - 0.2e1 * mrSges(4,2) * t70 + m(4) * (t70 ^ 2 + t71 ^ 2) * pkin(2)) * pkin(2); -t72 * t19 + t75 * t20 + t28 - t29 - t119 * t35 + m(7) * (-t1 * t72 + t2 * t75) + m(6) * t13 + m(5) * t39 + m(4) * t63 + t92; 0; m(4) + m(5) + t111; m(7) * (qJ(5) * t6 - t113 * t89) - t86 * t113 + m(6) * (-pkin(4) * t14 + qJ(5) * t16) + (-pkin(4) * t36 - qJ(5) * t35) * mrSges(6,1) + t81 - qJ(5) * t18; t104 - t103 + (-pkin(4) + t44) * mrSges(6,2) + m(7) * (-t41 * t93 + t98) + m(6) * (-pkin(4) * t44 + t98) + t84 + (-t41 + t113) * t94 + t118 * (qJ(5) + t42); 0; -0.2e1 * pkin(4) * mrSges(6,2) + qJ(5) * t90 - t113 * t85 + m(7) * (t113 ^ 2 * t99 + t78) + m(6) * (pkin(4) ^ 2 + t78) + t84; m(6) * t14 + m(7) * t89 + t36 * mrSges(6,1) + t86; m(6) * t44 + t41 * t55 + t88; 0; -m(6) * pkin(4) - m(7) * t93 + t88; t111; mrSges(7,1) * t1 - mrSges(7,2) * t2 + t96; t41 * t87 + t91; -t57; -t113 * t109 + t64 + (mrSges(7,2) * t113 - Ifges(7,6)) * t72; t87; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t3(1) t3(2) t3(4) t3(7) t3(11) t3(16); t3(2) t3(3) t3(5) t3(8) t3(12) t3(17); t3(4) t3(5) t3(6) t3(9) t3(13) t3(18); t3(7) t3(8) t3(9) t3(10) t3(14) t3(19); t3(11) t3(12) t3(13) t3(14) t3(15) t3(20); t3(16) t3(17) t3(18) t3(19) t3(20) t3(21);];
Mq  = res;
