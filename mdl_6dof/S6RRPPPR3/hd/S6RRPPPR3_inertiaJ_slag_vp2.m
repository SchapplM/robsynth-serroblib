% Calculate joint inertia matrix for
% S6RRPPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta5]';
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
% Datum: 2018-11-23 16:43
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRPPPR3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR3_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR3_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR3_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR3_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR3_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:43:25
% EndTime: 2018-11-23 16:43:26
% DurationCPUTime: 0.67s
% Computational Cost: add. (841->210), mult. (1396->274), div. (0->0), fcn. (1204->6), ass. (0->74)
t73 = sin(qJ(2));
t74 = cos(qJ(2));
t102 = t73 ^ 2 + t74 ^ 2;
t69 = sin(pkin(9));
t70 = cos(pkin(9));
t72 = sin(qJ(6));
t93 = cos(qJ(6));
t35 = t72 * t69 - t93 * t70;
t36 = t93 * t69 + t72 * t70;
t12 = -t35 * mrSges(7,1) - t36 * mrSges(7,2);
t101 = t35 ^ 2 + t36 ^ 2;
t85 = t69 ^ 2 + t70 ^ 2;
t40 = m(6) * t85;
t100 = m(7) * t101 + m(5) + t40;
t56 = t74 * qJ(4);
t50 = t74 * pkin(7) - t56;
t99 = t50 ^ 2;
t98 = -2 * mrSges(5,3);
t45 = -t74 * pkin(2) - t73 * qJ(3) - pkin(1);
t97 = -0.2e1 * t45;
t96 = -t70 / 0.2e1;
t75 = -pkin(2) - pkin(3);
t64 = -qJ(5) + t75;
t94 = pkin(8) - t64;
t92 = Ifges(6,4) * t69;
t91 = Ifges(6,4) * t70;
t89 = t69 * t74;
t88 = t70 * t74;
t87 = mrSges(4,2) - mrSges(5,3);
t71 = qJ(3) + pkin(4);
t32 = t74 * pkin(3) - t45;
t19 = t73 * pkin(4) + t74 * qJ(5) + t32;
t49 = (pkin(7) - qJ(4)) * t73;
t9 = t69 * t19 + t70 * t49;
t86 = t102 * pkin(7) ^ 2;
t25 = t36 * t74;
t26 = t35 * t74;
t84 = Ifges(7,5) * t26 + Ifges(7,6) * t25 + Ifges(7,3) * t73;
t83 = -m(4) * pkin(2) - mrSges(4,1);
t82 = t85 * mrSges(6,3);
t7 = -t25 * mrSges(7,1) + t26 * mrSges(7,2);
t46 = t70 * mrSges(6,1) - t69 * mrSges(6,2);
t8 = t70 * t19 - t69 * t49;
t3 = t73 * pkin(5) + pkin(8) * t88 + t8;
t4 = pkin(8) * t89 + t9;
t1 = t93 * t3 - t72 * t4;
t2 = t72 * t3 + t93 * t4;
t80 = t36 * t1 + t35 * t2;
t79 = -t69 * t8 + t70 * t9;
t41 = t94 * t69;
t42 = t94 * t70;
t10 = t93 * t41 + t72 * t42;
t11 = t72 * t41 - t93 * t42;
t78 = t36 * t10 + t35 * t11;
t28 = (-t69 * mrSges(6,1) - t70 * mrSges(6,2)) * t74;
t76 = qJ(3) ^ 2;
t57 = t73 * mrSges(5,1);
t52 = t70 * pkin(5) + t71;
t48 = -Ifges(6,1) * t69 - t91;
t47 = -Ifges(6,2) * t70 - t92;
t44 = t73 * mrSges(6,1) + mrSges(6,3) * t88;
t43 = -t73 * mrSges(6,2) + mrSges(6,3) * t89;
t31 = Ifges(7,5) * t36;
t30 = Ifges(7,6) * t35;
t27 = -t56 + (-pkin(5) * t69 + pkin(7)) * t74;
t24 = Ifges(6,5) * t73 + (-Ifges(6,1) * t70 + t92) * t74;
t23 = Ifges(6,6) * t73 + (Ifges(6,2) * t69 - t91) * t74;
t18 = t73 * mrSges(7,1) - t26 * mrSges(7,3);
t17 = -t73 * mrSges(7,2) + t25 * mrSges(7,3);
t14 = -Ifges(7,1) * t36 + Ifges(7,4) * t35;
t13 = -Ifges(7,4) * t36 + Ifges(7,2) * t35;
t6 = Ifges(7,1) * t26 + Ifges(7,4) * t25 + Ifges(7,5) * t73;
t5 = Ifges(7,4) * t26 + Ifges(7,2) * t25 + Ifges(7,6) * t73;
t15 = [0.2e1 * t1 * t18 + 0.2e1 * t2 * t17 + t25 * t5 + t26 * t6 + 0.2e1 * t27 * t7 + 0.2e1 * t50 * t28 + 0.2e1 * t32 * t57 + 0.2e1 * t9 * t43 + 0.2e1 * t8 * t44 + Ifges(2,3) + m(6) * (t8 ^ 2 + t9 ^ 2 + t99) + m(5) * (t32 ^ 2 + t49 ^ 2 + t99) + m(4) * (t45 ^ 2 + t86) + m(3) * (pkin(1) ^ 2 + t86) + m(7) * (t1 ^ 2 + t2 ^ 2 + t27 ^ 2) + (0.2e1 * pkin(1) * mrSges(3,1) + mrSges(4,1) * t97 - 0.2e1 * t32 * mrSges(5,2) + t50 * t98 + t69 * t23 - t70 * t24 + (Ifges(4,3) + Ifges(3,2) + Ifges(5,1)) * t74) * t74 + (-0.2e1 * pkin(1) * mrSges(3,2) + mrSges(4,3) * t97 + t49 * t98 + (Ifges(3,1) + Ifges(6,3) + Ifges(5,2) + Ifges(4,1)) * t73 + (-Ifges(6,5) * t70 + Ifges(6,6) * t69 + (2 * Ifges(3,4)) + (2 * Ifges(5,4)) - (2 * Ifges(4,5))) * t74 + t84) * t73 + 0.2e1 * (mrSges(4,2) + mrSges(3,3)) * pkin(7) * t102; t71 * t28 + t49 * mrSges(5,2) + t52 * t7 - t36 * t6 / 0.2e1 + t27 * t12 + t35 * t5 / 0.2e1 + t11 * t17 + t10 * t18 + t25 * t13 / 0.2e1 + t26 * t14 / 0.2e1 + (t46 + mrSges(5,1)) * t50 + t80 * mrSges(7,3) + (-t23 / 0.2e1 + t64 * t43 - t9 * mrSges(6,3)) * t70 + (-t24 / 0.2e1 - t64 * t44 + t8 * mrSges(6,3)) * t69 + m(6) * (t71 * t50 + t79 * t64) + m(5) * (qJ(3) * t50 + t75 * t49) + m(7) * (t10 * t1 + t11 * t2 + t52 * t27) + (-Ifges(6,5) * t69 / 0.2e1 + Ifges(6,6) * t96 - t31 / 0.2e1 + t30 / 0.2e1 + Ifges(5,6) + Ifges(4,4) + Ifges(3,5) - pkin(2) * mrSges(4,2) - t75 * mrSges(5,3) + (-mrSges(3,1) + t83) * pkin(7)) * t73 + (Ifges(5,5) - Ifges(4,6) + Ifges(3,6) + t69 * t47 / 0.2e1 + t48 * t96 + t87 * qJ(3) + (m(4) * qJ(3) - mrSges(3,2) + mrSges(4,3)) * pkin(7)) * t74; 0.2e1 * pkin(2) * mrSges(4,1) + 0.2e1 * t75 * mrSges(5,2) + 0.2e1 * t52 * t12 + t35 * t13 - t36 * t14 + 0.2e1 * t71 * t46 - t70 * t47 - t69 * t48 + Ifges(4,2) + Ifges(3,3) + Ifges(5,3) + m(7) * (t10 ^ 2 + t11 ^ 2 + t52 ^ 2) + m(6) * (t85 * t64 ^ 2 + t71 ^ 2) + m(5) * (t75 ^ 2 + t76) + m(4) * (pkin(2) ^ 2 + t76) + 0.2e1 * (mrSges(5,1) + mrSges(4,3)) * qJ(3) + 0.2e1 * t78 * mrSges(7,3) - 0.2e1 * t64 * t82; -t35 * t17 - t36 * t18 + t70 * t43 - t69 * t44 + (m(4) * pkin(7) + t87) * t73 - m(7) * t80 + m(6) * t79 + m(5) * t49; m(5) * t75 - m(7) * t78 - t101 * mrSges(7,3) + t64 * t40 + mrSges(5,2) - t82 + t83; m(4) + t100; -t74 * mrSges(5,2) + t36 * t17 - t35 * t18 + t69 * t43 + t70 * t44 + t57 + m(7) * (-t35 * t1 + t36 * t2) + m(6) * (t69 * t9 + t70 * t8) + m(5) * t32; m(7) * (-t35 * t10 + t36 * t11); 0; t100; m(6) * t50 + m(7) * t27 + t28 + t7; m(6) * t71 + m(7) * t52 + t12 + t46; 0; 0; m(6) + m(7); t1 * mrSges(7,1) - t2 * mrSges(7,2) + t84; t10 * mrSges(7,1) - t11 * mrSges(7,2) + t30 - t31; -t36 * mrSges(7,1) + t35 * mrSges(7,2); t12; 0; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t15(1) t15(2) t15(4) t15(7) t15(11) t15(16); t15(2) t15(3) t15(5) t15(8) t15(12) t15(17); t15(4) t15(5) t15(6) t15(9) t15(13) t15(18); t15(7) t15(8) t15(9) t15(10) t15(14) t15(19); t15(11) t15(12) t15(13) t15(14) t15(15) t15(20); t15(16) t15(17) t15(18) t15(19) t15(20) t15(21);];
Mq  = res;
