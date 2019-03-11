% Calculate joint inertia matrix for
% S6RPRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
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
% Datum: 2019-03-09 03:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRR5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR5_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR5_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR5_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR5_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR5_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:48:00
% EndTime: 2019-03-09 03:48:02
% DurationCPUTime: 0.69s
% Computational Cost: add. (1406->204), mult. (2582->271), div. (0->0), fcn. (2810->8), ass. (0->79)
t55 = sin(qJ(6));
t57 = cos(qJ(6));
t76 = t55 ^ 2 + t57 ^ 2;
t71 = t76 * mrSges(7,3);
t37 = -t57 * mrSges(7,1) + mrSges(7,2) * t55;
t78 = -t37 + mrSges(6,1);
t53 = sin(pkin(10));
t54 = cos(pkin(10));
t90 = sin(qJ(3));
t91 = cos(qJ(3));
t30 = t90 * t53 - t91 * t54;
t31 = t91 * t53 + t90 * t54;
t56 = sin(qJ(5));
t58 = cos(qJ(5));
t18 = -t58 * t30 + t31 * t56;
t19 = t30 * t56 + t31 * t58;
t87 = t19 * t55;
t10 = -mrSges(7,2) * t18 - mrSges(7,3) * t87;
t86 = t19 * t57;
t11 = mrSges(7,1) * t18 - mrSges(7,3) * t86;
t101 = t57 * t10 - t55 * t11;
t100 = m(7) * pkin(5) + t78;
t80 = pkin(7) + qJ(2);
t72 = t80 * t54;
t73 = t80 * t53;
t23 = t91 * t72 - t90 * t73;
t14 = pkin(8) * t30 + t23;
t21 = t90 * t72 + t91 * t73;
t62 = -t31 * pkin(8) + t21;
t6 = t14 * t56 - t58 * t62;
t99 = t6 ^ 2;
t98 = 0.2e1 * t6;
t48 = t54 ^ 2;
t59 = -pkin(3) - pkin(4);
t42 = -pkin(2) * t54 - pkin(1);
t65 = qJ(4) * t31 - t42;
t12 = t59 * t30 + t65;
t97 = 0.2e1 * t12;
t96 = -0.2e1 * t37;
t95 = 0.2e1 * t42;
t94 = t57 / 0.2e1;
t93 = t58 * t6;
t8 = t58 * t14 + t56 * t62;
t92 = t8 * mrSges(6,2);
t89 = Ifges(7,4) * t55;
t88 = Ifges(7,4) * t57;
t35 = -qJ(4) * t56 + t58 * t59;
t85 = t35 * mrSges(6,1);
t36 = t58 * qJ(4) + t56 * t59;
t84 = t36 * mrSges(6,2);
t81 = mrSges(5,2) + mrSges(4,3);
t79 = Ifges(7,5) * t86 + Ifges(7,3) * t18;
t77 = t53 ^ 2 + t48;
t38 = Ifges(7,5) * t55 + Ifges(7,6) * t57;
t75 = t38 / 0.2e1 - Ifges(6,6);
t74 = t21 ^ 2 + t23 ^ 2;
t34 = -pkin(9) + t36;
t70 = t76 * t34;
t69 = t76 * t56;
t68 = -t54 * mrSges(3,1) + t53 * mrSges(3,2);
t3 = pkin(5) * t18 - pkin(9) * t19 + t12;
t1 = t3 * t57 - t55 * t8;
t2 = t3 * t55 + t57 * t8;
t67 = -t1 * t55 + t2 * t57;
t66 = mrSges(7,1) * t55 + mrSges(7,2) * t57;
t39 = Ifges(7,2) * t57 + t89;
t40 = Ifges(7,1) * t55 + t88;
t64 = t57 * t39 + t55 * t40 + Ifges(6,3);
t63 = t40 * t94 - t55 * t39 / 0.2e1 + Ifges(6,5);
t52 = t58 ^ 2;
t50 = t56 ^ 2;
t33 = pkin(5) - t35;
t26 = t31 * mrSges(4,2);
t25 = t30 * mrSges(5,1);
t16 = pkin(3) * t30 - t65;
t9 = t66 * t19;
t5 = Ifges(7,5) * t18 + (Ifges(7,1) * t57 - t89) * t19;
t4 = Ifges(7,6) * t18 + (-Ifges(7,2) * t55 + t88) * t19;
t7 = [Ifges(3,2) * t48 + t26 * t95 - 0.2e1 * pkin(1) * t68 + 0.2e1 * t16 * t25 + t9 * t98 + 0.2e1 * t2 * t10 + 0.2e1 * t1 * t11 + Ifges(2,3) + (Ifges(3,1) * t53 + 0.2e1 * Ifges(3,4) * t54) * t53 + 0.2e1 * t77 * qJ(2) * mrSges(3,3) + (mrSges(6,1) * t97 - 0.2e1 * t8 * mrSges(6,3) + Ifges(6,2) * t18 + t79) * t18 + (mrSges(6,2) * t97 + mrSges(6,3) * t98 + Ifges(6,1) * t19 - t4 * t55 + t57 * t5 + (-Ifges(7,6) * t55 - (2 * Ifges(6,4))) * t18) * t19 + (-0.2e1 * t16 * mrSges(5,3) + (Ifges(5,1) + Ifges(4,1)) * t31 + 0.2e1 * t81 * t21) * t31 + (mrSges(4,1) * t95 + (Ifges(5,3) + Ifges(4,2)) * t30 + 0.2e1 * (-Ifges(4,4) + Ifges(5,5)) * t31 - 0.2e1 * t81 * t23) * t30 + m(3) * (t77 * qJ(2) ^ 2 + pkin(1) ^ 2) + m(4) * (t42 ^ 2 + t74) + m(5) * (t16 ^ 2 + t74) + m(6) * (t12 ^ 2 + t8 ^ 2 + t99) + m(7) * (t1 ^ 2 + t2 ^ 2 + t99); -m(3) * pkin(1) + t30 * mrSges(4,1) - t18 * mrSges(6,1) - t19 * mrSges(6,2) - t31 * mrSges(5,3) - t55 * t10 - t57 * t11 + t25 + t26 + m(7) * (-t1 * t57 - t2 * t55) - m(6) * t12 + m(5) * t16 + m(4) * t42 + t68; m(7) * t76 + m(3) + m(4) + m(5) + m(6); t92 + t33 * t9 + t78 * t6 + (-mrSges(4,2) + mrSges(5,3)) * t23 + (-mrSges(5,1) - mrSges(4,1)) * t21 + (-t4 / 0.2e1 - t2 * mrSges(7,3) + t34 * t10) * t57 + (-t5 / 0.2e1 + t1 * mrSges(7,3) - t34 * t11) * t55 + (-mrSges(5,2) * pkin(3) + Ifges(5,4) + Ifges(4,5)) * t31 + (-mrSges(5,2) * qJ(4) - Ifges(4,6) + Ifges(5,6)) * t30 + (-t36 * mrSges(6,3) - t75) * t18 + m(7) * (t33 * t6 + t67 * t34) + m(6) * (-t35 * t6 + t36 * t8) + m(5) * (-pkin(3) * t21 + qJ(4) * t23) + (-t35 * mrSges(6,3) - t63) * t19; 0; 0.2e1 * pkin(3) * mrSges(5,1) - 0.2e1 * t85 + 0.2e1 * t84 + 0.2e1 * qJ(4) * mrSges(5,3) + t33 * t96 + Ifges(5,2) + Ifges(4,3) - 0.2e1 * t34 * t71 + m(7) * (t76 * t34 ^ 2 + t33 ^ 2) + m(6) * (t35 ^ 2 + t36 ^ 2) + m(5) * (pkin(3) ^ 2 + qJ(4) ^ 2) + t64; t31 * mrSges(5,2) + (-t19 * mrSges(6,3) - t9) * t58 + (-t18 * mrSges(6,3) + t101) * t56 + m(7) * (t67 * t56 - t93) + m(6) * (t56 * t8 - t93) + m(5) * t21; 0; -m(5) * pkin(3) - mrSges(5,1) - t78 * t58 + (mrSges(6,2) - t71) * t56 + m(7) * (-t33 * t58 + t34 * t69) + m(6) * (t35 * t58 + t36 * t56); m(5) + m(6) * (t50 + t52) + m(7) * (t76 * t50 + t52); t55 * t5 / 0.2e1 + t4 * t94 - pkin(5) * t9 - t92 + t75 * t18 + t67 * mrSges(7,3) + t63 * t19 - t100 * t6 + (m(7) * t67 + t101) * pkin(9); 0; m(7) * (-pkin(5) * t33 + pkin(9) * t70) - t84 + t85 + (pkin(5) + t33) * t37 + (-t76 * pkin(9) + t70) * mrSges(7,3) - t64; -t56 * mrSges(6,2) + (m(7) * pkin(9) + mrSges(7,3)) * t69 + t100 * t58; pkin(5) * t96 + m(7) * (t76 * pkin(9) ^ 2 + pkin(5) ^ 2) + 0.2e1 * pkin(9) * t71 + t64; mrSges(7,1) * t1 - mrSges(7,2) * t2 - Ifges(7,6) * t87 + t79; t37; -t66 * t34 - t38; -t66 * t56; -t66 * pkin(9) + t38; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t7(1) t7(2) t7(4) t7(7) t7(11) t7(16); t7(2) t7(3) t7(5) t7(8) t7(12) t7(17); t7(4) t7(5) t7(6) t7(9) t7(13) t7(18); t7(7) t7(8) t7(9) t7(10) t7(14) t7(19); t7(11) t7(12) t7(13) t7(14) t7(15) t7(20); t7(16) t7(17) t7(18) t7(19) t7(20) t7(21);];
Mq  = res;
