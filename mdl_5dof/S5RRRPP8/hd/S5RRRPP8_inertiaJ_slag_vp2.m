% Calculate joint inertia matrix for
% S5RRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPP8_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP8_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP8_inertiaJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP8_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP8_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP8_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:07:31
% EndTime: 2019-12-31 21:07:33
% DurationCPUTime: 0.74s
% Computational Cost: add. (467->214), mult. (929->272), div. (0->0), fcn. (662->4), ass. (0->79)
t95 = -Ifges(6,4) - Ifges(5,5);
t62 = sin(qJ(3));
t64 = cos(qJ(3));
t94 = t62 ^ 2 + t64 ^ 2;
t93 = 2 * pkin(6);
t92 = 2 * mrSges(6,1);
t91 = pkin(4) + pkin(7);
t65 = cos(qJ(2));
t90 = pkin(6) * t65;
t89 = Ifges(4,4) * t62;
t88 = Ifges(4,4) * t64;
t87 = Ifges(5,4) * t65;
t86 = Ifges(4,6) * t65;
t85 = Ifges(5,6) * t62;
t84 = Ifges(5,6) * t64;
t83 = Ifges(6,6) * t62;
t82 = Ifges(6,6) * t64;
t63 = sin(qJ(2));
t81 = t62 * t63;
t80 = t63 * t64;
t79 = t64 * mrSges(6,1);
t78 = mrSges(5,1) + mrSges(6,1);
t61 = -pkin(3) - qJ(5);
t77 = pkin(3) * t81 + t63 * pkin(6);
t27 = -t65 * pkin(2) - t63 * pkin(7) - pkin(1);
t7 = t62 * t27 + t64 * t90;
t22 = t65 * mrSges(6,3) + t63 * t79;
t76 = t94 * pkin(7) ^ 2;
t75 = qJ(4) * t64;
t74 = Ifges(4,3) + Ifges(6,1) + Ifges(5,1);
t25 = mrSges(5,1) * t80 - t65 * mrSges(5,2);
t73 = -t62 * qJ(4) - pkin(2);
t46 = t62 * t90;
t6 = t64 * t27 - t46;
t72 = t95 * t81 + (-Ifges(4,5) - Ifges(6,5)) * t80;
t4 = t65 * qJ(4) - t7;
t70 = t62 * mrSges(4,1) + t64 * mrSges(4,2);
t69 = -t62 * mrSges(5,2) - t64 * mrSges(5,3);
t24 = -mrSges(6,1) * t81 - t65 * mrSges(6,2);
t68 = pkin(6) ^ 2;
t66 = qJ(4) ^ 2;
t60 = t65 ^ 2;
t58 = t63 ^ 2;
t56 = t65 * pkin(3);
t53 = t58 * t68;
t51 = Ifges(4,5) * t62;
t50 = Ifges(6,5) * t62;
t49 = Ifges(4,6) * t64;
t38 = t91 * t64;
t37 = t91 * t62;
t36 = Ifges(4,1) * t62 + t88;
t35 = Ifges(4,2) * t64 + t89;
t34 = -Ifges(5,2) * t62 - t84;
t33 = -Ifges(6,2) * t64 + t83;
t32 = -Ifges(5,3) * t64 - t85;
t31 = Ifges(6,3) * t62 - t82;
t30 = -t62 * mrSges(6,2) - t64 * mrSges(6,3);
t29 = -t64 * mrSges(4,1) + t62 * mrSges(4,2);
t28 = t64 * mrSges(5,2) - t62 * mrSges(5,3);
t26 = -t64 * pkin(3) + t73;
t23 = mrSges(5,1) * t81 + t65 * mrSges(5,3);
t21 = -t65 * mrSges(4,1) - mrSges(4,3) * t80;
t20 = t65 * mrSges(4,2) - mrSges(4,3) * t81;
t18 = t70 * t63;
t17 = t69 * t63;
t16 = (-mrSges(6,2) * t64 + mrSges(6,3) * t62) * t63;
t15 = t61 * t64 + t73;
t14 = -t63 * t75 + t77;
t13 = -t87 + (-Ifges(5,2) * t64 + t85) * t63;
t12 = -Ifges(6,4) * t65 + (Ifges(6,2) * t62 + t82) * t63;
t11 = -Ifges(5,5) * t65 + (Ifges(5,3) * t62 - t84) * t63;
t10 = -Ifges(6,5) * t65 + (Ifges(6,3) * t64 + t83) * t63;
t9 = -Ifges(4,5) * t65 + (Ifges(4,1) * t64 - t89) * t63;
t8 = -t86 + (-Ifges(4,2) * t62 + t88) * t63;
t5 = t56 - t6;
t3 = (qJ(5) * t62 - t75) * t63 + t77;
t2 = -pkin(4) * t81 - t4;
t1 = t65 * qJ(5) + t46 + t56 + (pkin(4) * t63 - t27) * t64;
t19 = [0.2e1 * t1 * t22 + 0.2e1 * t14 * t17 + 0.2e1 * t3 * t16 + 0.2e1 * t2 * t24 + 0.2e1 * t7 * t20 + 0.2e1 * t6 * t21 + 0.2e1 * t4 * t23 + 0.2e1 * t5 * t25 + Ifges(2,3) + (t58 + t60) * mrSges(3,3) * t93 + m(6) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(5) * (t14 ^ 2 + t4 ^ 2 + t5 ^ 2) + m(4) * (t6 ^ 2 + t7 ^ 2 + t53) + m(3) * (pkin(1) ^ 2 + t60 * t68 + t53) + (0.2e1 * pkin(1) * mrSges(3,1) + (Ifges(3,2) + t74) * t65 + t72) * t65 + (-0.2e1 * pkin(1) * mrSges(3,2) + Ifges(3,1) * t63 + 0.2e1 * Ifges(3,4) * t65 + t18 * t93 + (t10 - t13 + t9 + t87) * t64 + (t11 + t12 - t8 + t86) * t62) * t63; -pkin(2) * t18 + t15 * t16 + t26 * t17 + t37 * t22 + t38 * t24 + t3 * t30 + m(6) * (t37 * t1 + t15 * t3 + t38 * t2) + (-t50 / 0.2e1 - t51 / 0.2e1 - t49 / 0.2e1 + Ifges(3,6) - pkin(6) * mrSges(3,2)) * t65 + (t8 / 0.2e1 - t11 / 0.2e1 - t12 / 0.2e1 + t2 * mrSges(6,1) - t4 * mrSges(5,1) + t7 * mrSges(4,3) + (Ifges(5,5) / 0.2e1 + Ifges(6,4) / 0.2e1) * t65) * t64 + (t87 / 0.2e1 + t9 / 0.2e1 + t10 / 0.2e1 - t13 / 0.2e1 + t1 * mrSges(6,1) + t5 * mrSges(5,1) - t6 * mrSges(4,3)) * t62 + ((t20 - t23) * t64 + (-t21 + t25) * t62 + m(5) * (-t4 * t64 + t5 * t62) + m(4) * (-t6 * t62 + t7 * t64)) * pkin(7) + (Ifges(3,5) + (t31 / 0.2e1 - t34 / 0.2e1 + t36 / 0.2e1) * t64 + (t32 / 0.2e1 + t33 / 0.2e1 - t35 / 0.2e1) * t62 + (-m(4) * pkin(2) - mrSges(3,1) + t29) * pkin(6)) * t63 + (m(5) * t26 + t28) * t14; -0.2e1 * pkin(2) * t29 + 0.2e1 * t15 * t30 + 0.2e1 * t26 * t28 + Ifges(3,3) + m(5) * (t26 ^ 2 + t76) + m(6) * (t15 ^ 2 + t37 ^ 2 + t38 ^ 2) + m(4) * (pkin(2) ^ 2 + t76) + (t38 * t92 - t32 - t33 + t35) * t64 + (t37 * t92 + t31 - t34 + t36) * t62 + 0.2e1 * (mrSges(5,1) + mrSges(4,3)) * pkin(7) * t94; t6 * mrSges(4,1) - t7 * mrSges(4,2) + t5 * mrSges(5,2) + t2 * mrSges(6,2) - t4 * mrSges(5,3) - t1 * mrSges(6,3) - pkin(3) * t25 + t61 * t22 + (-Ifges(5,4) * t64 - Ifges(4,6) * t62) * t63 + (-t23 + t24) * qJ(4) + m(6) * (qJ(4) * t2 + t61 * t1) + m(5) * (-pkin(3) * t5 - qJ(4) * t4) - t74 * t65 - t72; m(6) * (qJ(4) * t38 + t61 * t37) + t50 + t51 + t49 + t38 * mrSges(6,2) - t37 * mrSges(6,3) + (-pkin(3) * mrSges(5,1) + t61 * mrSges(6,1) - Ifges(5,4)) * t62 + (t78 * qJ(4) + t95) * t64 + (m(5) * (-pkin(3) * t62 + t75) - t69 - t70) * pkin(7); -0.2e1 * pkin(3) * mrSges(5,2) - 0.2e1 * t61 * mrSges(6,3) + 0.2e1 * (mrSges(5,3) + mrSges(6,2)) * qJ(4) + m(5) * (pkin(3) ^ 2 + t66) + m(6) * (t61 ^ 2 + t66) + t74; m(5) * t5 + m(6) * t1 + t22 + t25; m(6) * t37 + (m(5) * pkin(7) + t78) * t62; -m(5) * pkin(3) + m(6) * t61 + mrSges(5,2) - mrSges(6,3); m(5) + m(6); m(6) * t2 + t24; m(6) * t38 + t79; m(6) * qJ(4) + mrSges(6,2); 0; m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t19(1), t19(2), t19(4), t19(7), t19(11); t19(2), t19(3), t19(5), t19(8), t19(12); t19(4), t19(5), t19(6), t19(9), t19(13); t19(7), t19(8), t19(9), t19(10), t19(14); t19(11), t19(12), t19(13), t19(14), t19(15);];
Mq = res;
