% Calculate joint inertia matrix for
% S5RRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:36
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR11_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR11_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR11_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR11_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR11_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR11_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:32:26
% EndTime: 2019-12-31 21:32:29
% DurationCPUTime: 0.84s
% Computational Cost: add. (715->230), mult. (1383->315), div. (0->0), fcn. (1180->6), ass. (0->92)
t73 = sin(qJ(3));
t76 = cos(qJ(3));
t109 = t73 ^ 2 + t76 ^ 2;
t108 = 2 * pkin(6);
t74 = sin(qJ(2));
t94 = t74 * t76;
t95 = t73 * t74;
t107 = -Ifges(5,6) * t95 + (-Ifges(5,4) - Ifges(4,5)) * t94;
t106 = pkin(3) + pkin(4);
t105 = pkin(7) - pkin(8);
t77 = cos(qJ(2));
t104 = pkin(6) * t77;
t103 = Ifges(4,4) * t73;
t102 = Ifges(4,4) * t76;
t101 = Ifges(5,5) * t73;
t100 = Ifges(5,5) * t76;
t99 = Ifges(4,6) * t77;
t98 = Ifges(5,6) * t77;
t72 = sin(qJ(5));
t75 = cos(qJ(5));
t41 = -t72 * qJ(4) - t106 * t75;
t97 = t41 * mrSges(6,1);
t42 = t75 * qJ(4) - t106 * t72;
t96 = t42 * mrSges(6,2);
t93 = Ifges(5,2) + Ifges(4,3);
t39 = t77 * mrSges(5,1) + mrSges(5,2) * t94;
t44 = -t77 * pkin(2) - t74 * pkin(7) - pkin(1);
t19 = t76 * t104 + t73 * t44;
t92 = t109 * pkin(7) ^ 2;
t91 = qJ(4) * t76;
t36 = -t72 * t76 + t75 * t73;
t27 = t36 * t74;
t83 = t72 * t73 + t75 * t76;
t28 = t83 * t74;
t90 = Ifges(6,5) * t28 + Ifges(6,6) * t27 + Ifges(6,3) * t77;
t89 = t73 * qJ(4) + pkin(2);
t57 = t73 * t104;
t18 = t76 * t44 - t57;
t14 = -t77 * qJ(4) + t19;
t87 = t73 * mrSges(4,1) + t76 * mrSges(4,2);
t86 = t73 * mrSges(5,1) - t76 * mrSges(5,3);
t85 = t75 * mrSges(6,1) - t72 * mrSges(6,2);
t84 = -pkin(3) * t73 + t91;
t51 = t105 * t73;
t52 = t105 * t76;
t12 = t75 * t51 - t72 * t52;
t13 = t72 * t51 + t75 * t52;
t32 = Ifges(6,6) * t83;
t33 = Ifges(6,5) * t36;
t82 = t12 * mrSges(6,1) - t13 * mrSges(6,2) - t32 + t33;
t67 = t77 * pkin(3);
t6 = t77 * pkin(4) + t57 + t67 + (-pkin(8) * t74 - t44) * t76;
t7 = pkin(8) * t95 + t14;
t1 = t75 * t6 - t72 * t7;
t2 = t72 * t6 + t75 * t7;
t81 = t1 * mrSges(6,1) - t2 * mrSges(6,2) + t90;
t80 = pkin(6) ^ 2;
t71 = t77 ^ 2;
t69 = t74 ^ 2;
t65 = t69 * t80;
t63 = Ifges(5,4) * t73;
t62 = Ifges(4,5) * t73;
t61 = Ifges(4,6) * t76;
t50 = Ifges(4,1) * t73 + t102;
t49 = Ifges(5,1) * t73 - t100;
t48 = Ifges(4,2) * t76 + t103;
t47 = -Ifges(5,3) * t76 + t101;
t46 = -t76 * mrSges(4,1) + t73 * mrSges(4,2);
t45 = -t76 * mrSges(5,1) - t73 * mrSges(5,3);
t43 = -t76 * pkin(3) - t89;
t40 = -mrSges(5,2) * t95 - t77 * mrSges(5,3);
t38 = -t77 * mrSges(4,1) - mrSges(4,3) * t94;
t37 = t77 * mrSges(4,2) - mrSges(4,3) * t95;
t31 = t106 * t76 + t89;
t30 = t87 * t74;
t29 = t86 * t74;
t26 = (pkin(6) - t84) * t74;
t25 = -Ifges(4,5) * t77 + (Ifges(4,1) * t76 - t103) * t74;
t24 = -Ifges(5,4) * t77 + (Ifges(5,1) * t76 + t101) * t74;
t23 = -t99 + (-Ifges(4,2) * t73 + t102) * t74;
t22 = -t98 + (Ifges(5,3) * t73 + t100) * t74;
t17 = t77 * mrSges(6,1) - t28 * mrSges(6,3);
t16 = -t77 * mrSges(6,2) + t27 * mrSges(6,3);
t15 = -t18 + t67;
t11 = (-t106 * t73 - pkin(6) + t91) * t74;
t10 = Ifges(6,1) * t36 - Ifges(6,4) * t83;
t9 = Ifges(6,4) * t36 - Ifges(6,2) * t83;
t8 = mrSges(6,1) * t83 + t36 * mrSges(6,2);
t5 = -t27 * mrSges(6,1) + t28 * mrSges(6,2);
t4 = Ifges(6,1) * t28 + Ifges(6,4) * t27 + Ifges(6,5) * t77;
t3 = Ifges(6,4) * t28 + Ifges(6,2) * t27 + Ifges(6,6) * t77;
t20 = [0.2e1 * t1 * t17 + 0.2e1 * t11 * t5 + 0.2e1 * t14 * t40 + 0.2e1 * t15 * t39 + 0.2e1 * t2 * t16 + 0.2e1 * t18 * t38 + 0.2e1 * t19 * t37 + 0.2e1 * t26 * t29 + t27 * t3 + t28 * t4 + Ifges(2,3) + (t69 + t71) * mrSges(3,3) * t108 + m(6) * (t1 ^ 2 + t11 ^ 2 + t2 ^ 2) + m(5) * (t14 ^ 2 + t15 ^ 2 + t26 ^ 2) + m(4) * (t18 ^ 2 + t19 ^ 2 + t65) + m(3) * (pkin(1) ^ 2 + t71 * t80 + t65) + (0.2e1 * pkin(1) * mrSges(3,1) + (Ifges(3,2) + t93) * t77 + t90 + t107) * t77 + (-0.2e1 * pkin(1) * mrSges(3,2) + Ifges(3,1) * t74 + 0.2e1 * Ifges(3,4) * t77 + t30 * t108 + (t24 + t25) * t76 + (t22 - t23 + t99) * t73) * t74; -t83 * t3 / 0.2e1 + t36 * t4 / 0.2e1 + t43 * t29 + t27 * t9 / 0.2e1 + t28 * t10 / 0.2e1 - pkin(2) * t30 + t31 * t5 + t11 * t8 + t13 * t16 + t12 * t17 + (-t1 * t36 - t2 * t83) * mrSges(6,3) + m(6) * (t12 * t1 + t31 * t11 + t13 * t2) + (-t63 / 0.2e1 - t62 / 0.2e1 - t61 / 0.2e1 + t33 / 0.2e1 - t32 / 0.2e1 + Ifges(3,6) - pkin(6) * mrSges(3,2)) * t77 + (t23 / 0.2e1 - t22 / 0.2e1 + t98 / 0.2e1 + t14 * mrSges(5,2) + t19 * mrSges(4,3)) * t76 + (t24 / 0.2e1 + t25 / 0.2e1 + t15 * mrSges(5,2) - t18 * mrSges(4,3)) * t73 + ((t37 + t40) * t76 + (-t38 + t39) * t73 + m(5) * (t14 * t76 + t15 * t73) + m(4) * (-t18 * t73 + t19 * t76)) * pkin(7) + (Ifges(3,5) + (t49 / 0.2e1 + t50 / 0.2e1) * t76 + (t47 / 0.2e1 - t48 / 0.2e1) * t73 + (-m(4) * pkin(2) - mrSges(3,1) + t46) * pkin(6)) * t74 + (m(5) * t43 + t45) * t26; -0.2e1 * pkin(2) * t46 + t36 * t10 + 0.2e1 * t31 * t8 - t83 * t9 + 0.2e1 * t43 * t45 + Ifges(3,3) + (-t47 + t48) * t76 + (t50 + t49) * t73 + 0.2e1 * (-t12 * t36 - t13 * t83) * mrSges(6,3) + m(6) * (t12 ^ 2 + t13 ^ 2 + t31 ^ 2) + m(5) * (t43 ^ 2 + t92) + m(4) * (pkin(2) ^ 2 + t92) + 0.2e1 * (mrSges(5,2) + mrSges(4,3)) * pkin(7) * t109; m(6) * (t41 * t1 + t42 * t2) + m(5) * (-pkin(3) * t15 + qJ(4) * t14) - t93 * t77 - pkin(3) * t39 + qJ(4) * t40 + t41 * t17 + t42 * t16 + t14 * mrSges(5,3) - t15 * mrSges(5,1) + t18 * mrSges(4,1) - t19 * mrSges(4,2) - t81 - Ifges(4,6) * t95 - t107; m(6) * (t41 * t12 + t42 * t13) + t63 - Ifges(5,6) * t76 + t62 + t61 + (-t41 * t36 - t42 * t83) * mrSges(6,3) + t84 * mrSges(5,2) + (m(5) * t84 - t86 - t87) * pkin(7) - t82; 0.2e1 * pkin(3) * mrSges(5,1) - 0.2e1 * t97 + 0.2e1 * t96 + 0.2e1 * qJ(4) * mrSges(5,3) + Ifges(6,3) + m(5) * (pkin(3) ^ 2 + qJ(4) ^ 2) + m(6) * (t41 ^ 2 + t42 ^ 2) + t93; t72 * t16 + t75 * t17 + m(6) * (t75 * t1 + t72 * t2) + m(5) * t15 + t39; m(6) * (t75 * t12 + t72 * t13) + (m(5) * pkin(7) + mrSges(5,2)) * t73 + (-t75 * t36 - t72 * t83) * mrSges(6,3); -mrSges(5,1) - m(5) * pkin(3) + m(6) * (t75 * t41 + t72 * t42) - t85; m(5) + m(6) * (t72 ^ 2 + t75 ^ 2); t81; t82; -Ifges(6,3) - t96 + t97; t85; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t20(1), t20(2), t20(4), t20(7), t20(11); t20(2), t20(3), t20(5), t20(8), t20(12); t20(4), t20(5), t20(6), t20(9), t20(13); t20(7), t20(8), t20(9), t20(10), t20(14); t20(11), t20(12), t20(13), t20(14), t20(15);];
Mq = res;
