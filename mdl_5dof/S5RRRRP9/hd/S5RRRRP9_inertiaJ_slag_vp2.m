% Calculate joint inertia matrix for
% S5RRRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 22:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRP9_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP9_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP9_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP9_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP9_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP9_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:03:49
% EndTime: 2019-12-31 22:03:51
% DurationCPUTime: 0.80s
% Computational Cost: add. (909->234), mult. (1816->323), div. (0->0), fcn. (1658->6), ass. (0->84)
t107 = 2 * pkin(6);
t106 = 2 * mrSges(6,3);
t105 = -pkin(8) - pkin(7);
t85 = cos(qJ(2));
t104 = pkin(6) * t85;
t82 = sin(qJ(2));
t75 = t82 * pkin(6);
t81 = sin(qJ(3));
t103 = Ifges(4,4) * t81;
t84 = cos(qJ(3));
t102 = Ifges(4,4) * t84;
t101 = t81 * t82;
t100 = t82 * t84;
t98 = -Ifges(5,3) - Ifges(6,2);
t59 = -pkin(2) * t85 - pkin(7) * t82 - pkin(1);
t52 = t84 * t59;
t16 = -pkin(8) * t100 + t52 + (-pkin(6) * t81 - pkin(3)) * t85;
t36 = t84 * t104 + t81 * t59;
t24 = -pkin(8) * t101 + t36;
t80 = sin(qJ(4));
t83 = cos(qJ(4));
t6 = t80 * t16 + t83 * t24;
t53 = t80 * t81 - t83 * t84;
t45 = t53 * t82;
t34 = t85 * mrSges(6,1) - t45 * mrSges(6,2);
t58 = pkin(3) * t101 + t75;
t97 = t81 ^ 2 + t84 ^ 2;
t63 = t105 * t84;
t94 = t105 * t81;
t28 = -t63 * t80 - t83 * t94;
t30 = -t83 * t63 + t80 * t94;
t96 = t28 ^ 2 + t30 ^ 2;
t95 = -Ifges(4,3) + t98;
t70 = -pkin(3) * t84 - pkin(2);
t54 = t80 * t84 + t81 * t83;
t44 = t54 * t82;
t93 = (Ifges(6,4) + Ifges(5,5)) * t45 + (Ifges(5,6) - Ifges(6,6)) * t44;
t92 = mrSges(4,1) * t81 + mrSges(4,2) * t84;
t5 = t16 * t83 - t24 * t80;
t91 = (mrSges(5,1) * t83 - mrSges(5,2) * t80) * pkin(3);
t2 = -qJ(5) * t85 + t6;
t3 = pkin(4) * t85 - t5;
t90 = t5 * mrSges(5,1) - t3 * mrSges(6,1) - t6 * mrSges(5,2) + t2 * mrSges(6,3) - t93;
t47 = Ifges(6,6) * t53;
t48 = Ifges(5,6) * t53;
t49 = Ifges(5,5) * t54;
t50 = Ifges(6,4) * t54;
t89 = t47 - t48 + t49 + t50 + (-mrSges(5,2) + mrSges(6,3)) * t30 + (-mrSges(5,1) - mrSges(6,1)) * t28;
t87 = pkin(6) ^ 2;
t79 = t85 ^ 2;
t77 = t82 ^ 2;
t74 = t77 * t87;
t73 = Ifges(4,5) * t81;
t72 = Ifges(4,6) * t84;
t69 = -pkin(3) * t83 - pkin(4);
t67 = pkin(3) * t80 + qJ(5);
t64 = Ifges(4,5) * t100;
t62 = Ifges(4,1) * t81 + t102;
t61 = Ifges(4,2) * t84 + t103;
t60 = -mrSges(4,1) * t84 + mrSges(4,2) * t81;
t57 = -mrSges(4,1) * t85 - mrSges(4,3) * t100;
t56 = mrSges(4,2) * t85 - mrSges(4,3) * t101;
t46 = t92 * t82;
t43 = -Ifges(4,5) * t85 + (Ifges(4,1) * t84 - t103) * t82;
t42 = -Ifges(4,6) * t85 + (-Ifges(4,2) * t81 + t102) * t82;
t35 = -t81 * t104 + t52;
t33 = -mrSges(5,1) * t85 + mrSges(5,3) * t45;
t32 = mrSges(5,2) * t85 - mrSges(5,3) * t44;
t31 = -mrSges(6,2) * t44 - mrSges(6,3) * t85;
t23 = Ifges(5,1) * t54 - Ifges(5,4) * t53;
t22 = Ifges(6,1) * t54 + Ifges(6,5) * t53;
t21 = Ifges(5,4) * t54 - Ifges(5,2) * t53;
t20 = Ifges(6,5) * t54 + Ifges(6,3) * t53;
t19 = mrSges(5,1) * t53 + mrSges(5,2) * t54;
t18 = mrSges(6,1) * t53 - mrSges(6,3) * t54;
t15 = pkin(4) * t53 - qJ(5) * t54 + t70;
t13 = mrSges(5,1) * t44 - mrSges(5,2) * t45;
t12 = mrSges(6,1) * t44 + mrSges(6,3) * t45;
t11 = -Ifges(5,1) * t45 - Ifges(5,4) * t44 - Ifges(5,5) * t85;
t10 = -Ifges(6,1) * t45 - Ifges(6,4) * t85 + Ifges(6,5) * t44;
t9 = -Ifges(5,4) * t45 - Ifges(5,2) * t44 - Ifges(5,6) * t85;
t8 = -Ifges(6,5) * t45 - Ifges(6,6) * t85 + Ifges(6,3) * t44;
t7 = pkin(4) * t44 + qJ(5) * t45 + t58;
t1 = [0.2e1 * t7 * t12 + 0.2e1 * t58 * t13 + 0.2e1 * t2 * t31 + 0.2e1 * t3 * t34 + 0.2e1 * t6 * t32 + 0.2e1 * t5 * t33 + 0.2e1 * t35 * t57 + 0.2e1 * t36 * t56 + Ifges(2,3) - (t10 + t11) * t45 + (t8 - t9) * t44 + (t77 + t79) * mrSges(3,3) * t107 + (-0.2e1 * mrSges(3,2) * pkin(1) + Ifges(3,1) * t82 + t46 * t107 - t42 * t81 + t43 * t84) * t82 + m(3) * (pkin(1) ^ 2 + t79 * t87 + t74) + m(4) * (t35 ^ 2 + t36 ^ 2 + t74) + m(5) * (t5 ^ 2 + t58 ^ 2 + t6 ^ 2) + m(6) * (t2 ^ 2 + t3 ^ 2 + t7 ^ 2) + (0.2e1 * pkin(1) * mrSges(3,1) - t64 + (Ifges(3,2) - t95) * t85 + (Ifges(4,6) * t81 + (2 * Ifges(3,4))) * t82 + t93) * t85; -pkin(2) * t46 + t15 * t12 + t70 * t13 + t7 * t18 + t58 * t19 + (-t73 / 0.2e1 - t72 / 0.2e1 - t49 / 0.2e1 + t48 / 0.2e1 - t50 / 0.2e1 - t47 / 0.2e1 + Ifges(3,6) - pkin(6) * mrSges(3,2)) * t85 - (t22 / 0.2e1 + t23 / 0.2e1) * t45 + (t20 / 0.2e1 - t21 / 0.2e1) * t44 + (t31 + t32) * t30 + (-t33 + t34) * t28 + (t42 / 0.2e1 + pkin(7) * t56 + t36 * mrSges(4,3)) * t84 + (t43 / 0.2e1 - pkin(7) * t57 - t35 * mrSges(4,3)) * t81 + (Ifges(3,5) + t84 * t62 / 0.2e1 - t81 * t61 / 0.2e1 + (t60 - mrSges(3,1)) * pkin(6)) * t82 + m(4) * (-pkin(2) * t75 + (-t35 * t81 + t36 * t84) * pkin(7)) + m(5) * (-t28 * t5 + t30 * t6 + t58 * t70) + m(6) * (t15 * t7 + t2 * t30 + t28 * t3) + (-t5 * mrSges(5,3) + t3 * mrSges(6,2) + t10 / 0.2e1 + t11 / 0.2e1) * t54 + (-t2 * mrSges(6,2) - t6 * mrSges(5,3) + t8 / 0.2e1 - t9 / 0.2e1) * t53; -0.2e1 * pkin(2) * t60 + 0.2e1 * t15 * t18 + 0.2e1 * t70 * t19 + t84 * t61 + t81 * t62 + Ifges(3,3) + m(6) * (t15 ^ 2 + t96) + m(5) * (t70 ^ 2 + t96) + m(4) * (t97 * pkin(7) ^ 2 + pkin(2) ^ 2) + (t22 + t23) * t54 + (t20 - t21) * t53 + 0.2e1 * t97 * pkin(7) * mrSges(4,3) + 0.2e1 * (t28 * t54 - t30 * t53) * (mrSges(6,2) + mrSges(5,3)); m(6) * (t2 * t67 + t3 * t69) + t95 * t85 + t64 - Ifges(4,6) * t101 + t67 * t31 + t69 * t34 + t35 * mrSges(4,1) - t36 * mrSges(4,2) + t90 + (t80 * t32 + t83 * t33 + m(5) * (t5 * t83 + t6 * t80)) * pkin(3); m(6) * (t28 * t69 + t30 * t67) + t73 + t72 - t92 * pkin(7) + (-t53 * t67 + t54 * t69) * mrSges(6,2) + (m(5) * (-t28 * t83 + t30 * t80) + (-t53 * t80 - t54 * t83) * mrSges(5,3)) * pkin(3) + t89; -0.2e1 * t69 * mrSges(6,1) + t67 * t106 + 0.2e1 * t91 + m(6) * (t67 ^ 2 + t69 ^ 2) + m(5) * (t80 ^ 2 + t83 ^ 2) * pkin(3) ^ 2 - t95; -pkin(4) * t34 + m(6) * (-pkin(4) * t3 + qJ(5) * t2) + qJ(5) * t31 + t98 * t85 + t90; m(6) * (-pkin(4) * t28 + qJ(5) * t30) + (-pkin(4) * t54 - qJ(5) * t53) * mrSges(6,2) + t89; m(6) * (-pkin(4) * t69 + qJ(5) * t67) + t91 + (qJ(5) + t67) * mrSges(6,3) + (pkin(4) - t69) * mrSges(6,1) - t98; 0.2e1 * pkin(4) * mrSges(6,1) + qJ(5) * t106 + m(6) * (pkin(4) ^ 2 + qJ(5) ^ 2) - t98; m(6) * t3 + t34; m(6) * t28 + t54 * mrSges(6,2); m(6) * t69 - mrSges(6,1); -m(6) * pkin(4) - mrSges(6,1); m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
