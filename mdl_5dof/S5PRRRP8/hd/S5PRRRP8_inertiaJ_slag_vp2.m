% Calculate joint inertia matrix for
% S5PRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 17:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRP8_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP8_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP8_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP8_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP8_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP8_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:58:16
% EndTime: 2019-12-05 16:58:18
% DurationCPUTime: 0.61s
% Computational Cost: add. (403->176), mult. (938->237), div. (0->0), fcn. (791->8), ass. (0->76)
t56 = sin(qJ(4));
t59 = cos(qJ(4));
t97 = t56 ^ 2 + t59 ^ 2;
t55 = cos(pkin(5));
t57 = sin(qJ(3));
t60 = cos(qJ(3));
t54 = sin(pkin(5));
t58 = sin(qJ(2));
t81 = t54 * t58;
t18 = t55 * t57 + t60 * t81;
t61 = cos(qJ(2));
t80 = t54 * t61;
t3 = t18 * t56 + t59 * t80;
t5 = t18 * t59 - t56 * t80;
t96 = t3 * t56 + t5 * t59;
t95 = 2 * pkin(7);
t94 = mrSges(5,3) + mrSges(6,2);
t93 = -m(6) * pkin(4) - mrSges(6,1);
t16 = -t55 * t60 + t57 * t81;
t15 = t16 ^ 2;
t92 = pkin(7) * t60;
t89 = Ifges(5,4) * t56;
t88 = Ifges(5,4) * t59;
t87 = Ifges(6,5) * t56;
t86 = Ifges(6,5) * t59;
t85 = Ifges(5,6) * t60;
t84 = Ifges(6,6) * t60;
t83 = t16 * t57;
t82 = t18 * t60;
t79 = t56 * t57;
t78 = t57 * t59;
t27 = -t60 * pkin(3) - t57 * pkin(8) - pkin(2);
t77 = t59 * t27;
t76 = Ifges(6,2) + Ifges(5,3);
t22 = t60 * mrSges(5,2) - mrSges(5,3) * t79;
t25 = -mrSges(6,2) * t79 - t60 * mrSges(6,3);
t75 = t22 + t25;
t23 = -t60 * mrSges(5,1) - mrSges(5,3) * t78;
t24 = t60 * mrSges(6,1) + mrSges(6,2) * t78;
t74 = -t23 + t24;
t29 = -t59 * mrSges(5,1) + t56 * mrSges(5,2);
t73 = t29 - mrSges(4,1);
t9 = t56 * t27 + t59 * t92;
t72 = t97 * pkin(8) ^ 2;
t71 = -Ifges(6,6) * t79 + (-Ifges(6,4) - Ifges(5,5)) * t78;
t69 = t96 * pkin(8);
t66 = t56 * mrSges(5,1) + t59 * mrSges(5,2);
t65 = t56 * mrSges(6,1) - t59 * mrSges(6,3);
t64 = -pkin(4) * t56 + qJ(5) * t59;
t63 = pkin(7) ^ 2;
t53 = t60 ^ 2;
t51 = t57 ^ 2;
t49 = t54 ^ 2;
t47 = t51 * t63;
t45 = Ifges(6,4) * t56;
t44 = Ifges(5,5) * t56;
t43 = Ifges(5,6) * t59;
t40 = t49 * t61 ^ 2;
t34 = Ifges(5,1) * t56 + t88;
t33 = Ifges(6,1) * t56 - t86;
t32 = Ifges(5,2) * t59 + t89;
t31 = -Ifges(6,3) * t59 + t87;
t30 = -t60 * mrSges(4,1) + t57 * mrSges(4,2);
t28 = -t59 * mrSges(6,1) - t56 * mrSges(6,3);
t26 = -t59 * pkin(4) - t56 * qJ(5) - pkin(3);
t20 = t66 * t57;
t19 = t65 * t57;
t14 = (pkin(7) - t64) * t57;
t13 = -Ifges(5,5) * t60 + (Ifges(5,1) * t59 - t89) * t57;
t12 = -Ifges(6,4) * t60 + (Ifges(6,1) * t59 + t87) * t57;
t11 = -t85 + (-Ifges(5,2) * t56 + t88) * t57;
t10 = -t84 + (Ifges(6,3) * t56 + t86) * t57;
t8 = -t56 * t92 + t77;
t7 = -t77 + (pkin(7) * t56 + pkin(4)) * t60;
t6 = -t60 * qJ(5) + t9;
t1 = [m(2) + m(4) * (t18 ^ 2 + t15 + t40) + m(3) * (t49 * t58 ^ 2 + t55 ^ 2 + t40) + (m(5) + m(6)) * (t3 ^ 2 + t5 ^ 2 + t15); mrSges(4,3) * t82 + t75 * t5 + t74 * t3 + (-t58 * mrSges(3,2) + (mrSges(3,1) - t30) * t61) * t54 + (t57 * mrSges(4,3) + t19 + t20) * t16 + m(5) * (pkin(7) * t83 - t8 * t3 + t9 * t5) + m(6) * (t14 * t16 + t7 * t3 + t6 * t5) + m(4) * (pkin(2) * t80 + (t82 + t83) * pkin(7)); -0.2e1 * pkin(2) * t30 + 0.2e1 * t14 * t19 + 0.2e1 * t9 * t22 + 0.2e1 * t8 * t23 + 0.2e1 * t7 * t24 + 0.2e1 * t6 * t25 + Ifges(3,3) + (t51 + t53) * mrSges(4,3) * t95 + m(5) * (t8 ^ 2 + t9 ^ 2 + t47) + m(6) * (t14 ^ 2 + t6 ^ 2 + t7 ^ 2) + m(4) * (pkin(2) ^ 2 + t53 * t63 + t47) + ((Ifges(4,2) + t76) * t60 + t71) * t60 + (Ifges(4,1) * t57 + 0.2e1 * Ifges(4,4) * t60 + t20 * t95 + (t12 + t13) * t59 + (t10 - t11 + t85) * t56) * t57; -t18 * mrSges(4,2) + (t28 + t73) * t16 + m(5) * (-pkin(3) * t16 + t69) + m(6) * (t26 * t16 + t69) + t94 * t96; -pkin(3) * t20 + t26 * t19 + (m(6) * t26 + t28) * t14 + (-pkin(7) * mrSges(4,2) + Ifges(4,6) - t44 / 0.2e1 - t43 / 0.2e1 - t45 / 0.2e1) * t60 + (t6 * mrSges(6,2) + t9 * mrSges(5,3) - t10 / 0.2e1 + t11 / 0.2e1 + t84 / 0.2e1) * t59 + (t7 * mrSges(6,2) - t8 * mrSges(5,3) + t12 / 0.2e1 + t13 / 0.2e1) * t56 + (t75 * t59 + t74 * t56 + m(6) * (t7 * t56 + t6 * t59) + m(5) * (-t8 * t56 + t9 * t59)) * pkin(8) + (Ifges(4,5) + (t33 / 0.2e1 + t34 / 0.2e1) * t59 + (t31 / 0.2e1 - t32 / 0.2e1) * t56 + (-m(5) * pkin(3) + t73) * pkin(7)) * t57; -0.2e1 * pkin(3) * t29 + 0.2e1 * t26 * t28 + Ifges(4,3) + (t32 - t31) * t59 + (t33 + t34) * t56 + m(6) * (t26 ^ 2 + t72) + m(5) * (pkin(3) ^ 2 + t72) + 0.2e1 * t94 * pkin(8) * t97; (m(6) * qJ(5) - mrSges(5,2) + mrSges(6,3)) * t5 + (-mrSges(5,1) + t93) * t3; -Ifges(5,6) * t79 - pkin(4) * t24 + m(6) * (-pkin(4) * t7 + qJ(5) * t6) + qJ(5) * t25 + t6 * mrSges(6,3) - t9 * mrSges(5,2) + t8 * mrSges(5,1) - t7 * mrSges(6,1) - t76 * t60 - t71; -Ifges(6,6) * t59 + t43 + t44 + t45 + t64 * mrSges(6,2) + (m(6) * t64 - t65 - t66) * pkin(8); 0.2e1 * pkin(4) * mrSges(6,1) + 0.2e1 * qJ(5) * mrSges(6,3) + m(6) * (pkin(4) ^ 2 + qJ(5) ^ 2) + t76; m(6) * t3; m(6) * t7 + t24; (m(6) * pkin(8) + mrSges(6,2)) * t56; t93; m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
