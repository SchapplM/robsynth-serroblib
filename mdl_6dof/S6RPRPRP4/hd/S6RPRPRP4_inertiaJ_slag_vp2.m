% Calculate joint inertia matrix for
% S6RPRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
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
% Datum: 2018-11-23 15:59
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRP4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP4_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP4_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP4_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP4_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP4_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:58:45
% EndTime: 2018-11-23 15:58:46
% DurationCPUTime: 0.64s
% Computational Cost: add. (601->194), mult. (1083->248), div. (0->0), fcn. (766->6), ass. (0->79)
t55 = sin(pkin(9));
t39 = t55 * pkin(1) + pkin(7);
t104 = pkin(4) + t39;
t57 = sin(qJ(5));
t59 = cos(qJ(5));
t81 = -t57 ^ 2 - t59 ^ 2;
t98 = (mrSges(7,2) + mrSges(6,3)) * t81;
t103 = Ifges(7,2) + Ifges(6,3);
t58 = sin(qJ(3));
t51 = t58 ^ 2;
t60 = cos(qJ(3));
t53 = t60 ^ 2;
t102 = t51 + t53;
t80 = m(6) / 0.2e1 + m(7) / 0.2e1;
t101 = 0.2e1 * t80;
t100 = 0.2e1 * t102;
t68 = pkin(5) * t59 + qJ(6) * t57;
t69 = t59 * mrSges(7,1) + t57 * mrSges(7,3);
t70 = t59 * mrSges(6,1) - t57 * mrSges(6,2);
t64 = m(7) * t68 + t69 + t70;
t88 = t60 * mrSges(6,3);
t23 = -t58 * mrSges(6,2) - t59 * t88;
t89 = t60 * mrSges(7,2);
t24 = t58 * mrSges(7,3) - t59 * t89;
t86 = t23 + t24;
t21 = t58 * mrSges(6,1) + t57 * t88;
t22 = -t58 * mrSges(7,1) - t57 * t89;
t87 = t21 - t22;
t97 = -t86 * t57 - t87 * t59;
t61 = -pkin(3) - pkin(8);
t96 = t60 * pkin(3);
t17 = t104 * t58;
t56 = cos(pkin(9));
t40 = -t56 * pkin(1) - pkin(2);
t43 = t58 * qJ(4);
t76 = t40 - t43;
t7 = t61 * t60 + t76;
t4 = t57 * t17 + t59 * t7;
t95 = Ifges(6,4) * t57;
t94 = Ifges(6,4) * t59;
t93 = Ifges(7,5) * t57;
t92 = Ifges(7,5) * t59;
t91 = Ifges(6,6) * t57;
t90 = t59 * t60;
t27 = t57 * mrSges(6,1) + t59 * mrSges(6,2);
t85 = t27 + mrSges(5,3);
t84 = t102 * t39 ^ 2;
t18 = t104 * t60;
t83 = t81 * t61 ^ 2;
t79 = Ifges(7,6) * t90 + t103 * t58;
t78 = -m(5) * pkin(3) + mrSges(5,2);
t77 = t81 * t60;
t73 = t61 * t77;
t1 = t58 * qJ(6) + t4;
t3 = t59 * t17 - t57 * t7;
t2 = -t58 * pkin(5) - t3;
t72 = t57 * t1 - t59 * t2;
t71 = t59 * t3 + t57 * t4;
t67 = -0.2e1 * t80 * t81;
t66 = -Ifges(6,6) * t59 + (-Ifges(7,4) - Ifges(6,5)) * t57;
t62 = qJ(4) ^ 2;
t48 = Ifges(7,4) * t59;
t47 = Ifges(6,5) * t59;
t45 = Ifges(7,6) * t57;
t31 = Ifges(6,1) * t59 - t95;
t30 = Ifges(7,1) * t59 + t93;
t29 = -Ifges(6,2) * t57 + t94;
t28 = Ifges(7,3) * t57 + t92;
t26 = t57 * mrSges(7,1) - t59 * mrSges(7,3);
t25 = t57 * pkin(5) - t59 * qJ(6) + qJ(4);
t16 = t70 * t60;
t15 = t69 * t60;
t14 = t76 - t96;
t11 = Ifges(6,5) * t58 + (-Ifges(6,1) * t57 - t94) * t60;
t10 = Ifges(7,4) * t58 + (-Ifges(7,1) * t57 + t92) * t60;
t9 = Ifges(6,6) * t58 + (-Ifges(6,2) * t59 - t95) * t60;
t8 = Ifges(7,6) * t58 + (Ifges(7,3) * t59 - t93) * t60;
t5 = t68 * t60 + t18;
t6 = [0.2e1 * t1 * t24 + 0.2e1 * t5 * t15 + 0.2e1 * t18 * t16 + 0.2e1 * t2 * t22 + 0.2e1 * t3 * t21 + 0.2e1 * t4 * t23 + Ifges(2,3) + Ifges(3,3) + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) + m(6) * (t18 ^ 2 + t3 ^ 2 + t4 ^ 2) + m(5) * (t14 ^ 2 + t84) + m(4) * (t40 ^ 2 + t84) + (0.2e1 * t40 * mrSges(4,2) - 0.2e1 * t14 * mrSges(5,3) + (Ifges(5,2) + Ifges(4,1)) * t58 + t79) * t58 + (-0.2e1 * t40 * mrSges(4,1) + 0.2e1 * t14 * mrSges(5,2) + (Ifges(5,3) + Ifges(4,2)) * t60 + (t8 - t9) * t59 + (-t10 - t11) * t57 + ((2 * Ifges(4,4)) + (2 * Ifges(5,6)) + t66) * t58) * t60 + (mrSges(5,1) + mrSges(4,3)) * t39 * t100 + (0.2e1 * t56 * mrSges(3,1) - 0.2e1 * t55 * mrSges(3,2) + m(3) * (t55 ^ 2 + t56 ^ 2) * pkin(1)) * pkin(1); (t15 + t16) * t58 + t97 * t60 + m(7) * (t5 * t58 - t72 * t60) + m(6) * (t18 * t58 - t71 * t60); m(3) + (m(5) / 0.2e1 + m(4) / 0.2e1) * t100 + (-t81 * t53 + t51) * t101; qJ(4) * t16 + t25 * t15 + t18 * t27 + t5 * t26 + (t10 / 0.2e1 + t11 / 0.2e1 + t2 * mrSges(7,2) - t3 * mrSges(6,3) + t87 * t61) * t59 + (t8 / 0.2e1 - t9 / 0.2e1 - t1 * mrSges(7,2) - t4 * mrSges(6,3) + t86 * t61) * t57 + m(6) * (qJ(4) * t18 + t71 * t61) + m(7) * (t25 * t5 + t72 * t61) + (-t91 / 0.2e1 + t47 / 0.2e1 + t48 / 0.2e1 + t45 / 0.2e1 - Ifges(5,4) + Ifges(4,5) - pkin(3) * mrSges(5,1) + (-mrSges(4,1) + t78) * t39) * t58 + (qJ(4) * mrSges(5,1) - Ifges(5,5) + Ifges(4,6) + (t28 / 0.2e1 - t29 / 0.2e1) * t59 + (-t30 / 0.2e1 - t31 / 0.2e1) * t57 + (m(5) * qJ(4) - mrSges(4,2) + mrSges(5,3)) * t39) * t60; (-mrSges(4,2) + t26 + t85) * t58 + (mrSges(4,1) - mrSges(5,2) - t98) * t60 + m(7) * (t25 * t58 + t73) + m(5) * (t43 + t96) + m(6) * (t43 + t73); -0.2e1 * pkin(3) * mrSges(5,2) + 0.2e1 * t25 * t26 + Ifges(5,1) + Ifges(4,3) + (t30 + t31) * t59 + (t28 - t29) * t57 + 0.2e1 * t85 * qJ(4) + m(7) * (t25 ^ 2 - t83) + m(6) * (t62 - t83) + m(5) * (pkin(3) ^ 2 + t62) + 0.2e1 * t61 * t98; (m(5) * t39 + mrSges(5,1)) * t58 + m(7) * t72 + m(6) * t71 - t97; -m(5) * t60 + t77 * t101; t61 * t67 + t78 + t98; m(5) + t67; -pkin(5) * t22 + qJ(6) * t24 + m(7) * (-pkin(5) * t2 + qJ(6) * t1) + t1 * mrSges(7,3) - t4 * mrSges(6,2) + t3 * mrSges(6,1) - t2 * mrSges(7,1) + t66 * t60 + t79; -t64 * t60; -t68 * mrSges(7,2) + t64 * t61 + t45 + t47 + t48 - t91; t64; 0.2e1 * pkin(5) * mrSges(7,1) + 0.2e1 * qJ(6) * mrSges(7,3) + m(7) * (pkin(5) ^ 2 + qJ(6) ^ 2) + t103; m(7) * t2 + t22; m(7) * t90; (-m(7) * t61 + mrSges(7,2)) * t59; -m(7) * t59; -m(7) * pkin(5) - mrSges(7,1); m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t6(1) t6(2) t6(4) t6(7) t6(11) t6(16); t6(2) t6(3) t6(5) t6(8) t6(12) t6(17); t6(4) t6(5) t6(6) t6(9) t6(13) t6(18); t6(7) t6(8) t6(9) t6(10) t6(14) t6(19); t6(11) t6(12) t6(13) t6(14) t6(15) t6(20); t6(16) t6(17) t6(18) t6(19) t6(20) t6(21);];
Mq  = res;
