% Calculate joint inertia matrix for
% S6RPPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
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
% Datum: 2018-11-23 15:46
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRP6_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP6_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP6_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP6_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP6_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP6_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:46:22
% EndTime: 2018-11-23 15:46:23
% DurationCPUTime: 0.51s
% Computational Cost: add. (463->183), mult. (828->229), div. (0->0), fcn. (555->4), ass. (0->72)
t93 = Ifges(7,2) + Ifges(6,3);
t52 = sin(qJ(5));
t54 = cos(qJ(5));
t71 = t52 ^ 2 + t54 ^ 2;
t49 = (qJ(2) - pkin(7));
t92 = -2 * t49;
t91 = m(7) + m(6);
t90 = (mrSges(7,2) + mrSges(6,3)) * t71;
t50 = (pkin(1) + qJ(3));
t89 = t50 ^ 2;
t88 = 2 * t50;
t86 = Ifges(6,4) * t52;
t85 = Ifges(6,4) * t54;
t84 = Ifges(7,5) * t52;
t83 = Ifges(7,5) * t54;
t53 = sin(qJ(4));
t82 = Ifges(6,6) * t53;
t81 = Ifges(7,6) * t53;
t80 = t52 * t53;
t55 = cos(qJ(4));
t79 = t52 * t55;
t18 = pkin(4) * t53 - pkin(8) * t55 + t50;
t78 = t54 * t18;
t77 = t54 * t55;
t21 = -t54 * mrSges(6,1) + mrSges(6,2) * t52;
t76 = mrSges(5,1) - t21;
t14 = -t53 * mrSges(6,2) - mrSges(6,3) * t79;
t17 = -mrSges(7,2) * t79 + t53 * mrSges(7,3);
t75 = t14 + t17;
t15 = mrSges(6,1) * t53 - mrSges(6,3) * t77;
t16 = -mrSges(7,1) * t53 + mrSges(7,2) * t77;
t74 = -t15 + t16;
t4 = t49 * t53 * t54 + t18 * t52;
t73 = t71 * pkin(8) * t53;
t72 = t71 * pkin(8) ^ 2;
t46 = t53 ^ 2;
t48 = t55 ^ 2;
t70 = t48 + t46;
t68 = t70 * mrSges(5,3);
t66 = Ifges(7,6) * t79 + (Ifges(7,4) + Ifges(6,5)) * t77 + t93 * t53;
t1 = qJ(6) * t53 + t4;
t2 = -t78 + (t49 * t52 - pkin(5)) * t53;
t65 = t1 * t54 + t2 * t52;
t3 = -t49 * t80 + t78;
t64 = -t3 * t52 + t4 * t54;
t63 = t52 * mrSges(6,1) + t54 * mrSges(6,2);
t20 = -t54 * mrSges(7,1) - t52 * mrSges(7,3);
t62 = t52 * mrSges(7,1) - t54 * mrSges(7,3);
t61 = -pkin(5) * t54 - qJ(6) * t52;
t60 = -pkin(5) * t52 + qJ(6) * t54;
t59 = t52 * t74 + t54 * t75;
t58 = m(7) * t60 - t62 - t63;
t56 = (qJ(2) ^ 2);
t44 = t49 ^ 2;
t41 = Ifges(7,4) * t52;
t40 = Ifges(6,5) * t52;
t38 = Ifges(6,6) * t54;
t35 = t48 * t49;
t34 = t48 * t44;
t25 = Ifges(6,1) * t52 + t85;
t24 = Ifges(7,1) * t52 - t83;
t23 = Ifges(6,2) * t54 + t86;
t22 = -Ifges(7,3) * t54 + t84;
t19 = -pkin(4) + t61;
t12 = t63 * t55;
t11 = t62 * t55;
t9 = Ifges(6,5) * t53 + (Ifges(6,1) * t54 - t86) * t55;
t8 = Ifges(7,4) * t53 + (Ifges(7,1) * t54 + t84) * t55;
t7 = t82 + (-Ifges(6,2) * t52 + t85) * t55;
t6 = t81 + (Ifges(7,3) * t52 + t83) * t55;
t5 = (-t49 - t60) * t55;
t10 = [-(2 * pkin(1) * mrSges(3,2)) + (mrSges(4,3) * t88) + 0.2e1 * t1 * t17 + 0.2e1 * t5 * t11 + 0.2e1 * t4 * t14 + 0.2e1 * t3 * t15 + 0.2e1 * t2 * t16 + Ifges(3,1) + Ifges(4,1) + Ifges(2,3) + (mrSges(5,1) * t88 + Ifges(5,2) * t53 + t66) * t53 + ((mrSges(5,2) * t88) + Ifges(5,1) * t55 - 0.2e1 * Ifges(5,4) * t53 + t12 * t92 + (t8 + t9) * t54 + (t6 - t7 - t82) * t52) * t55 + m(4) * (t56 + t89) + (m(3) * (pkin(1) ^ 2 + t56)) + m(5) * (t44 * t46 + t34 + t89) + m(6) * (t3 ^ 2 + t4 ^ 2 + t34) + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) + (2 * (mrSges(3,3) + mrSges(4,2)) * qJ(2)) + t68 * t92; -(m(3) * pkin(1)) - t53 * mrSges(5,1) - t55 * mrSges(5,2) + mrSges(3,2) - mrSges(4,3) + t74 * t54 - t75 * t52 + m(7) * (-t1 * t52 + t2 * t54) + m(6) * (-t3 * t54 - t4 * t52) + (-m(5) / 0.2e1 - m(4) / 0.2e1) * t88; t71 * t91 + m(3) + m(4) + m(5); m(4) * qJ(2) + mrSges(4,2) + (-t11 - t12) * t55 - t68 + t59 * t53 + m(7) * (-t55 * t5 + t53 * t65) + m(6) * (t53 * t64 + t35) + m(5) * (t46 * t49 + t35); 0; m(4) + m(5) * t70 + (t46 * t71 + t48) * t91; -pkin(4) * t12 + t5 * t20 + (m(7) * t5 + t11) * t19 + (t40 / 0.2e1 + t38 / 0.2e1 - Ifges(5,6) + t41 / 0.2e1 - (t49 * mrSges(5,2))) * t53 + (-t81 / 0.2e1 - t6 / 0.2e1 + t7 / 0.2e1 + t1 * mrSges(7,2) + t4 * mrSges(6,3)) * t54 + (t8 / 0.2e1 + t9 / 0.2e1 + t2 * mrSges(7,2) - t3 * mrSges(6,3)) * t52 + (m(6) * t64 + m(7) * t65 + t59) * pkin(8) + (Ifges(5,5) + (t24 / 0.2e1 + t25 / 0.2e1) * t54 + (t22 / 0.2e1 - t23 / 0.2e1) * t52 + (m(6) * pkin(4) + t76) * t49) * t55; 0; (-t20 + t76) * t55 + m(7) * (-t19 * t55 + t73) + m(6) * (pkin(4) * t55 + t73) + (-mrSges(5,2) + t90) * t53; -0.2e1 * pkin(4) * t21 + 0.2e1 * t19 * t20 + Ifges(5,3) + (-t22 + t23) * t54 + (t25 + t24) * t52 + m(7) * (t19 ^ 2 + t72) + m(6) * (pkin(4) ^ 2 + t72) + 0.2e1 * pkin(8) * t90; -Ifges(6,6) * t79 - pkin(5) * t16 + m(7) * (-pkin(5) * t2 + qJ(6) * t1) + qJ(6) * t17 + t1 * mrSges(7,3) - t4 * mrSges(6,2) + t3 * mrSges(6,1) - t2 * mrSges(7,1) + t66; m(7) * t61 + t20 + t21; t58 * t53; mrSges(7,2) * t60 - Ifges(7,6) * t54 + pkin(8) * t58 + t38 + t40 + t41; 0.2e1 * pkin(5) * mrSges(7,1) + 0.2e1 * qJ(6) * mrSges(7,3) + m(7) * (pkin(5) ^ 2 + qJ(6) ^ 2) + t93; m(7) * t2 + t16; m(7) * t54; m(7) * t80; (m(7) * pkin(8) + mrSges(7,2)) * t52; -m(7) * pkin(5) - mrSges(7,1); m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t10(1) t10(2) t10(4) t10(7) t10(11) t10(16); t10(2) t10(3) t10(5) t10(8) t10(12) t10(17); t10(4) t10(5) t10(6) t10(9) t10(13) t10(18); t10(7) t10(8) t10(9) t10(10) t10(14) t10(19); t10(11) t10(12) t10(13) t10(14) t10(15) t10(20); t10(16) t10(17) t10(18) t10(19) t10(20) t10(21);];
Mq  = res;
