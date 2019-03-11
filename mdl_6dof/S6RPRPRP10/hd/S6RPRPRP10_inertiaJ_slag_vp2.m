% Calculate joint inertia matrix for
% S6RPRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5]';
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
% Datum: 2019-03-09 03:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRP10_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP10_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPRP10_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP10_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP10_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP10_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:30:58
% EndTime: 2019-03-09 03:31:00
% DurationCPUTime: 0.68s
% Computational Cost: add. (580->199), mult. (1011->243), div. (0->0), fcn. (679->4), ass. (0->78)
t105 = Ifges(7,2) + Ifges(6,3);
t57 = sin(qJ(5));
t59 = cos(qJ(5));
t82 = -t57 ^ 2 - t59 ^ 2;
t58 = sin(qJ(3));
t52 = t58 ^ 2;
t60 = cos(qJ(3));
t54 = t60 ^ 2;
t81 = t54 + t52;
t100 = (-mrSges(6,3) - mrSges(7,2)) * t82;
t104 = (mrSges(5,1) + mrSges(4,3)) * t81;
t80 = m(6) / 0.2e1 + m(7) / 0.2e1;
t103 = 0.2e1 * t80;
t69 = pkin(5) * t59 + qJ(6) * t57;
t70 = -t59 * mrSges(7,1) - t57 * mrSges(7,3);
t71 = -t59 * mrSges(6,1) + t57 * mrSges(6,2);
t65 = m(7) * t69 - t70 - t71;
t88 = t58 * t59;
t21 = -t60 * mrSges(6,2) + mrSges(6,3) * t88;
t22 = mrSges(7,2) * t88 + t60 * mrSges(7,3);
t86 = t21 + t22;
t89 = t57 * t58;
t19 = t60 * mrSges(6,1) - mrSges(6,3) * t89;
t20 = -t60 * mrSges(7,1) + mrSges(7,2) * t89;
t87 = t19 - t20;
t99 = -t86 * t57 - t87 * t59;
t24 = t58 * pkin(3) - t60 * qJ(4) + qJ(2);
t98 = -0.2e1 * t24;
t97 = 0.2e1 * qJ(2);
t96 = m(7) * t59;
t95 = Ifges(6,4) * t57;
t94 = Ifges(6,4) * t59;
t93 = Ifges(7,5) * t57;
t92 = Ifges(7,5) * t59;
t91 = Ifges(6,6) * t57;
t90 = Ifges(7,6) * t60;
t12 = t58 * pkin(8) + t24;
t62 = -pkin(1) - pkin(7);
t26 = (pkin(4) - t62) * t60;
t4 = t59 * t12 + t57 * t26;
t28 = t57 * mrSges(6,1) + t59 * mrSges(6,2);
t85 = t28 + mrSges(5,3);
t61 = -pkin(3) - pkin(8);
t84 = t82 * t61 ^ 2;
t83 = t81 * t62 ^ 2;
t79 = -m(5) * pkin(3) + mrSges(5,2);
t78 = t82 * t60;
t75 = t61 * t78;
t74 = Ifges(6,6) * t88 + (Ifges(7,4) + Ifges(6,5)) * t89 + t105 * t60;
t1 = t60 * qJ(6) + t4;
t3 = -t57 * t12 + t59 * t26;
t2 = -t60 * pkin(5) - t3;
t73 = t57 * t1 - t59 * t2;
t72 = t59 * t3 + t57 * t4;
t68 = 0.2e1 * (m(5) / 0.2e1 + m(4) / 0.2e1) * t81;
t67 = -0.2e1 * t80 * t82;
t64 = qJ(2) ^ 2;
t63 = qJ(4) ^ 2;
t49 = Ifges(7,4) * t59;
t48 = Ifges(6,5) * t59;
t46 = Ifges(7,6) * t57;
t44 = qJ(4) * t58;
t41 = t58 * t62;
t32 = Ifges(6,1) * t59 - t95;
t31 = Ifges(7,1) * t59 + t93;
t30 = -Ifges(6,2) * t57 + t94;
t29 = Ifges(7,3) * t57 + t92;
t27 = t57 * mrSges(7,1) - t59 * mrSges(7,3);
t25 = -t58 * pkin(4) + t41;
t23 = t57 * pkin(5) - t59 * qJ(6) + qJ(4);
t14 = t71 * t58;
t13 = t70 * t58;
t9 = Ifges(6,5) * t60 + (Ifges(6,1) * t57 + t94) * t58;
t8 = Ifges(7,4) * t60 + (Ifges(7,1) * t57 - t92) * t58;
t7 = Ifges(6,6) * t60 + (Ifges(6,2) * t59 + t95) * t58;
t6 = t90 + (-Ifges(7,3) * t59 + t93) * t58;
t5 = t41 + (-pkin(4) - t69) * t58;
t10 = [-(2 * pkin(1) * mrSges(3,2)) + mrSges(3,3) * t97 + 0.2e1 * t1 * t22 + 0.2e1 * t5 * t13 + 0.2e1 * t25 * t14 + 0.2e1 * t3 * t19 + 0.2e1 * t2 * t20 + 0.2e1 * t4 * t21 + Ifges(3,1) + Ifges(2,3) + (mrSges(4,2) * t97 + mrSges(5,3) * t98 + (Ifges(4,1) + Ifges(5,2)) * t60 + t74) * t60 + (mrSges(4,1) * t97 + mrSges(5,2) * t98 + (Ifges(5,3) + Ifges(4,2)) * t58 + 0.2e1 * (-Ifges(4,4) - Ifges(5,6)) * t60 + (t8 + t9) * t57 + (-t6 + t7 - t90) * t59) * t58 + m(4) * (t64 + t83) + m(3) * ((pkin(1) ^ 2) + t64) + m(5) * (t24 ^ 2 + t83) + m(6) * (t25 ^ 2 + t3 ^ 2 + t4 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) - 0.2e1 * t62 * t104; -m(3) * pkin(1) + mrSges(3,2) + (t13 + t14) * t58 + t99 * t60 + m(7) * (t58 * t5 - t73 * t60) + m(6) * (t58 * t25 - t72 * t60) + t62 * t68 - t104; m(3) + t68 + (-t82 * t54 + t52) * t103; qJ(4) * t14 + t23 * t13 + t25 * t28 + t5 * t27 + (t8 / 0.2e1 + t9 / 0.2e1 + t2 * mrSges(7,2) - t3 * mrSges(6,3) + t87 * t61) * t59 + (-t1 * mrSges(7,2) - t4 * mrSges(6,3) + t6 / 0.2e1 - t7 / 0.2e1 + t86 * t61) * t57 + m(7) * (t23 * t5 + t73 * t61) + m(6) * (qJ(4) * t25 + t72 * t61) + (-t91 / 0.2e1 + t48 / 0.2e1 + t49 / 0.2e1 + t46 / 0.2e1 - Ifges(5,4) + Ifges(4,5) - pkin(3) * mrSges(5,1) + (mrSges(4,1) - t79) * t62) * t60 + (-qJ(4) * mrSges(5,1) + Ifges(5,5) - Ifges(4,6) + (-t29 / 0.2e1 + t30 / 0.2e1) * t59 + (t31 / 0.2e1 + t32 / 0.2e1) * t57 + (m(5) * qJ(4) - mrSges(4,2) + mrSges(5,3)) * t62) * t58; (-mrSges(4,2) + t27 + t85) * t58 + (mrSges(4,1) - mrSges(5,2) + t100) * t60 + m(6) * (t44 + t75) + m(7) * (t23 * t58 + t75) + m(5) * (pkin(3) * t60 + t44); -0.2e1 * pkin(3) * mrSges(5,2) + 0.2e1 * t23 * t27 + Ifges(5,1) + Ifges(4,3) + (t31 + t32) * t59 + (-t30 + t29) * t57 + 0.2e1 * t85 * qJ(4) + m(7) * (t23 ^ 2 - t84) + m(6) * (t63 - t84) + m(5) * (pkin(3) ^ 2 + t63) - 0.2e1 * t61 * t100; (-m(5) * t62 + mrSges(5,1)) * t60 + m(7) * t73 + m(6) * t72 - t99; -m(5) * t60 + t78 * t103; t61 * t67 - t100 + t79; m(5) + t67; -pkin(5) * t20 + m(7) * (-pkin(5) * t2 + qJ(6) * t1) + qJ(6) * t22 + t1 * mrSges(7,3) - t4 * mrSges(6,2) + t3 * mrSges(6,1) - Ifges(7,6) * t88 - t2 * mrSges(7,1) + t74; -t65 * t60; -t69 * mrSges(7,2) + t65 * t61 + t46 + t48 + t49 - t91; t65; 0.2e1 * pkin(5) * mrSges(7,1) + 0.2e1 * qJ(6) * mrSges(7,3) + m(7) * (pkin(5) ^ 2 + qJ(6) ^ 2) + t105; m(7) * t2 + t20; t60 * t96; (-m(7) * t61 + mrSges(7,2)) * t59; -t96; -m(7) * pkin(5) - mrSges(7,1); m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t10(1) t10(2) t10(4) t10(7) t10(11) t10(16); t10(2) t10(3) t10(5) t10(8) t10(12) t10(17); t10(4) t10(5) t10(6) t10(9) t10(13) t10(18); t10(7) t10(8) t10(9) t10(10) t10(14) t10(19); t10(11) t10(12) t10(13) t10(14) t10(15) t10(20); t10(16) t10(17) t10(18) t10(19) t10(20) t10(21);];
Mq  = res;
