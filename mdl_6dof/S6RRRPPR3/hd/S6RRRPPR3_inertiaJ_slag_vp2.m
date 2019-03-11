% Calculate joint inertia matrix for
% S6RRRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6]';
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
% Datum: 2019-03-09 15:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPPR3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR3_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPPR3_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR3_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR3_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR3_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:28:26
% EndTime: 2019-03-09 15:28:28
% DurationCPUTime: 0.67s
% Computational Cost: add. (1040->213), mult. (1766->265), div. (0->0), fcn. (1657->6), ass. (0->78)
t62 = sin(qJ(6));
t65 = cos(qJ(6));
t89 = t62 ^ 2 + t65 ^ 2;
t84 = t89 * mrSges(7,3);
t113 = (mrSges(6,1) + mrSges(5,3));
t103 = -pkin(8) - pkin(7);
t67 = cos(qJ(2));
t44 = t103 * t67;
t63 = sin(qJ(3));
t66 = cos(qJ(3));
t64 = sin(qJ(2));
t85 = t103 * t64;
t26 = -t66 * t44 + t63 * t85;
t37 = t63 * t64 - t66 * t67;
t13 = t37 * qJ(5) + t26;
t112 = t13 ^ 2;
t102 = pkin(2) * t63;
t49 = qJ(4) + t102;
t111 = t49 ^ 2;
t110 = 2 * mrSges(5,1);
t109 = 2 * mrSges(6,2);
t38 = t63 * t67 + t64 * t66;
t53 = -t67 * pkin(2) - pkin(1);
t82 = -t37 * pkin(3) + t38 * qJ(4) - t53;
t10 = -pkin(4) * t37 + t82;
t108 = 0.2e1 * t10;
t107 = 0.2e1 * t13;
t41 = mrSges(7,1) * t65 - mrSges(7,2) * t62;
t106 = 0.2e1 * t41;
t105 = 0.2e1 * t53;
t104 = -pkin(4) - pkin(9);
t101 = pkin(2) * t66;
t24 = -t44 * t63 - t66 * t85;
t12 = -qJ(5) * t38 + t24;
t4 = pkin(5) * t38 + t104 * t37 + t82;
t2 = -t12 * t62 + t4 * t65;
t100 = t2 * t62;
t40 = m(7) * t89;
t99 = m(6) + t40;
t98 = Ifges(7,4) * t62;
t97 = Ifges(7,4) * t65;
t96 = t37 * t62;
t95 = t37 * t65;
t94 = t65 * mrSges(7,3);
t92 = mrSges(5,2) - mrSges(6,3);
t91 = Ifges(7,5) * t95 + Ifges(7,3) * t38;
t90 = t64 ^ 2 + t67 ^ 2;
t88 = qJ(4) * t49;
t86 = t24 ^ 2 + t26 ^ 2;
t52 = -pkin(3) - t101;
t56 = -pkin(3) + t104;
t83 = t89 * t56;
t48 = -pkin(4) + t52;
t81 = 2 * t113;
t3 = t12 * t65 + t4 * t62;
t80 = t3 * t65 - t100;
t79 = mrSges(7,1) * t62 + t65 * mrSges(7,2);
t78 = -Ifges(7,5) * t62 - Ifges(7,6) * t65;
t17 = -mrSges(7,2) * t38 - mrSges(7,3) * t96;
t18 = mrSges(7,1) * t38 - t37 * t94;
t77 = t65 * t17 - t62 * t18;
t76 = -0.2e1 * t84;
t75 = -mrSges(5,1) + mrSges(6,2) - t84;
t74 = (mrSges(4,1) * t66 - mrSges(4,2) * t63) * pkin(2);
t42 = -Ifges(7,2) * t65 - t98;
t43 = -Ifges(7,1) * t62 - t97;
t73 = -t65 * t42 - t62 * t43 + Ifges(5,2) + Ifges(4,3) + Ifges(6,3);
t8 = Ifges(7,6) * t38 + (-Ifges(7,2) * t62 + t97) * t37;
t9 = Ifges(7,5) * t38 + (Ifges(7,1) * t65 - t98) * t37;
t72 = mrSges(7,3) * t100 - t3 * t94 + t12 * mrSges(6,2) - t42 * t96 / 0.2e1 + t43 * t95 / 0.2e1 - t62 * t9 / 0.2e1 - t65 * t8 / 0.2e1 + (Ifges(5,6) - Ifges(4,6)) * t37 + (-mrSges(4,2) + mrSges(5,3)) * t26 + (-mrSges(5,1) - mrSges(4,1)) * t24 + (mrSges(6,1) + t41) * t13 + (t78 / 0.2e1 + Ifges(4,5) + Ifges(5,4)) * t38;
t69 = qJ(4) ^ 2;
t68 = -pkin(3) - pkin(4);
t61 = qJ(4) + pkin(5);
t47 = pkin(5) + t49;
t46 = -pkin(9) + t48;
t30 = t37 * mrSges(6,2);
t15 = t79 * t37;
t1 = [t64 * (Ifges(3,1) * t64 + Ifges(3,4) * t67) - 0.2e1 * pkin(1) * (-t67 * mrSges(3,1) + t64 * mrSges(3,2)) + t67 * (Ifges(3,4) * t64 + Ifges(3,2) * t67) + t30 * t108 + t15 * t107 + 0.2e1 * t3 * t17 + 0.2e1 * t2 * t18 + Ifges(2,3) + m(3) * (t90 * pkin(7) ^ 2 + pkin(1) ^ 2) + m(4) * (t53 ^ 2 + t86) + m(5) * (t82 ^ 2 + t86) + m(7) * (t2 ^ 2 + t3 ^ 2 + t112) + m(6) * (t10 ^ 2 + t12 ^ 2 + t112) + (mrSges(6,1) * t108 + mrSges(4,2) * t105 + 0.2e1 * t82 * mrSges(5,3) - 0.2e1 * t12 * mrSges(6,3) + (Ifges(4,1) + Ifges(5,1) + Ifges(6,2)) * t38 + t91) * t38 + (mrSges(4,1) * t105 - t82 * t110 + mrSges(6,3) * t107 - t62 * t8 + t65 * t9 + (Ifges(4,2) + Ifges(5,3) + Ifges(6,1)) * t37 + (-Ifges(7,6) * t62 - (2 * Ifges(4,4)) - (2 * Ifges(6,4)) + (2 * Ifges(5,5))) * t38) * t37 + 0.2e1 * t90 * pkin(7) * mrSges(3,3) + 0.2e1 * (t24 * t38 - t26 * t37) * (mrSges(5,2) + mrSges(4,3)); (t52 * mrSges(5,2) - mrSges(4,3) * t101 - t48 * mrSges(6,3) + Ifges(6,6)) * t38 + (-t64 * mrSges(3,1) - t67 * mrSges(3,2)) * pkin(7) + Ifges(3,6) * t67 + Ifges(3,5) * t64 + t47 * t15 + t72 + m(4) * (-t24 * t66 + t26 * t63) * pkin(2) + m(6) * (t12 * t48 + t13 * t49) + m(5) * (t24 * t52 + t26 * t49) + t77 * t46 + (-mrSges(4,3) * t102 - t92 * t49 - Ifges(6,5)) * t37 + m(7) * (t13 * t47 + t80 * t46); -0.2e1 * t52 * mrSges(5,1) + t48 * t109 + t47 * t106 + Ifges(3,3) + t49 * t81 + 0.2e1 * t74 + t46 * t76 + m(6) * (t48 ^ 2 + t111) + m(7) * (t89 * t46 ^ 2 + t47 ^ 2) + m(5) * (t52 ^ 2 + t111) + m(4) * (t63 ^ 2 + t66 ^ 2) * pkin(2) ^ 2 + t73; (-t92 * qJ(4) - Ifges(6,5)) * t37 + t77 * t56 + m(6) * (qJ(4) * t13 + t12 * t68) + m(5) * (-pkin(3) * t24 + qJ(4) * t26) + (-mrSges(5,2) * pkin(3) - mrSges(6,3) * t68 + Ifges(6,6)) * t38 + t61 * t15 + t72 + m(7) * (t13 * t61 + t80 * t56); (t47 + t61) * t41 + t74 + (t48 + t68) * mrSges(6,2) + (-t52 + pkin(3)) * mrSges(5,1) + m(6) * (t48 * t68 + t88) + m(7) * (t46 * t83 + t47 * t61) + m(5) * (-pkin(3) * t52 + t88) + t73 + (-t46 - t56) * t84 + t113 * (t49 + qJ(4)); pkin(3) * t110 + t68 * t109 + t61 * t106 + qJ(4) * t81 + t56 * t76 + m(7) * (t89 * t56 ^ 2 + t61 ^ 2) + m(5) * (pkin(3) ^ 2 + t69) + m(6) * (t68 ^ 2 + t69) + t73; m(5) * t24 + m(6) * t12 + m(7) * t80 + t92 * t38 + t77; m(5) * t52 + m(6) * t48 + t46 * t40 + t75; -m(5) * pkin(3) + m(6) * t68 + m(7) * t83 + t75; m(5) + t99; t38 * mrSges(6,1) + t62 * t17 + t65 * t18 + t30 + m(7) * (t2 * t65 + t3 * t62) + m(6) * t10; 0; 0; 0; t99; mrSges(7,1) * t2 - mrSges(7,2) * t3 - Ifges(7,6) * t96 + t91; -t79 * t46 + t78; -t79 * t56 + t78; -t79; t41; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
