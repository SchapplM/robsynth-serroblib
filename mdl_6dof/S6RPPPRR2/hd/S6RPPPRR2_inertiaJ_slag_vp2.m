% Calculate joint inertia matrix for
% S6RPPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2,theta4]';
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
% Datum: 2019-03-09 01:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPPRR2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR2_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR2_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR2_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR2_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR2_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:31:39
% EndTime: 2019-03-09 01:31:40
% DurationCPUTime: 0.40s
% Computational Cost: add. (663->135), mult. (1154->191), div. (0->0), fcn. (1098->8), ass. (0->61)
t45 = sin(pkin(10));
t47 = cos(pkin(10));
t50 = sin(qJ(5));
t68 = t50 * t47;
t74 = cos(qJ(5));
t22 = t45 * t74 + t68;
t61 = t74 * t47;
t21 = t50 * t45 - t61;
t49 = sin(qJ(6));
t71 = t21 * t49;
t10 = -t22 * mrSges(7,2) + mrSges(7,3) * t71;
t51 = cos(qJ(6));
t70 = t21 * t51;
t11 = t22 * mrSges(7,1) + mrSges(7,3) * t70;
t54 = -t51 * t10 + t49 * t11;
t27 = -t51 * mrSges(7,1) + t49 * mrSges(7,2);
t86 = -m(7) * pkin(5) - mrSges(6,1) + t27;
t19 = t22 ^ 2;
t20 = t21 ^ 2;
t42 = t47 ^ 2;
t63 = t45 ^ 2 + t42;
t25 = m(5) * t63;
t85 = m(4) + t25 + m(6) * (t19 + t20);
t48 = cos(pkin(9));
t35 = -t48 * pkin(1) - pkin(2);
t31 = -qJ(4) + t35;
t75 = -pkin(7) + t31;
t16 = t75 * t45;
t6 = t50 * t16 - t61 * t75;
t84 = t6 ^ 2;
t46 = sin(pkin(9));
t32 = t46 * pkin(1) + qJ(3);
t83 = t32 ^ 2;
t82 = -2 * mrSges(6,3);
t80 = t49 / 0.2e1;
t79 = pkin(5) * t22;
t78 = t21 * t6;
t77 = t6 * t22;
t73 = Ifges(7,4) * t49;
t72 = Ifges(7,4) * t51;
t66 = -Ifges(7,5) * t70 + Ifges(7,3) * t22;
t65 = t45 * mrSges(5,1) + t47 * mrSges(5,2);
t64 = Ifges(7,5) * t49 + Ifges(7,6) * t51;
t62 = t49 ^ 2 + t51 ^ 2;
t60 = t63 * mrSges(5,3);
t59 = t62 * t21;
t17 = t21 * mrSges(6,2);
t57 = -t22 * mrSges(6,1) + t17;
t26 = t45 * pkin(4) + t32;
t5 = t21 * pkin(8) + t26 + t79;
t8 = t16 * t74 + t68 * t75;
t1 = -t49 * t8 + t51 * t5;
t2 = t49 * t5 + t51 * t8;
t56 = t1 * t49 - t2 * t51;
t55 = -mrSges(7,1) * t49 - mrSges(7,2) * t51;
t29 = Ifges(7,1) * t49 + t72;
t28 = Ifges(7,2) * t51 + t73;
t9 = mrSges(7,1) * t71 + mrSges(7,2) * t70;
t4 = Ifges(7,5) * t22 + (-Ifges(7,1) * t51 + t73) * t21;
t3 = Ifges(7,6) * t22 + (Ifges(7,2) * t49 - t72) * t21;
t7 = [Ifges(5,1) * t42 + 0.2e1 * t35 * mrSges(4,2) - 0.2e1 * t26 * t17 - 0.2e1 * t6 * t9 + 0.2e1 * t2 * t10 + 0.2e1 * t1 * t11 + Ifges(4,1) + Ifges(2,3) + Ifges(3,3) + (-0.2e1 * Ifges(5,4) * t47 + Ifges(5,2) * t45) * t45 + (0.2e1 * t26 * mrSges(6,1) + Ifges(6,2) * t22 + t8 * t82 + t66) * t22 + (t6 * t82 + Ifges(6,1) * t21 + t49 * t3 - t51 * t4 + (Ifges(7,6) * t49 + (2 * Ifges(6,4))) * t22) * t21 + m(7) * (t1 ^ 2 + t2 ^ 2 + t84) + m(6) * (t26 ^ 2 + t8 ^ 2 + t84) + m(4) * (t35 ^ 2 + t83) + m(5) * (t31 ^ 2 * t63 + t83) + m(3) * (t46 ^ 2 + t48 ^ 2) * pkin(1) ^ 2 + 0.2e1 * (mrSges(4,3) + t65) * t32 + 0.2e1 * (t48 * mrSges(3,1) - t46 * mrSges(3,2)) * pkin(1) - 0.2e1 * t31 * t60; -t22 * t9 + t54 * t21 + m(6) * (-t8 * t21 + t77) + m(7) * (t21 * t56 + t77); m(3) + m(7) * (t20 * t62 + t19) + t85; -t20 * mrSges(6,3) - t21 * t9 + mrSges(4,2) - t60 + (-mrSges(6,3) * t22 - t54) * t22 + m(7) * (-t22 * t56 + t78) + m(6) * (t22 * t8 + t78) + m(4) * t35 + t31 * t25; m(7) * (0.1e1 - t62) * t22 * t21; m(7) * (t19 * t62 + t20) + t85; t49 * t10 + t51 * t11 + m(7) * (t51 * t1 + t49 * t2) + m(6) * t26 + m(5) * t32 - t57 + t65; 0; 0; m(7) * t62 + m(5) + m(6); pkin(5) * t9 + t4 * t80 + t51 * t3 / 0.2e1 - t8 * mrSges(6,2) - t56 * mrSges(7,3) + (-t51 * t29 / 0.2e1 + t28 * t80 - Ifges(6,5)) * t21 + (t64 / 0.2e1 - Ifges(6,6)) * t22 + t86 * t6 + (-m(7) * t56 - t54) * pkin(8); t22 * t27 + m(7) * (-pkin(8) * t59 - t79) - mrSges(7,3) * t59 + t57; t86 * t21 + (-mrSges(6,2) + (m(7) * pkin(8) + mrSges(7,3)) * t62) * t22; 0; Ifges(6,3) + m(7) * (pkin(8) ^ 2 * t62 + pkin(5) ^ 2) - 0.2e1 * pkin(5) * t27 + t49 * t29 + t51 * t28 + 0.2e1 * t62 * pkin(8) * mrSges(7,3); t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,6) * t71 + t66; t9; t55 * t22; -t27; pkin(8) * t55 + t64; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t7(1) t7(2) t7(4) t7(7) t7(11) t7(16); t7(2) t7(3) t7(5) t7(8) t7(12) t7(17); t7(4) t7(5) t7(6) t7(9) t7(13) t7(18); t7(7) t7(8) t7(9) t7(10) t7(14) t7(19); t7(11) t7(12) t7(13) t7(14) t7(15) t7(20); t7(16) t7(17) t7(18) t7(19) t7(20) t7(21);];
Mq  = res;
