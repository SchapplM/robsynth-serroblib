% Calculate joint inertia matrix for
% S6RPPRRP5
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRP5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP5_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP5_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP5_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP5_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP5_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:07:52
% EndTime: 2019-03-09 02:07:53
% DurationCPUTime: 0.75s
% Computational Cost: add. (460->188), mult. (819->239), div. (0->0), fcn. (569->4), ass. (0->70)
t92 = Ifges(6,5) + Ifges(7,5);
t91 = Ifges(6,6) + Ifges(7,6);
t90 = Ifges(6,3) + Ifges(7,3);
t89 = -2 * mrSges(7,3);
t50 = (qJ(2) - pkin(7));
t88 = -2 * t50;
t87 = m(6) + m(7);
t53 = sin(qJ(5));
t55 = cos(qJ(5));
t86 = t92 * t53 + t91 * t55;
t51 = (pkin(1) + qJ(3));
t85 = t51 ^ 2;
t84 = 2 * t51;
t83 = m(7) * pkin(5);
t81 = Ifges(6,4) * t53;
t80 = Ifges(6,4) * t55;
t79 = Ifges(7,4) * t53;
t78 = Ifges(7,4) * t55;
t54 = sin(qJ(4));
t77 = t50 * t54;
t56 = cos(qJ(4));
t76 = t53 * t56;
t75 = t55 * t56;
t38 = t53 * mrSges(6,2);
t22 = -t55 * mrSges(6,1) + t38;
t74 = mrSges(5,1) - t22;
t72 = -qJ(6) - pkin(8);
t15 = -t54 * mrSges(7,2) - mrSges(7,3) * t76;
t16 = -t54 * mrSges(6,2) - mrSges(6,3) * t76;
t71 = t15 + t16;
t17 = t54 * mrSges(7,1) - mrSges(7,3) * t75;
t18 = t54 * mrSges(6,1) - mrSges(6,3) * t75;
t70 = -t17 - t18;
t19 = t54 * pkin(4) - t56 * pkin(8) + t51;
t4 = t53 * t19 + t55 * t77;
t10 = mrSges(7,1) * t76 + mrSges(7,2) * t75;
t67 = t53 ^ 2 + t55 ^ 2;
t47 = t54 ^ 2;
t49 = t56 ^ 2;
t66 = t49 + t47;
t65 = qJ(6) * t56;
t63 = t67 * mrSges(6,3);
t62 = t66 * mrSges(5,3);
t37 = t53 * mrSges(7,2);
t21 = -t55 * mrSges(7,1) + t37;
t61 = t90 * t54 + t92 * t75;
t60 = -mrSges(6,1) - mrSges(7,1) - t83;
t59 = mrSges(6,1) * t53 + mrSges(6,2) * t55;
t57 = (qJ(2) ^ 2);
t45 = t50 ^ 2;
t36 = -t55 * pkin(5) - pkin(4);
t35 = t49 * t50;
t34 = t49 * t45;
t27 = Ifges(6,1) * t53 + t80;
t26 = Ifges(7,1) * t53 + t78;
t25 = Ifges(6,2) * t55 + t81;
t24 = Ifges(7,2) * t55 + t79;
t23 = t72 * t55;
t20 = t72 * t53;
t14 = (pkin(5) * t53 - t50) * t56;
t13 = t55 * t19;
t11 = t59 * t56;
t8 = Ifges(6,5) * t54 + (Ifges(6,1) * t55 - t81) * t56;
t7 = Ifges(7,5) * t54 + (Ifges(7,1) * t55 - t79) * t56;
t6 = Ifges(6,6) * t54 + (-Ifges(6,2) * t53 + t80) * t56;
t5 = Ifges(7,6) * t54 + (-Ifges(7,2) * t53 + t78) * t56;
t3 = -t53 * t77 + t13;
t2 = -t53 * t65 + t4;
t1 = -t55 * t65 + t13 + (-t50 * t53 + pkin(5)) * t54;
t9 = [-(2 * pkin(1) * mrSges(3,2)) + (mrSges(4,3) * t84) + 0.2e1 * t1 * t17 + 0.2e1 * t14 * t10 + 0.2e1 * t2 * t15 + 0.2e1 * t4 * t16 + 0.2e1 * t3 * t18 + Ifges(3,1) + Ifges(4,1) + Ifges(2,3) + (mrSges(5,1) * t84 + Ifges(5,2) * t54 + t61) * t54 + ((mrSges(5,2) * t84) + Ifges(5,1) * t56 - 0.2e1 * Ifges(5,4) * t54 + t11 * t88 + (t7 + t8) * t55 + (-t54 * t91 - t5 - t6) * t53) * t56 + m(7) * (t1 ^ 2 + t14 ^ 2 + t2 ^ 2) + m(6) * (t3 ^ 2 + t4 ^ 2 + t34) + m(5) * (t47 * t45 + t34 + t85) + (m(3) * (pkin(1) ^ 2 + t57)) + m(4) * (t57 + t85) + (2 * (mrSges(3,3) + mrSges(4,2)) * qJ(2)) + t62 * t88; -(m(3) * pkin(1)) - t54 * mrSges(5,1) - t56 * mrSges(5,2) + mrSges(3,2) - mrSges(4,3) + t70 * t55 - t71 * t53 + m(7) * (-t55 * t1 - t53 * t2) + m(6) * (-t55 * t3 - t53 * t4) + (-m(5) / 0.2e1 - m(4) / 0.2e1) * t84; t67 * t87 + m(3) + m(4) + m(5); m(4) * qJ(2) + mrSges(4,2) + (-t10 - t11) * t56 - t62 + (t70 * t53 + t71 * t55) * t54 + m(7) * (-t56 * t14 + (-t1 * t53 + t2 * t55) * t54) + m(6) * (t35 + (-t3 * t53 + t4 * t55) * t54) + m(5) * (t47 * t50 + t35); 0; m(4) + m(5) * t66 + (t67 * t47 + t49) * t87; t36 * t10 + t20 * t17 + t14 * t21 - t23 * t15 - pkin(4) * t11 + m(7) * (t20 * t1 + t36 * t14 - t23 * t2) + (t5 / 0.2e1 + t6 / 0.2e1 + t2 * mrSges(7,3) + t4 * mrSges(6,3) + (m(6) * t4 + t16) * pkin(8)) * t55 + (t7 / 0.2e1 + t8 / 0.2e1 - t1 * mrSges(7,3) - t3 * mrSges(6,3) + (-m(6) * t3 - t18) * pkin(8)) * t53 + (Ifges(5,5) + (m(6) * pkin(4) + t74) * t50 + (t26 / 0.2e1 + t27 / 0.2e1) * t55 + (-t24 / 0.2e1 - t25 / 0.2e1) * t53) * t56 + (-Ifges(5,6) - (t50 * mrSges(5,2)) + t86 / 0.2e1) * t54; m(7) * (-t20 * t55 + t23 * t53); (-t21 + t74) * t56 + (t67 * mrSges(7,3) - mrSges(5,2) + t63) * t54 + m(6) * (t67 * t54 * pkin(8) + pkin(4) * t56) + m(7) * (-t36 * t56 + (-t20 * t53 - t23 * t55) * t54); -0.2e1 * pkin(4) * t22 + 0.2e1 * t36 * t21 + Ifges(5,3) + 0.2e1 * pkin(8) * t63 + m(7) * (t20 ^ 2 + t23 ^ 2 + t36 ^ 2) + m(6) * (t67 * pkin(8) ^ 2 + pkin(4) ^ 2) + (t23 * t89 + t24 + t25) * t55 + (t20 * t89 + t26 + t27) * t53; t3 * mrSges(6,1) + t1 * mrSges(7,1) - t4 * mrSges(6,2) - t2 * mrSges(7,2) - t91 * t76 + (m(7) * t1 + t17) * pkin(5) + t61; t60 * t55 + t37 + t38; ((-mrSges(6,2) - mrSges(7,2)) * t55 + t60 * t53) * t54; t20 * mrSges(7,1) + t23 * mrSges(7,2) - t59 * pkin(8) + (m(7) * t20 - t53 * mrSges(7,3)) * pkin(5) + t86; (0.2e1 * mrSges(7,1) + t83) * pkin(5) + t90; m(7) * t14 + t10; 0; -m(7) * t56; m(7) * t36 + t21; 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t9(1) t9(2) t9(4) t9(7) t9(11) t9(16); t9(2) t9(3) t9(5) t9(8) t9(12) t9(17); t9(4) t9(5) t9(6) t9(9) t9(13) t9(18); t9(7) t9(8) t9(9) t9(10) t9(14) t9(19); t9(11) t9(12) t9(13) t9(14) t9(15) t9(20); t9(16) t9(17) t9(18) t9(19) t9(20) t9(21);];
Mq  = res;
