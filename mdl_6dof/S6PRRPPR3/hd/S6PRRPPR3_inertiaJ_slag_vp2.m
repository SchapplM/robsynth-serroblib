% Calculate joint inertia matrix for
% S6PRRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
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
% Datum: 2018-11-23 15:09
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6PRRPPR3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR3_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR3_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPPR3_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:09:29
% EndTime: 2018-11-23 15:09:30
% DurationCPUTime: 0.55s
% Computational Cost: add. (474->185), mult. (950->231), div. (0->0), fcn. (757->8), ass. (0->68)
t52 = sin(qJ(3));
t55 = cos(qJ(3));
t91 = t52 ^ 2 + t55 ^ 2;
t49 = cos(pkin(6));
t48 = sin(pkin(6));
t53 = sin(qJ(2));
t79 = t48 * t53;
t13 = t49 * t52 + t55 * t79;
t10 = t13 ^ 2;
t87 = t55 * pkin(8);
t25 = -t55 * qJ(5) + t87;
t89 = t25 ^ 2;
t88 = -2 * mrSges(6,3);
t57 = -pkin(3) - pkin(4);
t51 = sin(qJ(6));
t54 = cos(qJ(6));
t69 = t51 ^ 2 + t54 ^ 2;
t17 = m(7) * t69;
t86 = m(6) + t17;
t85 = Ifges(7,4) * t51;
t84 = Ifges(7,4) * t54;
t83 = Ifges(7,5) * t51;
t82 = Ifges(7,6) * t54;
t11 = -t49 * t55 + t52 * t79;
t81 = t11 * t52;
t80 = t25 * t13;
t56 = cos(qJ(2));
t78 = t48 * t56;
t77 = t51 * t55;
t76 = t54 * t55;
t75 = mrSges(5,2) - mrSges(6,3);
t74 = mrSges(5,3) - mrSges(4,2);
t21 = t54 * mrSges(7,1) - t51 * mrSges(7,2);
t72 = t21 + mrSges(6,1);
t71 = Ifges(7,6) * t77 + Ifges(7,3) * t52;
t70 = t91 * pkin(8) ^ 2;
t68 = t13 * qJ(4);
t20 = -t55 * pkin(3) - t52 * qJ(4) - pkin(2);
t67 = mrSges(4,3) + t75;
t66 = -m(5) * pkin(3) - mrSges(5,1);
t65 = t69 * mrSges(7,3);
t27 = t52 * mrSges(6,1) - t55 * mrSges(6,2);
t16 = t55 * pkin(4) - t20;
t64 = pkin(8) * t81 + t13 * t87;
t22 = (pkin(8) - qJ(5)) * t52;
t5 = t52 * pkin(5) + t55 * pkin(9) + t16;
t1 = -t51 * t22 + t54 * t5;
t2 = t54 * t22 + t51 * t5;
t62 = -t51 * t1 + t54 * t2;
t3 = -t51 * t11 + t54 * t78;
t4 = t54 * t11 + t51 * t78;
t61 = t51 * t3 - t54 * t4;
t60 = -t51 * mrSges(7,1) - t54 * mrSges(7,2);
t58 = qJ(4) ^ 2;
t50 = qJ(4) + pkin(5);
t43 = -pkin(9) + t57;
t42 = t48 ^ 2;
t31 = t42 * t56 ^ 2;
t29 = -Ifges(7,1) * t51 - t84;
t28 = -Ifges(7,2) * t54 - t85;
t24 = -t55 * mrSges(4,1) + t52 * mrSges(4,2);
t23 = -t55 * mrSges(5,1) - t52 * mrSges(5,3);
t19 = t52 * mrSges(7,1) + mrSges(7,3) * t76;
t18 = -t52 * mrSges(7,2) + mrSges(7,3) * t77;
t15 = t60 * t55;
t8 = Ifges(7,5) * t52 + (-Ifges(7,1) * t54 + t85) * t55;
t7 = Ifges(7,6) * t52 + (Ifges(7,2) * t51 - t84) * t55;
t6 = [m(2) + m(7) * (t3 ^ 2 + t4 ^ 2 + t10) + m(3) * (t42 * t53 ^ 2 + t49 ^ 2 + t31) + (m(6) + m(4) + m(5)) * (t11 ^ 2 + t10 + t31); t4 * t18 + t3 * t19 + t67 * t81 + (t55 * t67 + t15) * t13 + (-t53 * mrSges(3,2) + (mrSges(3,1) - t23 - t24 + t27) * t56) * t48 + m(6) * (t22 * t11 + t16 * t78 + t80) + m(7) * (t1 * t3 + t2 * t4 + t80) + m(4) * (pkin(2) * t78 + t64) + m(5) * (-t20 * t78 + t64); -0.2e1 * pkin(2) * t24 + 0.2e1 * t1 * t19 + 0.2e1 * t25 * t15 + 0.2e1 * t16 * t27 + 0.2e1 * t2 * t18 + 0.2e1 * t20 * t23 + Ifges(3,3) + m(5) * (t20 ^ 2 + t70) + m(4) * (pkin(2) ^ 2 + t70) + m(7) * (t1 ^ 2 + t2 ^ 2 + t89) + m(6) * (t16 ^ 2 + t22 ^ 2 + t89) + (t22 * t88 + (Ifges(4,1) + Ifges(6,2) + Ifges(5,1)) * t52 + t71) * t52 + (t25 * t88 + t51 * t7 - t54 * t8 + (Ifges(4,2) + Ifges(6,1) + Ifges(5,3)) * t55 + (-Ifges(7,5) * t54 + (2 * Ifges(4,4)) + (2 * Ifges(6,4)) - (2 * Ifges(5,5))) * t52) * t55 + 0.2e1 * (mrSges(5,2) + mrSges(4,3)) * pkin(8) * t91; t61 * mrSges(7,3) + (-mrSges(4,1) - mrSges(5,1) + mrSges(6,2)) * t11 + (t72 + t74) * t13 + m(6) * (t57 * t11 + t68) + m(7) * (t50 * t13 - t43 * t61) + m(5) * (-pkin(3) * t11 + t68); t22 * mrSges(6,2) + t50 * t15 + t72 * t25 + (t43 * t18 - t2 * mrSges(7,3) - t7 / 0.2e1) * t54 + (-t43 * t19 + t1 * mrSges(7,3) - t8 / 0.2e1) * t51 + m(7) * (t50 * t25 + t43 * t62) + m(6) * (qJ(4) * t25 + t57 * t22) + (-t57 * mrSges(6,3) - pkin(3) * mrSges(5,2) - t83 / 0.2e1 - t82 / 0.2e1 + Ifges(4,5) + Ifges(6,6) + Ifges(5,4) + (-mrSges(4,1) + t66) * pkin(8)) * t52 + (-t54 * t29 / 0.2e1 + t51 * t28 / 0.2e1 + Ifges(4,6) + Ifges(6,5) - Ifges(5,6) + t75 * qJ(4) + (m(5) * qJ(4) + t74) * pkin(8)) * t55; 0.2e1 * pkin(3) * mrSges(5,1) + 0.2e1 * t57 * mrSges(6,2) + 0.2e1 * t50 * t21 - t54 * t28 - t51 * t29 + Ifges(5,2) + Ifges(4,3) + Ifges(6,3) + m(7) * (t43 ^ 2 * t69 + t50 ^ 2) + m(5) * (pkin(3) ^ 2 + t58) + m(6) * (t57 ^ 2 + t58) + 0.2e1 * (mrSges(5,3) + mrSges(6,1)) * qJ(4) - 0.2e1 * t43 * t65; -m(7) * t61 + 0.2e1 * (m(6) / 0.2e1 + m(5) / 0.2e1) * t11; t54 * t18 - t51 * t19 + m(7) * t62 + m(6) * t22 + (m(5) * pkin(8) + t75) * t52; m(6) * t57 + t17 * t43 + mrSges(6,2) - t65 + t66; m(5) + t86; m(6) * t78 + m(7) * (t54 * t3 + t51 * t4); t51 * t18 + t54 * t19 + m(7) * (t54 * t1 + t51 * t2) + m(6) * t16 + t27; 0; 0; t86; t3 * mrSges(7,1) - t4 * mrSges(7,2); t1 * mrSges(7,1) - t2 * mrSges(7,2) - Ifges(7,5) * t76 + t71; t43 * t60 - t82 - t83; t60; t21; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t6(1) t6(2) t6(4) t6(7) t6(11) t6(16); t6(2) t6(3) t6(5) t6(8) t6(12) t6(17); t6(4) t6(5) t6(6) t6(9) t6(13) t6(18); t6(7) t6(8) t6(9) t6(10) t6(14) t6(19); t6(11) t6(12) t6(13) t6(14) t6(15) t6(20); t6(16) t6(17) t6(18) t6(19) t6(20) t6(21);];
Mq  = res;
