% Calculate joint inertia matrix for
% S6RPRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4,theta5]';
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
% Datum: 2018-11-23 15:55
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RPRPPR6_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR6_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR6_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR6_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR6_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR6_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:54:42
% EndTime: 2018-11-23 15:54:43
% DurationCPUTime: 0.72s
% Computational Cost: add. (1307->217), mult. (2343->318), div. (0->0), fcn. (2481->8), ass. (0->81)
t78 = cos(qJ(3));
t109 = t78 ^ 2;
t76 = sin(qJ(3));
t79 = -pkin(1) - pkin(7);
t91 = -qJ(4) + t79;
t54 = t91 * t76;
t72 = sin(pkin(9));
t74 = cos(pkin(9));
t86 = t91 * t78;
t34 = t54 * t72 - t74 * t86;
t108 = t34 ^ 2;
t49 = t72 * t76 - t74 * t78;
t107 = t49 ^ 2;
t106 = -2 * mrSges(5,3);
t71 = sin(pkin(10));
t103 = t71 / 0.2e1;
t60 = pkin(3) * t72 + qJ(5);
t102 = pkin(8) + t60;
t73 = cos(pkin(10));
t96 = t49 * t73;
t97 = t49 * t71;
t26 = -mrSges(6,1) * t97 - mrSges(6,2) * t96;
t75 = sin(qJ(6));
t77 = cos(qJ(6));
t53 = t71 * t77 + t73 * t75;
t19 = t53 * t49;
t82 = t71 * t75 - t73 * t77;
t21 = t82 * t49;
t7 = -t19 * mrSges(7,1) + t21 * mrSges(7,2);
t101 = t26 + t7;
t100 = Ifges(6,4) * t71;
t99 = Ifges(6,4) * t73;
t98 = t34 * t49;
t50 = t72 * t78 + t74 * t76;
t63 = t76 * pkin(3) + qJ(2);
t27 = pkin(4) * t50 + qJ(5) * t49 + t63;
t36 = t74 * t54 + t72 * t86;
t9 = t71 * t27 + t73 * t36;
t95 = Ifges(7,5) * t53 - Ifges(7,6) * t82;
t56 = -t73 * mrSges(6,1) + t71 * mrSges(6,2);
t94 = t56 - mrSges(5,1);
t93 = t71 ^ 2 + t73 ^ 2;
t92 = t76 ^ 2 + t109;
t90 = Ifges(7,5) * t21 + Ifges(7,6) * t19 + Ifges(7,3) * t50;
t62 = -pkin(3) * t74 - pkin(4);
t89 = m(4) * t92;
t88 = t93 * t60;
t87 = t92 * mrSges(4,3);
t8 = t73 * t27 - t36 * t71;
t85 = -t8 * t71 + t9 * t73;
t31 = mrSges(7,1) * t82 + t53 * mrSges(7,2);
t28 = -mrSges(6,2) * t50 + mrSges(6,3) * t97;
t29 = mrSges(6,1) * t50 + mrSges(6,3) * t96;
t84 = t73 * t28 - t71 * t29;
t83 = t74 * t49 - t72 * t50;
t80 = qJ(2) ^ 2;
t58 = Ifges(6,1) * t71 + t99;
t57 = Ifges(6,2) * t73 + t100;
t55 = -pkin(5) * t73 + t62;
t47 = t50 ^ 2;
t41 = t49 * mrSges(5,2);
t40 = t102 * t73;
t39 = t102 * t71;
t33 = Ifges(7,1) * t53 - Ifges(7,4) * t82;
t32 = Ifges(7,4) * t53 - Ifges(7,2) * t82;
t25 = -t39 * t75 + t40 * t77;
t24 = -t39 * t77 - t40 * t75;
t20 = t82 * t50;
t18 = t53 * t50;
t14 = -pkin(5) * t97 + t34;
t13 = Ifges(6,5) * t50 + (-Ifges(6,1) * t73 + t100) * t49;
t12 = Ifges(6,6) * t50 + (Ifges(6,2) * t71 - t99) * t49;
t11 = mrSges(7,1) * t50 - mrSges(7,3) * t21;
t10 = -mrSges(7,2) * t50 + mrSges(7,3) * t19;
t6 = pkin(8) * t97 + t9;
t5 = Ifges(7,1) * t21 + Ifges(7,4) * t19 + Ifges(7,5) * t50;
t4 = Ifges(7,4) * t21 + Ifges(7,2) * t19 + Ifges(7,6) * t50;
t3 = pkin(5) * t50 + pkin(8) * t96 + t8;
t2 = t3 * t75 + t6 * t77;
t1 = t3 * t77 - t6 * t75;
t15 = [Ifges(4,1) * t109 - 0.2e1 * t63 * t41 + 0.2e1 * t9 * t28 + 0.2e1 * t8 * t29 + 0.2e1 * t34 * t26 + t19 * t4 + t21 * t5 + 0.2e1 * t2 * t10 + 0.2e1 * t1 * t11 + 0.2e1 * t14 * t7 + Ifges(3,1) + Ifges(2,3) - (2 * pkin(1) * mrSges(3,2)) - 0.2e1 * t79 * t87 + (0.2e1 * t63 * mrSges(5,1) + t36 * t106 + (Ifges(6,3) + Ifges(5,2)) * t50 + t90) * t50 + (t34 * t106 + Ifges(5,1) * t49 + t71 * t12 - t73 * t13 + (-Ifges(6,5) * t73 + Ifges(6,6) * t71 + (2 * Ifges(5,4))) * t50) * t49 + m(4) * (t79 ^ 2 * t92 + t80) + m(3) * ((pkin(1) ^ 2) + t80) + m(5) * (t36 ^ 2 + t63 ^ 2 + t108) + m(6) * (t8 ^ 2 + t9 ^ 2 + t108) + m(7) * (t1 ^ 2 + t14 ^ 2 + t2 ^ 2) + (-0.2e1 * Ifges(4,4) * t78 + Ifges(4,2) * t76) * t76 + 0.2e1 * (mrSges(4,1) * t76 + mrSges(4,2) * t78 + mrSges(3,3)) * qJ(2); -m(3) * pkin(1) - t107 * mrSges(5,3) - t20 * t10 - t18 * t11 + mrSges(3,2) + t101 * t49 - t87 + (-mrSges(5,3) * t50 + t84) * t50 + m(7) * (-t1 * t18 + t14 * t49 - t2 * t20) + m(6) * (t50 * t85 + t98) + m(5) * (t36 * t50 + t98) + t79 * t89; m(3) + m(7) * (t18 ^ 2 + t20 ^ 2 + t107) + m(5) * (t47 + t107) + m(6) * (t47 * t93 + t107) + t89; m(7) * (t1 * t24 + t14 * t55 + t2 * t25) + (t79 * mrSges(4,1) + Ifges(4,5)) * t78 + (-t79 * mrSges(4,2) - Ifges(4,6)) * t76 + (-t1 * t53 - t2 * t82) * mrSges(7,3) + (-t73 * t58 / 0.2e1 - Ifges(5,5) + t57 * t103) * t49 + t94 * t34 + t84 * t60 + t85 * mrSges(6,3) + m(6) * (t34 * t62 + t60 * t85) + (m(5) * (-t34 * t74 + t36 * t72) + t83 * mrSges(5,3)) * pkin(3) + t73 * t12 / 0.2e1 + t53 * t5 / 0.2e1 + t55 * t7 + t62 * t26 - Ifges(5,6) * t50 - t82 * t4 / 0.2e1 + t14 * t31 + t19 * t32 / 0.2e1 + t21 * t33 / 0.2e1 - t36 * mrSges(5,2) + t24 * t11 + t25 * t10 + t13 * t103 + (Ifges(6,5) * t71 + Ifges(6,6) * t73 + t95) * t50 / 0.2e1; t78 * mrSges(4,1) - t76 * mrSges(4,2) + (t18 * t53 + t20 * t82) * mrSges(7,3) + (mrSges(6,3) * t93 - mrSges(5,2)) * t50 + (t31 + t94) * t49 + m(7) * (-t18 * t24 - t20 * t25 + t49 * t55) + m(6) * (t49 * t62 + t50 * t88) - m(5) * t83 * pkin(3); 0.2e1 * t55 * t31 - t82 * t32 + t53 * t33 + 0.2e1 * t62 * t56 + t73 * t57 + t71 * t58 + Ifges(4,3) + Ifges(5,3) + m(7) * (t24 ^ 2 + t25 ^ 2 + t55 ^ 2) + m(6) * (t60 ^ 2 * t93 + t62 ^ 2) + m(5) * (t72 ^ 2 + t74 ^ 2) * pkin(3) ^ 2 + 0.2e1 * (mrSges(5,1) * t74 - mrSges(5,2) * t72) * pkin(3) + 0.2e1 * (-t24 * t53 - t25 * t82) * mrSges(7,3) + 0.2e1 * mrSges(6,3) * t88; t50 * mrSges(5,1) + t53 * t10 - t82 * t11 + t71 * t28 + t73 * t29 - t41 + m(7) * (-t1 * t82 + t2 * t53) + m(6) * (t71 * t9 + t73 * t8) + m(5) * t63; m(7) * (t18 * t82 - t20 * t53); m(7) * (-t24 * t82 + t25 * t53); m(5) + m(6) * t93 + m(7) * (t53 ^ 2 + t82 ^ 2); m(6) * t34 + m(7) * t14 + t101; 0.2e1 * (m(7) / 0.2e1 + m(6) / 0.2e1) * t49; m(6) * t62 + m(7) * t55 + t31 + t56; 0; m(6) + m(7); mrSges(7,1) * t1 - mrSges(7,2) * t2 + t90; -mrSges(7,1) * t18 + mrSges(7,2) * t20; mrSges(7,1) * t24 - t25 * mrSges(7,2) + t95; -t31; 0; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t15(1) t15(2) t15(4) t15(7) t15(11) t15(16); t15(2) t15(3) t15(5) t15(8) t15(12) t15(17); t15(4) t15(5) t15(6) t15(9) t15(13) t15(18); t15(7) t15(8) t15(9) t15(10) t15(14) t15(19); t15(11) t15(12) t15(13) t15(14) t15(15) t15(20); t15(16) t15(17) t15(18) t15(19) t15(20) t15(21);];
Mq  = res;
