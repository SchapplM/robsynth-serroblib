% Calculate joint inertia matrix for
% S6RPRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta5]';
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
% Datum: 2018-11-23 15:54
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RPRPPR5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR5_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR5_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR5_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR5_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR5_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:54:13
% EndTime: 2018-11-23 15:54:14
% DurationCPUTime: 0.62s
% Computational Cost: add. (1290->204), mult. (2375->279), div. (0->0), fcn. (2578->8), ass. (0->83)
t76 = cos(pkin(9));
t71 = t76 ^ 2;
t104 = cos(qJ(3));
t74 = sin(pkin(9));
t79 = sin(qJ(3));
t48 = -t104 * t76 + t74 * t79;
t51 = t104 * t74 + t79 * t76;
t64 = -pkin(2) * t76 - pkin(1);
t84 = -qJ(4) * t51 + t64;
t27 = pkin(3) * t48 + t84;
t108 = -0.2e1 * t27;
t107 = 0.2e1 * t64;
t75 = cos(pkin(10));
t106 = t75 / 0.2e1;
t97 = pkin(3) + qJ(5);
t105 = -pkin(8) - t97;
t73 = sin(pkin(10));
t103 = Ifges(6,4) * t73;
t102 = Ifges(6,4) * t75;
t101 = t48 * t73;
t100 = t48 * t75;
t98 = mrSges(5,2) - mrSges(4,1);
t96 = pkin(7) + qJ(2);
t16 = t97 * t48 + t84;
t55 = t96 * t74;
t57 = t96 * t76;
t35 = t104 * t55 + t57 * t79;
t19 = pkin(4) * t51 + t35;
t6 = t75 * t16 + t73 * t19;
t56 = t73 * mrSges(6,1) + t75 * mrSges(6,2);
t95 = t74 ^ 2 + t71;
t94 = t73 ^ 2 + t75 ^ 2;
t78 = sin(qJ(6));
t80 = cos(qJ(6));
t50 = -t73 * t78 + t75 * t80;
t85 = t80 * t73 + t78 * t75;
t93 = t50 ^ 2 + t85 ^ 2;
t24 = t50 * t48;
t25 = t85 * t48;
t92 = Ifges(7,5) * t25 + Ifges(7,6) * t24 + Ifges(7,3) * t51;
t37 = t104 * t57 - t79 * t55;
t91 = t35 ^ 2 + t37 ^ 2;
t52 = m(6) * t94;
t90 = t94 * mrSges(6,3);
t89 = -t76 * mrSges(3,1) + t74 * mrSges(3,2);
t9 = -t24 * mrSges(7,1) + t25 * mrSges(7,2);
t26 = -mrSges(6,1) * t100 + mrSges(6,2) * t101;
t18 = t75 * t19;
t3 = pkin(5) * t51 + t18 + (-pkin(8) * t48 - t16) * t73;
t4 = pkin(8) * t100 + t6;
t1 = t3 * t80 - t4 * t78;
t2 = t3 * t78 + t4 * t80;
t88 = t1 * t50 + t2 * t85;
t5 = -t16 * t73 + t18;
t87 = t5 * t75 + t6 * t73;
t32 = mrSges(7,1) * t85 + t50 * mrSges(7,2);
t53 = t105 * t73;
t54 = t105 * t75;
t30 = -t53 * t78 + t54 * t80;
t31 = t53 * t80 + t54 * t78;
t86 = t30 * t50 + t31 * t85;
t83 = m(7) * t93 + m(5) + t52;
t81 = qJ(4) ^ 2;
t60 = pkin(5) * t73 + qJ(4);
t59 = Ifges(6,1) * t75 - t103;
t58 = -Ifges(6,2) * t73 + t102;
t44 = Ifges(7,5) * t50;
t43 = Ifges(7,6) * t85;
t41 = t51 * mrSges(5,3);
t40 = t51 * mrSges(4,2);
t34 = Ifges(7,1) * t50 - Ifges(7,4) * t85;
t33 = Ifges(7,4) * t50 - Ifges(7,2) * t85;
t29 = -mrSges(6,2) * t51 + mrSges(6,3) * t100;
t28 = mrSges(6,1) * t51 - mrSges(6,3) * t101;
t20 = -t48 * pkin(4) + t37;
t15 = Ifges(6,5) * t51 + (Ifges(6,1) * t73 + t102) * t48;
t14 = Ifges(6,6) * t51 + (Ifges(6,2) * t75 + t103) * t48;
t12 = mrSges(7,1) * t51 - mrSges(7,3) * t25;
t11 = -mrSges(7,2) * t51 + mrSges(7,3) * t24;
t10 = (-pkin(5) * t75 - pkin(4)) * t48 + t37;
t8 = Ifges(7,1) * t25 + Ifges(7,4) * t24 + Ifges(7,5) * t51;
t7 = Ifges(7,4) * t25 + Ifges(7,2) * t24 + Ifges(7,6) * t51;
t13 = [Ifges(3,2) * t71 - 0.2e1 * pkin(1) * t89 + t40 * t107 + t41 * t108 + t24 * t7 + t25 * t8 + 0.2e1 * t20 * t26 + 0.2e1 * t5 * t28 + 0.2e1 * t6 * t29 + 0.2e1 * t10 * t9 + 0.2e1 * t2 * t11 + 0.2e1 * t1 * t12 + Ifges(2,3) + (Ifges(3,1) * t74 + 0.2e1 * Ifges(3,4) * t76) * t74 + (mrSges(4,1) * t107 + mrSges(5,2) * t108 + t75 * t14 + t73 * t15 + (Ifges(4,2) + Ifges(5,3)) * t48) * t48 + m(3) * (t95 * qJ(2) ^ 2 + pkin(1) ^ 2) + m(4) * (t64 ^ 2 + t91) + m(7) * (t1 ^ 2 + t10 ^ 2 + t2 ^ 2) + m(6) * (t20 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(5) * (t27 ^ 2 + t91) + ((Ifges(6,3) + Ifges(4,1) + Ifges(5,2)) * t51 + (Ifges(6,5) * t73 + Ifges(6,6) * t75 - (2 * Ifges(4,4)) - (2 * Ifges(5,6))) * t48 + t92) * t51 + 0.2e1 * t95 * qJ(2) * mrSges(3,3) + 0.2e1 * (t35 * t51 - t37 * t48) * (mrSges(5,1) + mrSges(4,3)); -m(3) * pkin(1) + t50 * t11 - t85 * t12 - t73 * t28 + t75 * t29 + t40 - t41 - t98 * t48 + m(7) * (-t1 * t85 + t2 * t50) + m(6) * (-t5 * t73 + t6 * t75) + m(5) * t27 + m(4) * t64 + t89; m(3) + m(4) + t83; t20 * t56 + t60 * t9 + t50 * t8 / 0.2e1 - t85 * t7 / 0.2e1 + qJ(4) * t26 + t30 * t12 + t31 * t11 + t10 * t32 + t24 * t33 / 0.2e1 + t25 * t34 / 0.2e1 + (-mrSges(4,2) + mrSges(5,3)) * t37 + t98 * t35 - t88 * mrSges(7,3) + (-t5 * mrSges(6,3) - t97 * t28 + t15 / 0.2e1) * t75 + (-t97 * t29 - t6 * mrSges(6,3) - t14 / 0.2e1) * t73 + (-pkin(3) * mrSges(5,1) + Ifges(6,5) * t106 - Ifges(6,6) * t73 / 0.2e1 + t44 / 0.2e1 - t43 / 0.2e1 + Ifges(4,5) - Ifges(5,4)) * t51 + m(6) * (qJ(4) * t20 - t87 * t97) + m(7) * (t1 * t30 + t10 * t60 + t2 * t31) + m(5) * (-pkin(3) * t35 + qJ(4) * t37) + (t58 * t106 + t73 * t59 / 0.2e1 - qJ(4) * mrSges(5,1) + Ifges(5,5) - Ifges(4,6)) * t48; m(7) * (-t30 * t85 + t31 * t50); -0.2e1 * pkin(3) * mrSges(5,2) + 0.2e1 * t60 * t32 - t85 * t33 + t50 * t34 - t73 * t58 + t75 * t59 + Ifges(5,1) + Ifges(4,3) + m(7) * (t30 ^ 2 + t31 ^ 2 + t60 ^ 2) + m(6) * (t94 * t97 ^ 2 + t81) + m(5) * (pkin(3) ^ 2 + t81) + 0.2e1 * (mrSges(5,3) + t56) * qJ(4) - 0.2e1 * t86 * mrSges(7,3) + 0.2e1 * t97 * t90; m(5) * t35 + m(6) * t87 + m(7) * t88 + t51 * mrSges(5,1) + t11 * t85 + t50 * t12 + t75 * t28 + t73 * t29; 0; -m(5) * pkin(3) + m(7) * t86 - t93 * mrSges(7,3) - t52 * t97 + mrSges(5,2) - t90; t83; m(6) * t20 + m(7) * t10 + t26 + t9; 0; m(6) * qJ(4) + m(7) * t60 + t32 + t56; 0; m(6) + m(7); mrSges(7,1) * t1 - mrSges(7,2) * t2 + t92; -t32; mrSges(7,1) * t30 - t31 * mrSges(7,2) - t43 + t44; mrSges(7,1) * t50 - mrSges(7,2) * t85; 0; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t13(1) t13(2) t13(4) t13(7) t13(11) t13(16); t13(2) t13(3) t13(5) t13(8) t13(12) t13(17); t13(4) t13(5) t13(6) t13(9) t13(13) t13(18); t13(7) t13(8) t13(9) t13(10) t13(14) t13(19); t13(11) t13(12) t13(13) t13(14) t13(15) t13(20); t13(16) t13(17) t13(18) t13(19) t13(20) t13(21);];
Mq  = res;
