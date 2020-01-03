% Calculate joint inertia matrix for
% S5RRPPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR10_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR10_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR10_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR10_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR10_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR10_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:43:03
% EndTime: 2019-12-31 19:43:05
% DurationCPUTime: 0.72s
% Computational Cost: add. (563->199), mult. (1149->269), div. (0->0), fcn. (989->6), ass. (0->82)
t70 = sin(pkin(8));
t71 = cos(pkin(8));
t98 = t70 ^ 2 + t71 ^ 2;
t97 = 2 * pkin(6);
t96 = m(4) * pkin(2);
t95 = pkin(3) + pkin(4);
t82 = t70 * qJ(4) + pkin(2);
t42 = -t71 * pkin(3) - t82;
t94 = m(5) * t42;
t75 = cos(qJ(2));
t93 = pkin(6) * t75;
t92 = Ifges(4,4) * t70;
t91 = Ifges(4,4) * t71;
t90 = Ifges(5,5) * t70;
t89 = Ifges(5,5) * t71;
t88 = t70 * mrSges(5,3);
t73 = sin(qJ(2));
t87 = t70 * t73;
t86 = t71 * t73;
t85 = -pkin(7) + qJ(3);
t40 = t75 * mrSges(5,1) + mrSges(5,2) * t86;
t32 = mrSges(4,1) * t87 + mrSges(4,2) * t86;
t43 = -t75 * pkin(2) - t73 * qJ(3) - pkin(1);
t20 = t70 * t43 + t71 * t93;
t84 = t98 * qJ(3) ^ 2;
t72 = sin(qJ(5));
t74 = cos(qJ(5));
t37 = t74 * t70 - t72 * t71;
t28 = t37 * t73;
t78 = t72 * t70 + t74 * t71;
t29 = t78 * t73;
t83 = Ifges(6,5) * t29 + Ifges(6,6) * t28 + Ifges(6,3) * t75;
t57 = t70 * t93;
t19 = t71 * t43 - t57;
t15 = -t75 * qJ(4) + t20;
t5 = -t28 * mrSges(6,1) + t29 * mrSges(6,2);
t8 = mrSges(6,1) * t78 + t37 * mrSges(6,2);
t65 = t75 * pkin(3);
t16 = -t19 + t65;
t80 = t15 * t71 + t16 * t70;
t79 = -t19 * t70 + t20 * t71;
t77 = pkin(6) ^ 2;
t69 = t75 ^ 2;
t68 = t73 ^ 2;
t64 = t68 * t77;
t61 = t70 * mrSges(4,2);
t53 = mrSges(5,1) * t87;
t52 = qJ(4) * t86;
t51 = Ifges(4,1) * t70 + t91;
t50 = Ifges(5,1) * t70 - t89;
t49 = Ifges(4,2) * t71 + t92;
t48 = -Ifges(5,3) * t71 + t90;
t47 = t85 * t71;
t46 = -t71 * mrSges(4,1) + t61;
t45 = -t71 * mrSges(5,1) - t88;
t44 = t85 * t70;
t41 = -mrSges(5,2) * t87 - t75 * mrSges(5,3);
t39 = -t75 * mrSges(4,1) - mrSges(4,3) * t86;
t38 = t75 * mrSges(4,2) - mrSges(4,3) * t87;
t34 = Ifges(6,5) * t37;
t33 = Ifges(6,6) * t78;
t31 = -mrSges(5,3) * t86 + t53;
t30 = t95 * t71 + t82;
t27 = -t52 + (pkin(3) * t70 + pkin(6)) * t73;
t26 = -Ifges(4,5) * t75 + (Ifges(4,1) * t71 - t92) * t73;
t25 = -Ifges(5,4) * t75 + (Ifges(5,1) * t71 + t90) * t73;
t24 = -Ifges(4,6) * t75 + (-Ifges(4,2) * t70 + t91) * t73;
t23 = -Ifges(5,6) * t75 + (Ifges(5,3) * t70 + t89) * t73;
t18 = t75 * mrSges(6,1) - t29 * mrSges(6,3);
t17 = -t75 * mrSges(6,2) + t28 * mrSges(6,3);
t13 = -t52 - (-t95 * t70 - pkin(6)) * t73;
t12 = t72 * t44 + t74 * t47;
t11 = t74 * t44 - t72 * t47;
t10 = Ifges(6,1) * t37 - Ifges(6,4) * t78;
t9 = Ifges(6,4) * t37 - Ifges(6,2) * t78;
t7 = pkin(7) * t87 + t15;
t6 = t75 * pkin(4) + t57 + t65 + (-pkin(7) * t73 - t43) * t71;
t4 = Ifges(6,1) * t29 + Ifges(6,4) * t28 + Ifges(6,5) * t75;
t3 = Ifges(6,4) * t29 + Ifges(6,2) * t28 + Ifges(6,6) * t75;
t2 = t72 * t6 + t74 * t7;
t1 = t74 * t6 - t72 * t7;
t14 = [0.2e1 * t1 * t18 - 0.2e1 * t13 * t5 + 0.2e1 * t15 * t41 + 0.2e1 * t16 * t40 + 0.2e1 * t2 * t17 + 0.2e1 * t19 * t39 + 0.2e1 * t20 * t38 + 0.2e1 * t27 * t31 + t28 * t3 + t29 * t4 + Ifges(2,3) + (t68 + t69) * mrSges(3,3) * t97 + m(3) * (pkin(1) ^ 2 + t69 * t77 + t64) + m(4) * (t19 ^ 2 + t20 ^ 2 + t64) + m(5) * (t15 ^ 2 + t16 ^ 2 + t27 ^ 2) + m(6) * (t1 ^ 2 + t13 ^ 2 + t2 ^ 2) + (0.2e1 * pkin(1) * mrSges(3,1) + (Ifges(3,2) + Ifges(5,2) + Ifges(4,3)) * t75 + t83) * t75 + (-0.2e1 * pkin(1) * mrSges(3,2) + Ifges(3,1) * t73 + t32 * t97 + (t25 + t26) * t71 + (t23 - t24) * t70 + ((2 * Ifges(3,4)) + (-Ifges(5,4) - Ifges(4,5)) * t71 + (Ifges(4,6) - Ifges(5,6)) * t70) * t75) * t73; -t78 * t3 / 0.2e1 + t37 * t4 / 0.2e1 + t42 * t31 + t28 * t9 / 0.2e1 + t29 * t10 / 0.2e1 + t30 * t5 - pkin(2) * t32 - t13 * t8 + t12 * t17 + t11 * t18 + (-t23 / 0.2e1 + t24 / 0.2e1) * t71 + (t25 / 0.2e1 + t26 / 0.2e1) * t70 + (-pkin(6) * mrSges(3,2) + t34 / 0.2e1 - t33 / 0.2e1 + Ifges(3,6) + (Ifges(5,6) / 0.2e1 - Ifges(4,6) / 0.2e1) * t71 + (-Ifges(5,4) / 0.2e1 - Ifges(4,5) / 0.2e1) * t70) * t75 + (-t1 * t37 - t2 * t78) * mrSges(6,3) + t79 * mrSges(4,3) + t80 * mrSges(5,2) + m(6) * (t11 * t1 + t12 * t2 - t30 * t13) + ((t38 + t41) * t71 + (-t39 + t40) * t70 + m(4) * t79 + m(5) * t80) * qJ(3) + (Ifges(3,5) + (t50 / 0.2e1 + t51 / 0.2e1) * t71 + (t48 / 0.2e1 - t49 / 0.2e1) * t70 + (-mrSges(3,1) + t46 - t96) * pkin(6)) * t73 + (t45 + t94) * t27; -0.2e1 * pkin(2) * t46 + t37 * t10 + 0.2e1 * t30 * t8 - t78 * t9 + 0.2e1 * t42 * t45 + Ifges(3,3) + (-t48 + t49) * t71 + (t50 + t51) * t70 + 0.2e1 * (-t11 * t37 - t12 * t78) * mrSges(6,3) + m(6) * (t11 ^ 2 + t12 ^ 2 + t30 ^ 2) + m(4) * (pkin(2) ^ 2 + t84) + m(5) * (t42 ^ 2 + t84) + 0.2e1 * (mrSges(5,2) + mrSges(4,3)) * qJ(3) * t98; t53 + (m(4) * pkin(6) - t71 * mrSges(5,3)) * t73 + m(5) * t27 + m(6) * t13 - t5 + t32; -t96 - t88 + t61 + (-mrSges(5,1) - mrSges(4,1)) * t71 + t94 - m(6) * t30 - t8; m(4) + m(5) + m(6); t72 * t17 + t74 * t18 + m(6) * (t74 * t1 + t72 * t2) + m(5) * t16 + t40; m(6) * (t74 * t11 + t72 * t12) + (m(5) * qJ(3) + mrSges(5,2)) * t70 + (-t74 * t37 - t72 * t78) * mrSges(6,3); 0; m(5) + m(6) * (t72 ^ 2 + t74 ^ 2); t1 * mrSges(6,1) - t2 * mrSges(6,2) + t83; t11 * mrSges(6,1) - t12 * mrSges(6,2) - t33 + t34; 0; t74 * mrSges(6,1) - t72 * mrSges(6,2); Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t14(1), t14(2), t14(4), t14(7), t14(11); t14(2), t14(3), t14(5), t14(8), t14(12); t14(4), t14(5), t14(6), t14(9), t14(13); t14(7), t14(8), t14(9), t14(10), t14(14); t14(11), t14(12), t14(13), t14(14), t14(15);];
Mq = res;
