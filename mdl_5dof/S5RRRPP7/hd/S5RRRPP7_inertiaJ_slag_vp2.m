% Calculate joint inertia matrix for
% S5RRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 21:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPP7_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP7_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP7_inertiaJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP7_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP7_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP7_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:03:58
% EndTime: 2019-12-31 21:04:00
% DurationCPUTime: 0.78s
% Computational Cost: add. (467->213), mult. (932->268), div. (0->0), fcn. (664->4), ass. (0->78)
t61 = sin(qJ(3));
t63 = cos(qJ(3));
t95 = t61 ^ 2 + t63 ^ 2;
t94 = 2 * pkin(6);
t93 = -2 * mrSges(6,3);
t92 = pkin(3) + pkin(4);
t64 = cos(qJ(2));
t91 = pkin(6) * t64;
t90 = Ifges(4,4) * t61;
t89 = Ifges(4,4) * t63;
t88 = Ifges(6,4) * t61;
t87 = Ifges(6,4) * t63;
t86 = Ifges(5,5) * t61;
t85 = Ifges(5,5) * t63;
t84 = Ifges(6,5) * t64;
t62 = sin(qJ(2));
t83 = t61 * t62;
t82 = t62 * t63;
t81 = mrSges(5,2) - mrSges(6,3);
t80 = Ifges(4,6) + Ifges(6,6);
t79 = pkin(7) - qJ(5);
t24 = t64 * mrSges(5,1) + mrSges(5,2) * t82;
t27 = -t64 * pkin(2) - t62 * pkin(7) - pkin(1);
t7 = t61 * t27 + t63 * t91;
t30 = t63 * mrSges(6,1) + t61 * mrSges(6,2);
t78 = t95 * pkin(7) ^ 2;
t77 = qJ(4) * t63;
t76 = qJ(5) * t62;
t75 = Ifges(4,3) + Ifges(6,3) + Ifges(5,2);
t74 = -Ifges(5,6) * t83 + (-Ifges(5,4) - Ifges(4,5)) * t82;
t73 = t61 * qJ(4) + pkin(2);
t44 = t61 * t91;
t6 = t63 * t27 - t44;
t16 = -mrSges(6,1) * t83 + mrSges(6,2) * t82;
t22 = t64 * mrSges(6,1) - mrSges(6,3) * t82;
t4 = -t64 * qJ(4) + t7;
t71 = t61 * mrSges(4,1) + t63 * mrSges(4,2);
t70 = t61 * mrSges(5,1) - t63 * mrSges(5,3);
t69 = -pkin(3) * t61 + t77;
t68 = pkin(6) ^ 2;
t66 = qJ(4) ^ 2;
t60 = t64 ^ 2;
t58 = t62 ^ 2;
t56 = t64 * pkin(3);
t54 = t58 * t68;
t52 = Ifges(5,4) * t61;
t51 = Ifges(4,5) * t61;
t50 = Ifges(4,6) * t63;
t38 = Ifges(4,1) * t61 + t89;
t37 = Ifges(5,1) * t61 - t85;
t36 = Ifges(6,1) * t61 - t87;
t35 = Ifges(4,2) * t63 + t90;
t34 = -Ifges(6,2) * t63 + t88;
t33 = -Ifges(5,3) * t63 + t86;
t32 = t79 * t63;
t31 = -t63 * mrSges(4,1) + t61 * mrSges(4,2);
t29 = -t63 * mrSges(5,1) - t61 * mrSges(5,3);
t28 = t79 * t61;
t26 = -t63 * pkin(3) - t73;
t25 = -mrSges(5,2) * t83 - t64 * mrSges(5,3);
t23 = -t64 * mrSges(4,1) - mrSges(4,3) * t82;
t21 = t64 * mrSges(4,2) - mrSges(4,3) * t83;
t20 = -t64 * mrSges(6,2) + mrSges(6,3) * t83;
t18 = t92 * t63 + t73;
t17 = t71 * t62;
t15 = t70 * t62;
t14 = (pkin(6) - t69) * t62;
t13 = -Ifges(4,5) * t64 + (Ifges(4,1) * t63 - t90) * t62;
t12 = -Ifges(5,4) * t64 + (Ifges(5,1) * t63 + t86) * t62;
t11 = t84 + (Ifges(6,1) * t63 + t88) * t62;
t10 = -Ifges(4,6) * t64 + (-Ifges(4,2) * t61 + t89) * t62;
t9 = Ifges(6,6) * t64 + (Ifges(6,2) * t61 + t87) * t62;
t8 = -Ifges(5,6) * t64 + (Ifges(5,3) * t61 + t85) * t62;
t5 = t56 - t6;
t3 = (-t92 * t61 - pkin(6) + t77) * t62;
t2 = t61 * t76 + t4;
t1 = t64 * pkin(4) + t44 + t56 + (-t27 - t76) * t63;
t19 = [0.2e1 * t1 * t22 + 0.2e1 * t14 * t15 + 0.2e1 * t3 * t16 + 0.2e1 * t2 * t20 + 0.2e1 * t7 * t21 + 0.2e1 * t6 * t23 + 0.2e1 * t5 * t24 + 0.2e1 * t4 * t25 + Ifges(2,3) + (t58 + t60) * mrSges(3,3) * t94 + m(6) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(5) * (t14 ^ 2 + t4 ^ 2 + t5 ^ 2) + m(4) * (t6 ^ 2 + t7 ^ 2 + t54) + m(3) * (pkin(1) ^ 2 + t60 * t68 + t54) + (0.2e1 * pkin(1) * mrSges(3,1) + (Ifges(3,2) + t75) * t64 + t74) * t64 + (-0.2e1 * pkin(1) * mrSges(3,2) + Ifges(3,1) * t62 + 0.2e1 * Ifges(3,4) * t64 + t17 * t94 + (t11 + t12 + t13 + t84) * t63 + (t80 * t64 - t10 + t8 + t9) * t61) * t62; -pkin(2) * t17 + t26 * t15 + t18 * t16 + t32 * t20 + t28 * t22 + t3 * t30 + m(6) * (t28 * t1 + t18 * t3 + t32 * t2) + (-t52 / 0.2e1 - t51 / 0.2e1 - t50 / 0.2e1 + Ifges(3,6) - pkin(6) * mrSges(3,2)) * t64 + (-t8 / 0.2e1 - t9 / 0.2e1 + t10 / 0.2e1 - t2 * mrSges(6,3) + t4 * mrSges(5,2) + t7 * mrSges(4,3) + (Ifges(5,6) / 0.2e1 - Ifges(6,6) / 0.2e1) * t64) * t63 + (t84 / 0.2e1 + t11 / 0.2e1 + t12 / 0.2e1 + t13 / 0.2e1 - t1 * mrSges(6,3) + t5 * mrSges(5,2) - t6 * mrSges(4,3)) * t61 + ((t21 + t25) * t63 + (-t23 + t24) * t61 + m(5) * (t4 * t63 + t5 * t61) + m(4) * (-t6 * t61 + t7 * t63)) * pkin(7) + (Ifges(3,5) + (t36 / 0.2e1 + t37 / 0.2e1 + t38 / 0.2e1) * t63 + (t33 / 0.2e1 + t34 / 0.2e1 - t35 / 0.2e1) * t61 + (-m(4) * pkin(2) - mrSges(3,1) + t31) * pkin(6)) * t62 + (m(5) * t26 + t29) * t14; -0.2e1 * pkin(2) * t31 + 0.2e1 * t18 * t30 + 0.2e1 * t26 * t29 + Ifges(3,3) + m(6) * (t18 ^ 2 + t28 ^ 2 + t32 ^ 2) + m(5) * (t26 ^ 2 + t78) + m(4) * (pkin(2) ^ 2 + t78) + (t32 * t93 - t33 - t34 + t35) * t63 + (t28 * t93 + t36 + t37 + t38) * t61 + 0.2e1 * (mrSges(5,2) + mrSges(4,3)) * pkin(7) * t95; t6 * mrSges(4,1) - t5 * mrSges(5,1) - t1 * mrSges(6,1) - t7 * mrSges(4,2) + t2 * mrSges(6,2) + t4 * mrSges(5,3) - pkin(3) * t24 - t92 * t22 + (t20 + t25) * qJ(4) + m(5) * (-pkin(3) * t5 + qJ(4) * t4) + m(6) * (qJ(4) * t2 - t1 * t92) - t75 * t64 + (-Ifges(6,5) * t63 - t80 * t61) * t62 - t74; m(6) * (qJ(4) * t32 - t28 * t92) - t28 * mrSges(6,1) + t32 * mrSges(6,2) + t52 + t51 + t50 + (-pkin(3) * mrSges(5,2) + mrSges(6,3) * t92 - Ifges(6,5)) * t61 + (t81 * qJ(4) - Ifges(5,6) + Ifges(6,6)) * t63 + (m(5) * t69 - t70 - t71) * pkin(7); 0.2e1 * pkin(3) * mrSges(5,1) + 0.2e1 * t92 * mrSges(6,1) + 0.2e1 * (mrSges(5,3) + mrSges(6,2)) * qJ(4) + m(5) * (pkin(3) ^ 2 + t66) + m(6) * (t92 ^ 2 + t66) + t75; m(5) * t5 + m(6) * t1 + t22 + t24; m(6) * t28 + (m(5) * pkin(7) + t81) * t61; -m(5) * pkin(3) - m(6) * t92 - mrSges(5,1) - mrSges(6,1); m(5) + m(6); m(6) * t3 + t16; m(6) * t18 + t30; 0; 0; m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t19(1), t19(2), t19(4), t19(7), t19(11); t19(2), t19(3), t19(5), t19(8), t19(12); t19(4), t19(5), t19(6), t19(9), t19(13); t19(7), t19(8), t19(9), t19(10), t19(14); t19(11), t19(12), t19(13), t19(14), t19(15);];
Mq = res;
