% Calculate joint inertia matrix for
% S5RRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 21:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRP3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP3_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP3_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP3_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP3_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP3_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:49:08
% EndTime: 2019-12-31 21:49:08
% DurationCPUTime: 0.32s
% Computational Cost: add. (454->113), mult. (845->147), div. (0->0), fcn. (517->6), ass. (0->72)
t63 = cos(qJ(4));
t59 = t63 ^ 2;
t60 = sin(qJ(4));
t58 = t60 ^ 2;
t99 = t58 + t59;
t61 = sin(qJ(3));
t88 = t61 * pkin(2);
t43 = pkin(8) + t88;
t82 = t43 * t59;
t83 = t43 * t58;
t98 = t82 + t83;
t97 = -m(6) * pkin(4) - mrSges(6,1);
t34 = -t63 * mrSges(6,1) - t60 * mrSges(6,3);
t65 = cos(qJ(2));
t45 = t65 * pkin(1) + pkin(2);
t64 = cos(qJ(3));
t62 = sin(qJ(2));
t92 = pkin(1) * t62;
t20 = t64 * t45 - t61 * t92;
t27 = -t63 * pkin(4) - t60 * qJ(5) - pkin(3);
t5 = -t20 + t27;
t1 = t5 * t34;
t21 = t61 * t45 + t64 * t92;
t19 = pkin(8) + t21;
t85 = t19 * t59;
t10 = mrSges(6,2) * t85;
t16 = t20 * mrSges(4,1);
t18 = -pkin(3) - t20;
t35 = -t63 * mrSges(5,1) + t60 * mrSges(5,2);
t2 = t18 * t35;
t86 = t19 * t58;
t7 = mrSges(5,3) * t86;
t8 = mrSges(6,2) * t86;
t9 = mrSges(5,3) * t85;
t96 = t1 + t10 + t16 + t2 + t7 + t8 + t9;
t87 = t64 * pkin(2);
t44 = -pkin(3) - t87;
t22 = t44 * t35;
t28 = mrSges(5,3) * t83;
t29 = mrSges(6,2) * t83;
t30 = mrSges(5,3) * t82;
t31 = mrSges(6,2) * t82;
t51 = mrSges(4,1) * t87;
t23 = t27 - t87;
t6 = t23 * t34;
t95 = t22 + t28 + t29 + t30 + t31 + t51 + t6;
t94 = t98 * t19;
t93 = m(6) * t60;
t91 = pkin(3) * t35;
t90 = pkin(8) * t58;
t89 = pkin(8) * t59;
t84 = t21 * mrSges(4,2);
t52 = t60 * mrSges(6,2);
t81 = (t85 + t86) * pkin(8);
t80 = t99 * t19 ^ 2;
t79 = t98 * pkin(8);
t78 = t99 * t43 ^ 2;
t77 = t99 * pkin(8) ^ 2;
t76 = qJ(5) * t63;
t75 = mrSges(4,2) * t88;
t74 = Ifges(4,3) + (Ifges(6,3) + Ifges(5,2)) * t59 + ((Ifges(6,1) + Ifges(5,1)) * t60 + 0.2e1 * (Ifges(5,4) - Ifges(6,5)) * t63) * t60;
t73 = (t65 * mrSges(3,1) - t62 * mrSges(3,2)) * pkin(1);
t72 = Ifges(3,3) + t74;
t71 = mrSges(6,2) * t76 - pkin(4) * t52 + (Ifges(5,6) - Ifges(6,6)) * t63 + (Ifges(6,4) + Ifges(5,5)) * t60;
t15 = t27 * t34;
t47 = mrSges(5,3) * t90;
t48 = mrSges(6,2) * t90;
t49 = mrSges(5,3) * t89;
t50 = mrSges(6,2) * t89;
t70 = t15 + t47 + t48 + t49 + t50 + t74 - t91;
t69 = m(6) * t76 + (mrSges(6,3) - mrSges(5,2)) * t63 + (-mrSges(5,1) + t97) * t60;
t3 = [t72 + 0.2e1 * t73 + m(3) * (t62 ^ 2 + t65 ^ 2) * pkin(1) ^ 2 + m(5) * (t18 ^ 2 + t80) + m(6) * (t5 ^ 2 + t80) + m(4) * (t20 ^ 2 + t21 ^ 2) + 0.2e1 * t16 + 0.2e1 * t9 + 0.2e1 * t10 + 0.2e1 * t7 + 0.2e1 * t8 + 0.2e1 * t2 + 0.2e1 * t1 - 0.2e1 * t84 + Ifges(2,3); t95 + t72 + t73 + (-t21 - t88) * mrSges(4,2) + m(5) * (t44 * t18 + t94) + m(6) * (t23 * t5 + t94) + m(4) * (t20 * t64 + t21 * t61) * pkin(2) + t96; -0.2e1 * t75 + 0.2e1 * t22 + 0.2e1 * t28 + 0.2e1 * t29 + 0.2e1 * t30 + 0.2e1 * t31 + 0.2e1 * t51 + 0.2e1 * t6 + m(6) * (t23 ^ 2 + t78) + m(5) * (t44 ^ 2 + t78) + m(4) * (t61 ^ 2 + t64 ^ 2) * pkin(2) ^ 2 + t72; t70 + m(5) * (-pkin(3) * t18 + t81) + m(6) * (t27 * t5 + t81) - t84 + t96; t70 + m(6) * (t27 * t23 + t79) + m(5) * (-pkin(3) * t44 + t79) - t75 + t95; -0.2e1 * t91 + 0.2e1 * t15 + 0.2e1 * t47 + 0.2e1 * t48 + 0.2e1 * t49 + 0.2e1 * t50 + m(6) * (t27 ^ 2 + t77) + m(5) * (pkin(3) ^ 2 + t77) + t74; t69 * t19 + t71; t69 * t43 + t71; t69 * pkin(8) + t71; Ifges(6,2) + Ifges(5,3) + 0.2e1 * pkin(4) * mrSges(6,1) + 0.2e1 * qJ(5) * mrSges(6,3) + m(6) * (pkin(4) ^ 2 + qJ(5) ^ 2); t19 * t93 + t52; t43 * t93 + t52; pkin(8) * t93 + t52; t97; m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t3(1), t3(2), t3(4), t3(7), t3(11); t3(2), t3(3), t3(5), t3(8), t3(12); t3(4), t3(5), t3(6), t3(9), t3(13); t3(7), t3(8), t3(9), t3(10), t3(14); t3(11), t3(12), t3(13), t3(14), t3(15);];
Mq = res;
