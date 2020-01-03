% Calculate joint inertia matrix for
% S5RRRRP4
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
% Datum: 2019-12-31 21:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRP4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP4_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP4_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP4_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP4_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP4_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:50:41
% EndTime: 2019-12-31 21:50:42
% DurationCPUTime: 0.48s
% Computational Cost: add. (639->144), mult. (1178->182), div. (0->0), fcn. (1013->6), ass. (0->57)
t95 = -mrSges(5,1) - mrSges(6,1);
t94 = -mrSges(5,2) + mrSges(6,3);
t55 = cos(qJ(3));
t93 = t55 ^ 2;
t91 = (mrSges(5,3) + mrSges(6,2));
t92 = 2 * t91;
t53 = sin(qJ(2));
t42 = t53 * pkin(1) + pkin(7);
t48 = t55 * pkin(8);
t30 = t55 * t42 + t48;
t51 = sin(qJ(4));
t54 = cos(qJ(4));
t52 = sin(qJ(3));
t73 = (-pkin(8) - t42) * t52;
t10 = t51 * t30 - t54 * t73;
t12 = t54 * t30 + t51 * t73;
t90 = t95 * t10 + t94 * t12;
t38 = t55 * pkin(7) + t48;
t74 = (-pkin(8) - pkin(7)) * t52;
t19 = t51 * t38 - t54 * t74;
t21 = t54 * t38 + t51 * t74;
t89 = t95 * t19 + t94 * t21;
t88 = 2 * mrSges(6,3);
t33 = t51 * t52 - t54 * t55;
t34 = t51 * t55 + t54 * t52;
t14 = t33 * mrSges(6,1) - t34 * mrSges(6,3);
t87 = 0.2e1 * t14;
t15 = t33 * mrSges(5,1) + t34 * mrSges(5,2);
t86 = 0.2e1 * t15;
t85 = t51 * pkin(3);
t56 = cos(qJ(2));
t84 = t56 * pkin(1);
t25 = t34 * mrSges(6,2);
t79 = Ifges(6,2) + Ifges(5,3);
t78 = t52 ^ 2 + t93;
t77 = t10 ^ 2 + t12 ^ 2;
t76 = t19 ^ 2 + t21 ^ 2;
t75 = t54 * t34 * mrSges(5,3);
t45 = -t55 * pkin(3) - pkin(2);
t72 = t19 * t10 + t21 * t12;
t71 = t78 * t42;
t70 = (Ifges(6,4) + Ifges(5,5)) * t34 + (-Ifges(5,6) + Ifges(6,6)) * t33;
t69 = -t52 * mrSges(4,1) - t55 * mrSges(4,2);
t68 = 0.2e1 * mrSges(4,3) * t78;
t67 = Ifges(4,2) * t93 + Ifges(3,3) + (Ifges(4,1) * t52 + 0.2e1 * Ifges(4,4) * t55) * t52 + (Ifges(6,1) + Ifges(5,1)) * t34 ^ 2 + ((Ifges(6,3) + Ifges(5,2)) * t33 + 0.2e1 * (-Ifges(5,4) + Ifges(6,5)) * t34) * t33;
t66 = (t56 * mrSges(3,1) - t53 * mrSges(3,2)) * pkin(1);
t65 = (t54 * mrSges(5,1) - t51 * mrSges(5,2)) * pkin(3);
t13 = t33 * pkin(4) - t34 * qJ(5) + t45;
t61 = (-pkin(4) * t34 - qJ(5) * t33) * mrSges(6,2) + t70;
t40 = qJ(5) + t85;
t43 = -t54 * pkin(3) - pkin(4);
t60 = Ifges(4,5) * t52 + Ifges(4,6) * t55 + t43 * t25 + t70 + (-mrSges(6,2) * t40 - mrSges(5,3) * t85) * t33;
t44 = -pkin(2) - t84;
t37 = -t55 * mrSges(4,1) + t52 * mrSges(4,2);
t36 = t45 - t84;
t6 = t13 - t84;
t1 = [t67 + m(4) * (t42 ^ 2 * t78 + t44 ^ 2) + m(3) * (t53 ^ 2 + t56 ^ 2) * pkin(1) ^ 2 + t42 * t68 + m(6) * (t6 ^ 2 + t77) + m(5) * (t36 ^ 2 + t77) + 0.2e1 * t44 * t37 + t36 * t86 + t6 * t87 + Ifges(2,3) + 0.2e1 * t66 + (t10 * t34 - t12 * t33) * t92; t67 + t66 + m(4) * (-pkin(2) * t44 + pkin(7) * t71) + (pkin(7) * t78 + t71) * mrSges(4,3) + m(6) * (t13 * t6 + t72) + m(5) * (t45 * t36 + t72) + (t44 - pkin(2)) * t37 + (t45 + t36) * t15 + (t13 + t6) * t14 + t91 * ((t10 + t19) * t34 + (-t12 - t21) * t33); -0.2e1 * pkin(2) * t37 + t13 * t87 + t45 * t86 + pkin(7) * t68 + m(6) * (t13 ^ 2 + t76) + m(5) * (t45 ^ 2 + t76) + m(4) * (pkin(7) ^ 2 * t78 + pkin(2) ^ 2) + t67 + (t19 * t34 - t21 * t33) * t92; t60 + (-t75 + m(5) * (-t10 * t54 + t12 * t51)) * pkin(3) + m(6) * (t43 * t10 + t40 * t12) + t69 * t42 + t90; t60 + (-t75 + m(5) * (-t19 * t54 + t21 * t51)) * pkin(3) + m(6) * (t43 * t19 + t40 * t21) + t69 * pkin(7) + t89; -0.2e1 * t43 * mrSges(6,1) + t40 * t88 + Ifges(4,3) + 0.2e1 * t65 + m(6) * (t40 ^ 2 + t43 ^ 2) + m(5) * (t51 ^ 2 + t54 ^ 2) * pkin(3) ^ 2 + t79; m(6) * (-pkin(4) * t10 + qJ(5) * t12) + t61 + t90; m(6) * (-pkin(4) * t19 + qJ(5) * t21) + t61 + t89; m(6) * (-pkin(4) * t43 + qJ(5) * t40) + t65 + (t40 + qJ(5)) * mrSges(6,3) + (-t43 + pkin(4)) * mrSges(6,1) + t79; 0.2e1 * pkin(4) * mrSges(6,1) + qJ(5) * t88 + m(6) * (pkin(4) ^ 2 + qJ(5) ^ 2) + t79; m(6) * t10 + t25; m(6) * t19 + t25; m(6) * t43 - mrSges(6,1); -m(6) * pkin(4) - mrSges(6,1); m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
