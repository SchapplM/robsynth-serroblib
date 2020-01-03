% Calculate joint inertia matrix for
% S5RRRPP3
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
% Datum: 2019-12-31 20:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPP3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP3_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP3_inertiaJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP3_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP3_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP3_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:53:08
% EndTime: 2019-12-31 20:53:10
% DurationCPUTime: 0.36s
% Computational Cost: add. (303->128), mult. (561->140), div. (0->0), fcn. (305->4), ass. (0->47)
t43 = cos(qJ(3));
t39 = t43 ^ 2;
t41 = sin(qJ(3));
t72 = t41 ^ 2 + t39;
t71 = -m(5) * pkin(3) + mrSges(5,2);
t70 = (mrSges(4,3) + mrSges(5,1)) * t72;
t69 = 2 * mrSges(6,1);
t12 = t43 * mrSges(5,2) - t41 * mrSges(5,3);
t68 = 0.2e1 * t12;
t14 = -t41 * mrSges(6,2) - t43 * mrSges(6,3);
t67 = 0.2e1 * t14;
t65 = t41 * pkin(7);
t44 = cos(qJ(2));
t64 = t44 * pkin(1);
t28 = t41 * mrSges(6,1);
t42 = sin(qJ(2));
t24 = t42 * pkin(1) + pkin(7);
t63 = t41 * t24;
t40 = -pkin(3) - qJ(5);
t62 = t72 * pkin(7) * t24;
t61 = (Ifges(6,2) + Ifges(5,3)) * t43 + (Ifges(5,6) - Ifges(6,6)) * t41;
t60 = t72 * t24 ^ 2;
t59 = t41 * mrSges(5,1) + t28;
t58 = t72 * pkin(7) ^ 2;
t56 = qJ(4) * t43;
t55 = Ifges(4,2) * t39 + Ifges(3,3) + ((Ifges(6,3) + Ifges(4,1)) * t41 + ((2 * Ifges(4,4)) - Ifges(6,6)) * t43) * t41;
t54 = -t41 * qJ(4) - pkin(2);
t9 = -t43 * pkin(3) + t54;
t52 = (t44 * mrSges(3,1) - t42 * mrSges(3,2)) * pkin(1);
t2 = t40 * t43 + t54;
t50 = 0.2e1 * t70;
t49 = t40 * t28 + (mrSges(6,1) + mrSges(5,1)) * t56 + (Ifges(4,6) - Ifges(6,4) - Ifges(5,5)) * t43 + (-pkin(3) * mrSges(5,1) - Ifges(5,4) + Ifges(4,5) + Ifges(6,5)) * t41;
t48 = m(5) * t56 + (mrSges(5,3) - mrSges(4,2)) * t43 + (-mrSges(4,1) + t71) * t41;
t45 = qJ(4) ^ 2;
t37 = t43 * pkin(4);
t35 = t41 * pkin(4);
t30 = t43 * mrSges(6,1);
t25 = -pkin(2) - t64;
t19 = t43 * pkin(7) + t37;
t18 = t35 + t65;
t17 = -Ifges(5,2) * t41 - Ifges(5,6) * t43;
t13 = -t43 * mrSges(4,1) + t41 * mrSges(4,2);
t5 = t43 * t24 + t37;
t4 = t35 + t63;
t3 = t9 - t64;
t1 = t2 - t64;
t6 = [t1 * t67 + t3 * t68 + 0.2e1 * t25 * t13 + Ifges(2,3) + (t4 * t69 - t17) * t41 + 0.2e1 * t52 + (t5 * t69 + t61) * t43 + t50 * t24 + m(6) * (t1 ^ 2 + t4 ^ 2 + t5 ^ 2) + m(5) * (t3 ^ 2 + t60) + m(4) * (t25 ^ 2 + t60) + m(3) * (t42 ^ 2 + t44 ^ 2) * pkin(1) ^ 2 + t55; -t41 * t17 + t61 * t43 + (t2 + t1) * t14 + (-pkin(2) + t25) * t13 + (t9 + t3) * t12 + t52 + m(5) * (t9 * t3 + t62) + m(6) * (t2 * t1 + t18 * t4 + t19 * t5) + m(4) * (-pkin(2) * t25 + t62) + ((t19 + t5) * t43 + (t18 + t4) * t41) * mrSges(6,1) + t55 + (pkin(7) + t24) * t70; -0.2e1 * pkin(2) * t13 + t9 * t68 + t2 * t67 + (t18 * t69 - t17) * t41 + (t19 * t69 + t61) * t43 + m(5) * (t9 ^ 2 + t58) + m(6) * (t18 ^ 2 + t19 ^ 2 + t2 ^ 2) + m(4) * (pkin(2) ^ 2 + t58) + t50 * pkin(7) + t55; m(6) * (qJ(4) * t5 + t40 * t4) - t4 * mrSges(6,3) + t5 * mrSges(6,2) + t48 * t24 + t49; m(6) * (qJ(4) * t19 + t40 * t18) + t19 * mrSges(6,2) - t18 * mrSges(6,3) + t48 * pkin(7) + t49; -0.2e1 * pkin(3) * mrSges(5,2) - 0.2e1 * t40 * mrSges(6,3) + Ifges(5,1) + Ifges(6,1) + Ifges(4,3) + 0.2e1 * (mrSges(5,3) + mrSges(6,2)) * qJ(4) + m(5) * (pkin(3) ^ 2 + t45) + m(6) * (t40 ^ 2 + t45); m(5) * t63 + m(6) * t4 + t59; m(5) * t65 + m(6) * t18 + t59; m(6) * t40 - mrSges(6,3) + t71; m(5) + m(6); m(6) * t5 + t30; m(6) * t19 + t30; m(6) * qJ(4) + mrSges(6,2); 0; m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t6(1), t6(2), t6(4), t6(7), t6(11); t6(2), t6(3), t6(5), t6(8), t6(12); t6(4), t6(5), t6(6), t6(9), t6(13); t6(7), t6(8), t6(9), t6(10), t6(14); t6(11), t6(12), t6(13), t6(14), t6(15);];
Mq = res;
