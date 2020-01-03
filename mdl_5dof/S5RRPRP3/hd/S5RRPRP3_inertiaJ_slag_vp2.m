% Calculate joint inertia matrix for
% S5RRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 19:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP3_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP3_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP3_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP3_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP3_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:50:56
% EndTime: 2019-12-31 19:50:57
% DurationCPUTime: 0.34s
% Computational Cost: add. (523->111), mult. (969->133), div. (0->0), fcn. (855->6), ass. (0->49)
t45 = cos(pkin(8));
t77 = t45 ^ 2;
t75 = (mrSges(5,3) + mrSges(6,2));
t76 = 2 * t75;
t74 = -m(6) * pkin(4) - mrSges(6,1);
t73 = -mrSges(5,1) + t74;
t72 = m(6) * qJ(5) - mrSges(5,2) + mrSges(6,3);
t44 = sin(pkin(8));
t46 = sin(qJ(4));
t66 = t46 * t44;
t67 = cos(qJ(4));
t30 = -t67 * t45 + t66;
t61 = t67 * t44;
t31 = t46 * t45 + t61;
t12 = t30 * mrSges(6,1) - t31 * mrSges(6,3);
t71 = 0.2e1 * t12;
t13 = t30 * mrSges(5,1) + t31 * mrSges(5,2);
t70 = 0.2e1 * t13;
t48 = cos(qJ(2));
t69 = t48 * pkin(1);
t47 = sin(qJ(2));
t37 = t47 * pkin(1) + qJ(3);
t68 = -pkin(7) - t37;
t65 = t44 ^ 2 + t77;
t41 = t45 * pkin(7);
t19 = t45 * t37 + t41;
t10 = t67 * t19 + t68 * t66;
t8 = t46 * t19 - t68 * t61;
t64 = t10 ^ 2 + t8 ^ 2;
t34 = t45 * qJ(3) + t41;
t60 = (-pkin(7) - qJ(3)) * t44;
t15 = t46 * t34 - t67 * t60;
t17 = t67 * t34 + t46 * t60;
t63 = t15 ^ 2 + t17 ^ 2;
t38 = -t45 * pkin(3) - pkin(2);
t62 = t17 * t10 + t15 * t8;
t33 = -t45 * mrSges(4,1) + t44 * mrSges(4,2);
t59 = t65 * qJ(3);
t58 = 0.2e1 * t65 * mrSges(4,3);
t57 = Ifges(4,2) * t77 + Ifges(3,3) + (Ifges(4,1) * t44 + 0.2e1 * Ifges(4,4) * t45) * t44 + (Ifges(6,1) + Ifges(5,1)) * t31 ^ 2 + ((Ifges(6,3) + Ifges(5,2)) * t30 + 0.2e1 * (-Ifges(5,4) + Ifges(6,5)) * t31) * t30;
t56 = (t48 * mrSges(3,1) - t47 * mrSges(3,2)) * pkin(1);
t11 = t30 * pkin(4) - t31 * qJ(5) + t38;
t52 = t12 + t13 + t33;
t51 = (-mrSges(6,2) * pkin(4) + Ifges(6,4) + Ifges(5,5)) * t31 + (-mrSges(6,2) * qJ(5) - Ifges(5,6) + Ifges(6,6)) * t30;
t39 = -pkin(2) - t69;
t32 = t38 - t69;
t21 = t31 * mrSges(6,2);
t6 = t11 - t69;
t1 = [t57 + m(3) * (t47 ^ 2 + t48 ^ 2) * pkin(1) ^ 2 + m(4) * (t65 * t37 ^ 2 + t39 ^ 2) + t37 * t58 + m(6) * (t6 ^ 2 + t64) + m(5) * (t32 ^ 2 + t64) + t32 * t70 + 0.2e1 * t39 * t33 + t6 * t71 + Ifges(2,3) + 0.2e1 * t56 + (-t10 * t30 + t8 * t31) * t76; t57 + t56 + m(4) * (-pkin(2) * t39 + t37 * t59) + (t65 * t37 + t59) * mrSges(4,3) + m(6) * (t11 * t6 + t62) + m(5) * (t38 * t32 + t62) + (-pkin(2) + t39) * t33 + (t32 + t38) * t13 + (t11 + t6) * t12 + t75 * ((t15 + t8) * t31 + (-t10 - t17) * t30); -0.2e1 * pkin(2) * t33 + t11 * t71 + t38 * t70 + qJ(3) * t58 + m(6) * (t11 ^ 2 + t63) + m(5) * (t38 ^ 2 + t63) + m(4) * (t65 * qJ(3) ^ 2 + pkin(2) ^ 2) + t57 + (t15 * t31 - t17 * t30) * t76; m(4) * t39 + m(5) * t32 + m(6) * t6 + t52; -m(4) * pkin(2) + m(5) * t38 + m(6) * t11 + t52; m(4) + m(5) + m(6); t10 * t72 + t73 * t8 + t51; t15 * t73 + t17 * t72 + t51; 0; Ifges(6,2) + Ifges(5,3) + 0.2e1 * pkin(4) * mrSges(6,1) + 0.2e1 * qJ(5) * mrSges(6,3) + m(6) * (pkin(4) ^ 2 + qJ(5) ^ 2); m(6) * t8 + t21; m(6) * t15 + t21; 0; t74; m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
