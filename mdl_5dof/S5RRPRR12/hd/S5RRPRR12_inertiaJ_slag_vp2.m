% Calculate joint inertia matrix for
% S5RRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR12_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR12_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR12_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR12_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR12_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR12_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:29:04
% EndTime: 2019-12-31 20:29:06
% DurationCPUTime: 0.64s
% Computational Cost: add. (637->168), mult. (1166->225), div. (0->0), fcn. (1015->6), ass. (0->65)
t45 = sin(qJ(2));
t48 = cos(qJ(2));
t86 = t45 ^ 2 + t48 ^ 2;
t43 = sin(qJ(5));
t46 = cos(qJ(5));
t62 = t43 ^ 2 + t46 ^ 2;
t59 = t62 * mrSges(6,3);
t24 = -mrSges(6,1) * t46 + mrSges(6,2) * t43;
t64 = -t24 + mrSges(5,1);
t44 = sin(qJ(4));
t47 = cos(qJ(4));
t15 = t44 * t45 + t47 * t48;
t16 = -t44 * t48 + t45 * t47;
t69 = t16 * t43;
t7 = -mrSges(6,2) * t15 - mrSges(6,3) * t69;
t68 = t16 * t46;
t8 = mrSges(6,1) * t15 - mrSges(6,3) * t68;
t85 = -t43 * t8 + t46 * t7;
t84 = -m(4) * pkin(2) - mrSges(4,1);
t83 = m(6) * pkin(4) + t64;
t76 = pkin(6) - pkin(7);
t28 = t76 * t48;
t60 = t76 * t45;
t9 = t44 * t28 - t47 * t60;
t82 = t9 ^ 2;
t81 = 0.2e1 * t9;
t23 = -t48 * pkin(2) - t45 * qJ(3) - pkin(1);
t13 = t48 * pkin(3) - t23;
t80 = 0.2e1 * t13;
t79 = -0.2e1 * t23;
t78 = -0.2e1 * t24;
t77 = t46 / 0.2e1;
t73 = t47 * t9;
t72 = Ifges(6,4) * t43;
t71 = Ifges(6,4) * t46;
t11 = t47 * t28 + t44 * t60;
t70 = t11 * mrSges(5,2);
t49 = -pkin(2) - pkin(3);
t21 = -qJ(3) * t44 + t47 * t49;
t67 = t21 * mrSges(5,1);
t22 = t47 * qJ(3) + t44 * t49;
t66 = t22 * mrSges(5,2);
t65 = Ifges(6,5) * t68 + Ifges(6,3) * t15;
t63 = t86 * pkin(6) ^ 2;
t25 = Ifges(6,5) * t43 + Ifges(6,6) * t46;
t61 = t25 / 0.2e1 - Ifges(5,6);
t20 = -pkin(8) + t22;
t58 = t62 * t20;
t57 = t62 * t44;
t5 = pkin(4) * t15 - pkin(8) * t16 + t13;
t1 = -t11 * t43 + t46 * t5;
t2 = t11 * t46 + t43 * t5;
t55 = -t1 * t43 + t2 * t46;
t54 = mrSges(6,1) * t43 + mrSges(6,2) * t46;
t26 = Ifges(6,2) * t46 + t72;
t27 = Ifges(6,1) * t43 + t71;
t53 = t46 * t26 + t43 * t27 + Ifges(5,3);
t52 = t27 * t77 - t43 * t26 / 0.2e1 + Ifges(5,5);
t41 = t47 ^ 2;
t38 = t44 ^ 2;
t19 = pkin(4) - t21;
t6 = t54 * t16;
t4 = Ifges(6,5) * t15 + (Ifges(6,1) * t46 - t72) * t16;
t3 = Ifges(6,6) * t15 + (-Ifges(6,2) * t43 + t71) * t16;
t10 = [0.2e1 * t1 * t8 + 0.2e1 * t2 * t7 + t6 * t81 + Ifges(2,3) + (0.2e1 * pkin(1) * mrSges(3,1) + mrSges(4,1) * t79 + (Ifges(3,2) + Ifges(4,3)) * t48) * t48 + (-0.2e1 * pkin(1) * mrSges(3,2) + mrSges(4,3) * t79 + (Ifges(4,1) + Ifges(3,1)) * t45 + 0.2e1 * (Ifges(3,4) - Ifges(4,5)) * t48) * t45 + (mrSges(5,1) * t80 - 0.2e1 * t11 * mrSges(5,3) + Ifges(5,2) * t15 + t65) * t15 + (mrSges(5,2) * t80 + mrSges(5,3) * t81 + Ifges(5,1) * t16 - t43 * t3 + t46 * t4 + (-Ifges(6,6) * t43 - (2 * Ifges(5,4))) * t15) * t16 + m(6) * (t1 ^ 2 + t2 ^ 2 + t82) + m(5) * (t11 ^ 2 + t13 ^ 2 + t82) + m(3) * (pkin(1) ^ 2 + t63) + m(4) * (t23 ^ 2 + t63) + 0.2e1 * (mrSges(4,2) + mrSges(3,3)) * pkin(6) * t86; t70 + t19 * t6 + t64 * t9 + (-t3 / 0.2e1 - t2 * mrSges(6,3) + t20 * t7) * t46 + (-t4 / 0.2e1 + t1 * mrSges(6,3) - t20 * t8) * t43 + m(5) * (t11 * t22 - t21 * t9) + m(6) * (t19 * t9 + t20 * t55) + (qJ(3) * mrSges(4,2) + Ifges(3,6) - Ifges(4,6)) * t48 + (-pkin(2) * mrSges(4,2) + Ifges(4,4) + Ifges(3,5)) * t45 + (-t22 * mrSges(5,3) - t61) * t15 + (-t21 * mrSges(5,3) - t52) * t16 + ((m(4) * qJ(3) - mrSges(3,2) + mrSges(4,3)) * t48 + (-mrSges(3,1) + t84) * t45) * pkin(6); 0.2e1 * pkin(2) * mrSges(4,1) - 0.2e1 * t67 + 0.2e1 * t66 + 0.2e1 * qJ(3) * mrSges(4,3) + t19 * t78 + Ifges(4,2) + Ifges(3,3) - 0.2e1 * t20 * t59 + m(6) * (t20 ^ 2 * t62 + t19 ^ 2) + m(5) * (t21 ^ 2 + t22 ^ 2) + m(4) * (pkin(2) ^ 2 + qJ(3) ^ 2) + t53; (-t16 * mrSges(5,3) - t6) * t47 + (m(4) * pkin(6) + mrSges(4,2)) * t45 + (-t15 * mrSges(5,3) + t85) * t44 + m(6) * (t44 * t55 - t73) + m(5) * (t11 * t44 - t73); -t64 * t47 + (mrSges(5,2) - t59) * t44 + m(6) * (-t47 * t19 + t20 * t57) + m(5) * (t21 * t47 + t22 * t44) + t84; m(4) + m(5) * (t38 + t41) + m(6) * (t38 * t62 + t41); t43 * t4 / 0.2e1 + t3 * t77 - pkin(4) * t6 - t70 + t61 * t15 + t55 * mrSges(6,3) + t52 * t16 - t83 * t9 + (m(6) * t55 + t85) * pkin(8); m(6) * (-pkin(4) * t19 + pkin(8) * t58) - t66 + t67 + (t19 + pkin(4)) * t24 + (-pkin(8) * t62 + t58) * mrSges(6,3) - t53; -t44 * mrSges(5,2) + (m(6) * pkin(8) + mrSges(6,3)) * t57 + t83 * t47; m(6) * (pkin(8) ^ 2 * t62 + pkin(4) ^ 2) + pkin(4) * t78 + 0.2e1 * pkin(8) * t59 + t53; mrSges(6,1) * t1 - mrSges(6,2) * t2 - Ifges(6,6) * t69 + t65; -t20 * t54 - t25; -t54 * t44; -pkin(8) * t54 + t25; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t10(1), t10(2), t10(4), t10(7), t10(11); t10(2), t10(3), t10(5), t10(8), t10(12); t10(4), t10(5), t10(6), t10(9), t10(13); t10(7), t10(8), t10(9), t10(10), t10(14); t10(11), t10(12), t10(13), t10(14), t10(15);];
Mq = res;
