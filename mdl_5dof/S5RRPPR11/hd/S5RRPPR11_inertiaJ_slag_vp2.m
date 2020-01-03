% Calculate joint inertia matrix for
% S5RRPPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
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
% Datum: 2019-12-31 19:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR11_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR11_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR11_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR11_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR11_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR11_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:46:26
% EndTime: 2019-12-31 19:46:28
% DurationCPUTime: 0.56s
% Computational Cost: add. (570->167), mult. (1052->231), div. (0->0), fcn. (914->6), ass. (0->67)
t92 = pkin(3) + pkin(6);
t67 = sin(qJ(2));
t69 = cos(qJ(2));
t91 = t67 ^ 2 + t69 ^ 2;
t90 = -m(4) * pkin(2) + mrSges(4,2);
t63 = sin(pkin(8));
t89 = -t63 / 0.2e1;
t65 = -pkin(2) - qJ(4);
t88 = -pkin(7) + t65;
t87 = Ifges(5,4) * t63;
t64 = cos(pkin(8));
t86 = Ifges(5,4) * t64;
t85 = t63 * t69;
t84 = t64 * t69;
t77 = -t67 * qJ(3) - pkin(1);
t30 = t65 * t69 + t77;
t45 = t92 * t67;
t9 = t64 * t30 + t63 * t45;
t42 = t63 * mrSges(5,1) + t64 * mrSges(5,2);
t83 = t91 * pkin(6) ^ 2;
t46 = t92 * t69;
t82 = t63 ^ 2 + t64 ^ 2;
t66 = sin(qJ(5));
t68 = cos(qJ(5));
t32 = -t63 * t68 - t64 * t66;
t72 = t63 * t66 - t64 * t68;
t81 = t32 ^ 2 + t72 ^ 2;
t22 = t72 * t69;
t23 = t32 * t69;
t80 = Ifges(6,5) * t23 + Ifges(6,6) * t22 + Ifges(6,3) * t67;
t79 = m(5) * t82;
t78 = t82 * mrSges(5,3);
t6 = -t22 * mrSges(6,1) + mrSges(6,2) * t23;
t12 = -mrSges(6,1) * t32 - mrSges(6,2) * t72;
t26 = mrSges(5,1) * t84 - mrSges(5,2) * t85;
t36 = t64 * t45;
t5 = t67 * pkin(4) + t36 + (pkin(7) * t69 - t30) * t63;
t7 = -pkin(7) * t84 + t9;
t1 = t5 * t68 - t66 * t7;
t2 = t5 * t66 + t68 * t7;
t75 = -t1 * t72 - t2 * t32;
t8 = -t30 * t63 + t36;
t74 = t63 * t9 + t64 * t8;
t39 = t88 * t63;
t40 = t88 * t64;
t10 = -t39 * t66 + t40 * t68;
t11 = t39 * t68 + t40 * t66;
t73 = -t10 * t72 - t11 * t32;
t70 = qJ(3) ^ 2;
t48 = pkin(4) * t63 + qJ(3);
t44 = Ifges(5,1) * t64 - t87;
t43 = -Ifges(5,2) * t63 + t86;
t41 = -pkin(2) * t69 + t77;
t38 = -mrSges(5,2) * t67 - mrSges(5,3) * t84;
t37 = mrSges(5,1) * t67 + mrSges(5,3) * t85;
t29 = Ifges(6,5) * t72;
t28 = Ifges(6,6) * t32;
t25 = pkin(4) * t84 + t46;
t21 = Ifges(5,5) * t67 + (-Ifges(5,1) * t63 - t86) * t69;
t20 = Ifges(5,6) * t67 + (-Ifges(5,2) * t64 - t87) * t69;
t16 = mrSges(6,1) * t67 - mrSges(6,3) * t23;
t15 = -mrSges(6,2) * t67 + mrSges(6,3) * t22;
t14 = -Ifges(6,1) * t72 + Ifges(6,4) * t32;
t13 = -Ifges(6,4) * t72 + Ifges(6,2) * t32;
t4 = Ifges(6,1) * t23 + Ifges(6,4) * t22 + Ifges(6,5) * t67;
t3 = Ifges(6,4) * t23 + Ifges(6,2) * t22 + Ifges(6,6) * t67;
t17 = [0.2e1 * t1 * t16 + 0.2e1 * t2 * t15 + t22 * t3 + t23 * t4 + 0.2e1 * t25 * t6 + 0.2e1 * t46 * t26 + 0.2e1 * t8 * t37 + 0.2e1 * t9 * t38 + Ifges(2,3) + (0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * t41 * mrSges(4,2) - t64 * t20 - t63 * t21 + (Ifges(3,2) + Ifges(4,3)) * t69) * t69 + (-0.2e1 * pkin(1) * mrSges(3,2) - 0.2e1 * t41 * mrSges(4,3) + (Ifges(4,2) + Ifges(3,1) + Ifges(5,3)) * t67 + (-Ifges(5,5) * t63 - Ifges(5,6) * t64 + (2 * Ifges(3,4)) + (2 * Ifges(4,6))) * t69 + t80) * t67 + m(6) * (t1 ^ 2 + t2 ^ 2 + t25 ^ 2) + m(5) * (t46 ^ 2 + t8 ^ 2 + t9 ^ 2) + m(3) * (pkin(1) ^ 2 + t83) + m(4) * (t41 ^ 2 + t83) + 0.2e1 * (mrSges(4,1) + mrSges(3,3)) * pkin(6) * t91; -t72 * t4 / 0.2e1 + t46 * t42 + t48 * t6 + t23 * t14 / 0.2e1 + t25 * t12 + qJ(3) * t26 + t32 * t3 / 0.2e1 + t11 * t15 + t10 * t16 + t22 * t13 / 0.2e1 - t75 * mrSges(6,3) + (t21 / 0.2e1 - t8 * mrSges(5,3) + t65 * t37) * t64 + (-t20 / 0.2e1 - t9 * mrSges(5,3) + t65 * t38) * t63 + m(6) * (t1 * t10 + t11 * t2 + t25 * t48) + m(5) * (qJ(3) * t46 + t65 * t74) + (Ifges(5,5) * t64 / 0.2e1 + Ifges(5,6) * t89 - t29 / 0.2e1 + t28 / 0.2e1 - Ifges(4,4) + Ifges(3,5) - pkin(2) * mrSges(4,1)) * t67 + (-Ifges(4,5) + Ifges(3,6) + t44 * t89 - t64 * t43 / 0.2e1 + qJ(3) * mrSges(4,1)) * t69 + ((m(4) * qJ(3) - mrSges(3,2) + mrSges(4,3)) * t69 + (-mrSges(3,1) + t90) * t67) * pkin(6); -0.2e1 * pkin(2) * mrSges(4,2) + 0.2e1 * t48 * t12 + t32 * t13 - t72 * t14 - t63 * t43 + t64 * t44 + Ifges(4,1) + Ifges(3,3) + m(6) * (t10 ^ 2 + t11 ^ 2 + t48 ^ 2) + m(5) * (t65 ^ 2 * t82 + t70) + m(4) * (pkin(2) ^ 2 + t70) + 0.2e1 * (mrSges(4,3) + t42) * qJ(3) - 0.2e1 * t73 * mrSges(6,3) - 0.2e1 * t65 * t78; -t32 * t15 - t72 * t16 + t64 * t37 + t63 * t38 + (m(4) * pkin(6) + mrSges(4,1)) * t67 + m(6) * t75 + m(5) * t74; m(6) * t73 - mrSges(6,3) * t81 + t65 * t79 - t78 + t90; m(6) * t81 + m(4) + t79; m(5) * t46 + m(6) * t25 + t26 + t6; m(5) * qJ(3) + m(6) * t48 + t12 + t42; 0; m(5) + m(6); mrSges(6,1) * t1 - mrSges(6,2) * t2 + t80; mrSges(6,1) * t10 - mrSges(6,2) * t11 + t28 - t29; -mrSges(6,1) * t72 + mrSges(6,2) * t32; 0; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t17(1), t17(2), t17(4), t17(7), t17(11); t17(2), t17(3), t17(5), t17(8), t17(12); t17(4), t17(5), t17(6), t17(9), t17(13); t17(7), t17(8), t17(9), t17(10), t17(14); t17(11), t17(12), t17(13), t17(14), t17(15);];
Mq = res;
