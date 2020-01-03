% Calculate joint inertia matrix for
% S5RPRPR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-31 18:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR15_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR15_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR15_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR15_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR15_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR15_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:36:38
% EndTime: 2019-12-31 18:36:40
% DurationCPUTime: 0.45s
% Computational Cost: add. (517->165), mult. (1019->235), div. (0->0), fcn. (898->6), ass. (0->73)
t62 = (-pkin(1) - pkin(6));
t85 = -2 * t62;
t84 = 2 * qJ(2);
t83 = m(5) * pkin(3);
t56 = sin(pkin(8));
t82 = t56 / 0.2e1;
t57 = cos(pkin(8));
t81 = t57 / 0.2e1;
t80 = m(5) + m(6);
t61 = cos(qJ(3));
t75 = t57 * t61;
t76 = t56 * t61;
t26 = mrSges(5,1) * t76 + mrSges(5,2) * t75;
t58 = sin(qJ(5));
t60 = cos(qJ(5));
t34 = t60 * t56 + t58 * t57;
t23 = t34 * t61;
t33 = -t58 * t56 + t60 * t57;
t25 = t33 * t61;
t5 = t23 * mrSges(6,1) + t25 * mrSges(6,2);
t79 = -t26 - t5;
t78 = Ifges(5,4) * t56;
t77 = Ifges(5,4) * t57;
t59 = sin(qJ(3));
t74 = t59 * t62;
t39 = -t57 * mrSges(5,1) + t56 * mrSges(5,2);
t73 = mrSges(4,1) - t39;
t72 = pkin(7) + qJ(4);
t37 = t59 * pkin(3) - t61 * qJ(4) + qJ(2);
t16 = t56 * t37 + t57 * t74;
t71 = t56 ^ 2 + t57 ^ 2;
t53 = t59 ^ 2;
t54 = t61 ^ 2;
t70 = t54 + t53;
t69 = Ifges(6,5) * t25 - Ifges(6,6) * t23 + Ifges(6,3) * t59;
t68 = t70 * mrSges(4,3);
t8 = -t33 * mrSges(6,1) + t34 * mrSges(6,2);
t67 = qJ(4) * t71;
t31 = t57 * t37;
t15 = -t56 * t74 + t31;
t66 = -t15 * t56 + t16 * t57;
t35 = -t59 * mrSges(5,2) - mrSges(5,3) * t76;
t36 = t59 * mrSges(5,1) - mrSges(5,3) * t75;
t65 = t57 * t35 - t56 * t36;
t64 = qJ(2) ^ 2;
t55 = t62 ^ 2;
t48 = t54 * t62;
t47 = t54 * t55;
t46 = -t57 * pkin(4) - pkin(3);
t42 = Ifges(5,1) * t56 + t77;
t41 = Ifges(5,2) * t57 + t78;
t40 = t72 * t57;
t38 = t72 * t56;
t32 = (pkin(4) * t56 - t62) * t61;
t29 = Ifges(6,5) * t34;
t28 = Ifges(6,6) * t33;
t24 = t33 * t59;
t22 = t34 * t59;
t21 = Ifges(5,5) * t59 + (Ifges(5,1) * t57 - t78) * t61;
t20 = Ifges(5,6) * t59 + (-Ifges(5,2) * t56 + t77) * t61;
t14 = t59 * mrSges(6,1) - t25 * mrSges(6,3);
t13 = -t59 * mrSges(6,2) - t23 * mrSges(6,3);
t12 = -t58 * t38 + t60 * t40;
t11 = -t60 * t38 - t58 * t40;
t10 = Ifges(6,1) * t34 + Ifges(6,4) * t33;
t9 = Ifges(6,4) * t34 + Ifges(6,2) * t33;
t7 = -pkin(7) * t76 + t16;
t6 = -pkin(7) * t75 + t31 + (-t56 * t62 + pkin(4)) * t59;
t4 = Ifges(6,1) * t25 - Ifges(6,4) * t23 + Ifges(6,5) * t59;
t3 = Ifges(6,4) * t25 - Ifges(6,2) * t23 + Ifges(6,6) * t59;
t2 = t58 * t6 + t60 * t7;
t1 = -t58 * t7 + t60 * t6;
t17 = [-(2 * pkin(1) * mrSges(3,2)) + (mrSges(3,3) * t84) + 0.2e1 * t1 * t14 + 0.2e1 * t2 * t13 + 0.2e1 * t15 * t36 + 0.2e1 * t16 * t35 - t23 * t3 + t25 * t4 + 0.2e1 * t32 * t5 + Ifges(3,1) + Ifges(2,3) + t68 * t85 + (mrSges(4,1) * t84 + (Ifges(4,2) + Ifges(5,3)) * t59 + t69) * t59 + ((mrSges(4,2) * t84) + Ifges(4,1) * t61 - t56 * t20 + t57 * t21 + t26 * t85 + (Ifges(5,5) * t57 - Ifges(5,6) * t56 - (2 * Ifges(4,4))) * t59) * t61 + m(6) * (t1 ^ 2 + t2 ^ 2 + t32 ^ 2) + m(5) * (t15 ^ 2 + t16 ^ 2 + t47) + (m(3) * (pkin(1) ^ 2 + t64)) + m(4) * (t53 * t55 + t47 + t64); -(m(3) * pkin(1)) + t24 * t13 - t22 * t14 + mrSges(3,2) + t79 * t61 + t65 * t59 - t68 + m(6) * (-t22 * t1 + t24 * t2 - t61 * t32) + m(5) * (t66 * t59 + t48) + m(4) * (t53 * t62 + t48); m(3) + m(4) * t70 + m(5) * (t71 * t53 + t54) + m(6) * (t22 ^ 2 + t24 ^ 2 + t54); t20 * t81 + t21 * t82 + t33 * t3 / 0.2e1 + t34 * t4 / 0.2e1 + t46 * t5 - t23 * t9 / 0.2e1 + t25 * t10 / 0.2e1 - pkin(3) * t26 + t32 * t8 + t12 * t13 + t11 * t14 + m(6) * (t11 * t1 + t12 * t2 + t46 * t32) + (Ifges(5,5) * t82 + Ifges(5,6) * t81 + t29 / 0.2e1 + t28 / 0.2e1 - Ifges(4,6) - (t62 * mrSges(4,2))) * t59 + (-t1 * t34 + t2 * t33) * mrSges(6,3) + t66 * mrSges(5,3) + (m(5) * t66 + t65) * qJ(4) + (Ifges(4,5) - t56 * t41 / 0.2e1 + t42 * t81 + (t73 + t83) * t62) * t61; (t22 * t34 + t24 * t33) * mrSges(6,3) + (-t8 + t73) * t61 + (t71 * mrSges(5,3) - mrSges(4,2)) * t59 + m(5) * (pkin(3) * t61 + t59 * t67) + m(6) * (-t11 * t22 + t12 * t24 - t46 * t61); -0.2e1 * pkin(3) * t39 + t34 * t10 + t33 * t9 + t57 * t41 + t56 * t42 + 0.2e1 * t46 * t8 + Ifges(4,3) + m(6) * (t11 ^ 2 + t12 ^ 2 + t46 ^ 2) + m(5) * (t71 * qJ(4) ^ 2 + pkin(3) ^ 2) + 0.2e1 * (-t11 * t34 + t12 * t33) * mrSges(6,3) + 0.2e1 * mrSges(5,3) * t67; -m(5) * t61 * t62 + m(6) * t32 - t79; -t80 * t61; m(6) * t46 + t39 + t8 - t83; t80; t1 * mrSges(6,1) - t2 * mrSges(6,2) + t69; -t22 * mrSges(6,1) - t24 * mrSges(6,2); t11 * mrSges(6,1) - t12 * mrSges(6,2) + t28 + t29; 0; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t17(1), t17(2), t17(4), t17(7), t17(11); t17(2), t17(3), t17(5), t17(8), t17(12); t17(4), t17(5), t17(6), t17(9), t17(13); t17(7), t17(8), t17(9), t17(10), t17(14); t17(11), t17(12), t17(13), t17(14), t17(15);];
Mq = res;
