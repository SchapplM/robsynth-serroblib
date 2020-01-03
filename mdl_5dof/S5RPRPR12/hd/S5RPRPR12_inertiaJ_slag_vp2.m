% Calculate joint inertia matrix for
% S5RPRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-31 18:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR12_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR12_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR12_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR12_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR12_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR12_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:29:27
% EndTime: 2019-12-31 18:29:28
% DurationCPUTime: 0.45s
% Computational Cost: add. (854->159), mult. (1674->238), div. (0->0), fcn. (1795->8), ass. (0->67)
t64 = cos(pkin(8));
t78 = pkin(6) + qJ(2);
t49 = t78 * t64;
t66 = sin(qJ(3));
t62 = sin(pkin(8));
t72 = t78 * t62;
t83 = cos(qJ(3));
t31 = t66 * t49 + t83 * t72;
t88 = t31 ^ 2;
t60 = t64 ^ 2;
t87 = 0.2e1 * t31;
t54 = -t64 * pkin(2) - pkin(1);
t86 = 0.2e1 * t54;
t63 = cos(pkin(9));
t84 = t63 / 0.2e1;
t61 = sin(pkin(9));
t82 = Ifges(5,4) * t61;
t81 = Ifges(5,4) * t63;
t45 = t83 * t62 + t66 * t64;
t80 = t45 * t61;
t79 = t45 * t63;
t77 = pkin(7) + qJ(4);
t43 = t66 * t62 - t83 * t64;
t23 = t43 * pkin(3) - t45 * qJ(4) + t54;
t34 = t83 * t49 - t66 * t72;
t9 = t61 * t23 + t63 * t34;
t22 = mrSges(5,1) * t80 + mrSges(5,2) * t79;
t65 = sin(qJ(5));
t67 = cos(qJ(5));
t42 = -t65 * t61 + t67 * t63;
t44 = t67 * t61 + t65 * t63;
t76 = Ifges(6,5) * t44 + Ifges(6,6) * t42;
t75 = t61 ^ 2 + t63 ^ 2;
t74 = t62 ^ 2 + t60;
t18 = t44 * t45;
t19 = t42 * t45;
t73 = Ifges(6,5) * t19 - Ifges(6,6) * t18 + Ifges(6,3) * t43;
t71 = -t64 * mrSges(3,1) + t62 * mrSges(3,2);
t47 = -t63 * mrSges(5,1) + t61 * mrSges(5,2);
t7 = t18 * mrSges(6,1) + t19 * mrSges(6,2);
t8 = t63 * t23 - t61 * t34;
t70 = -t8 * t61 + t9 * t63;
t27 = -t42 * mrSges(6,1) + t44 * mrSges(6,2);
t53 = -t63 * pkin(4) - pkin(3);
t51 = Ifges(5,1) * t61 + t81;
t50 = Ifges(5,2) * t63 + t82;
t48 = t77 * t63;
t46 = t77 * t61;
t37 = t45 * mrSges(4,2);
t33 = -t65 * t46 + t67 * t48;
t30 = -t67 * t46 - t65 * t48;
t29 = Ifges(6,1) * t44 + Ifges(6,4) * t42;
t28 = Ifges(6,4) * t44 + Ifges(6,2) * t42;
t25 = t43 * mrSges(5,1) - mrSges(5,3) * t79;
t24 = -t43 * mrSges(5,2) - mrSges(5,3) * t80;
t14 = pkin(4) * t80 + t31;
t13 = Ifges(5,5) * t43 + (Ifges(5,1) * t63 - t82) * t45;
t12 = Ifges(5,6) * t43 + (-Ifges(5,2) * t61 + t81) * t45;
t11 = t43 * mrSges(6,1) - t19 * mrSges(6,3);
t10 = -t43 * mrSges(6,2) - t18 * mrSges(6,3);
t6 = -pkin(7) * t80 + t9;
t5 = Ifges(6,1) * t19 - Ifges(6,4) * t18 + Ifges(6,5) * t43;
t4 = Ifges(6,4) * t19 - Ifges(6,2) * t18 + Ifges(6,6) * t43;
t3 = t43 * pkin(4) - pkin(7) * t79 + t8;
t2 = t65 * t3 + t67 * t6;
t1 = t67 * t3 - t65 * t6;
t15 = [Ifges(3,2) * t60 - 0.2e1 * pkin(1) * t71 + t37 * t86 + t22 * t87 - t18 * t4 + t19 * t5 + 0.2e1 * t9 * t24 + 0.2e1 * t8 * t25 + 0.2e1 * t2 * t10 + 0.2e1 * t1 * t11 + 0.2e1 * t14 * t7 + Ifges(2,3) + (Ifges(3,1) * t62 + 0.2e1 * Ifges(3,4) * t64) * t62 + 0.2e1 * t74 * qJ(2) * mrSges(3,3) + (mrSges(4,1) * t86 - 0.2e1 * t34 * mrSges(4,3) + (Ifges(5,3) + Ifges(4,2)) * t43 + t73) * t43 + m(6) * (t1 ^ 2 + t14 ^ 2 + t2 ^ 2) + m(5) * (t8 ^ 2 + t9 ^ 2 + t88) + m(4) * (t34 ^ 2 + t54 ^ 2 + t88) + m(3) * (t74 * qJ(2) ^ 2 + pkin(1) ^ 2) + (mrSges(4,3) * t87 + Ifges(4,1) * t45 - t61 * t12 + t63 * t13 + (Ifges(5,5) * t63 - Ifges(5,6) * t61 - (2 * Ifges(4,4))) * t43) * t45; -m(3) * pkin(1) + t43 * mrSges(4,1) + t44 * t10 + t42 * t11 + t61 * t24 + t63 * t25 + t37 + m(6) * (t42 * t1 + t44 * t2) + m(5) * (t61 * t9 + t63 * t8) + m(4) * t54 + t71; m(3) + m(4) + m(5) * t75 + m(6) * (t42 ^ 2 + t44 ^ 2); t61 * t13 / 0.2e1 + t12 * t84 + t53 * t7 - Ifges(4,6) * t43 + t44 * t5 / 0.2e1 + t33 * t10 - t34 * mrSges(4,2) + t42 * t4 / 0.2e1 - pkin(3) * t22 + t14 * t27 - t18 * t28 / 0.2e1 + t19 * t29 / 0.2e1 + t30 * t11 + (t47 - mrSges(4,1)) * t31 + (t63 * t24 - t61 * t25) * qJ(4) + (-t1 * t44 + t2 * t42) * mrSges(6,3) + t70 * mrSges(5,3) + m(6) * (t30 * t1 + t53 * t14 + t33 * t2) + m(5) * (-pkin(3) * t31 + qJ(4) * t70) + (Ifges(4,5) - t61 * t50 / 0.2e1 + t51 * t84) * t45 + (Ifges(5,5) * t61 + Ifges(5,6) * t63 + t76) * t43 / 0.2e1; m(6) * (t30 * t42 + t33 * t44); -0.2e1 * pkin(3) * t47 + 0.2e1 * t53 * t27 + t42 * t28 + t44 * t29 + t63 * t50 + t61 * t51 + Ifges(4,3) + m(6) * (t30 ^ 2 + t33 ^ 2 + t53 ^ 2) + m(5) * (t75 * qJ(4) ^ 2 + pkin(3) ^ 2) + 0.2e1 * (-t30 * t44 + t33 * t42) * mrSges(6,3) + 0.2e1 * t75 * qJ(4) * mrSges(5,3); m(5) * t31 + m(6) * t14 + t22 + t7; 0; -m(5) * pkin(3) + m(6) * t53 + t27 + t47; m(5) + m(6); t1 * mrSges(6,1) - t2 * mrSges(6,2) + t73; -t27; t30 * mrSges(6,1) - t33 * mrSges(6,2) + t76; 0; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t15(1), t15(2), t15(4), t15(7), t15(11); t15(2), t15(3), t15(5), t15(8), t15(12); t15(4), t15(5), t15(6), t15(9), t15(13); t15(7), t15(8), t15(9), t15(10), t15(14); t15(11), t15(12), t15(13), t15(14), t15(15);];
Mq = res;
