% Calculate time derivative of joint inertia matrix for
% S5RPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR5_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR5_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR5_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR5_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:16:08
% EndTime: 2019-12-05 18:16:10
% DurationCPUTime: 0.59s
% Computational Cost: add. (1145->138), mult. (2503->206), div. (0->0), fcn. (1944->8), ass. (0->78)
t58 = sin(qJ(5));
t59 = sin(qJ(4));
t61 = cos(qJ(5));
t62 = cos(qJ(4));
t40 = t58 * t62 + t59 * t61;
t99 = qJD(4) + qJD(5);
t21 = t99 * t40;
t39 = -t58 * t59 + t61 * t62;
t103 = t39 * t21;
t20 = t99 * t39;
t102 = t40 * t20;
t51 = cos(pkin(9)) * pkin(1) + pkin(2);
t60 = sin(qJ(3));
t63 = cos(qJ(3));
t95 = pkin(1) * sin(pkin(9));
t69 = t63 * t51 - t60 * t95;
t29 = t69 * qJD(3);
t84 = t60 * t51 + t63 * t95;
t34 = pkin(7) + t84;
t92 = -pkin(8) - t34;
t68 = qJD(4) * t92;
t14 = t62 * t29 + t59 * t68;
t15 = -t59 * t29 + t62 * t68;
t23 = t92 * t59;
t54 = t62 * pkin(8);
t24 = t34 * t62 + t54;
t8 = t23 * t61 - t24 * t58;
t2 = qJD(5) * t8 + t61 * t14 + t58 * t15;
t9 = t23 * t58 + t24 * t61;
t3 = -qJD(5) * t9 - t58 * t14 + t61 * t15;
t101 = t3 * mrSges(6,1) - t2 * mrSges(6,2);
t96 = -pkin(8) - pkin(7);
t47 = t96 * t59;
t48 = pkin(7) * t62 + t54;
t26 = t47 * t61 - t48 * t58;
t73 = qJD(4) * t96;
t42 = t59 * t73;
t43 = t62 * t73;
t12 = qJD(5) * t26 + t61 * t42 + t58 * t43;
t27 = t47 * t58 + t48 * t61;
t13 = -qJD(5) * t27 - t58 * t42 + t61 * t43;
t100 = t13 * mrSges(6,1) - t12 * mrSges(6,2);
t67 = mrSges(5,1) * t59 + mrSges(5,2) * t62;
t41 = t67 * qJD(4);
t98 = 2 * m(6);
t6 = mrSges(6,1) * t21 + t20 * mrSges(6,2);
t97 = 0.2e1 * t6;
t93 = t62 * pkin(4);
t91 = Ifges(5,4) * t59;
t89 = Ifges(5,6) * t59;
t22 = -mrSges(6,1) * t39 + mrSges(6,2) * t40;
t30 = t84 * qJD(3);
t81 = qJD(4) * t59;
t74 = pkin(4) * t81;
t25 = t30 + t74;
t87 = t25 * t22;
t85 = Ifges(6,5) * t20 - Ifges(6,6) * t21;
t83 = t59 ^ 2 + t62 ^ 2;
t82 = pkin(4) * qJD(5);
t80 = qJD(4) * t62;
t79 = qJD(5) * t58;
t78 = qJD(5) * t61;
t77 = 2 * mrSges(6,3);
t76 = t61 * t20 * mrSges(6,3);
t75 = mrSges(6,3) * t82;
t44 = -mrSges(5,1) * t62 + mrSges(5,2) * t59;
t72 = (-mrSges(4,1) + t44) * t30;
t71 = t83 * mrSges(5,3);
t70 = t83 * t29;
t33 = -pkin(3) - t69;
t66 = t61 * t39 * t75 + Ifges(5,5) * t80 + t85 + (-mrSges(6,3) * pkin(4) * t21 + t40 * t75) * t58;
t65 = 0.2e1 * Ifges(6,1) * t102 - 0.2e1 * Ifges(6,2) * t103 + (Ifges(5,1) * t62 - t91) * t81 + (0.2e1 * Ifges(5,4) * t62 + (Ifges(5,1) - Ifges(5,2)) * t59) * t80 + 0.2e1 * (t20 * t39 - t21 * t40) * Ifges(6,4);
t45 = Ifges(5,2) * t62 + t91;
t64 = -t45 * t81 + t65;
t52 = -pkin(3) - t93;
t38 = (-mrSges(6,1) * t58 - mrSges(6,2) * t61) * t82;
t28 = t33 - t93;
t1 = [(t2 * t39 - t8 * t20 - t9 * t21 - t3 * t40) * t77 - 0.2e1 * t29 * mrSges(4,2) + 0.2e1 * m(5) * (t33 * t30 + t34 * t70) + 0.2e1 * m(4) * (t29 * t84 - t30 * t69) + (t2 * t9 + t25 * t28 + t3 * t8) * t98 + t64 + 0.2e1 * t33 * t41 + 0.2e1 * t87 + t28 * t97 + 0.2e1 * t71 * t29 + 0.2e1 * t72; m(6) * (t2 * t40 + t20 * t9 - t21 * t8 + t3 * t39); (t102 - t103) * t98; t87 + (t52 + t28) * t6 + (-pkin(3) + t33) * t41 + t72 + (pkin(4) * t22 - t45) * t81 + (-mrSges(4,2) + t71) * t29 + m(6) * (t12 * t9 + t13 * t8 + t2 * t27 + t25 * t52 + t26 * t3 + t28 * t74) + m(5) * (-pkin(3) * t30 + pkin(7) * t70) + ((-t13 - t3) * t40 + (t12 + t2) * t39 - (t27 + t9) * t21 + (-t26 - t8) * t20) * mrSges(6,3) + t65; m(6) * (t12 * t40 + t13 * t39 + t20 * t27 - t21 * t26); 0.2e1 * t22 * t74 + t52 * t97 - 0.2e1 * pkin(3) * t41 + (t12 * t27 + t13 * t26 + t52 * t74) * t98 + (t12 * t39 - t13 * t40 - t26 * t20 - t27 * t21) * t77 + t64; -t67 * t29 + (t34 * t44 - t89) * qJD(4) + (-t76 + m(6) * (t2 * t58 + t3 * t61 + t78 * t9 - t79 * t8)) * pkin(4) + t66 + t101; -t41 + m(6) * (t20 * t58 - t21 * t61 + (-t39 * t58 + t40 * t61) * qJD(5)) * pkin(4) - t6; (pkin(7) * t44 - t89) * qJD(4) + (-t76 + m(6) * (t12 * t58 + t13 * t61 - t26 * t79 + t27 * t78)) * pkin(4) + t66 + t100; 0.2e1 * t38; t85 + t101; -t6; t85 + t100; t38; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
