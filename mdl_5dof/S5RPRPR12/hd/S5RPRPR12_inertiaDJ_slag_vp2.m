% Calculate time derivative of joint inertia matrix for
% S5RPRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR12_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR12_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR12_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR12_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR12_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR12_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR12_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:29:27
% EndTime: 2019-12-31 18:29:30
% DurationCPUTime: 1.20s
% Computational Cost: add. (1949->211), mult. (4539->336), div. (0->0), fcn. (4352->8), ass. (0->91)
t101 = pkin(6) + qJ(2);
t85 = sin(pkin(8));
t75 = t101 * t85;
t87 = cos(pkin(8));
t77 = t101 * t87;
t89 = sin(qJ(3));
t91 = cos(qJ(3));
t116 = -t91 * t75 - t77 * t89;
t84 = sin(pkin(9));
t86 = cos(pkin(9));
t88 = sin(qJ(5));
t90 = cos(qJ(5));
t92 = t84 * t88 - t86 * t90;
t65 = t92 * qJD(5);
t115 = 2 * m(5);
t114 = 2 * m(6);
t82 = t86 ^ 2;
t112 = Ifges(5,4) * t84;
t111 = Ifges(5,4) * t86;
t110 = Ifges(5,2) * t84;
t58 = -t89 * t75 + t77 * t91;
t73 = t85 * t91 + t89 * t87;
t37 = t73 * qJD(2) + t58 * qJD(3);
t109 = t37 * t116;
t68 = t73 * qJD(3);
t108 = t68 * Ifges(5,5);
t107 = t68 * Ifges(5,6);
t106 = t73 * t84;
t105 = t73 * t86;
t71 = t85 * t89 - t91 * t87;
t67 = t71 * qJD(3);
t103 = t84 * t67;
t102 = t86 * t67;
t100 = pkin(7) + qJ(4);
t32 = pkin(3) * t68 + qJ(4) * t67 - qJD(4) * t73;
t36 = -t71 * qJD(2) + t116 * qJD(3);
t13 = t84 * t32 + t86 * t36;
t96 = -pkin(2) * t87 - pkin(1);
t49 = pkin(3) * t71 - qJ(4) * t73 + t96;
t24 = t84 * t49 + t86 * t58;
t39 = -mrSges(5,1) * t103 - mrSges(5,2) * t102;
t72 = t84 * t90 + t86 * t88;
t66 = t72 * qJD(5);
t99 = -Ifges(6,5) * t65 - Ifges(6,6) * t66;
t21 = -t73 * t66 + t92 * t67;
t22 = t73 * t65 + t72 * t67;
t97 = Ifges(6,5) * t21 + Ifges(6,6) * t22 + Ifges(6,3) * t68;
t95 = t68 * mrSges(4,1) - t67 * mrSges(4,2);
t7 = -t22 * mrSges(6,1) + t21 * mrSges(6,2);
t12 = t86 * t32 - t36 * t84;
t23 = t86 * t49 - t58 * t84;
t93 = Ifges(5,5) * t86 - Ifges(5,6) * t84;
t14 = pkin(4) * t71 - pkin(7) * t105 + t23;
t17 = -pkin(7) * t106 + t24;
t3 = t14 * t90 - t17 * t88;
t4 = t14 * t88 + t17 * t90;
t74 = t100 * t84;
t76 = t100 * t86;
t55 = -t74 * t90 - t76 * t88;
t57 = -t74 * t88 + t76 * t90;
t79 = -pkin(4) * t86 - pkin(3);
t54 = Ifges(6,1) * t72 - Ifges(6,4) * t92;
t53 = Ifges(6,4) * t72 - Ifges(6,2) * t92;
t51 = mrSges(5,1) * t71 - mrSges(5,3) * t105;
t50 = -mrSges(5,2) * t71 - mrSges(5,3) * t106;
t48 = -Ifges(6,1) * t65 - Ifges(6,4) * t66;
t47 = -Ifges(6,4) * t65 - Ifges(6,2) * t66;
t46 = mrSges(6,1) * t66 - mrSges(6,2) * t65;
t45 = mrSges(5,1) * t68 + mrSges(5,3) * t102;
t44 = -mrSges(5,2) * t68 + mrSges(5,3) * t103;
t41 = t92 * t73;
t40 = t72 * t73;
t38 = pkin(4) * t106 - t116;
t35 = -t72 * qJD(4) - t57 * qJD(5);
t34 = -t92 * qJD(4) + t55 * qJD(5);
t29 = mrSges(6,1) * t71 + mrSges(6,3) * t41;
t28 = -mrSges(6,2) * t71 - mrSges(6,3) * t40;
t27 = t108 - (Ifges(5,1) * t86 - t112) * t67;
t26 = t107 - (-t110 + t111) * t67;
t25 = -pkin(4) * t103 + t37;
t16 = -Ifges(6,1) * t41 - Ifges(6,4) * t40 + Ifges(6,5) * t71;
t15 = -Ifges(6,4) * t41 - Ifges(6,2) * t40 + Ifges(6,6) * t71;
t11 = -mrSges(6,2) * t68 + mrSges(6,3) * t22;
t10 = mrSges(6,1) * t68 - mrSges(6,3) * t21;
t9 = pkin(7) * t103 + t13;
t8 = pkin(4) * t68 + pkin(7) * t102 + t12;
t6 = Ifges(6,1) * t21 + Ifges(6,4) * t22 + t68 * Ifges(6,5);
t5 = Ifges(6,4) * t21 + Ifges(6,2) * t22 + t68 * Ifges(6,6);
t2 = -t4 * qJD(5) + t8 * t90 - t88 * t9;
t1 = t3 * qJD(5) + t8 * t88 + t9 * t90;
t18 = [-0.2e1 * t116 * t39 + 0.2e1 * t96 * t95 - t40 * t5 - t41 * t6 + 0.2e1 * t24 * t44 + 0.2e1 * t23 * t45 + 0.2e1 * t13 * t50 + 0.2e1 * t12 * t51 + 0.2e1 * t38 * t7 + t22 * t15 + 0.2e1 * t1 * t28 + 0.2e1 * t2 * t29 + 0.2e1 * t4 * t11 + t21 * t16 + 0.2e1 * t3 * t10 + (t1 * t4 + t2 * t3 + t25 * t38) * t114 + (t12 * t23 + t13 * t24 - t109) * t115 + (-0.2e1 * t36 * mrSges(4,3) + ((2 * Ifges(4,2)) + (2 * Ifges(5,3)) + Ifges(6,3)) * t68 - 0.2e1 * (-Ifges(4,4) + t93) * t67 + t97) * t71 + 0.2e1 * (t116 * t67 - t58 * t68) * mrSges(4,3) + (-t84 * t26 + t86 * t27 + (-0.2e1 * Ifges(4,4) + t93) * t68 - (Ifges(5,1) * t82 + (2 * Ifges(4,1)) + (t110 - 0.2e1 * t111) * t84) * t67 + 0.2e1 * (mrSges(5,1) * t84 + mrSges(5,2) * t86 + mrSges(4,3)) * t37) * t73 + 0.2e1 * t25 * (t40 * mrSges(6,1) - t41 * mrSges(6,2)) + t68 * (-Ifges(6,5) * t41 - Ifges(6,6) * t40) + 0.2e1 * m(4) * (t36 * t58 - t109) + 0.2e1 * (m(3) * qJ(2) + mrSges(3,3)) * (t85 ^ 2 + t87 ^ 2) * qJD(2); -t92 * t10 + t72 * t11 - t65 * t28 - t66 * t29 + t84 * t44 + t86 * t45 + m(6) * (t1 * t72 - t2 * t92 - t3 * t66 - t4 * t65) + m(5) * (t12 * t86 + t13 * t84) + t95; (-t65 * t72 + t66 * t92) * t114; m(6) * (t1 * t57 + t2 * t55 + t25 * t79 + t3 * t35 + t34 * t4) + t79 * t7 + t72 * t6 / 0.2e1 - t65 * t16 / 0.2e1 - t66 * t15 / 0.2e1 - Ifges(4,6) * t68 + t22 * t53 / 0.2e1 + t21 * t54 / 0.2e1 + t55 * t10 + t57 * t11 + t38 * t46 - t40 * t47 / 0.2e1 - t41 * t48 / 0.2e1 + t34 * t28 + t35 * t29 - t36 * mrSges(4,2) - t37 * mrSges(4,1) - pkin(3) * t39 - t92 * t5 / 0.2e1 + t68 * (Ifges(6,5) * t72 - Ifges(6,6) * t92) / 0.2e1 + t25 * (mrSges(6,1) * t92 + mrSges(6,2) * t72) + (-t1 * t92 - t2 * t72 + t3 * t65 - t4 * t66) * mrSges(6,3) + t71 * t99 / 0.2e1 + (qJ(4) * t44 + qJD(4) * t50 + t13 * mrSges(5,3) + t107 / 0.2e1 - t37 * mrSges(5,1) + t26 / 0.2e1) * t86 + (t27 / 0.2e1 + t108 / 0.2e1 + t37 * mrSges(5,2) - qJ(4) * t45 - qJD(4) * t51 - t12 * mrSges(5,3)) * t84 - (t86 * (Ifges(5,1) * t84 + t111) / 0.2e1 - t84 * (Ifges(5,2) * t86 + t112) / 0.2e1 + Ifges(4,5)) * t67 + m(5) * (-pkin(3) * t37 + (-t23 * t84 + t24 * t86) * qJD(4) + (-t12 * t84 + t13 * t86) * qJ(4)); m(6) * (t34 * t72 - t35 * t92 - t55 * t66 - t57 * t65); (t34 * t57 + t35 * t55) * t114 - t65 * t54 + t72 * t48 - t66 * t53 - t92 * t47 + 0.2e1 * t79 * t46 + 0.2e1 * (-t34 * t92 - t35 * t72 + t55 * t65 - t57 * t66) * mrSges(6,3) + (qJ(4) * t115 + 0.2e1 * mrSges(5,3)) * qJD(4) * (t84 ^ 2 + t82); m(5) * t37 + m(6) * t25 + t39 + t7; 0; t46; 0; mrSges(6,1) * t2 - mrSges(6,2) * t1 + t97; -t46; mrSges(6,1) * t35 - mrSges(6,2) * t34 + t99; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t18(1), t18(2), t18(4), t18(7), t18(11); t18(2), t18(3), t18(5), t18(8), t18(12); t18(4), t18(5), t18(6), t18(9), t18(13); t18(7), t18(8), t18(9), t18(10), t18(14); t18(11), t18(12), t18(13), t18(14), t18(15);];
Mq = res;
