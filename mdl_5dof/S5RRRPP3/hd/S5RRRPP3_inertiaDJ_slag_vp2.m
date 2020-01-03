% Calculate time derivative of joint inertia matrix for
% S5RRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPP3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP3_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP3_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP3_inertiaDJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP3_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP3_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP3_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:53:09
% EndTime: 2019-12-31 20:53:11
% DurationCPUTime: 0.78s
% Computational Cost: add. (448->167), mult. (1094->198), div. (0->0), fcn. (553->4), ass. (0->75)
t59 = sin(qJ(3));
t61 = cos(qJ(3));
t111 = t59 ^ 2 + t61 ^ 2;
t110 = 0.2e1 * t61;
t109 = Ifges(4,1) + Ifges(6,3);
t108 = -Ifges(6,2) - Ifges(5,3);
t62 = cos(qJ(2));
t89 = pkin(1) * qJD(2);
t80 = t62 * t89;
t107 = t111 * t80;
t106 = (-Ifges(4,4) - Ifges(5,6)) * t59;
t105 = (-mrSges(3,2) + (mrSges(5,1) + mrSges(4,3)) * t111) * t80;
t104 = 2 * m(6);
t21 = (-mrSges(6,2) * t61 + mrSges(6,3) * t59) * qJD(3);
t103 = 0.2e1 * t21;
t30 = t61 * mrSges(5,2) - t59 * mrSges(5,3);
t102 = 0.2e1 * t30;
t32 = -mrSges(6,2) * t59 - mrSges(6,3) * t61;
t101 = 0.2e1 * t32;
t100 = pkin(1) * t62;
t94 = Ifges(6,6) * t61;
t58 = -pkin(3) - qJ(5);
t60 = sin(qJ(2));
t41 = pkin(1) * t60 + pkin(7);
t93 = t107 * t41;
t92 = t107 * pkin(7);
t86 = qJD(3) * t61;
t45 = mrSges(6,1) * t86;
t91 = mrSges(5,1) * t86 + t45;
t88 = t59 * qJ(4);
t87 = qJD(3) * t59;
t85 = qJD(4) * t61;
t84 = qJ(4) * qJD(3);
t83 = qJ(4) * qJD(4);
t82 = 0.2e1 * mrSges(6,1);
t81 = m(5) * t85;
t79 = pkin(7) * t86;
t53 = t60 * t89;
t78 = mrSges(6,1) * t87;
t77 = ((-Ifges(4,2) + Ifges(5,2)) * t61 + t106) * qJD(3);
t76 = pkin(3) * t87 - qJD(4) * t59;
t72 = (-t94 + Ifges(5,6) * t110 + (Ifges(5,2) + t108) * t59) * qJD(3);
t31 = -t61 * mrSges(4,1) + t59 * mrSges(4,2);
t71 = t59 * mrSges(4,1) + t61 * mrSges(4,2);
t70 = -t59 * mrSges(5,2) - t61 * mrSges(5,3);
t69 = -t61 * pkin(3) - t88;
t68 = (-mrSges(3,1) + t31) * t53;
t29 = -pkin(2) + t69;
t67 = (-t94 + Ifges(4,4) * t110 + (-Ifges(4,2) + t109) * t59) * t86 + (0.2e1 * Ifges(6,6) * t59 + (t108 + t109) * t61 + t106) * t87;
t11 = t58 * t61 - pkin(2) - t88;
t7 = -t61 * t84 + t76;
t66 = t41 * t86 + t59 * t80;
t65 = mrSges(5,1) * t85 + t58 * t45 + (Ifges(6,4) + Ifges(5,5)) * t87 + (Ifges(4,5) + Ifges(6,5)) * t86 + (-qJD(5) * t59 + t85) * mrSges(6,1);
t2 = qJ(5) * t87 + (-qJD(5) - t84) * t61 + t76;
t64 = (-Ifges(4,6) + (-mrSges(5,1) - mrSges(6,1)) * qJ(4)) * t59 + (-pkin(3) * mrSges(5,1) - Ifges(5,4)) * t61;
t63 = m(5) * t69 + t30 + t31;
t55 = t61 * pkin(4);
t54 = t59 * pkin(4);
t52 = pkin(4) * t86;
t42 = -pkin(2) - t100;
t36 = pkin(7) * t61 + t55;
t35 = pkin(7) * t59 + t54;
t28 = t52 + t79;
t27 = (-pkin(4) - pkin(7)) * t87;
t23 = t71 * qJD(3);
t22 = t70 * qJD(3);
t18 = t41 * t61 + t55;
t17 = t41 * t59 + t54;
t12 = t29 - t100;
t6 = t11 - t100;
t5 = t53 + t7;
t4 = t52 + t66;
t3 = t61 * t80 + (-pkin(4) - t41) * t87;
t1 = t53 + t2;
t8 = [0.2e1 * t105 + 0.2e1 * m(4) * (t42 * t53 + t93) + (t1 * t6 + t17 * t4 + t18 * t3) * t104 + 0.2e1 * m(5) * (t12 * t5 + t93) + t67 + t5 * t102 + t1 * t101 + 0.2e1 * t42 * t23 + t6 * t103 + 0.2e1 * t12 * t22 + 0.2e1 * t68 + ((qJD(3) * t17 + t3) * t82 + t72) * t61 + ((-qJD(3) * t18 + t4) * t82 + t77) * t59; t68 + t105 + t72 * t61 + t77 * t59 + m(4) * (-pkin(2) * t53 + t92) + m(6) * (t1 * t11 + t17 * t28 + t18 * t27 + t2 * t6 + t3 * t36 + t35 * t4) + m(5) * (t12 * t7 + t29 * t5 + t92) + (t2 + t1) * t32 + (t7 + t5) * t30 + (t42 - pkin(2)) * t23 + (t29 + t12) * t22 + (t11 + t6) * t21 + t67 + ((t27 + t3) * t61 + (t28 + t4) * t59 + ((t17 + t35) * t61 + (-t18 - t36) * t59) * qJD(3)) * mrSges(6,1); t7 * t102 - 0.2e1 * pkin(2) * t23 + (t11 * t2 + t27 * t36 + t28 * t35) * t104 + t2 * t101 + t11 * t103 + 0.2e1 * (m(5) * t7 + t22) * t29 + ((-qJD(3) * t36 + t28) * t82 + t77) * t59 + ((qJD(3) * t35 + t27) * t82 + t72) * t61 + t67; t3 * mrSges(6,2) - t4 * mrSges(6,3) + t41 * t81 + m(6) * (qJ(4) * t3 + qJD(4) * t18 - qJD(5) * t17 + t4 * t58) + (m(5) * (-pkin(3) * t59 + qJ(4) * t61) - t70 - t71) * t80 + (t63 * t41 + t64) * qJD(3) + t65; t27 * mrSges(6,2) - t28 * mrSges(6,3) + pkin(7) * t81 + m(6) * (qJ(4) * t27 + qJD(4) * t36 - qJD(5) * t35 + t28 * t58) + (t63 * pkin(7) + t64) * qJD(3) + t65; 0.2e1 * m(5) * t83 + 0.2e1 * qJD(5) * mrSges(6,3) + 0.2e1 * m(6) * (-qJD(5) * t58 + t83) + 0.2e1 * (mrSges(5,3) + mrSges(6,2)) * qJD(4); t66 * m(5) + m(6) * t4 + t91; m(5) * t79 + m(6) * t28 + t91; -m(6) * qJD(5); 0; m(6) * t3 - t78; m(6) * t27 - t78; m(6) * qJD(4); 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t8(1), t8(2), t8(4), t8(7), t8(11); t8(2), t8(3), t8(5), t8(8), t8(12); t8(4), t8(5), t8(6), t8(9), t8(13); t8(7), t8(8), t8(9), t8(10), t8(14); t8(11), t8(12), t8(13), t8(14), t8(15);];
Mq = res;
