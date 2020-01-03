% Calculate time derivative of joint inertia matrix for
% S5RPRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR14_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR14_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR14_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR14_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR14_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR14_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR14_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:34:32
% EndTime: 2019-12-31 18:34:35
% DurationCPUTime: 1.16s
% Computational Cost: add. (1128->177), mult. (2315->278), div. (0->0), fcn. (1969->6), ass. (0->89)
t106 = cos(qJ(3));
t49 = sin(pkin(8));
t51 = sin(qJ(3));
t85 = cos(pkin(8));
t30 = t106 * t49 + t51 * t85;
t27 = t30 * qJD(3);
t64 = t85 * t106;
t29 = t49 * t51 - t64;
t52 = cos(qJ(5));
t50 = sin(qJ(5));
t83 = qJD(5) * t50;
t57 = -t52 * t27 + t29 * t83;
t95 = t27 * t29;
t118 = -Ifges(4,1) + Ifges(4,2);
t82 = qJD(5) * t52;
t58 = t50 * t27 + t29 * t82;
t68 = qJD(3) * t106;
t84 = qJD(3) * t51;
t117 = -mrSges(4,1) * t84 - mrSges(4,2) * t68;
t43 = t51 * pkin(3) + qJ(2);
t16 = pkin(4) * t30 + pkin(7) * t29 + t43;
t53 = -pkin(1) - pkin(6);
t86 = qJ(4) - t53;
t34 = t86 * t51;
t56 = t86 * t106;
t20 = -t34 * t85 - t49 * t56;
t6 = t16 * t52 - t20 * t50;
t7 = t16 * t50 + t20 * t52;
t116 = -t6 * t82 - t7 * t83;
t26 = -qJD(3) * t64 + t49 * t84;
t115 = t26 * t49 + t85 * t27;
t8 = -mrSges(6,1) * t26 - mrSges(6,3) * t57;
t9 = mrSges(6,2) * t26 + mrSges(6,3) * t58;
t114 = -t50 * t8 + t52 * t9;
t71 = t26 * (-t50 ^ 2 - t52 ^ 2);
t35 = -mrSges(6,1) * t52 + mrSges(6,2) * t50;
t42 = -pkin(3) * t85 - pkin(4);
t113 = m(6) * t42 - mrSges(5,1) + t35;
t112 = -2 * mrSges(5,3);
t111 = 2 * qJD(2);
t110 = m(5) * pkin(3);
t104 = Ifges(6,4) * t50;
t103 = Ifges(6,4) * t52;
t102 = Ifges(6,6) * t50;
t23 = -qJD(3) * t56 - t51 * qJD(4);
t55 = -qJD(4) * t106 + t84 * t86;
t12 = t23 * t49 - t55 * t85;
t19 = -t34 * t49 + t56 * t85;
t101 = t12 * t19;
t100 = t19 * t27;
t97 = t26 * t50;
t96 = t26 * t52;
t94 = t29 * t50;
t93 = t29 * t52;
t92 = t50 * mrSges(6,3);
t36 = Ifges(6,2) * t52 + t104;
t90 = t50 * t36;
t89 = t52 * mrSges(6,3);
t37 = Ifges(6,1) * t50 + t103;
t87 = t52 * t37;
t38 = pkin(3) * t68 + qJD(2);
t81 = 0.2e1 * t26 * mrSges(5,3);
t80 = t29 * t112;
t69 = -t26 * mrSges(5,1) - t27 * mrSges(5,2);
t13 = t23 * t85 + t49 * t55;
t14 = -pkin(4) * t26 + pkin(7) * t27 + t38;
t1 = t6 * qJD(5) + t13 * t52 + t14 * t50;
t2 = -t7 * qJD(5) - t13 * t50 + t14 * t52;
t65 = t1 * t52 - t2 * t50;
t63 = mrSges(6,1) * t50 + mrSges(6,2) * t52;
t62 = Ifges(6,1) * t52 - t104;
t61 = -Ifges(6,2) * t50 + t103;
t60 = t12 * t29 + t100;
t17 = -mrSges(6,2) * t30 + t29 * t92;
t18 = mrSges(6,1) * t30 + t29 * t89;
t59 = t52 * t17 - t50 * t18;
t54 = t57 * Ifges(6,5) + Ifges(6,6) * t58 - Ifges(6,3) * t26;
t44 = Ifges(6,5) * t82;
t41 = pkin(3) * t49 + pkin(7);
t33 = t62 * qJD(5);
t32 = t61 * qJD(5);
t31 = t63 * qJD(5);
t15 = t63 * t29;
t11 = t30 * Ifges(6,5) - t29 * t62;
t10 = t30 * Ifges(6,6) - t29 * t61;
t5 = -mrSges(6,1) * t58 + mrSges(6,2) * t57;
t4 = Ifges(6,1) * t57 + Ifges(6,4) * t58 - t26 * Ifges(6,5);
t3 = Ifges(6,4) * t57 + Ifges(6,2) * t58 - t26 * Ifges(6,6);
t21 = [(t80 - 0.2e1 * t15) * t12 + t3 * t94 + t20 * t81 + 0.2e1 * t38 * (t30 * mrSges(5,1) - t29 * mrSges(5,2)) + t57 * t11 + t58 * t10 + (-0.2e1 * Ifges(4,4) * t106 + t118 * t51) * t68 + (0.2e1 * Ifges(4,4) * t51 + t118 * t106) * t84 + 0.2e1 * Ifges(5,1) * t95 + 0.2e1 * t1 * t17 + 0.2e1 * t2 * t18 + 0.2e1 * t19 * t5 + 0.2e1 * t6 * t8 + 0.2e1 * t7 * t9 + (t13 * t30 + t100) * t112 + (t51 * mrSges(4,1) + mrSges(4,2) * t106 + mrSges(3,3)) * t111 + (-(-Ifges(6,5) * t52 + t102) * t29 + (-(2 * Ifges(5,2)) - Ifges(6,3)) * t30) * t26 + 0.2e1 * m(5) * (t13 * t20 + t38 * t43 + t101) + 0.2e1 * m(6) * (t1 * t7 + t2 * t6 + t101) - t4 * t93 + 0.2e1 * t43 * t69 + t30 * t54 + 0.2e1 * (-t26 * t29 + t27 * t30) * Ifges(5,4) + (0.2e1 * (mrSges(4,1) * t106 - mrSges(4,2) * t51) * qJD(3) + ((m(4) + m(3)) * t111)) * qJ(2); t29 * t5 + (-t15 + t80) * t27 - t59 * t26 + m(6) * (t6 * t97 - t7 * t96 + t60) + m(5) * (-t20 * t26 + t60) + (t81 + (-t50 * t17 - t52 * t18) * qJD(5) + m(6) * (t65 + t116) + m(5) * t13 + t114) * t30; 0.2e1 * m(5) * (-t30 * t26 + t95) + 0.2e1 * m(6) * (t30 * t71 + t95); t1 * t89 + t50 * t4 / 0.2e1 - t26 * (Ifges(6,5) * t50 + Ifges(6,6) * t52) / 0.2e1 + t52 * t3 / 0.2e1 + t42 * t5 + t19 * t31 + Ifges(5,6) * t26 - t33 * t93 / 0.2e1 - t2 * t92 - Ifges(4,5) * t84 - t10 * t83 / 0.2e1 + t30 * (-Ifges(6,6) * t83 + t44) / 0.2e1 - Ifges(4,6) * t68 + (qJD(5) * t37 + t32) * t94 / 0.2e1 + (t29 * t36 + t11) * t82 / 0.2e1 + t117 * t53 + (t49 * t110 - mrSges(5,2)) * t13 + t115 * mrSges(5,3) * pkin(3) + t116 * mrSges(6,3) + (-Ifges(5,5) + t90 / 0.2e1 - t87 / 0.2e1) * t27 + (-t85 * t110 + t113) * t12 + (m(6) * ((-t50 * t7 - t52 * t6) * qJD(5) + t65) - t18 * t82 - t17 * t83 + t114) * t41; t26 * mrSges(5,2) - t115 * t110 + t113 * t27 + t29 * t31 + t117 + (m(6) * t41 + mrSges(6,3)) * t71; 0.2e1 * t31 * t42 + t32 * t52 + t33 * t50 + (t87 - t90) * qJD(5); t50 * t9 + t52 * t8 + t59 * qJD(5) + m(6) * (t1 * t50 + t2 * t52 + (-t50 * t6 + t52 * t7) * qJD(5)) + m(5) * t38 + t69; 0; 0; 0; mrSges(6,1) * t2 - mrSges(6,2) * t1 + t54; (t30 * t83 + t96) * mrSges(6,2) + (-t30 * t82 + t97) * mrSges(6,1); t44 + (t35 * t41 - t102) * qJD(5); -t31; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t21(1), t21(2), t21(4), t21(7), t21(11); t21(2), t21(3), t21(5), t21(8), t21(12); t21(4), t21(5), t21(6), t21(9), t21(13); t21(7), t21(8), t21(9), t21(10), t21(14); t21(11), t21(12), t21(13), t21(14), t21(15);];
Mq = res;
