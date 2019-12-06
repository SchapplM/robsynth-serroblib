% Calculate time derivative of joint inertia matrix for
% S5RPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-05 18:00
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP1_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_inertiaDJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP1_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP1_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP1_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:59:15
% EndTime: 2019-12-05 17:59:18
% DurationCPUTime: 0.88s
% Computational Cost: add. (840->111), mult. (1702->173), div. (0->0), fcn. (1348->4), ass. (0->56)
t47 = cos(qJ(4));
t48 = cos(qJ(3));
t79 = sin(qJ(4));
t80 = sin(qJ(3));
t34 = -t47 * t80 - t48 * t79;
t93 = qJD(3) + qJD(4);
t20 = t93 * t34;
t44 = pkin(3) * t47 + pkin(4);
t35 = t47 * t48 - t79 * t80;
t58 = qJD(4) * t79;
t55 = t35 * t58;
t19 = t35 * t93;
t63 = pkin(3) * t79;
t66 = qJD(4) * t47;
t65 = pkin(3) * t66;
t70 = t19 * t63 - t34 * t65;
t92 = pkin(3) * t55 - t70;
t108 = t44 * t20 - t92;
t107 = mrSges(5,1) + mrSges(6,1);
t49 = -pkin(1) - pkin(6);
t103 = -pkin(7) + t49;
t36 = t103 * t80;
t37 = t103 * t48;
t21 = -t36 * t79 + t47 * t37;
t22 = t47 * t36 + t79 * t37;
t32 = t36 * qJD(3);
t33 = qJD(3) * t37;
t6 = -t79 * t32 + t47 * t33 - t36 * t58 + t37 * t66;
t7 = -t22 * qJD(4) - t47 * t32 - t79 * t33;
t106 = -t19 * t22 - t20 * t21 + t6 * t34 - t7 * t35;
t2 = -t19 * qJ(5) + t34 * qJD(5) + t6;
t3 = -t20 * qJ(5) - t35 * qJD(5) + t7;
t9 = qJ(5) * t34 + t22;
t105 = -t19 * t9 + t2 * t34 - t3 * t35;
t71 = -mrSges(5,2) - mrSges(6,2);
t104 = t107 * t20 + t19 * t71;
t102 = t19 * t34;
t100 = t20 * t35;
t98 = -Ifges(4,1) + Ifges(4,2);
t59 = qJD(3) * t80;
t67 = qJD(3) * t48;
t97 = -mrSges(4,1) * t59 - mrSges(4,2) * t67;
t87 = 2 * qJD(2);
t86 = m(6) * pkin(4);
t78 = mrSges(6,3) * t20;
t74 = t20 * t47;
t68 = pkin(3) * qJD(4);
t41 = t80 * pkin(3) + qJ(2);
t38 = pkin(3) * t67 + qJD(2);
t60 = t19 * mrSges(6,1) + t20 * mrSges(6,2);
t53 = -t100 + t102;
t50 = t7 * mrSges(5,1) + t3 * mrSges(6,1) - t6 * mrSges(5,2) - t2 * mrSges(6,2) + (Ifges(5,5) + Ifges(6,5)) * t20 + (-Ifges(5,6) - Ifges(6,6)) * t19;
t23 = -pkin(4) * t34 + t41;
t10 = pkin(4) * t19 + t38;
t8 = -t35 * qJ(5) + t21;
t1 = [0.2e1 * m(6) * (t10 * t23 + t2 * t9 + t3 * t8) + 0.2e1 * m(5) * (t21 * t7 + t22 * t6 + t38 * t41) - 0.2e1 * t8 * t78 + 0.2e1 * t23 * t60 + 0.2e1 * t41 * (t19 * mrSges(5,1) + t20 * mrSges(5,2)) + 0.2e1 * t10 * (-t34 * mrSges(6,1) + t35 * mrSges(6,2)) + 0.2e1 * t38 * (-t34 * mrSges(5,1) + t35 * mrSges(5,2)) + (-0.2e1 * Ifges(4,4) * t48 + t98 * t80) * t67 + (0.2e1 * Ifges(4,4) * t80 + t98 * t48) * t59 + (mrSges(4,1) * t80 + t48 * mrSges(4,2) + mrSges(3,3)) * t87 - 0.2e1 * (Ifges(6,2) + Ifges(5,2)) * t102 + 0.2e1 * (Ifges(6,1) + Ifges(5,1)) * t100 + 0.2e1 * (Ifges(6,4) + Ifges(5,4)) * (-t19 * t35 + t20 * t34) + (0.2e1 * (mrSges(4,1) * t48 - mrSges(4,2) * t80) * qJD(3) + ((m(4) + m(3)) * t87)) * qJ(2) + 0.2e1 * t106 * mrSges(5,3) + 0.2e1 * t105 * mrSges(6,3); m(6) * (t20 * t8 - t105) - m(5) * t106 + 0.2e1 * (mrSges(6,3) + mrSges(5,3)) * t53; -0.4e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * t53; -Ifges(4,6) * t67 + t50 + m(6) * (t44 * t3 + (t79 * t2 + (t47 * t9 - t79 * t8) * qJD(4)) * pkin(3)) - Ifges(4,5) * t59 + m(5) * (t79 * t6 + t47 * t7 + (-t21 * t79 + t22 * t47) * qJD(4)) * pkin(3) + t97 * t49 - t108 * mrSges(6,3) + (-pkin(3) * t74 + t92) * mrSges(5,3); m(5) * ((-t55 + t74) * pkin(3) + t70) + m(6) * t108 + t97 + t104; 0.2e1 * m(6) * (-t44 * t79 + t47 * t63) * t68 + 0.2e1 * t71 * t65 - 0.2e1 * t107 * pkin(3) * t58; (m(6) * t3 - t78) * pkin(4) + t50; t20 * t86 + t104; (t47 * t71 + (-t86 - t107) * t79) * t68; 0; m(6) * t10 + t60; 0; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
