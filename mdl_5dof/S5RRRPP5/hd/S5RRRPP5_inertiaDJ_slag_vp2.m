% Calculate time derivative of joint inertia matrix for
% S5RRRPP5
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
% Datum: 2019-12-31 20:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPP5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP5_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP5_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP5_inertiaDJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP5_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP5_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP5_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:57:36
% EndTime: 2019-12-31 20:57:38
% DurationCPUTime: 0.80s
% Computational Cost: add. (889->144), mult. (2082->203), div. (0->0), fcn. (1583->4), ass. (0->63)
t47 = sin(qJ(3));
t48 = sin(qJ(2));
t49 = cos(qJ(3));
t50 = cos(qJ(2));
t29 = t47 * t48 - t49 * t50;
t51 = -pkin(3) - pkin(4);
t30 = t47 * t50 + t48 * t49;
t44 = -pkin(2) * t50 - pkin(1);
t55 = qJ(4) * t30 - t44;
t11 = t29 * t51 + t55;
t86 = -0.2e1 * t11;
t72 = (mrSges(5,2) + mrSges(4,3));
t85 = 2 * t72;
t84 = -mrSges(4,1) - mrSges(5,1);
t83 = mrSges(6,2) + mrSges(5,3);
t82 = ((-mrSges(6,1) + t84) * t47 - mrSges(4,2) * t49) * pkin(2) * qJD(3);
t81 = 2 * mrSges(6,3);
t75 = -pkin(7) - pkin(6);
t36 = t75 * t48;
t37 = t75 * t50;
t23 = t47 * t36 - t49 * t37;
t80 = qJD(2) + qJD(3);
t79 = 2 * m(5);
t78 = 2 * m(6);
t77 = 0.2e1 * t44;
t76 = m(5) + m(6);
t22 = -t49 * t36 - t37 * t47;
t73 = t22 * t47;
t71 = mrSges(5,2) - mrSges(6,3);
t66 = qJD(3) * t49;
t39 = pkin(2) * t66 + qJD(4);
t41 = pkin(2) * t47 + qJ(4);
t70 = qJ(4) * t39 + qJD(4) * t41;
t69 = t83 * qJD(4);
t67 = qJD(3) * t47;
t65 = 0.2e1 * t50;
t64 = pkin(2) * t67;
t43 = -pkin(2) * t49 - pkin(3);
t62 = qJD(2) * t75;
t34 = t48 * t62;
t59 = t50 * t62;
t10 = t23 * qJD(3) + t34 * t47 - t49 * t59;
t9 = t49 * t34 + t36 * t66 + t37 * t67 + t47 * t59;
t63 = t10 * t22 + t23 * t9;
t61 = t83 * t39;
t21 = t80 * t30;
t58 = t21 * t41 + t29 * t39;
t57 = -2 * Ifges(4,4) + 2 * Ifges(6,4) + 2 * Ifges(5,5);
t56 = -qJ(4) * t21 - qJD(4) * t29;
t20 = t80 * t29;
t53 = -qJD(2) * t48 * pkin(2) - qJ(4) * t20 + qJD(4) * t30;
t3 = qJ(5) * t21 + qJD(5) * t29 + t9;
t4 = qJ(5) * t20 - qJD(5) * t30 + t10;
t52 = -t4 * mrSges(6,1) + t3 * mrSges(6,2) + (-mrSges(4,2) + mrSges(5,3)) * t9 + t84 * t10 + (-Ifges(4,6) + Ifges(5,6) - Ifges(6,6)) * t21 - (Ifges(5,4) + Ifges(4,5) - Ifges(6,5)) * t20;
t40 = -pkin(4) + t43;
t27 = t41 * t39;
t15 = t20 * mrSges(6,2);
t14 = pkin(3) * t29 - t55;
t13 = qJ(5) * t29 + t23;
t12 = -qJ(5) * t30 + t22;
t6 = pkin(3) * t21 - t53;
t1 = t21 * t51 + t53;
t2 = [0.2e1 * t6 * (t29 * mrSges(5,1) - t30 * mrSges(5,3)) + 0.2e1 * t1 * (-t29 * mrSges(6,1) + t30 * mrSges(6,2)) + t15 * t86 + (t3 * t29 - t4 * t30) * t81 + 0.2e1 * m(4) * t63 + (t1 * t11 + t12 * t4 + t13 * t3) * t78 + (t14 * t6 + t63) * t79 + ((-pkin(1) * mrSges(3,2) + Ifges(3,4) * t50) * t65 + (m(4) * pkin(2) * t77 + 0.2e1 * pkin(2) * (mrSges(4,1) * t29 + mrSges(4,2) * t30) - 0.2e1 * pkin(1) * mrSges(3,1) - 0.2e1 * Ifges(3,4) * t48 + (Ifges(3,1) - Ifges(3,2)) * t65) * t48) * qJD(2) + (mrSges(4,1) * t77 + 0.2e1 * t14 * mrSges(5,1) + mrSges(6,1) * t86 + t13 * t81 + t30 * t57 - 0.2e1 * t72 * t23 + 0.2e1 * (Ifges(4,2) + Ifges(6,2) + Ifges(5,3)) * t29) * t21 - (mrSges(4,2) * t77 - 0.2e1 * mrSges(5,3) * t14 - 0.2e1 * mrSges(6,3) * t12 + t29 * t57 + t22 * t85 + 0.2e1 * (Ifges(4,1) + Ifges(5,1) + Ifges(6,1)) * t30) * t20 + (t10 * t30 - t9 * t29) * t85; (Ifges(3,5) * t50 - Ifges(3,6) * t48 + (-mrSges(3,1) * t50 + mrSges(3,2) * t48) * pkin(6)) * qJD(2) + m(6) * (t13 * t39 + t3 * t41 + t4 * t40) + m(5) * (t10 * t43 + t23 * t39 + t41 * t9) + (t20 * t40 + t58) * mrSges(6,3) + (-t20 * t43 - t58) * mrSges(5,2) + t52 + (m(4) * (-t10 * t49 + t47 * t9) + (t20 * t49 - t21 * t47) * mrSges(4,3) + (-t49 * t29 * mrSges(4,3) + m(5) * t73 + m(4) * (t23 * t49 + t73) + ((mrSges(4,3) + t71) * t30 + m(6) * t12) * t47) * qJD(3)) * pkin(2); (t43 * t64 + t27) * t79 + (t40 * t64 + t27) * t78 + 0.2e1 * t61 + 0.2e1 * t82; t52 + m(5) * (-pkin(3) * t10 + qJ(4) * t9 + qJD(4) * t23) + m(6) * (qJ(4) * t3 + qJD(4) * t13 + t4 * t51) + (t20 * t51 - t56) * mrSges(6,3) + (pkin(3) * t20 + t56) * mrSges(5,2); t61 + t82 + m(5) * (-pkin(3) * t64 + t70) + m(6) * (t51 * t64 + t70) + t69; 0.2e1 * qJ(4) * qJD(4) * t76 + 0.2e1 * t69; m(5) * t10 + m(6) * t4 - t20 * t71; t76 * t64; 0; 0; m(6) * t1 - t21 * mrSges(6,1) - t15; 0; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
