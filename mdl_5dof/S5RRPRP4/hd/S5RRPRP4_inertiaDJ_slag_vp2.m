% Calculate time derivative of joint inertia matrix for
% S5RRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% Datum: 2019-12-31 19:53
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP4_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP4_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP4_inertiaDJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP4_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP4_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP4_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:52:37
% EndTime: 2019-12-31 19:52:39
% DurationCPUTime: 0.42s
% Computational Cost: add. (262->94), mult. (623->121), div. (0->0), fcn. (282->4), ass. (0->50)
t79 = -Ifges(5,1) - Ifges(6,1);
t32 = sin(qJ(4));
t78 = (Ifges(5,4) - Ifges(6,5)) * t32;
t77 = 2 * qJD(3);
t34 = cos(qJ(4));
t56 = t32 ^ 2 + t34 ^ 2;
t33 = sin(qJ(2));
t55 = pkin(1) * qJD(2);
t49 = t33 * t55;
t76 = t49 * t56;
t75 = -(mrSges(5,3) + mrSges(6,2)) * t56 + mrSges(4,2) - mrSges(3,1);
t74 = 2 * m(4);
t18 = t32 * mrSges(5,1) + t34 * mrSges(5,2);
t72 = mrSges(4,3) + t18;
t70 = -t32 * pkin(4) + qJ(5) * t34;
t53 = qJD(4) * t32;
t43 = t34 * mrSges(5,1) - t32 * mrSges(5,2);
t12 = t43 * qJD(4);
t69 = 0.2e1 * t12;
t17 = t32 * mrSges(6,1) - t34 * mrSges(6,3);
t68 = 0.2e1 * t17;
t36 = -pkin(2) - pkin(7);
t35 = cos(qJ(2));
t65 = t35 * pkin(1);
t23 = -t65 + t36;
t67 = t23 * t76;
t66 = pkin(1) * t33;
t63 = Ifges(5,4) * t34;
t48 = t35 * t55;
t22 = qJD(3) + t48;
t25 = qJ(3) + t66;
t60 = t25 * t22;
t59 = t35 * mrSges(3,2);
t58 = t36 * t76;
t47 = mrSges(6,2) * t53;
t52 = qJD(4) * t34;
t57 = Ifges(6,6) * t52 + pkin(4) * t47;
t51 = qJD(5) * t32;
t14 = qJ(3) - t70;
t42 = t34 * mrSges(6,1) + t32 * mrSges(6,3);
t41 = qJ(3) * t22 + qJD(3) * t25;
t2 = pkin(4) * t52 + qJ(5) * t53 - qJD(5) * t34 + qJD(3);
t40 = (-mrSges(6,2) * qJ(5) - Ifges(5,6)) * t34 + (-Ifges(6,4) - Ifges(5,5)) * t32;
t39 = (t78 * t32 + (-t63 + (Ifges(5,2) + t79) * t32) * t34) * qJD(4) + ((Ifges(5,2) + Ifges(6,3)) * t34 + t78) * t53 + (0.2e1 * Ifges(6,5) * t34 - t63 + (Ifges(6,3) + t79) * t32) * t52;
t38 = m(6) * t70 - t17 - t18;
t37 = m(6) * t51 + t38 * qJD(4);
t11 = t42 * qJD(4);
t7 = t14 + t66;
t1 = t2 + t48;
t3 = [t1 * t68 + 0.2e1 * t7 * t11 + t25 * t69 + 0.2e1 * t72 * t22 + 0.2e1 * m(6) * (t1 * t7 + t67) + 0.2e1 * m(5) * (t60 + t67) + t60 * t74 + 0.2e1 * (-pkin(1) * t59 + ((-pkin(2) / 0.2e1 - t65 / 0.2e1) * t74 + t75) * t66) * qJD(2) + t39; (t2 + t1) * t17 + (qJ(3) + t25) * t12 + (t14 + t7) * t11 + m(4) * t41 + m(5) * (t41 + t58) + m(6) * (t1 * t14 + t2 * t7 + t58) + (-t59 + (-m(4) * pkin(2) + t75) * t33) * t55 + t39 + t72 * (t22 + qJD(3)); t2 * t68 + t72 * t77 + t39 + 0.2e1 * (m(6) * t2 + t11) * t14 + (t69 + (m(4) + m(5)) * t77) * qJ(3); (m(4) + 0.2e1 * (m(6) / 0.2e1 + m(5) / 0.2e1) * t56) * t49; 0; 0; (m(6) * t23 - mrSges(6,2)) * t51 + (m(6) * (pkin(4) * t34 + qJ(5) * t32) + t42 + t43) * t49 + (t38 * t23 + t40) * qJD(4) + t57; -mrSges(6,2) * t51 + t40 * qJD(4) + t37 * t36 + t57; t37; 0.2e1 * (m(6) * qJ(5) + mrSges(6,3)) * qJD(5); -t47 + (t23 * t53 - t34 * t49) * m(6); (m(6) * t36 - mrSges(6,2)) * t53; m(6) * t53; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t3(1), t3(2), t3(4), t3(7), t3(11); t3(2), t3(3), t3(5), t3(8), t3(12); t3(4), t3(5), t3(6), t3(9), t3(13); t3(7), t3(8), t3(9), t3(10), t3(14); t3(11), t3(12), t3(13), t3(14), t3(15);];
Mq = res;
