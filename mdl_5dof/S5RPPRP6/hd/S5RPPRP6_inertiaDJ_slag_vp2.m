% Calculate time derivative of joint inertia matrix for
% S5RPPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
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
% Datum: 2019-12-31 17:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRP6_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP6_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP6_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP6_inertiaDJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP6_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP6_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP6_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:54:58
% EndTime: 2019-12-31 17:54:59
% DurationCPUTime: 0.35s
% Computational Cost: add. (409->76), mult. (866->109), div. (0->0), fcn. (684->4), ass. (0->37)
t29 = sin(pkin(7));
t30 = cos(pkin(7));
t49 = sin(qJ(4));
t37 = qJD(4) * t49;
t50 = cos(qJ(4));
t38 = qJD(4) * t50;
t13 = -t29 * t38 - t30 * t37;
t14 = -t29 * t37 + t30 * t38;
t16 = t49 * t29 - t50 * t30;
t17 = t50 * t29 + t49 * t30;
t31 = -pkin(1) - qJ(3);
t51 = -pkin(6) + t31;
t18 = t51 * t29;
t19 = t51 * t30;
t8 = t49 * t18 - t50 * t19;
t4 = -t17 * qJD(3) - t8 * qJD(4);
t9 = t50 * t18 + t49 * t19;
t5 = -qJD(3) * t16 + t9 * qJD(4);
t61 = t13 * t8 - t14 * t9 - t16 * t5 - t17 * t4;
t60 = m(5) + m(6);
t58 = mrSges(5,1) + mrSges(6,1);
t57 = mrSges(6,2) + mrSges(5,3);
t56 = -t16 * t13 + t17 * t14;
t36 = (t29 ^ 2 + t30 ^ 2) * qJD(3);
t55 = m(6) * qJ(5) + mrSges(6,3);
t24 = t29 * pkin(3) + qJ(2);
t53 = 0.2e1 * t24;
t46 = Ifges(5,4) - Ifges(6,5);
t44 = 0.2e1 * t16;
t43 = 0.2e1 * t17;
t42 = t9 * t4 + t8 * t5;
t32 = pkin(4) * t13 + qJ(5) * t14 + qJD(5) * t17;
t12 = t14 * mrSges(6,1);
t11 = t13 * mrSges(5,2);
t7 = t17 * pkin(4) + t16 * qJ(5) + t24;
t3 = t14 * pkin(4) - t13 * qJ(5) + t16 * qJD(5) + qJD(2);
t1 = [t11 * t53 + 0.2e1 * t3 * (t17 * mrSges(6,1) + t16 * mrSges(6,3)) + 0.2e1 * t7 * t12 + 0.2e1 * m(5) * (t24 * qJD(2) + t42) + 0.2e1 * m(4) * (qJ(2) * qJD(2) - t31 * t36) + 0.2e1 * m(6) * (t7 * t3 + t42) + (mrSges(5,1) * t53 + (Ifges(5,2) + Ifges(6,3)) * t43 + t46 * t44) * t14 + (-0.2e1 * t7 * mrSges(6,3) - t46 * t43 + (-Ifges(5,1) - Ifges(6,1)) * t44) * t13 + 0.2e1 * (m(3) * qJ(2) + t29 * mrSges(4,1) + t17 * mrSges(5,1) + t30 * mrSges(4,2) - t16 * mrSges(5,2) + mrSges(3,3)) * qJD(2) + 0.2e1 * mrSges(4,3) * t36 + 0.2e1 * t61 * t57; -m(4) * t36 - 0.2e1 * t57 * t56 - t60 * t61; 0.2e1 * t60 * t56; t11 + t14 * mrSges(5,1) + m(6) * t3 - t13 * mrSges(6,3) + t12 + (m(5) + m(4)) * qJD(2); 0; 0; m(6) * qJD(5) * t9 + (Ifges(6,6) - Ifges(5,6)) * t14 + (Ifges(6,4) + Ifges(5,5)) * t13 - t32 * mrSges(6,2) + (-m(6) * pkin(4) - t58) * t5 + (-mrSges(5,2) + t55) * t4; m(6) * t32 + (-mrSges(5,2) + mrSges(6,3)) * t14 + t58 * t13; 0; 0.2e1 * t55 * qJD(5); m(6) * t5 + t13 * mrSges(6,2); -m(6) * t13; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
