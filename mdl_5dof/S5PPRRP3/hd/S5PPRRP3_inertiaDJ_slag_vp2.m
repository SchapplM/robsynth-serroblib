% Calculate time derivative of joint inertia matrix for
% S5PPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
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
% Datum: 2019-12-05 15:11
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRRP3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP3_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP3_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP3_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP3_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRP3_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRP3_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:10:34
% EndTime: 2019-12-05 15:10:35
% DurationCPUTime: 0.36s
% Computational Cost: add. (212->83), mult. (726->134), div. (0->0), fcn. (514->6), ass. (0->49)
t64 = m(5) / 0.2e1 + m(6) / 0.2e1;
t68 = 0.2e1 * t64;
t28 = cos(qJ(4));
t22 = t28 ^ 2;
t26 = sin(qJ(4));
t54 = t26 ^ 2 + t22;
t67 = m(5) + m(6);
t66 = mrSges(6,2) + mrSges(5,3);
t16 = -t28 * mrSges(6,1) - t26 * mrSges(6,3);
t65 = -t28 * mrSges(5,1) + t26 * mrSges(5,2) + t16;
t63 = m(6) * qJD(5);
t62 = m(6) * qJ(5) + mrSges(6,3);
t25 = cos(pkin(8));
t24 = sin(pkin(8));
t27 = sin(qJ(3));
t52 = qJD(3) * t27;
t45 = t24 * t52;
t29 = cos(qJ(3));
t57 = t24 * t29;
t48 = t28 * t57;
t4 = -qJD(4) * t48 + (qJD(4) * t25 + t45) * t26;
t59 = t26 * t4;
t9 = t25 * t28 + t26 * t57;
t5 = -qJD(4) * t9 - t28 * t45;
t58 = t5 * t28;
t56 = Ifges(5,4) - Ifges(6,5);
t51 = qJD(3) * t29;
t55 = t54 * pkin(6) * t51;
t53 = qJD(3) * t24;
t50 = qJD(4) * t26;
t49 = qJD(4) * t28;
t47 = -mrSges(4,1) + t65;
t44 = t24 * t51;
t43 = t27 * t49;
t42 = t27 * t51;
t37 = t26 * mrSges(6,1) - t28 * mrSges(6,3);
t11 = t37 * qJD(4);
t35 = pkin(4) * t26 - qJ(5) * t28;
t8 = t35 * qJD(4) - t26 * qJD(5);
t41 = m(6) * t8 + t11;
t39 = (m(6) * pkin(6) + mrSges(6,2)) * t28;
t38 = t26 * mrSges(5,1) + t28 * mrSges(5,2);
t36 = -t28 * pkin(4) - t26 * qJ(5);
t30 = m(6) * t36 + t65;
t13 = -pkin(3) + t36;
t12 = t38 * qJD(4);
t10 = -t25 * t26 + t48;
t3 = pkin(6) * t58;
t1 = [0.2e1 * t67 * (t24 ^ 2 * t42 + t10 * t5 - t9 * t4); t26 * (t9 * t51 + (-qJD(4) * t10 - t4) * t27) * t68 + t67 * (t28 * t10 * t51 - t29 ^ 2 * t53 + t9 * t43 + (t53 * t27 + t58) * t27); 0.4e1 * t64 * (-0.1e1 + t54) * t42; ((t11 + t12) * t27 + (mrSges(4,2) * t27 + t47 * t29) * qJD(3)) * t24 + m(5) * (-pkin(3) * t44 + t3) + m(6) * (t8 * t27 * t24 + t13 * t44 + t3) + pkin(6) * (-t10 * t50 + t9 * t49 - t59) * t68 + t66 * (-t59 + t58 + (-t10 * t26 + t28 * t9) * qJD(4)); t47 * t52 + m(6) * (t13 * t52 + t55) + m(5) * (-pkin(3) * t52 + t55) + (-t12 + (t66 * t54 - mrSges(4,2)) * qJD(3) - t41) * t29; -0.2e1 * pkin(3) * t12 + 0.2e1 * t8 * t16 + 0.2e1 * t41 * t13 + 0.2e1 * t56 * qJD(4) * t22 + 0.2e1 * (-t56 * t26 + (Ifges(5,1) + Ifges(6,1) - Ifges(5,2) - Ifges(6,3)) * t28) * t50; t10 * t63 + (-mrSges(5,2) + t62) * t5 + (m(6) * pkin(4) + mrSges(5,1) + mrSges(6,1)) * t4; (-m(6) * t35 - t37 - t38) * t51 + (t30 * qJD(4) + t28 * t63) * t27; qJD(5) * t39 + ((-pkin(4) * mrSges(6,2) + Ifges(6,4) + Ifges(5,5)) * t28 + (-qJ(5) * mrSges(6,2) - Ifges(5,6) + Ifges(6,6)) * t26 + t30 * pkin(6)) * qJD(4); 0.2e1 * t62 * qJD(5); -m(6) * t4; m(6) * (t26 * t51 + t43); qJD(4) * t39; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
