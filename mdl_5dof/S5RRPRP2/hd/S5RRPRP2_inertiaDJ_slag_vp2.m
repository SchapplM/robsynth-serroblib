% Calculate time derivative of joint inertia matrix for
% S5RRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 19:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP2_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP2_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP2_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP2_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP2_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP2_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:49:27
% EndTime: 2019-12-31 19:49:28
% DurationCPUTime: 0.37s
% Computational Cost: add. (304->86), mult. (814->118), div. (0->0), fcn. (468->6), ass. (0->52)
t41 = sin(qJ(4));
t43 = cos(qJ(4));
t81 = t41 ^ 2 + t43 ^ 2;
t80 = 2 * Ifges(5,4) - 2 * Ifges(6,5);
t40 = cos(pkin(8));
t44 = cos(qJ(2));
t61 = pkin(1) * qJD(2);
t39 = sin(pkin(8));
t42 = sin(qJ(2));
t65 = t39 * t42;
t10 = (t40 * t44 - t65) * t61;
t79 = t81 * t10;
t78 = Ifges(5,1) + Ifges(6,1) - Ifges(5,2) - Ifges(6,3);
t59 = qJD(4) * t43;
t77 = (-mrSges(4,2) + (mrSges(6,2) + mrSges(5,3)) * t81) * t10;
t51 = t41 * mrSges(6,1) - t43 * mrSges(6,3);
t19 = t51 * qJD(4);
t76 = 0.2e1 * t19;
t52 = t41 * mrSges(5,1) + t43 * mrSges(5,2);
t20 = t52 * qJD(4);
t75 = 0.2e1 * t20;
t30 = pkin(1) * t44 + pkin(2);
t64 = t40 * t42;
t63 = pkin(1) * t64 + t39 * t30;
t8 = pkin(7) + t63;
t74 = t79 * t8;
t28 = pkin(2) * t39 + pkin(7);
t73 = t79 * t28;
t72 = pkin(2) * t40;
t60 = qJD(4) * t41;
t58 = qJD(5) * t43;
t57 = m(6) * t58;
t23 = -t43 * mrSges(5,1) + t41 * mrSges(5,2);
t9 = (t39 * t44 + t64) * t61;
t56 = (t23 - mrSges(4,1)) * t9;
t55 = -pkin(1) * t65 + t30 * t40;
t54 = mrSges(6,2) * t58 + Ifges(6,6) * t60 + (Ifges(6,4) + Ifges(5,5)) * t59;
t22 = -t43 * mrSges(6,1) - t41 * mrSges(6,3);
t50 = -pkin(4) * t43 - qJ(5) * t41;
t49 = -pkin(3) + t50;
t48 = (-mrSges(3,1) * t42 - mrSges(3,2) * t44) * t61;
t11 = -pkin(4) * t60 + qJ(5) * t59 + t41 * qJD(5);
t47 = t50 * mrSges(6,2) - Ifges(5,6) * t41;
t46 = (-t80 * t41 + t78 * t43) * t60 + (t78 * t41 + t80 * t43) * t59;
t45 = m(6) * t50 + t22 + t23;
t32 = mrSges(6,2) * t59;
t29 = -pkin(3) - t72;
t15 = t49 - t72;
t7 = -pkin(3) - t55;
t6 = t49 - t55;
t3 = t9 - t11;
t1 = [t46 + 0.2e1 * m(6) * (t3 * t6 + t74) + 0.2e1 * m(5) * (t7 * t9 + t74) + 0.2e1 * m(4) * (t63 * t10 - t55 * t9) + t6 * t76 + t7 * t75 + 0.2e1 * t3 * t22 + 0.2e1 * t77 + 0.2e1 * t56 + 0.2e1 * t48; t46 + t56 + m(4) * (t10 * t39 - t40 * t9) * pkin(2) + t48 + m(6) * (-t11 * t6 + t15 * t3 + t73) + m(5) * (t29 * t9 + t73) + (-t11 + t3) * t22 + (t29 + t7) * t20 + (t15 + t6) * t19 + t77; t15 * t76 + t29 * t75 + 0.2e1 * (-m(6) * t15 - t22) * t11 + t46; 0; 0; 0; t8 * t57 + (m(6) * (-pkin(4) * t41 + qJ(5) * t43) - t51 - t52) * t10 + (t45 * t8 + t47) * qJD(4) + t54; t47 * qJD(4) + (t45 * qJD(4) + t57) * t28 + t54; m(6) * t11 + ((-mrSges(5,2) + mrSges(6,3)) * t43 + (-mrSges(5,1) - mrSges(6,1)) * t41) * qJD(4); 0.2e1 * (m(6) * qJ(5) + mrSges(6,3)) * qJD(5); t32 + (t10 * t41 + t8 * t59) * m(6); m(6) * t28 * t59 + t32; m(6) * t60; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
