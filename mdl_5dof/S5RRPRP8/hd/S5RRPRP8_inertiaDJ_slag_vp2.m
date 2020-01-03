% Calculate time derivative of joint inertia matrix for
% S5RRPRP8
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
% Datum: 2019-12-31 20:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP8_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP8_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_inertiaDJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP8_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP8_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP8_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:03:21
% EndTime: 2019-12-31 20:03:24
% DurationCPUTime: 0.68s
% Computational Cost: add. (777->144), mult. (1709->204), div. (0->0), fcn. (1268->4), ass. (0->63)
t87 = 2 * mrSges(5,3);
t86 = 2 * mrSges(6,3);
t46 = sin(qJ(4));
t48 = cos(qJ(4));
t50 = -pkin(2) - pkin(3);
t35 = -qJ(3) * t46 + t48 * t50;
t18 = t48 * qJD(3) + t35 * qJD(4);
t72 = mrSges(5,2) + mrSges(6,2);
t85 = t72 * t18;
t84 = t72 * t48;
t83 = m(4) * pkin(6) + mrSges(4,2);
t47 = sin(qJ(2));
t74 = pkin(6) - pkin(7);
t38 = t74 * t47;
t49 = cos(qJ(2));
t39 = t74 * t49;
t14 = t46 * t38 + t48 * t39;
t82 = -t49 * pkin(2) - t47 * qJ(3);
t81 = qJD(2) - qJD(4);
t55 = t46 * t49 - t47 * t48;
t10 = t81 * t55;
t66 = qJD(2) * t47;
t33 = t74 * t66;
t34 = qJD(2) * t39;
t62 = qJD(4) * t48;
t63 = qJD(4) * t46;
t3 = -t48 * t33 + t46 * t34 + t38 * t62 - t39 * t63;
t54 = t46 * t47 + t48 * t49;
t1 = -qJ(5) * t10 - qJD(5) * t54 + t3;
t4 = -t14 * qJD(4) + t33 * t46 + t48 * t34;
t80 = t4 * mrSges(5,1) - t3 * mrSges(5,2) - t1 * mrSges(6,2) + (-Ifges(6,6) - Ifges(5,6)) * t10;
t79 = 2 * m(5);
t78 = 2 * m(6);
t77 = -2 * pkin(1);
t37 = -pkin(1) + t82;
t76 = 0.2e1 * t37;
t73 = mrSges(5,1) + mrSges(6,1);
t71 = -Ifges(4,5) + Ifges(3,4);
t70 = -Ifges(6,5) - Ifges(5,5);
t65 = qJD(2) * t49;
t68 = qJ(3) * t65 + t47 * qJD(3);
t61 = 0.2e1 * t49;
t36 = t48 * qJ(3) + t46 * t50;
t19 = -t46 * qJD(3) - t36 * qJD(4);
t60 = t46 * t18 + t48 * t19 + t36 * t62;
t59 = (m(6) * pkin(4)) + mrSges(6,1);
t13 = t48 * t38 - t39 * t46;
t25 = t49 * pkin(3) - t37;
t58 = mrSges(5,1) + t59;
t57 = 2 * Ifges(5,4) + 2 * Ifges(6,4);
t56 = -t49 * mrSges(4,1) - t47 * mrSges(4,3);
t15 = t50 * t66 + t68;
t52 = -t10 * t36 - t18 * t54 + t19 * t55;
t32 = -pkin(4) + t35;
t12 = pkin(4) * t54 + t25;
t11 = t81 * t54;
t9 = t36 * t18;
t8 = t11 * mrSges(6,2);
t7 = -qJ(5) * t54 + t14;
t6 = qJ(5) * t55 + t13;
t5 = pkin(4) * t10 + t15;
t2 = -qJ(5) * t11 + qJD(5) * t55 + t4;
t16 = [0.2e1 * t15 * (mrSges(5,1) * t54 - mrSges(5,2) * t55) + 0.2e1 * t5 * (mrSges(6,1) * t54 - mrSges(6,2) * t55) + 0.2e1 * t12 * t8 + (t1 * t7 + t12 * t5 + t2 * t6) * t78 + (t13 * t4 + t14 * t3 + t25 * t15) * t79 + (m(4) * t76 + 0.2e1 * t56) * (pkin(2) * t66 - t68) + (-t1 * t54 + t2 * t55) * t86 + (-t3 * t54 + t4 * t55) * t87 + (0.2e1 * mrSges(5,2) * t25 - 0.2e1 * mrSges(5,3) * t13 - 0.2e1 * mrSges(6,3) * t6 - t54 * t57 - 0.2e1 * (Ifges(5,1) + Ifges(6,1)) * t55) * t11 - (-0.2e1 * t25 * mrSges(5,1) - 0.2e1 * t12 * mrSges(6,1) + t14 * t87 - t55 * t57 + t7 * t86 - 0.2e1 * (Ifges(5,2) + Ifges(6,2)) * t54) * t10 + (((mrSges(3,2) * t77) - 0.2e1 * t37 * mrSges(4,3) + t71 * t61) * t49 + ((mrSges(3,1) * t77) + mrSges(4,1) * t76 + (Ifges(3,1) + Ifges(4,1) - Ifges(3,2) - Ifges(4,3)) * t61 - 0.2e1 * t71 * t47) * t47) * qJD(2); -t2 * mrSges(6,1) + t70 * t11 + m(6) * (t1 * t36 + t18 * t7 + t19 * t6 + t2 * t32) + m(5) * (t13 * t19 + t14 * t18 + t3 * t36 + t35 * t4) + (-t11 * t32 + t52) * mrSges(6,3) + (-t11 * t35 + t52) * mrSges(5,3) + ((-pkin(2) * mrSges(4,2) + Ifges(4,4) + Ifges(3,5)) * t49 + (-qJ(3) * mrSges(4,2) - Ifges(3,6) + Ifges(4,6)) * t47 + (m(4) * t82 - t49 * mrSges(3,1) + t47 * mrSges(3,2) + t56) * pkin(6)) * qJD(2) + t83 * qJD(3) * t49 - t80; (t19 * t35 + t9) * t79 + (t19 * t32 + t9) * t78 - 0.2e1 * t73 * t19 + 0.2e1 * t85 + 0.2e1 * (m(4) * qJ(3) + mrSges(4,3)) * qJD(3); t83 * t65 + m(6) * (t1 * t46 + t2 * t48 + (-t46 * t6 + t48 * t7) * qJD(4)) + m(5) * (t3 * t46 + t4 * t48 + (-t13 * t46 + t14 * t48) * qJD(4)) + (mrSges(6,3) + mrSges(5,3)) * (-t46 * t10 - t48 * t11 + (-t46 * t55 - t48 * t54) * qJD(4)); (t73 * t46 + t84) * qJD(4) + m(5) * (-t35 * t63 + t60) + m(6) * (-t32 * t63 + t60); 0; t59 * t2 + (-(mrSges(6,3) * pkin(4)) - t70) * t11 + t80; t58 * t19 - t85; (-t58 * t46 - t84) * qJD(4); 0; m(6) * t5 + t10 * mrSges(6,1) + t8; 0; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t16(1), t16(2), t16(4), t16(7), t16(11); t16(2), t16(3), t16(5), t16(8), t16(12); t16(4), t16(5), t16(6), t16(9), t16(13); t16(7), t16(8), t16(9), t16(10), t16(14); t16(11), t16(12), t16(13), t16(14), t16(15);];
Mq = res;
