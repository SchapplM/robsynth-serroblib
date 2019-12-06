% Calculate time derivative of joint inertia matrix for
% S5RRPRP1
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
% Datum: 2019-12-05 18:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP1_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP1_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP1_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP1_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:22:00
% EndTime: 2019-12-05 18:22:03
% DurationCPUTime: 0.55s
% Computational Cost: add. (435->122), mult. (1047->165), div. (0->0), fcn. (650->6), ass. (0->69)
t91 = Ifges(6,4) + Ifges(5,4);
t90 = -Ifges(5,2) - Ifges(6,2);
t46 = cos(qJ(4));
t81 = 0.2e1 * t46;
t89 = Ifges(5,1) + Ifges(6,1);
t44 = sin(qJ(4));
t60 = qJD(4) * t46;
t61 = qJD(4) * t44;
t87 = t91 * t44;
t49 = (t89 * t46 - t87) * t61 + ((t89 + t90) * t44 + t91 * t81) * t60;
t68 = t90 * t46 - t87;
t88 = t68 * t61 + t49;
t85 = 2 * m(6);
t21 = mrSges(6,1) * t61 + mrSges(6,2) * t60;
t84 = 0.2e1 * t21;
t76 = mrSges(5,2) * t46;
t50 = mrSges(5,1) * t44 + t76;
t22 = t50 * qJD(4);
t83 = 0.2e1 * t22;
t82 = -0.2e1 * t44;
t80 = m(6) * pkin(4);
t79 = mrSges(6,3) * pkin(4);
t78 = pkin(4) * t46;
t24 = -mrSges(6,1) * t46 + mrSges(6,2) * t44;
t42 = sin(pkin(8));
t47 = cos(qJ(2));
t64 = pkin(1) * qJD(2);
t43 = cos(pkin(8));
t45 = sin(qJ(2));
t70 = t43 * t45;
t11 = (t42 * t47 + t70) * t64;
t58 = pkin(4) * t61;
t5 = t11 + t58;
t77 = t5 * t24;
t71 = t42 * t45;
t69 = -Ifges(5,6) - Ifges(6,6);
t33 = pkin(1) * t47 + pkin(2);
t67 = pkin(1) * t70 + t42 * t33;
t66 = (Ifges(5,5) + Ifges(6,5)) * t60;
t65 = t44 ^ 2 + t46 ^ 2;
t10 = pkin(7) + t67;
t63 = -qJ(5) - t10;
t31 = pkin(2) * t42 + pkin(7);
t62 = -qJ(5) - t31;
t59 = 0.2e1 * qJD(4);
t57 = mrSges(6,1) + t80;
t32 = -pkin(2) * t43 - pkin(3);
t56 = (-mrSges(5,1) * t46 + mrSges(5,2) * t44 - mrSges(4,1)) * t11;
t55 = t65 * mrSges(5,3);
t12 = (t43 * t47 - t71) * t64;
t54 = t65 * t12;
t53 = -pkin(1) * t71 + t33 * t43;
t52 = qJD(4) * t63;
t51 = qJD(4) * t62;
t9 = -pkin(3) - t53;
t48 = (-mrSges(3,1) * t45 - mrSges(3,2) * t47) * t64;
t39 = t46 * qJ(5);
t38 = t46 * qJD(5);
t23 = t32 - t78;
t18 = t31 * t46 + t39;
t17 = t62 * t44;
t8 = -qJD(5) * t44 + t46 * t51;
t7 = t44 * t51 + t38;
t6 = t9 - t78;
t4 = t10 * t46 + t39;
t3 = t63 * t44;
t2 = (-qJD(5) - t12) * t44 + t46 * t52;
t1 = t12 * t46 + t44 * t52 + t38;
t13 = [t6 * t84 + t9 * t83 + 0.2e1 * t77 - 0.2e1 * t12 * mrSges(4,2) + 0.2e1 * m(5) * (t10 * t54 + t11 * t9) + 0.2e1 * m(4) * (-t53 * t11 + t67 * t12) + (t1 * t4 + t2 * t3 + t5 * t6) * t85 + (t1 * t81 + t2 * t82 + (-t3 * t46 - t4 * t44) * t59) * mrSges(6,3) + 0.2e1 * t55 * t12 + 0.2e1 * t48 + 0.2e1 * t56 + t88; t77 + (t9 + t32) * t22 + (t6 + t23) * t21 + t56 + t48 + (-mrSges(4,2) + t55) * t12 + (pkin(4) * t24 + t68) * t61 + m(6) * (t1 * t18 + t17 * t2 + t23 * t5 + t3 * t8 + t4 * t7 + t6 * t58) + m(5) * (t11 * t32 + t31 * t54) + m(4) * (-t11 * t43 + t12 * t42) * pkin(2) + ((t1 + t7) * t46 + (-t2 - t8) * t44 + ((-t17 - t3) * t46 + (-t18 - t4) * t44) * qJD(4)) * mrSges(6,3) + t49; 0.2e1 * t24 * t58 + t23 * t84 + t32 * t83 + (t17 * t8 + t18 * t7 + t23 * t58) * t85 + (t8 * t82 + t7 * t81 + (-t17 * t46 - t18 * t44) * t59) * mrSges(6,3) + t88; m(6) * (t1 * t44 + t2 * t46 + (-t3 * t44 + t4 * t46) * qJD(4)); m(6) * (t44 * t7 + t46 * t8 + (-t17 * t44 + t18 * t46) * qJD(4)); 0; -mrSges(6,2) * t1 + t57 * t2 - t50 * t12 + ((-mrSges(5,1) * t10 - t79) * t46 + (mrSges(5,2) * t10 + t69) * t44) * qJD(4) + t66; -mrSges(6,2) * t7 + t57 * t8 + ((-mrSges(5,1) * t31 - t79) * t46 + (mrSges(5,2) * t31 + t69) * t44) * qJD(4) + t66; (-t76 + (-mrSges(5,1) - t80) * t44) * qJD(4) - t21; 0; m(6) * t5 + t21; m(6) * t58 + t21; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t13(1), t13(2), t13(4), t13(7), t13(11); t13(2), t13(3), t13(5), t13(8), t13(12); t13(4), t13(5), t13(6), t13(9), t13(13); t13(7), t13(8), t13(9), t13(10), t13(14); t13(11), t13(12), t13(13), t13(14), t13(15);];
Mq = res;
