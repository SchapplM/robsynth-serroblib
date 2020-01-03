% Calculate time derivative of joint inertia matrix for
% S5RRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2020-01-03 11:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR2_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR2_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR2_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR2_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:57:21
% EndTime: 2020-01-03 11:57:23
% DurationCPUTime: 0.55s
% Computational Cost: add. (693->140), mult. (1591->205), div. (0->0), fcn. (1158->8), ass. (0->82)
t95 = 2 * m(5);
t94 = 2 * m(6);
t40 = sin(pkin(9));
t36 = t40 ^ 2;
t44 = sin(qJ(5));
t46 = cos(qJ(5));
t70 = qJD(5) * t40;
t18 = (mrSges(6,1) * t46 - mrSges(6,2) * t44) * t70;
t93 = 0.2e1 * t18;
t41 = sin(pkin(8));
t43 = cos(pkin(8));
t47 = cos(qJ(2));
t73 = pkin(1) * qJD(2);
t63 = t47 * t73;
t45 = sin(qJ(2));
t64 = t45 * t73;
t22 = -t41 * t64 + t43 * t63;
t16 = qJD(4) + t22;
t77 = t43 * t45;
t21 = (t41 * t47 + t77) * t73;
t42 = cos(pkin(9));
t50 = -t42 * pkin(4) - t40 * pkin(7) - pkin(3);
t35 = t47 * pkin(1) + pkin(2);
t59 = -t41 * t45 * pkin(1) + t43 * t35;
t11 = t50 - t59;
t74 = pkin(1) * t77 + t41 * t35;
t17 = qJ(4) + t74;
t79 = t42 * t44;
t3 = t46 * t11 - t17 * t79;
t78 = t42 * t46;
t1 = qJD(5) * t3 + t16 * t78 + t44 * t21;
t87 = mrSges(6,3) * t40;
t26 = t42 * mrSges(6,2) - t44 * t87;
t92 = t1 * t26;
t4 = t44 * t11 + t17 * t78;
t2 = -qJD(5) * t4 - t16 * t79 + t46 * t21;
t27 = -t42 * mrSges(6,1) - t46 * t87;
t91 = t2 * t27;
t90 = t43 * pkin(2);
t71 = qJD(4) * t42;
t25 = t50 - t90;
t34 = t41 * pkin(2) + qJ(4);
t8 = t46 * t25 - t34 * t79;
t5 = qJD(5) * t8 + t46 * t71;
t89 = t5 * t26;
t9 = t44 * t25 + t34 * t78;
t6 = -qJD(5) * t9 - t44 * t71;
t88 = t6 * t27;
t86 = Ifges(6,4) * t44;
t85 = Ifges(6,4) * t46;
t53 = mrSges(6,1) * t44 + mrSges(6,2) * t46;
t24 = t53 * t40;
t84 = t16 * t24;
t83 = t22 * mrSges(4,2);
t82 = t36 * t16;
t37 = t42 ^ 2;
t81 = t37 * t16;
t76 = t44 * (-Ifges(6,2) * t46 - t86) * t70;
t67 = t36 * qJD(4);
t75 = t17 * t67 + t34 * t82;
t72 = qJD(4) * t24;
t69 = qJD(5) * t44;
t68 = qJD(5) * t46;
t66 = t37 * qJD(4);
t65 = 0.2e1 * mrSges(6,3);
t62 = t40 * t69;
t61 = t40 * t68;
t19 = (-Ifges(6,5) * t44 - Ifges(6,6) * t46) * t70;
t60 = -t42 * t19 + t36 * (-Ifges(6,1) * t44 - t85) * t68;
t56 = t26 * t68 - t27 * t69;
t55 = t3 * t44 - t4 * t46;
t54 = t44 * t8 - t46 * t9;
t14 = -Ifges(6,6) * t42 + (-Ifges(6,2) * t44 + t85) * t40;
t15 = -Ifges(6,5) * t42 + (Ifges(6,1) * t46 - t86) * t40;
t52 = -t14 * t46 - t15 * t44;
t51 = 0.2e1 * mrSges(5,3) * (t36 + t37);
t49 = -t26 * t69 - t27 * t68;
t48 = -t42 * t18 + (-t44 ^ 2 - t46 ^ 2) * t36 * qJD(5) * mrSges(6,3);
t29 = -t42 * mrSges(5,1) + t40 * mrSges(5,2);
t28 = t34 * t67;
t7 = t17 * t82;
t10 = [-0.2e1 * t83 + 0.2e1 * t92 + 0.2e1 * t91 + t16 * t51 + (t17 * t81 + t7 + (-pkin(3) - t59) * t21) * t95 + (t4 * t1 + t3 * t2 + t7) * t94 + 0.2e1 * m(4) * (-t59 * t21 + t74 * t22) + (0.2e1 * t84 + t17 * t93 - t76 + (t55 * t65 + t52) * qJD(5)) * t40 + t60 + 0.2e1 * (-mrSges(4,1) + t29) * t21 + 0.2e1 * (-mrSges(3,1) * t45 - mrSges(3,2) * t47) * t73; m(4) * (-t43 * t21 + t22 * t41) * pkin(2) + t91 + t92 + t21 * t29 - t21 * mrSges(4,1) + t89 + t88 - t83 + m(5) * ((-pkin(3) - t90) * t21 + (qJD(4) * t17 + t34 * t16) * t37 + t75) + m(6) * (t9 * t1 + t8 * t2 + t6 * t3 + t5 * t4 + t75) + t60 - mrSges(3,2) * t63 - mrSges(3,1) * t64 - t14 * t61 - t15 * t62 + (t81 + t82 + t66 + t67) * mrSges(5,3) + ((t17 + t34) * t18 + t72 + t84 - t76) * t40 + ((t3 + t8) * t62 + (-t4 - t9) * t61) * mrSges(6,3); 0.2e1 * t89 + 0.2e1 * t88 + qJD(4) * t51 + (t9 * t5 + t8 * t6 + t28) * t94 + (t34 * t66 + t28) * t95 + (0.2e1 * t72 + t34 * t93 - t76 + (t54 * t65 + t52) * qJD(5)) * t40 + t60; (m(6) * (t1 * t46 - t16 * t42 - t2 * t44 - t3 * t68 - t4 * t69) + t49) * t40 + t48; (m(6) * (-t44 * t6 + t46 * t5 - t8 * t68 - t9 * t69 - t71) + t49) * t40 + t48; 0; m(5) * t21 + m(6) * (-t55 * qJD(5) + t44 * t1 + t46 * t2) + t56; m(6) * (-t54 * qJD(5) + t44 * t5 + t46 * t6) + t56; 0; 0; t2 * mrSges(6,1) - t1 * mrSges(6,2) + t19; t6 * mrSges(6,1) - t5 * mrSges(6,2) + t19; -t18; -t53 * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t10(1), t10(2), t10(4), t10(7), t10(11); t10(2), t10(3), t10(5), t10(8), t10(12); t10(4), t10(5), t10(6), t10(9), t10(13); t10(7), t10(8), t10(9), t10(10), t10(14); t10(11), t10(12), t10(13), t10(14), t10(15);];
Mq = res;
