% Calculate time derivative of joint inertia matrix for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-12-05 17:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR4_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR4_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR4_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR4_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:44:06
% EndTime: 2019-12-05 17:44:11
% DurationCPUTime: 0.90s
% Computational Cost: add. (1596->145), mult. (3733->250), div. (0->0), fcn. (3615->8), ass. (0->75)
t66 = sin(pkin(8));
t68 = cos(pkin(8));
t58 = -pkin(2) * t68 - qJ(3) * t66 - pkin(1);
t67 = cos(pkin(9));
t53 = t67 * t58;
t65 = sin(pkin(9));
t86 = t66 * t67;
t29 = -pkin(6) * t86 + t53 + (-qJ(2) * t65 - pkin(3)) * t68;
t82 = qJ(2) * t68;
t83 = t65 * t58 + t67 * t82;
t87 = t65 * t66;
t34 = -pkin(6) * t87 + t83;
t70 = sin(qJ(4));
t72 = cos(qJ(4));
t20 = t72 * t29 - t34 * t70;
t63 = t66 ^ 2;
t21 = t70 * t29 + t72 * t34;
t74 = t70 * t65 - t67 * t72;
t50 = t74 * qJD(4);
t93 = -2 * mrSges(5,1);
t56 = t65 * t72 + t70 * t67;
t46 = t56 * t66;
t47 = t74 * t66;
t69 = sin(qJ(5));
t71 = cos(qJ(5));
t24 = -t46 * t71 + t47 * t69;
t92 = 0.2e1 * t24;
t25 = -t46 * t69 - t47 * t71;
t91 = 0.2e1 * t25;
t90 = -0.2e1 * t46;
t89 = -0.2e1 * t47;
t51 = t56 * qJD(4);
t41 = t66 * t51;
t42 = t66 * t50;
t13 = t24 * qJD(5) - t41 * t71 + t42 * t69;
t14 = -t25 * qJD(5) + t41 * t69 + t42 * t71;
t85 = Ifges(6,5) * t13 + Ifges(6,6) * t14;
t84 = -Ifges(5,5) * t41 + Ifges(5,6) * t42;
t57 = pkin(3) * t87 + t66 * qJ(2);
t81 = qJD(2) * t68;
t80 = qJD(3) * t66;
t79 = qJD(5) * t69;
t78 = qJD(5) * t71;
t77 = t66 * qJD(2);
t76 = qJ(2) * qJD(2);
t32 = -t56 * t69 - t71 * t74;
t18 = t32 * qJD(5) - t50 * t71 - t51 * t69;
t33 = t56 * t71 - t69 * t74;
t19 = -t33 * qJD(5) + t50 * t69 - t51 * t71;
t75 = t19 * mrSges(6,1) - t18 * mrSges(6,2);
t15 = -pkin(4) * t68 + t47 * pkin(7) + t20;
t16 = -pkin(7) * t46 + t21;
t4 = t15 * t71 - t16 * t69;
t5 = t15 * t69 + t16 * t71;
t48 = -t65 * t81 - t67 * t80;
t49 = -t65 * t80 + t67 * t81;
t8 = t20 * qJD(4) + t70 * t48 + t72 * t49;
t6 = t42 * pkin(7) + t8;
t9 = -t21 * qJD(4) + t72 * t48 - t70 * t49;
t7 = t41 * pkin(7) + t9;
t2 = qJD(5) * t4 + t6 * t71 + t69 * t7;
t3 = -qJD(5) * t5 - t6 * t69 + t7 * t71;
t73 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + t85;
t64 = t68 ^ 2;
t61 = t63 * t76;
t54 = (-mrSges(6,1) * t69 - mrSges(6,2) * t71) * qJD(5) * pkin(4);
t38 = t41 * mrSges(5,2);
t37 = -mrSges(5,1) * t68 + t47 * mrSges(5,3);
t36 = mrSges(5,2) * t68 - t46 * mrSges(5,3);
t35 = -t42 * pkin(4) + t77;
t31 = t46 * pkin(4) + t57;
t23 = -mrSges(6,1) * t68 - t25 * mrSges(6,3);
t22 = mrSges(6,2) * t68 + t24 * mrSges(6,3);
t10 = t13 * mrSges(6,2);
t1 = [0.2e1 * t49 * (t68 * mrSges(4,2) - mrSges(4,3) * t87) + 0.2e1 * t48 * (-t68 * mrSges(4,1) - mrSges(4,3) * t86) - t68 * t84 - t68 * t85 - 0.2e1 * t57 * t38 + 0.2e1 * t35 * (-t24 * mrSges(6,1) + t25 * mrSges(6,2)) + 0.2e1 * t8 * t36 + 0.2e1 * t9 * t37 + 0.2e1 * t3 * t23 + 0.2e1 * t31 * t10 + 0.2e1 * t2 * t22 + (0.2e1 * t21 * mrSges(5,3) + Ifges(5,4) * t89 + Ifges(5,2) * t90 - Ifges(5,6) * t68 + t57 * t93) * t42 - (-0.2e1 * t20 * mrSges(5,3) + Ifges(5,1) * t89 + Ifges(5,4) * t90 - Ifges(5,5) * t68) * t41 + (-0.2e1 * t31 * mrSges(6,1) + 0.2e1 * t5 * mrSges(6,3) + Ifges(6,4) * t91 + Ifges(6,2) * t92 - Ifges(6,6) * t68) * t14 + (-0.2e1 * t4 * mrSges(6,3) + Ifges(6,1) * t91 + Ifges(6,4) * t92 - Ifges(6,5) * t68) * t13 + ((mrSges(5,2) * t89 - t46 * t93) * t66 + 0.2e1 * (mrSges(4,1) * t65 + mrSges(4,2) * t67) * t63 + 0.2e1 * (t63 + t64) * mrSges(3,3)) * qJD(2) + 0.2e1 * m(5) * (t20 * t9 + t21 * t8 + t57 * t77) + 0.2e1 * m(4) * (t83 * t49 + (-t65 * t82 + t53) * t48 + t61) + 0.2e1 * m(3) * (t64 * t76 + t61) + 0.2e1 * m(6) * (t2 * t5 + t3 * t4 + t31 * t35); t18 * t22 + t19 * t23 - t50 * t36 - t51 * t37 + (-t13 * t32 + t14 * t33) * mrSges(6,3) + (-t41 * t74 + t42 * t56) * mrSges(5,3) + m(6) * (t18 * t5 + t19 * t4 + t2 * t33 + t3 * t32) + m(5) * (-t20 * t51 - t21 * t50 + t56 * t8 - t74 * t9) + m(4) * (t48 * t67 + t49 * t65); 0.2e1 * m(5) * (-t50 * t56 + t51 * t74) + 0.2e1 * m(6) * (t18 * t33 + t19 * t32); t10 - t14 * mrSges(6,1) + m(6) * t35 - t38 - t42 * mrSges(5,1) + (m(4) + m(5)) * t77; 0; 0; t9 * mrSges(5,1) - t8 * mrSges(5,2) + (m(6) * (t2 * t69 + t3 * t71 - t4 * t79 + t5 * t78) + t22 * t78 - t23 * t79 + (-t13 * t71 + t14 * t69) * mrSges(6,3)) * pkin(4) + t73 + t84; -t51 * mrSges(5,1) + t50 * mrSges(5,2) + m(6) * (t18 * t69 + t19 * t71 + (-t32 * t69 + t33 * t71) * qJD(5)) * pkin(4) + t75; 0; 0.2e1 * t54; t73; t75; 0; t54; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
