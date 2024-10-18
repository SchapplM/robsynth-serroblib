% Calculate time derivative of joint inertia matrix for
% S5PRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d2,d3,d4,d5,theta1]';
% m [6x1]
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
% Datum: 2024-09-28 18:09
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRR12_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR12_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR12_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR12_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR12_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR12_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR12_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-28 18:07:21
% EndTime: 2024-09-28 18:07:22
% DurationCPUTime: 1.02s
% Computational Cost: add. (1817->202), mult. (5503->317), div. (0->0), fcn. (5129->12), ass. (0->101)
t68 = sin(qJ(4));
t94 = qJD(4) * t68;
t118 = pkin(3) * t94;
t63 = sin(pkin(6));
t117 = t63 * t118;
t64 = sin(pkin(5));
t69 = sin(qJ(3));
t70 = sin(qJ(2));
t73 = cos(qJ(3));
t74 = cos(qJ(2));
t41 = (-t69 * t70 + t73 * t74) * t64;
t42 = (t69 * t74 + t70 * t73) * t64;
t67 = sin(qJ(5));
t101 = t63 * t67;
t65 = cos(pkin(6));
t71 = cos(qJ(5));
t98 = t65 * t71;
t50 = pkin(4) * t98 - pkin(10) * t101;
t44 = t50 * qJD(5);
t100 = t63 * t71;
t53 = -mrSges(6,2) * t65 + mrSges(6,3) * t100;
t30 = t44 * t53;
t99 = t65 * t67;
t51 = pkin(4) * t99 + pkin(10) * t100;
t45 = t51 * qJD(5);
t52 = mrSges(6,1) * t65 - mrSges(6,3) * t101;
t31 = t45 * t52;
t116 = t30 - t31;
t115 = qJD(2) + qJD(3);
t61 = t63 * pkin(10);
t54 = pkin(3) * t68 + t61;
t72 = cos(qJ(4));
t59 = pkin(3) * t72 + pkin(4);
t33 = -t54 * t67 + t59 * t98;
t95 = pkin(3) * qJD(4);
t20 = t33 * qJD(5) + (-t68 * t99 + t71 * t72) * t95;
t15 = t20 * t53;
t34 = t54 * t71 + t59 * t99;
t21 = -t34 * qJD(5) + (-t67 * t72 - t68 * t98) * t95;
t16 = t21 * t52;
t49 = (-mrSges(6,1) * t71 + mrSges(6,2) * t67) * t63;
t37 = t49 * t117;
t114 = t15 + t16 + t37;
t23 = t41 * t68 + t42 * t72;
t22 = t41 * t72 - t42 * t68;
t66 = cos(pkin(5));
t82 = t22 * t65 + t63 * t66;
t13 = t23 * t71 + t67 * t82;
t113 = 2 * m(6);
t62 = t63 ^ 2;
t43 = (mrSges(6,1) * t67 + mrSges(6,2) * t71) * t63 * qJD(5);
t112 = -0.2e1 * t43;
t111 = pkin(4) * t62;
t110 = Ifges(6,4) * t67;
t109 = Ifges(6,4) * t71;
t17 = -t22 * t63 + t65 * t66;
t108 = t17 * t63;
t60 = pkin(2) * t73 + pkin(3);
t93 = qJD(4) * t72;
t97 = t68 * t69;
t28 = t60 * t93 + (-t69 * t94 + (t72 * t73 - t97) * qJD(3)) * pkin(2);
t106 = t28 * mrSges(5,2);
t96 = t69 * t72;
t29 = -t60 * t94 + (-t69 * t93 + (-t68 * t73 - t96) * qJD(3)) * pkin(2);
t105 = t29 * t49;
t104 = (Ifges(6,6) * t65 + (Ifges(6,2) * t71 + t110) * t63) * t67;
t47 = -pkin(2) * t97 + t72 * t60;
t46 = pkin(4) + t47;
t103 = t46 * t62;
t102 = t59 * t62;
t48 = pkin(2) * t96 + t68 * t60;
t92 = qJD(5) * t67;
t91 = qJD(5) * t71;
t90 = 0.2e1 * mrSges(6,3);
t89 = t62 * t94;
t87 = t63 * t91;
t79 = -Ifges(6,6) * t63 * t92 + Ifges(6,5) * t87;
t86 = (Ifges(6,5) * t65 + (Ifges(6,1) * t67 + t109) * t63) * t87 + t65 * t79 + ((Ifges(6,1) * t71 - t110) * t92 + (-Ifges(6,2) * t67 + t109) * t91) * t62;
t85 = pkin(3) * t89;
t84 = mrSges(5,1) * t118;
t83 = pkin(3) * mrSges(5,2) * t93;
t39 = t61 + t48;
t18 = -t39 * t67 + t46 * t98;
t19 = t39 * t71 + t46 * t99;
t78 = (-mrSges(4,1) * t69 - mrSges(4,2) * t73) * qJD(3) * pkin(2);
t27 = t29 * mrSges(5,1);
t10 = qJD(5) * t18 + t28 * t71 + t29 * t99;
t5 = t10 * t53;
t11 = -qJD(5) * t19 - t28 * t67 + t29 * t98;
t6 = t11 * t52;
t77 = t27 + t5 + t6 + t86 - t106;
t12 = -t23 * t67 + t71 * t82;
t25 = t115 * t41;
t26 = t115 * t42;
t8 = qJD(4) * t22 + t25 * t72 - t26 * t68;
t9 = -qJD(4) * t23 - t25 * t68 - t26 * t72;
t3 = qJD(5) * t12 + t71 * t8 + t9 * t99;
t4 = -qJD(5) * t13 - t67 * t8 + t9 * t98;
t76 = t3 * t53 + t17 * t43 + t4 * t52 + t9 * mrSges(5,1) + (-t9 * t49 + (-t12 * t71 - t13 * t67) * qJD(5) * mrSges(6,3)) * t63 - t8 * mrSges(5,2);
t75 = -t26 * mrSges(4,1) - t25 * mrSges(4,2) + t76;
t1 = [0.2e1 * m(6) * (-t108 * t9 + t12 * t4 + t13 * t3) + 0.2e1 * m(5) * (t22 * t9 + t23 * t8) + 0.2e1 * m(4) * (t25 * t42 - t26 * t41); (-mrSges(3,1) * t70 - mrSges(3,2) * t74) * t64 * qJD(2) + m(6) * (t10 * t13 + t103 * t9 - t108 * t29 + t11 * t12 + t18 * t4 + t19 * t3) + m(5) * (t22 * t29 + t23 * t28 + t47 * t9 + t48 * t8) + m(4) * (t25 * t69 - t26 * t73 + (-t41 * t69 + t42 * t73) * qJD(3)) * pkin(2) + t75; -0.2e1 * t106 + 0.2e1 * t27 + 0.2e1 * t5 + 0.2e1 * t6 + 0.2e1 * t78 + (t10 * t19 + t103 * t29 + t11 * t18) * t113 + 0.2e1 * m(5) * (t28 * t48 + t29 * t47) + (-0.2e1 * t105 + t46 * t112 + (-t104 + (-t18 * t71 - t19 * t67) * t90) * qJD(5)) * t63 + t86; m(5) * (-t22 * t94 + t23 * t93 + t68 * t8 + t72 * t9) * pkin(3) + t75 + (t102 * t9 + t117 * t17 + t12 * t21 + t13 * t20 + t3 * t34 + t33 * t4) * m(6); m(6) * (t10 * t34 + t102 * t29 + t11 * t33 + t18 * t21 + t19 * t20) + t78 + ((-t68 * mrSges(5,1) - t72 * mrSges(5,2)) * qJD(4) + m(5) * (t28 * t68 + t29 * t72 - t47 * t94 + t48 * t93) - m(6) * t46 * t89) * pkin(3) + (-t105 + (-t46 - t59) * t43 + (-t104 + ((-t18 - t33) * t71 + (-t19 - t34) * t67) * mrSges(6,3)) * qJD(5)) * t63 + t77 + t114; 0.2e1 * t15 + 0.2e1 * t16 + 0.2e1 * t37 - 0.2e1 * t83 - 0.2e1 * t84 + (t20 * t34 + t21 * t33 - t59 * t85) * t113 + (t59 * t112 + (-t104 + (-t33 * t71 - t34 * t67) * t90) * qJD(5)) * t63 + t86; m(6) * (t111 * t9 - t12 * t45 + t13 * t44 + t3 * t51 + t4 * t50) + t76; m(6) * (t10 * t51 + t11 * t50 + t111 * t29 - t18 * t45 + t19 * t44) + (-t105 + (-pkin(4) - t46) * t43 + (-t104 + ((-t18 - t50) * t71 + (-t19 - t51) * t67) * mrSges(6,3)) * qJD(5)) * t63 + t77 + t116; m(6) * (-pkin(4) * t85 + t20 * t51 + t21 * t50 - t33 * t45 + t34 * t44) - t84 - t83 + ((-pkin(4) - t59) * t43 + (-t104 + ((-t33 - t50) * t71 + (-t34 - t51) * t67) * mrSges(6,3)) * qJD(5)) * t63 + t86 + t114 + t116; 0.2e1 * t30 - 0.2e1 * t31 + (t44 * t51 - t45 * t50) * t113 + (pkin(4) * t112 + (-t104 + (-t50 * t71 - t51 * t67) * t90) * qJD(5)) * t63 + t86; mrSges(6,1) * t4 - mrSges(6,2) * t3; mrSges(6,1) * t11 - mrSges(6,2) * t10 + t79; mrSges(6,1) * t21 - mrSges(6,2) * t20 + t79; -mrSges(6,1) * t45 - mrSges(6,2) * t44 + t79; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
