% Calculate time derivative of joint inertia matrix for
% S5RRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-05 18:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR5_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR5_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR5_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR5_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:58:10
% EndTime: 2019-12-05 18:58:13
% DurationCPUTime: 1.02s
% Computational Cost: add. (2049->201), mult. (4613->310), div. (0->0), fcn. (3463->8), ass. (0->129)
t94 = sin(qJ(4));
t91 = t94 ^ 2;
t98 = cos(qJ(4));
t92 = t98 ^ 2;
t156 = t91 + t92;
t155 = Ifges(5,1) - Ifges(5,2);
t96 = sin(qJ(2));
t99 = cos(qJ(3));
t126 = t96 * t99;
t100 = cos(qJ(2));
t86 = pkin(1) * t100 + pkin(2);
t95 = sin(qJ(3));
t60 = pkin(1) * t126 + t95 * t86;
t58 = pkin(8) + t60;
t140 = -pkin(9) - t58;
t109 = qJD(4) * t140;
t121 = qJD(3) * t99;
t122 = qJD(3) * t95;
t127 = t95 * t96;
t35 = t86 * t121 + (-t96 * t122 + (t100 * t99 - t127) * qJD(2)) * pkin(1);
t17 = t109 * t94 + t35 * t98;
t18 = t109 * t98 - t35 * t94;
t47 = t140 * t94;
t90 = t98 * pkin(9);
t48 = t58 * t98 + t90;
t93 = sin(qJ(5));
t97 = cos(qJ(5));
t26 = t47 * t97 - t48 * t93;
t3 = qJD(5) * t26 + t17 * t97 + t18 * t93;
t27 = t47 * t93 + t48 * t97;
t4 = -qJD(5) * t27 - t17 * t93 + t18 * t97;
t154 = t4 * mrSges(6,1) - t3 * mrSges(6,2);
t84 = pkin(2) * t95 + pkin(8);
t139 = -pkin(9) - t84;
t66 = t139 * t94;
t67 = t84 * t98 + t90;
t37 = t66 * t97 - t67 * t93;
t108 = qJD(4) * t139;
t112 = pkin(2) * t121;
t51 = t94 * t108 + t98 * t112;
t52 = t98 * t108 - t94 * t112;
t7 = qJD(5) * t37 + t51 * t97 + t52 * t93;
t38 = t66 * t93 + t67 * t97;
t8 = -qJD(5) * t38 - t51 * t93 + t52 * t97;
t153 = t8 * mrSges(6,1) - t7 * mrSges(6,2);
t146 = -pkin(9) - pkin(8);
t79 = t146 * t94;
t80 = pkin(8) * t98 + t90;
t49 = t79 * t97 - t80 * t93;
t111 = qJD(4) * t146;
t71 = t94 * t111;
t72 = t98 * t111;
t24 = qJD(5) * t49 + t71 * t97 + t72 * t93;
t50 = t79 * t93 + t80 * t97;
t25 = -qJD(5) * t50 - t71 * t93 + t72 * t97;
t152 = t25 * mrSges(6,1) - t24 * mrSges(6,2);
t151 = t156 * t99;
t150 = qJD(4) + qJD(5);
t68 = -t93 * t94 + t97 * t98;
t137 = mrSges(6,3) * t68;
t14 = t24 * t137;
t105 = mrSges(5,1) * t94 + mrSges(5,2) * t98;
t70 = t105 * qJD(4);
t144 = pkin(3) * t70;
t44 = t150 * t68;
t69 = t93 * t98 + t94 * t97;
t45 = t150 * t69;
t20 = mrSges(6,1) * t45 + mrSges(6,2) * t44;
t143 = pkin(4) * t98;
t87 = -pkin(3) - t143;
t16 = t87 * t20;
t138 = mrSges(6,3) * t45;
t28 = t50 * t138;
t46 = -mrSges(6,1) * t68 + mrSges(6,2) * t69;
t120 = qJD(4) * t94;
t89 = pkin(4) * t120;
t40 = t46 * t89;
t149 = t14 + t16 - t28 + t40 - t144;
t145 = pkin(2) * t99;
t74 = t87 - t145;
t13 = t74 * t20;
t19 = t38 * t138;
t113 = pkin(2) * t122;
t73 = t89 + t113;
t30 = t73 * t46;
t5 = t7 * t137;
t85 = -pkin(3) - t145;
t56 = t85 * t70;
t75 = -mrSges(5,1) * t98 + mrSges(5,2) * t94;
t63 = t75 * t113;
t107 = mrSges(5,3) * t112;
t77 = t91 * t107;
t78 = t92 * t107;
t148 = t13 - t19 + t30 + t5 + t56 + t63 + t77 + t78;
t147 = 2 * m(6);
t134 = Ifges(5,6) * t94;
t132 = t35 * mrSges(4,2);
t131 = t35 * t91;
t130 = t35 * t92;
t125 = Ifges(6,5) * t44 - Ifges(6,6) * t45;
t124 = pkin(2) * qJD(3);
t123 = pkin(4) * qJD(5);
t119 = qJD(4) * t98;
t118 = qJD(5) * t93;
t117 = qJD(5) * t97;
t116 = 0.2e1 * mrSges(6,3);
t115 = t97 * t44 * mrSges(6,3);
t110 = t156 * t35;
t59 = -pkin(1) * t127 + t86 * t99;
t57 = -pkin(3) - t59;
t106 = -t95 * mrSges(4,1) - t99 * mrSges(4,2);
t104 = -t93 * pkin(4) * t138 + Ifges(5,5) * t119 + t125 + (t68 * t97 + t69 * t93) * mrSges(6,3) * t123;
t103 = (-mrSges(3,1) * t96 - mrSges(3,2) * t100) * qJD(2) * pkin(1);
t102 = 0.2e1 * t69 * t44 * Ifges(6,1) - 0.2e1 * t68 * Ifges(6,2) * t45 + (-0.2e1 * Ifges(5,4) * t94 + t155 * t98) * t120 + (0.2e1 * Ifges(5,4) * t98 + t155 * t94) * t119 + 0.2e1 * (t68 * t44 - t69 * t45) * Ifges(6,4);
t36 = t86 * t122 + (t96 * t121 + (t100 * t95 + t126) * qJD(2)) * pkin(1);
t1 = t3 * t137;
t53 = t57 - t143;
t10 = t53 * t20;
t31 = t89 + t36;
t15 = t31 * t46;
t29 = t36 * t75;
t32 = mrSges(5,3) * t131;
t33 = mrSges(5,3) * t130;
t34 = t36 * mrSges(4,1);
t43 = t57 * t70;
t9 = t27 * t138;
t101 = t1 + t10 + t102 + t15 + t29 + t32 + t33 - t34 + t43 - t9 - t132;
t65 = (-mrSges(6,1) * t93 - mrSges(6,2) * t97) * t123;
t2 = [(t26 * t4 + t27 * t3 + t31 * t53) * t147 + 0.2e1 * m(4) * (t35 * t60 - t36 * t59) + t102 + 0.2e1 * t43 - 0.2e1 * t34 + 0.2e1 * t33 + 0.2e1 * t32 + 0.2e1 * t29 + 0.2e1 * t15 + 0.2e1 * t10 - 0.2e1 * t9 - 0.2e1 * t132 + 0.2e1 * t1 + 0.2e1 * t103 + 0.2e1 * m(5) * (t110 * t58 + t36 * t57) + (-t26 * t44 - t4 * t69) * t116; t101 + (m(4) * (t35 * t95 - t36 * t99) + (m(5) * (t151 * t58 + t57 * t95) + m(4) * (-t59 * t95 + t60 * t99) + t106) * qJD(3)) * pkin(2) + ((-t4 - t8) * t69 + (-t26 - t37) * t44) * mrSges(6,3) + m(5) * (t36 * t85 + (t130 + t131) * t84) + m(6) * (t26 * t8 + t27 * t7 + t3 * t38 + t31 * t74 + t37 * t4 + t53 * t73) + t103 + t148; (t37 * t8 + t38 * t7 + t73 * t74) * t147 + t102 + 0.2e1 * t77 + 0.2e1 * t78 + 0.2e1 * t63 + 0.2e1 * t56 + 0.2e1 * t30 - 0.2e1 * t19 + 0.2e1 * t13 + 0.2e1 * t5 + 0.2e1 * (m(5) * (t151 * t84 + t85 * t95) + t106) * t124 + (-t37 * t44 - t8 * t69) * t116; t101 + m(5) * (-pkin(3) * t36 + pkin(8) * t110) + m(6) * (t24 * t27 + t25 * t26 + t3 * t50 + t31 * t87 + t4 * t49 + t53 * t89) + ((-t25 - t4) * t69 + (-t26 - t49) * t44) * mrSges(6,3) + t149; t102 + m(6) * (t24 * t38 + t25 * t37 + t49 * t8 + t50 * t7 + t73 * t87 + t74 * t89) + t148 + ((-t25 - t8) * t69 + (-t37 - t49) * t44) * mrSges(6,3) + (m(5) * (-pkin(3) * t95 + t151 * pkin(8)) + t106) * t124 + t149; 0.2e1 * t14 - 0.2e1 * t28 + 0.2e1 * t40 + 0.2e1 * t16 - 0.2e1 * t144 + (t24 * t50 + t25 * t49 + t87 * t89) * t147 + (-t25 * t69 - t44 * t49) * t116 + t102; -t105 * t35 + (t58 * t75 - t134) * qJD(4) + (-t115 + m(6) * (t117 * t27 - t118 * t26 + t3 * t93 + t4 * t97)) * pkin(4) + t104 + t154; -t105 * t112 + (t75 * t84 - t134) * qJD(4) + (-t115 + m(6) * (t117 * t38 - t118 * t37 + t7 * t93 + t8 * t97)) * pkin(4) + t104 + t153; (pkin(8) * t75 - t134) * qJD(4) + (-t115 + m(6) * (t117 * t50 - t118 * t49 + t24 * t93 + t25 * t97)) * pkin(4) + t104 + t152; 0.2e1 * t65; t125 + t154; t125 + t153; t125 + t152; t65; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
