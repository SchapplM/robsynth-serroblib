% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPP3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP3_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP3_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP3_coriolisvecJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP3_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPP3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:12:20
% EndTime: 2019-12-31 18:12:26
% DurationCPUTime: 2.49s
% Computational Cost: add. (1219->250), mult. (3446->289), div. (0->0), fcn. (2216->4), ass. (0->115)
t146 = Ifges(4,1) + Ifges(6,3);
t145 = Ifges(6,4) + Ifges(5,5);
t144 = Ifges(4,5) + Ifges(6,5);
t140 = Ifges(6,2) + Ifges(5,3);
t110 = mrSges(5,2) - mrSges(4,1);
t108 = pkin(6) + qJ(2);
t78 = cos(pkin(7));
t68 = t108 * t78;
t65 = qJD(1) * t68;
t80 = sin(qJ(3));
t116 = t65 * t80;
t124 = cos(qJ(3));
t77 = sin(pkin(7));
t67 = t108 * t77;
t64 = qJD(1) * t67;
t36 = t124 * t64 + t116;
t63 = t124 * t77 + t80 * t78;
t57 = t63 * qJD(1);
t86 = pkin(4) * t57 + t36;
t143 = t86 + qJD(4);
t115 = t77 * t80;
t96 = t124 * t78;
t91 = qJD(1) * t96;
t56 = qJD(1) * t115 - t91;
t142 = -t56 * pkin(4) + qJD(5);
t139 = Ifges(5,6) - Ifges(6,6);
t120 = Ifges(5,6) * t57;
t51 = Ifges(6,6) * t57;
t138 = t145 * qJD(3) + t140 * t56 - t120 + t51;
t119 = Ifges(6,6) * t56;
t53 = Ifges(4,4) * t56;
t137 = t144 * qJD(3) + t146 * t57 + t119 - t53;
t136 = -qJD(4) - t36;
t37 = t124 * t65 - t80 * t64;
t135 = -t37 - t142;
t134 = Ifges(5,4) - t144;
t133 = -Ifges(4,6) + t145;
t132 = (m(3) * qJ(2) + mrSges(3,3)) * (t77 ^ 2 + t78 ^ 2);
t71 = qJD(2) * t96;
t94 = qJD(3) * t124;
t20 = (qJD(2) * t77 + qJD(3) * t68) * t80 + t67 * t94 - t71;
t104 = qJD(1) * t71 - t64 * t94;
t100 = qJD(1) * qJD(2);
t95 = t77 * t100;
t10 = qJD(3) * (-qJD(4) + t116) + t80 * t95 - t104;
t131 = -t56 / 0.2e1;
t130 = t56 / 0.2e1;
t129 = -t57 / 0.2e1;
t128 = t57 / 0.2e1;
t123 = mrSges(4,3) * t56;
t122 = mrSges(4,3) * t57;
t121 = Ifges(4,4) * t57;
t84 = t63 * qJD(2);
t83 = qJD(1) * t84;
t14 = qJD(3) * t37 + t83;
t38 = t124 * t67 + t68 * t80;
t118 = t14 * t38;
t117 = t14 * t63;
t113 = qJD(3) / 0.2e1;
t112 = -mrSges(5,1) - mrSges(6,1);
t111 = -mrSges(5,1) - mrSges(4,3);
t109 = pkin(3) + qJ(5);
t40 = -qJD(3) * mrSges(4,2) - t123;
t43 = mrSges(5,1) * t56 - qJD(3) * mrSges(5,3);
t107 = t40 - t43;
t106 = -t57 * mrSges(5,1) - t110 * qJD(3) - t122;
t44 = -mrSges(6,1) * t56 + qJD(3) * mrSges(6,2);
t105 = t43 - t44;
t101 = qJ(4) * t56;
t99 = -Ifges(4,4) - t139;
t98 = qJD(3) * t115;
t97 = -pkin(2) * t78 - pkin(1);
t49 = qJD(1) * t98 - qJD(3) * t91;
t90 = qJ(4) * t49 - qJD(4) * t57;
t58 = -t78 * t94 + t98;
t89 = qJ(4) * t58 - qJD(4) * t63;
t87 = -qJ(4) * t63 + t97;
t39 = t124 * t68 - t80 * t67;
t66 = qJD(1) * t97 + qJD(2);
t59 = t63 * qJD(3);
t31 = -qJD(3) * qJ(4) - t37;
t82 = -qJ(4) * t57 + t66;
t21 = qJD(3) * t39 + t84;
t62 = -t96 + t115;
t52 = Ifges(5,6) * t56;
t50 = qJD(1) * t59;
t48 = t50 * mrSges(6,3);
t47 = t49 * mrSges(4,2);
t46 = t49 * mrSges(5,3);
t42 = mrSges(6,1) * t57 - qJD(3) * mrSges(6,3);
t35 = pkin(3) * t62 + t87;
t34 = -mrSges(5,2) * t56 - mrSges(5,3) * t57;
t33 = pkin(3) * t57 + t101;
t32 = -mrSges(6,2) * t57 + mrSges(6,3) * t56;
t30 = -qJD(3) * pkin(3) - t136;
t28 = -Ifges(4,2) * t56 + Ifges(4,6) * qJD(3) + t121;
t27 = Ifges(5,4) * qJD(3) - Ifges(5,2) * t57 + t52;
t23 = -t62 * pkin(4) + t39;
t22 = pkin(4) * t63 + t38;
t19 = pkin(3) * t56 + t82;
t16 = t109 * t62 + t87;
t15 = pkin(3) * t59 + t89;
t13 = (-qJD(3) * t65 - t95) * t80 + t104;
t12 = t109 * t57 + t101;
t11 = -t31 + t142;
t9 = -qJD(3) * t109 + t143;
t8 = pkin(3) * t50 + t90;
t7 = -t58 * pkin(4) + t21;
t6 = -pkin(4) * t59 - t20;
t5 = t109 * t56 + t82;
t4 = -t49 * pkin(4) + t83 + (-qJD(5) + t37) * qJD(3);
t3 = -pkin(4) * t50 - t10;
t2 = qJD(5) * t62 + t109 * t59 + t89;
t1 = qJD(5) * t56 + t109 * t50 + t90;
t17 = [m(6) * (t1 * t16 + t11 * t6 + t2 * t5 + t22 * t4 + t23 * t3 + t7 * t9) + m(4) * (t13 * t39 - t20 * t37 + t21 * t36 + t118) + m(5) * (-t10 * t39 + t15 * t19 + t20 * t31 + t21 * t30 + t35 * t8 + t118) - t106 * t21 - t107 * t20 - t97 * t47 + (t10 * t62 + t117) * mrSges(5,1) + (-mrSges(6,1) * t22 + mrSges(6,2) * t16 + t111 * t38 + (-Ifges(5,2) - t146) * t63 - t99 * t62) * t49 + t16 * t48 + (t97 * mrSges(4,1) - t35 * mrSges(5,2) - t23 * mrSges(6,1) + t111 * t39 + t99 * t63 + (Ifges(4,2) + t140) * t62) * t50 + t35 * t46 + 0.2e1 * t132 * t100 + (-t3 * t62 + t4 * t63) * mrSges(6,1) + t1 * (-mrSges(6,2) * t63 + mrSges(6,3) * t62) + t8 * (-mrSges(5,2) * t62 - mrSges(5,3) * t63) + t7 * t42 + t6 * t44 + t2 * t32 + t15 * t34 + (-t13 * t62 + t117) * mrSges(4,3) + (-t37 * mrSges(4,3) + t31 * mrSges(5,1) + (-Ifges(4,4) + Ifges(6,6)) * t128 + t138 / 0.2e1 + t140 * t130 + t133 * t113 + t66 * mrSges(4,1) + t5 * mrSges(6,3) - t28 / 0.2e1 - t19 * mrSges(5,2) + Ifges(5,6) * t129 - Ifges(4,2) * t131 - t11 * mrSges(6,1)) * t59 + (-t36 * mrSges(4,3) - t30 * mrSges(5,1) - t146 * t128 - t137 / 0.2e1 + t139 * t130 + t134 * t113 - t66 * mrSges(4,2) + t5 * mrSges(6,2) + t27 / 0.2e1 + t19 * mrSges(5,3) + Ifges(5,2) * t129 - Ifges(4,4) * t131 - t9 * mrSges(6,1)) * t58; t49 * mrSges(6,2) + t46 - t47 + t48 - t110 * t50 + (-t42 + t106) * t57 + (t40 - t105) * t56 - m(4) * (t36 * t57 - t37 * t56) - t132 * qJD(1) ^ 2 + (t11 * t56 - t57 * t9 + t1) * m(6) + (-t30 * t57 - t31 * t56 + t8) * m(5); (-t146 * t56 - t121 + t138 + t51) * t129 + (Ifges(5,2) * t56 + t120 + t28) * t128 + (t106 + t122) * t37 + (t107 + t123) * t36 - t105 * qJD(4) + t110 * t14 + (-Ifges(4,2) * t57 + t137 - t53) * t130 + (t140 * t57 - t119 + t27 + t52) * t131 + (qJ(4) * t112 + t133) * t50 - (t133 * t57 + t134 * t56) * qJD(3) / 0.2e1 + t135 * t42 + (-pkin(3) * t14 - qJ(4) * t10 + t136 * t31 - t19 * t33 - t30 * t37) * m(5) + t86 * t44 + (mrSges(5,1) * pkin(3) + mrSges(6,1) * t109 + t134) * t49 - t66 * (mrSges(4,1) * t57 - mrSges(4,2) * t56) - t19 * (-t57 * mrSges(5,2) + t56 * mrSges(5,3)) - t5 * (t56 * mrSges(6,2) + t57 * mrSges(6,3)) - t12 * t32 - t33 * t34 - t13 * mrSges(4,2) + t3 * mrSges(6,2) - t4 * mrSges(6,3) - t10 * mrSges(5,3) + (qJ(4) * t3 - t4 * t109 + t143 * t11 - t12 * t5 + t135 * t9) * m(6) + (t30 * t56 - t31 * t57) * mrSges(5,1) + (t11 * t57 + t56 * t9) * mrSges(6,1); (t32 + t34) * t57 + t112 * t49 + t105 * qJD(3) + (-t11 * qJD(3) + t5 * t57 + t4) * m(6) + (qJD(3) * t31 + t19 * t57 + t14) * m(5); -t50 * mrSges(6,1) + qJD(3) * t42 - t56 * t32 + 0.2e1 * (t3 / 0.2e1 + t9 * t113 + t5 * t131) * m(6);];
tauc = t17(:);
