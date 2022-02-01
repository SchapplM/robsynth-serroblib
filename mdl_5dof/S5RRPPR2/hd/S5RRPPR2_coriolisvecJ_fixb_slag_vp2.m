% Calculate vector of centrifugal and Coriolis load on the joints for
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:06
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPPR2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR2_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR2_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:05:27
% EndTime: 2022-01-20 10:05:31
% DurationCPUTime: 1.27s
% Computational Cost: add. (1475->168), mult. (2930->260), div. (0->0), fcn. (1659->8), ass. (0->116)
t77 = sin(pkin(9));
t74 = t77 ^ 2;
t79 = cos(pkin(9));
t75 = t79 ^ 2;
t167 = mrSges(5,3) * (t74 + t75);
t76 = qJD(1) + qJD(2);
t135 = pkin(1) * qJD(1);
t82 = sin(qJ(2));
t122 = t82 * t135;
t84 = cos(qJ(2));
t121 = t84 * t135;
t62 = t76 * pkin(2) + t121;
t78 = sin(pkin(8));
t80 = cos(pkin(8));
t46 = t80 * t122 + t78 * t62;
t38 = qJ(4) * t76 + t46;
t22 = -t79 * qJD(3) + t38 * t77;
t149 = t22 * t77;
t23 = qJD(3) * t77 + t38 * t79;
t98 = t23 * t79 + t149;
t158 = -m(5) * t98 - t167 * t76;
t81 = sin(qJ(5));
t83 = cos(qJ(5));
t100 = mrSges(6,1) * t81 + mrSges(6,2) * t83;
t144 = t76 * t77;
t49 = t100 * t144;
t166 = -t77 * t49 + t158;
t165 = -mrSges(5,1) * t79 + mrSges(5,2) * t77 - mrSges(4,1);
t133 = qJD(4) * t76;
t134 = pkin(1) * qJD(2);
t117 = qJD(1) * t134;
t109 = t84 * t117;
t110 = t82 * t117;
t52 = t80 * t109 - t78 * t110;
t43 = t52 + t133;
t148 = t43 * t74;
t127 = t79 * qJD(4);
t139 = t79 * t81;
t138 = t79 * t83;
t154 = pkin(2) * t80;
t96 = -pkin(4) * t79 - pkin(7) * t77 - pkin(3);
t60 = t96 - t154;
t70 = pkin(2) * t78 + qJ(4);
t32 = t70 * t138 + t60 * t81;
t137 = t80 * t82;
t95 = pkin(1) * (t78 * t84 + t137);
t55 = qJD(1) * t95;
t66 = t78 * t122;
t57 = t80 * t121 - t66;
t163 = -t32 * qJD(5) - t81 * t127 + t57 * t139 - t55 * t83;
t31 = -t70 * t139 + t60 * t83;
t162 = t31 * qJD(5) + t83 * t127 - t57 * t138 - t55 * t81;
t118 = t74 * qJD(5) * t76;
t155 = pkin(1) * t84;
t156 = pkin(1) * t82;
t160 = (-mrSges(3,1) * t156 - mrSges(3,2) * t155) * t76;
t159 = (qJD(4) - t57) * t77;
t45 = t62 * t80 - t66;
t105 = qJD(4) - t45;
t157 = m(4) * t45 - m(5) * t105 + (m(5) * pkin(3) - t165) * t76;
t152 = mrSges(6,3) * t77;
t151 = Ifges(6,4) * t81;
t150 = Ifges(6,4) * t83;
t147 = t43 * t75;
t145 = t76 * mrSges(4,2);
t143 = t76 * t79;
t101 = mrSges(6,1) * t83 - mrSges(6,2) * t81;
t131 = qJD(5) * t77;
t89 = t101 * t131;
t44 = t76 * t89;
t142 = t77 * t44;
t71 = pkin(2) + t155;
t136 = pkin(1) * t137 + t78 * t71;
t130 = qJD(5) * t81;
t129 = qJD(5) * t83;
t68 = t78 * t156;
t126 = t81 * t152;
t125 = t83 * t152;
t120 = mrSges(6,3) * t131;
t115 = t71 * t80 - t68;
t112 = t81 * t120;
t111 = t83 * t120;
t14 = t96 * t76 + t105;
t8 = t14 * t83 - t23 * t81;
t9 = t14 * t81 + t23 * t83;
t106 = -t8 * t81 + t83 * t9;
t104 = t76 * t112;
t103 = t76 * t111;
t99 = -Ifges(6,5) * t81 - Ifges(6,6) * t83;
t65 = qJD(5) - t143;
t47 = -mrSges(6,2) * t65 - t76 * t126;
t48 = mrSges(6,1) * t65 - t76 * t125;
t97 = t83 * t47 - t81 * t48;
t58 = t80 * t84 * t134 - qJD(2) * t68;
t42 = -t115 + t96;
t54 = qJ(4) + t136;
t12 = t54 * t138 + t42 * t81;
t11 = -t54 * t139 + t42 * t83;
t94 = t81 * (-Ifges(6,2) * t83 - t151);
t93 = t83 * (-Ifges(6,1) * t81 - t150);
t92 = (t83 * Ifges(6,1) - t151) * t77;
t91 = (-t81 * Ifges(6,2) + t150) * t77;
t90 = t99 * qJD(5);
t56 = qJD(2) * t95;
t25 = Ifges(6,6) * t65 + t76 * t91;
t26 = Ifges(6,5) * t65 + t76 * t92;
t51 = qJD(1) * t56;
t3 = qJD(5) * t8 + t43 * t138 + t51 * t81;
t4 = -qJD(5) * t9 - t43 * t139 + t51 * t83;
t85 = -mrSges(3,1) * t110 - mrSges(3,2) * t109 - t52 * mrSges(4,2) - t9 * t111 + t3 * (mrSges(6,2) * t79 - t126) + t22 * t89 + t100 * t148 + t4 * (-mrSges(6,1) * t79 - t125) + t8 * t112 + (-t143 / 0.2e1 + t65 / 0.2e1) * t77 * t90 + t165 * t51 + (-t94 + t93) * t118 + (t148 + t147) * mrSges(5,3) - ((t26 + t76 * (-Ifges(6,5) * t79 + t92)) * t81 + (t25 + t76 * (-Ifges(6,6) * t79 + t91)) * t83) * t131 / 0.2e1;
t53 = qJD(4) + t58;
t24 = t70 * t148;
t13 = t54 * t148;
t6 = -t12 * qJD(5) - t53 * t139 + t56 * t83;
t5 = t11 * qJD(5) + t53 * t138 + t56 * t81;
t1 = [t85 + m(4) * (-t51 * t115 + t52 * t136 + t46 * t58) + m(6) * (t11 * t4 + t12 * t3 + t5 * t9 + t6 * t8 + t13) + m(5) * (t54 * t147 + t13 + t51 * (-pkin(3) - t115)) - t58 * t145 + t54 * t142 + t6 * t48 + t5 * t47 + t11 * t104 - t12 * t103 - t157 * t56 + t160 * qJD(2) + (m(6) * t149 - t166) * t53; -t32 * t103 + t70 * t142 + t85 + m(5) * (t70 * t147 + t24 + t51 * (-pkin(3) - t154) + t98 * qJD(4)) + t31 * t104 + m(4) * (-t51 * t80 + t52 * t78) * pkin(2) + t159 * t49 + t163 * t48 + t162 * t47 + t133 * t167 - t160 * qJD(1) + (t159 * t22 + t162 * t9 + t163 * t8 + t3 * t32 + t31 * t4 + t24) * m(6) + t157 * t55 + (-m(4) * t46 + t145 + t158) * t57; -t79 * t44 + (-t81 ^ 2 - t83 ^ 2) * mrSges(6,3) * t118 + (m(6) * (-t8 * t129 - t9 * t130 + t3 * t83 - t4 * t81 - t43 * t79) - t47 * t130 - t48 * t129) * t77; t97 * qJD(5) + m(5) * t51 + m(6) * (t106 * qJD(5) + t3 * t81 + t4 * t83) + (-t97 * t79 - m(6) * (t9 * t138 - t8 * t139 + t149) + t166) * t76; t4 * mrSges(6,1) - t3 * mrSges(6,2) - t8 * t47 + t9 * t48 + (-t22 * t101 + t81 * t26 / 0.2e1 + t83 * t25 / 0.2e1 - t65 * t99 / 0.2e1 + (-t93 / 0.2e1 + t94 / 0.2e1) * t144 + t90 + t106 * mrSges(6,3)) * t144;];
tauc = t1(:);
