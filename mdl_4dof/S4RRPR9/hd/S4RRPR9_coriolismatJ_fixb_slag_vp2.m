% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RRPR9_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR9_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR9_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR9_coriolismatJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR9_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR9_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR9_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:09:23
% EndTime: 2019-12-31 17:09:26
% DurationCPUTime: 1.30s
% Computational Cost: add. (2981->228), mult. (6595->340), div. (0->0), fcn. (6387->6), ass. (0->115)
t126 = sin(pkin(7));
t127 = cos(pkin(7));
t129 = cos(qJ(4));
t167 = sin(qJ(4));
t100 = -t129 * t126 - t127 * t167;
t128 = sin(qJ(2));
t83 = t100 * t128;
t188 = t83 * mrSges(5,3);
t99 = -t167 * t126 + t129 * t127;
t96 = Ifges(5,4) * t99;
t60 = -Ifges(5,1) * t100 + t96;
t187 = Ifges(5,2) * t100 + t60 + t96;
t130 = cos(qJ(2));
t112 = t128 * pkin(2) - t130 * qJ(3);
t147 = t126 * t128;
t75 = pkin(5) * t147 + t127 * t112;
t145 = t127 * t128;
t76 = -pkin(5) * t145 + t126 * t112;
t186 = -t75 * t126 + t76 * t127;
t183 = t127 ^ 2;
t182 = m(4) * pkin(5);
t163 = pkin(6) + qJ(3);
t110 = t163 * t126;
t111 = t163 * t127;
t64 = -t110 * t167 + t129 * t111;
t181 = -t64 / 0.2e1;
t180 = t83 / 0.2e1;
t84 = t100 * t130;
t178 = t84 / 0.2e1;
t85 = t99 * t128;
t177 = t85 / 0.2e1;
t86 = t99 * t130;
t176 = t86 / 0.2e1;
t174 = t99 / 0.2e1;
t172 = -t100 / 0.2e1;
t171 = t100 / 0.2e1;
t170 = -t126 / 0.2e1;
t169 = t127 / 0.2e1;
t168 = t128 / 0.2e1;
t166 = Ifges(5,4) * t85;
t165 = t84 * mrSges(5,1);
t164 = t86 * mrSges(5,2);
t162 = Ifges(5,5) * t83 - Ifges(5,6) * t85;
t161 = Ifges(5,5) * t99 + Ifges(5,6) * t100;
t160 = Ifges(4,4) * t126;
t159 = Ifges(4,4) * t127;
t158 = Ifges(5,4) * t100;
t157 = Ifges(4,5) * t127;
t156 = Ifges(4,6) * t126;
t103 = t130 * mrSges(4,2) - mrSges(4,3) * t147;
t146 = t126 * t130;
t104 = -t128 * mrSges(4,2) - mrSges(4,3) * t146;
t105 = -t130 * mrSges(4,1) - mrSges(4,3) * t145;
t144 = t127 * t130;
t106 = t128 * mrSges(4,1) - mrSges(4,3) * t144;
t142 = pkin(3) * t126 + pkin(5);
t107 = t142 * t128;
t108 = t142 * t130;
t135 = -Ifges(5,5) * t86 / 0.2e1 - Ifges(5,6) * t84 / 0.2e1;
t152 = t127 * mrSges(4,2);
t153 = t126 * Ifges(4,2);
t109 = -t130 * pkin(2) - t128 * qJ(3) - pkin(1);
t98 = t127 * t109;
t56 = -pkin(6) * t145 + t98 + (-pkin(5) * t126 - pkin(3)) * t130;
t70 = pkin(5) * t144 + t126 * t109;
t61 = -pkin(6) * t147 + t70;
t18 = t129 * t56 - t167 * t61;
t19 = t129 * t61 + t167 * t56;
t57 = t128 * pkin(3) - pkin(6) * t144 + t75;
t62 = -pkin(6) * t146 + t76;
t20 = t129 * t57 - t167 * t62;
t21 = t129 * t62 + t167 * t57;
t40 = Ifges(5,2) * t83 - Ifges(5,6) * t130 + t166;
t41 = Ifges(5,4) * t86 + Ifges(5,2) * t84 + Ifges(5,6) * t128;
t80 = Ifges(5,4) * t83;
t42 = Ifges(5,1) * t85 - Ifges(5,5) * t130 + t80;
t43 = Ifges(5,1) * t86 + Ifges(5,4) * t84 + Ifges(5,5) * t128;
t47 = t164 - t165;
t65 = t130 * mrSges(5,2) + t188;
t66 = -t128 * mrSges(5,2) + t84 * mrSges(5,3);
t67 = -t130 * mrSges(5,1) - t85 * mrSges(5,3);
t68 = t128 * mrSges(5,1) - t86 * mrSges(5,3);
t69 = -pkin(5) * t146 + t98;
t81 = Ifges(4,6) * t128 + (-t153 + t159) * t130;
t82 = Ifges(4,5) * t128 + (t127 * Ifges(4,1) - t160) * t130;
t154 = t126 * mrSges(4,1);
t90 = (t152 + t154) * t130;
t1 = t21 * t65 + t19 * t66 + t20 * t67 + t18 * t68 + t41 * t180 + t40 * t178 + t43 * t177 + t42 * t176 + t76 * t103 + t70 * t104 + t75 * t105 + t69 * t106 + t107 * t47 + t108 * (-t83 * mrSges(5,1) + t85 * mrSges(5,2)) + m(5) * (t107 * t108 + t18 * t20 + t19 * t21) + m(4) * (t69 * t75 + t70 * t76) + (pkin(5) * t90 + t81 * t170 + t82 * t169 - pkin(1) * mrSges(3,1) + Ifges(5,5) * t177 + Ifges(5,6) * t180 + (-Ifges(3,4) + t157 / 0.2e1 - t156 / 0.2e1) * t128) * t128 + (-pkin(1) * mrSges(3,2) + (t183 * Ifges(4,1) / 0.2e1 - Ifges(3,2) + Ifges(3,1) - Ifges(5,3) - Ifges(4,3) + (t152 + t182) * pkin(5) + (pkin(5) * mrSges(4,1) - t159 + t153 / 0.2e1) * t126) * t128 + t135 + (Ifges(3,4) + t156 - t157) * t130) * t130;
t155 = t1 * qJD(1);
t46 = t85 * mrSges(5,1) + t83 * mrSges(5,2);
t48 = -Ifges(5,2) * t85 + t80;
t49 = Ifges(5,1) * t83 - t166;
t4 = t18 * t65 - t19 * t67 - t130 * t162 / 0.2e1 + t107 * t46 + (-t19 * mrSges(5,3) + t49 / 0.2e1 - t40 / 0.2e1) * t85 + (-t18 * mrSges(5,3) + t42 / 0.2e1 + t48 / 0.2e1) * t83;
t151 = t4 * qJD(1);
t9 = m(5) * (-t18 * t85 + t19 * t83) + t83 * t65 - t85 * t67 + (-t126 * t103 - t127 * t105 + m(4) * (-t126 * t70 - t127 * t69)) * t128;
t148 = t9 * qJD(1);
t139 = Ifges(5,1) * t99 + t158;
t121 = -t127 * pkin(3) - pkin(2);
t58 = -t100 * mrSges(5,1) + t99 * mrSges(5,2);
t59 = Ifges(5,2) * t99 - t158;
t63 = -t129 * t110 - t111 * t167;
t131 = t107 * t58 / 0.2e1 + t121 * t46 / 0.2e1 - t130 * t161 / 0.4e1 + t67 * t181 + (t65 / 0.2e1 - t188 / 0.2e1) * t63 + (t181 * mrSges(5,3) - t59 / 0.4e1 + t139 / 0.4e1) * t85 + t187 * t83 / 0.4e1 + (t42 + t48) * t99 / 0.4e1 + (t40 / 0.4e1 - t49 / 0.4e1) * t100;
t133 = Ifges(5,3) * t168 + t20 * mrSges(5,1) / 0.2e1 - t21 * mrSges(5,2) / 0.2e1 - t135;
t2 = t131 - t133;
t7 = -t121 * t58 + t139 * t171 + t172 * t59 - t187 * t99 / 0.2e1;
t138 = t2 * qJD(1) - t7 * qJD(2);
t11 = -m(5) * (t63 * t100 + t64 * t99) + (-t100 ^ 2 - t99 ^ 2) * mrSges(5,3) + (m(4) * qJ(3) + mrSges(4,3)) * (-t126 ^ 2 - t183);
t132 = (t172 * t85 + t174 * t83) * mrSges(5,3) + m(4) * (-t126 * t69 + t127 * t70) / 0.2e1 + m(5) * (t100 * t18 + t99 * t19 - t63 * t85 + t64 * t83) / 0.2e1 + t67 * t171 + t105 * t170 + t103 * t169 + t65 * t174;
t134 = -m(5) * t108 / 0.2e1 + t165 / 0.2e1 - t164 / 0.2e1;
t6 = (-t152 / 0.2e1 - t154 / 0.2e1 - t182 / 0.2e1) * t130 + t132 + t134;
t137 = t6 * qJD(1) - t11 * qJD(2);
t136 = t46 * qJD(1) + t58 * qJD(2);
t5 = mrSges(4,2) * t144 / 0.2e1 + mrSges(4,1) * t146 / 0.2e1 + t130 * t182 / 0.2e1 + t132 - t134;
t3 = t131 + t133;
t8 = [t1 * qJD(2) + t9 * qJD(3) + t4 * qJD(4), t5 * qJD(3) + t3 * qJD(4) + t155 + (m(5) * (t121 * t108 + t63 * t20 + t64 * t21) + t64 * t66 + t63 * t68 + t59 * t178 + t60 * t176 - pkin(2) * t90 + t41 * t174 + t43 * t172 + t108 * (-t99 * mrSges(5,1) - t100 * mrSges(5,2)) + t121 * t47 + t126 * t82 / 0.2e1 + t81 * t169 + (m(4) * t186 + t127 * t104 - t126 * t106) * qJ(3) + ((Ifges(4,1) * t126 + t159) * t169 + (Ifges(4,2) * t127 + t160) * t170 + Ifges(3,5) + (-m(4) * pkin(2) - t127 * mrSges(4,1) + t126 * mrSges(4,2) - mrSges(3,1)) * pkin(5)) * t130 + (Ifges(4,5) * t126 - Ifges(5,5) * t100 + Ifges(4,6) * t127 + Ifges(5,6) * t99) * t168 + (pkin(5) * mrSges(3,2) - Ifges(3,6)) * t128 + (t20 * t100 + t21 * t99) * mrSges(5,3) + t186 * mrSges(4,3)) * qJD(2), t5 * qJD(2) + t148, t151 + t3 * qJD(2) + (-t19 * mrSges(5,1) - t18 * mrSges(5,2) + t162) * qJD(4); t6 * qJD(3) + t2 * qJD(4) - t155, -t11 * qJD(3) - t7 * qJD(4), t137, (-t64 * mrSges(5,1) - t63 * mrSges(5,2) + t161) * qJD(4) + t138; -t6 * qJD(2) + t46 * qJD(4) - t148, t58 * qJD(4) - t137, 0, t136; -t2 * qJD(2) - t46 * qJD(3) - t151, -t58 * qJD(3) - t138, -t136, 0;];
Cq = t8;
