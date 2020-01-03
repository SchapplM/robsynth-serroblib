% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRRP6_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP6_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP6_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP6_coriolisvecJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP6_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP6_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP6_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:18:10
% EndTime: 2019-12-31 17:18:17
% DurationCPUTime: 3.36s
% Computational Cost: add. (1345->297), mult. (3612->411), div. (0->0), fcn. (1963->4), ass. (0->150)
t206 = qJD(2) / 0.2e1;
t204 = Ifges(4,4) + Ifges(5,4);
t129 = Ifges(3,5) * t206;
t205 = Ifges(4,1) + Ifges(5,1);
t195 = Ifges(4,5) + Ifges(5,5);
t203 = Ifges(4,2) + Ifges(5,2);
t194 = Ifges(5,6) + Ifges(4,6);
t98 = sin(qJ(2));
t146 = qJD(1) * t98;
t97 = sin(qJ(3));
t135 = t97 * t146;
t99 = cos(qJ(3));
t143 = qJD(2) * t99;
t74 = -t135 + t143;
t202 = t204 * t74;
t134 = t99 * t146;
t144 = qJD(2) * t97;
t75 = t134 + t144;
t201 = t204 * t75;
t200 = t204 * t99;
t199 = t204 * t97;
t100 = cos(qJ(2));
t79 = -pkin(2) * t100 - pkin(6) * t98 - pkin(1);
t66 = t79 * qJD(1);
t139 = qJD(1) * t100;
t96 = pkin(5) * t139;
t86 = qJD(2) * pkin(6) + t96;
t32 = t66 * t99 - t86 * t97;
t33 = t66 * t97 + t86 * t99;
t107 = t32 * t99 + t33 * t97;
t16 = -qJ(4) * t75 + t32;
t91 = qJD(3) - t139;
t11 = pkin(3) * t91 + t16;
t119 = mrSges(5,1) * t97 + mrSges(5,2) * t99;
t121 = mrSges(4,1) * t97 + mrSges(4,2) * t99;
t159 = Ifges(5,6) * t97;
t160 = Ifges(4,6) * t97;
t161 = Ifges(5,5) * t99;
t162 = Ifges(4,5) * t99;
t17 = qJ(4) * t74 + t33;
t171 = t99 / 0.2e1;
t174 = -t97 / 0.2e1;
t177 = t75 / 0.2e1;
t186 = t195 * t91 + t205 * t75 + t202;
t187 = t194 * t91 + t203 * t74 + t201;
t188 = t91 / 0.2e1;
t189 = t74 / 0.2e1;
t190 = t205 * t99 - t199;
t192 = -t203 * t97 + t200;
t85 = -qJD(2) * pkin(2) + pkin(5) * t146;
t42 = -pkin(3) * t74 + qJD(4) + t85;
t184 = t107 * mrSges(4,3) + (t11 * t99 + t17 * t97) * mrSges(5,3) - t42 * t119 - t85 * t121 - t192 * t189 - t190 * t177 - (-t159 + t161 - t160 + t162) * t188 - t187 * t174 - t186 * t171;
t94 = Ifges(3,4) * t139;
t198 = -t184 + Ifges(3,1) * t146 / 0.2e1 + t129 + t94 / 0.2e1;
t137 = qJD(1) * qJD(2);
t130 = t98 * t137;
t136 = qJD(2) * qJD(3);
t138 = qJD(2) * t100;
t142 = qJD(3) * t97;
t43 = t99 * t136 + (t138 * t99 - t142 * t98) * qJD(1);
t141 = qJD(3) * t99;
t105 = t138 * t97 + t141 * t98;
t44 = -qJD(1) * t105 - t136 * t97;
t197 = t130 * t194 + t203 * t44 + t204 * t43;
t196 = t130 * t195 + t204 * t44 + t205 * t43;
t193 = t203 * t99 + t199;
t191 = t205 * t97 + t200;
t128 = -Ifges(3,6) * qJD(2) / 0.2e1;
t158 = pkin(5) * t100;
t92 = t99 * t158;
t53 = t79 * t97 + t92;
t185 = t194 * t99 + t195 * t97;
t131 = Ifges(5,3) / 0.2e1 + Ifges(4,3) / 0.2e1;
t132 = Ifges(5,6) / 0.2e1 + Ifges(4,6) / 0.2e1;
t133 = Ifges(5,5) / 0.2e1 + Ifges(4,5) / 0.2e1;
t167 = Ifges(3,4) * t98;
t183 = -t131 * t91 - t132 * t74 - t133 * t75 - t11 * mrSges(5,1) - t32 * mrSges(4,1) - t128 + (Ifges(3,2) * t100 + t167) * qJD(1) / 0.2e1 + t17 * mrSges(5,2) + t33 * mrSges(4,2) - t194 * t189 - (Ifges(5,3) + Ifges(4,3)) * t188 - t195 * t177;
t182 = t43 / 0.2e1;
t181 = t44 / 0.2e1;
t180 = -t74 / 0.2e1;
t178 = -t75 / 0.2e1;
t176 = -t91 / 0.2e1;
t170 = pkin(1) * mrSges(3,1);
t169 = pkin(1) * mrSges(3,2);
t168 = pkin(3) * t97;
t155 = t98 * t99;
t153 = -qJ(4) - pkin(6);
t123 = pkin(2) * t98 - pkin(6) * t100;
t77 = t123 * qJD(2);
t152 = t141 * t79 + t77 * t97;
t151 = pkin(5) * t144 * t98 + t77 * t99;
t76 = t123 * qJD(1);
t45 = pkin(5) * t135 + t76 * t99;
t150 = qJ(4) * t98;
t147 = qJ(4) * t100;
t145 = qJD(2) * mrSges(3,2);
t140 = qJD(4) * t99;
t12 = -t44 * mrSges(5,1) + mrSges(5,2) * t43;
t127 = qJD(3) * t153;
t126 = pkin(5) * t130;
t125 = m(4) * t85 - qJD(2) * mrSges(3,1) - mrSges(4,1) * t74 + mrSges(4,2) * t75 + mrSges(3,3) * t146;
t67 = qJD(1) * t77;
t5 = -t126 * t99 + t141 * t66 - t142 * t86 + t67 * t97;
t6 = -qJD(3) * t33 + t126 * t97 + t67 * t99;
t124 = t5 * t99 - t6 * t97;
t122 = mrSges(4,1) * t99 - mrSges(4,2) * t97;
t120 = mrSges(5,1) * t99 - mrSges(5,2) * t97;
t106 = pkin(3) * t98 - t147 * t99;
t1 = pkin(3) * t130 - qJ(4) * t43 - qJD(4) * t75 + t6;
t2 = qJ(4) * t44 + qJD(4) * t74 + t5;
t104 = -t6 * mrSges(4,1) - t1 * mrSges(5,1) + t5 * mrSges(4,2) + t2 * mrSges(5,2);
t93 = -pkin(3) * t99 - pkin(2);
t90 = Ifges(4,3) * t130;
t89 = Ifges(5,3) * t130;
t84 = t153 * t99;
t83 = t153 * t97;
t81 = mrSges(3,3) * t139 - t145;
t78 = (pkin(5) + t168) * t98;
t73 = t99 * t79;
t69 = t139 * t168 + t96;
t62 = t97 * t76;
t59 = -qJD(4) * t97 + t127 * t99;
t58 = t127 * t97 + t140;
t52 = -t158 * t97 + t73;
t51 = pkin(3) * t105 + pkin(5) * t138;
t50 = mrSges(4,1) * t91 - mrSges(4,3) * t75;
t49 = mrSges(5,1) * t91 - mrSges(5,3) * t75;
t48 = -mrSges(4,2) * t91 + mrSges(4,3) * t74;
t47 = -mrSges(5,2) * t91 + mrSges(5,3) * t74;
t46 = -pkin(5) * t134 + t62;
t41 = Ifges(4,5) * t43;
t40 = Ifges(5,5) * t43;
t39 = Ifges(4,6) * t44;
t38 = Ifges(5,6) * t44;
t36 = -t150 * t97 + t53;
t34 = -mrSges(5,1) * t74 + mrSges(5,2) * t75;
t31 = -t99 * t150 + t73 + (-pkin(5) * t97 - pkin(3)) * t100;
t30 = t62 + (-pkin(5) * t155 - t147 * t97) * qJD(1);
t29 = -pkin(3) * t44 + qJD(2) * t96;
t28 = -mrSges(4,2) * t130 + mrSges(4,3) * t44;
t27 = -mrSges(5,2) * t130 + mrSges(5,3) * t44;
t26 = mrSges(4,1) * t130 - mrSges(4,3) * t43;
t25 = mrSges(5,1) * t130 - mrSges(5,3) * t43;
t18 = qJD(1) * t106 + t45;
t15 = -qJD(3) * t53 + t151;
t14 = (-t100 * t142 - t143 * t98) * pkin(5) + t152;
t13 = -mrSges(4,1) * t44 + mrSges(4,2) * t43;
t4 = (-pkin(5) * qJD(2) - qJ(4) * qJD(3)) * t155 + (-qJD(4) * t98 + (-pkin(5) * qJD(3) - qJ(4) * qJD(2)) * t100) * t97 + t152;
t3 = -t98 * t140 + t106 * qJD(2) + (-t92 + (-t79 + t150) * t97) * qJD(3) + t151;
t7 = [t78 * t12 + t14 * t48 + t15 * t50 + t31 * t25 + t52 * t26 + t36 * t27 + t53 * t28 + t3 * t49 + t51 * t34 + t4 * t47 + m(4) * (t14 * t33 + t15 * t32 + t5 * t53 + t52 * t6) + m(5) * (t1 * t31 + t11 * t3 + t17 * t4 + t2 * t36 + t29 * t78 + t42 * t51) + (-t89 / 0.2e1 - t90 / 0.2e1 - t41 / 0.2e1 - t39 / 0.2e1 - t40 / 0.2e1 - t38 / 0.2e1 - t132 * t44 - t133 * t43 + ((-0.2e1 * t169 + 0.3e1 / 0.2e1 * Ifges(3,4) * t100) * qJD(1) + t125 * pkin(5) + t129 + t198) * qJD(2) + t104) * t100 + (t29 * t119 + pkin(5) * t13 + (-t1 * t99 - t2 * t97) * mrSges(5,3) + (-t5 * t97 - t6 * t99) * mrSges(4,3) + (-pkin(5) * t81 + t128 - t183) * qJD(2) + (t85 * t122 + t42 * t120 + (t11 * t97 - t17 * t99) * mrSges(5,3) + (t32 * t97 - t33 * t99) * mrSges(4,3) + t193 * t180 + t191 * t178 + t185 * t176 - t187 * t99 / 0.2e1) * qJD(3) + ((t161 / 0.2e1 - t159 / 0.2e1 + t162 / 0.2e1 - t160 / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(3,4)) * t98 - 0.2e1 * t170 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2) + (m(4) * pkin(5) + t121) * pkin(5) - t131) * t100) * t137 + t190 * t182 + t192 * t181 + (qJD(3) * t186 + t197) * t174 + t196 * t171) * t98; ((m(5) * t42 + t34) * t168 - t184) * qJD(3) + (t59 - t18) * t49 - m(5) * (t11 * t18 + t17 * t30 + t42 * t69) + m(5) * (t1 * t83 + t11 * t59 + t17 * t58 - t2 * t84 + t29 * t93) + (t58 - t30) * t47 + t124 * mrSges(4,3) - t29 * t120 + t93 * t12 + t83 * t25 - t84 * t27 - t69 * t34 - t46 * t48 - t45 * t50 - pkin(2) * t13 + (-t1 * t97 + t2 * t99) * mrSges(5,3) - m(4) * (t32 * t45 + t33 * t46) + ((t128 + (t170 + t167 / 0.2e1) * qJD(1) + (t81 + t145) * pkin(5) + t185 * t206 + t183) * t98 + ((t169 + (Ifges(3,2) / 0.2e1 - Ifges(3,1) / 0.2e1) * t98) * qJD(1) + ((-m(4) * pkin(2) - mrSges(3,1) - t122) * qJD(2) - t125) * pkin(5) - t94 / 0.2e1 + t129 - t198) * t100) * qJD(1) + t191 * t182 + t193 * t181 + t196 * t97 / 0.2e1 + t197 * t171 + ((-m(4) * t107 - t48 * t97 - t50 * t99) * qJD(3) + m(4) * t124 - t97 * t26 + t99 * t28) * pkin(6); (-(-t11 + t16) * t17 + (-t42 * t75 + t1) * pkin(3)) * m(5) - t104 + (-t75 * t34 + t25) * pkin(3) + (t11 * t74 + t17 * t75) * mrSges(5,3) + (t32 * t74 + t33 * t75) * mrSges(4,3) + t89 + t90 + t41 + t39 + t40 + t38 - t85 * (mrSges(4,1) * t75 + mrSges(4,2) * t74) - t42 * (mrSges(5,1) * t75 + mrSges(5,2) * t74) - t16 * t47 - t32 * t48 + t17 * t49 + t33 * t50 + (t205 * t74 - t201) * t178 + t187 * t177 + (-t194 * t75 + t195 * t74) * t176 + (-t203 * t75 + t186 + t202) * t180; -t74 * t47 + t75 * t49 + 0.2e1 * (t29 / 0.2e1 + t11 * t177 + t17 * t180) * m(5) + t12;];
tauc = t7(:);
