% Calculate matrix of centrifugal and coriolis load on the joints for
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPPRR4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR4_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR4_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR4_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR4_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:44:10
% EndTime: 2019-12-05 17:44:18
% DurationCPUTime: 2.36s
% Computational Cost: add. (7285->232), mult. (15891->328), div. (0->0), fcn. (17125->8), ass. (0->130)
t148 = sin(pkin(9));
t150 = cos(pkin(9));
t153 = sin(qJ(4));
t155 = cos(qJ(4));
t132 = t148 * t155 + t153 * t150;
t149 = sin(pkin(8));
t123 = t132 * t149;
t131 = -t153 * t148 + t150 * t155;
t125 = t131 * t149;
t152 = sin(qJ(5));
t154 = cos(qJ(5));
t96 = -t123 * t152 + t125 * t154;
t200 = t96 * mrSges(6,1);
t236 = t200 / 0.2e1;
t112 = t131 * t152 + t132 * t154;
t166 = t154 * t131 - t132 * t152;
t24 = -t112 * mrSges(6,1) - t166 * mrSges(6,2);
t235 = t24 * qJD(5);
t151 = cos(pkin(8));
t124 = t132 * t151;
t126 = t131 * t151;
t220 = m(6) * pkin(4);
t172 = t220 / 0.2e1;
t94 = -t124 * t154 - t126 * t152;
t97 = -t124 * t152 + t126 * t154;
t197 = t94 * mrSges(6,1) / 0.2e1 - t97 * mrSges(6,2) / 0.2e1;
t234 = t124 * mrSges(5,1) / 0.2e1 + t126 * mrSges(5,2) / 0.2e1 - (t152 * t97 + t154 * t94) * t172 - t197;
t167 = -t154 * t123 - t125 * t152;
t228 = Ifges(6,5) * t167;
t232 = Ifges(6,6) * t96;
t196 = t228 - t232;
t164 = -t228 / 0.2e1 + t232 / 0.2e1;
t214 = -t112 / 0.2e1;
t216 = -t166 / 0.2e1;
t74 = -mrSges(6,1) * t151 - t96 * mrSges(6,3);
t233 = t74 * t214 + (t167 * t216 + t214 * t96) * mrSges(6,3);
t88 = t167 * mrSges(6,2);
t170 = t88 / 0.2e1;
t229 = t131 / 0.2e1;
t137 = -pkin(2) * t151 - qJ(3) * t149 - pkin(1);
t130 = t150 * t137;
t179 = t149 * t150;
t107 = -pkin(6) * t179 + t130 + (-qJ(2) * t148 - pkin(3)) * t151;
t182 = qJ(2) * t151;
t117 = t148 * t137 + t150 * t182;
t180 = t148 * t149;
t113 = -pkin(6) * t180 + t117;
t64 = t153 * t107 + t113 * t155;
t53 = -t123 * pkin(7) + t64;
t187 = t154 * t53;
t63 = t155 * t107 - t113 * t153;
t52 = -pkin(7) * t125 + t63;
t49 = -pkin(4) * t151 + t52;
t29 = t152 * t49 + t187;
t225 = -t29 * mrSges(6,3) - Ifges(6,4) * t96;
t224 = t148 ^ 2;
t223 = t150 ^ 2;
t222 = m(5) / 0.2e1;
t221 = m(6) / 0.2e1;
t73 = mrSges(6,2) * t151 + mrSges(6,3) * t167;
t219 = t73 / 0.2e1;
t218 = -t167 / 0.2e1;
t215 = t166 / 0.2e1;
t213 = -t132 / 0.2e1;
t212 = t149 / 0.2e1;
t211 = -t151 / 0.2e1;
t210 = -t152 / 0.2e1;
t207 = t125 * pkin(4);
t189 = t152 * t53;
t28 = t154 * t49 - t189;
t206 = t28 * mrSges(6,2);
t205 = t28 * mrSges(6,3);
t204 = t29 * mrSges(6,1);
t32 = -t152 * t52 - t187;
t202 = t32 * mrSges(6,1);
t33 = t154 * t52 - t189;
t201 = t33 * mrSges(6,2);
t173 = pkin(3) * t180 + t149 * qJ(2);
t108 = pkin(4) * t123 + t173;
t114 = mrSges(5,2) * t151 - t123 * mrSges(5,3);
t190 = t125 * mrSges(5,3);
t115 = -mrSges(5,1) * t151 - t190;
t162 = -Ifges(5,6) * t125 + t196;
t168 = t125 * mrSges(5,1) - t123 * mrSges(5,2);
t169 = t88 + t200;
t50 = -mrSges(6,1) * t167 + mrSges(6,2) * t96;
t91 = Ifges(6,4) * t167;
t1 = -t50 * t207 - t108 * t169 - t33 * t73 - t32 * t74 - m(6) * (t108 * t207 + t28 * t32 + t29 * t33) - t173 * t168 + t125 ^ 2 * Ifges(5,4) - t63 * t114 + t162 * t151 - (t63 * mrSges(5,3) + Ifges(5,4) * t123 + Ifges(5,5) * t151 + (-Ifges(5,1) + Ifges(5,2)) * t125) * t123 + (-Ifges(6,1) * t96 - t91 / 0.2e1 + t205 + t218 * Ifges(6,4)) * t167 - (-Ifges(6,2) * t167 + t225) * t96 + (t115 + t190) * t64;
t195 = t1 * qJD(1);
t188 = t152 * t167;
t186 = t154 * t96;
t2 = t196 * t211 + t28 * t73 - t29 * t74 + (t108 * mrSges(6,2) + Ifges(6,5) * t211 - t205 + t91) * t167 + (t108 * mrSges(6,1) + Ifges(6,6) * t151 / 0.2e1 + (Ifges(6,1) - Ifges(6,2)) * t167 + t225) * t96;
t185 = t2 * qJD(1);
t159 = t166 * t219 + t233;
t7 = t159 - t197;
t184 = t7 * qJD(1);
t116 = -t148 * t182 + t130;
t134 = mrSges(4,2) * t151 - mrSges(4,3) * t180;
t135 = -mrSges(4,1) * t151 - mrSges(4,3) * t179;
t146 = t149 ^ 2;
t144 = t146 * qJ(2);
t147 = t151 ^ 2;
t9 = t126 * t114 - t124 * t115 + t97 * t73 + t94 * t74 + (t150 * t134 - t148 * t135) * t151 + (t147 + t146) * mrSges(3,3) + (t123 * mrSges(5,1) + t125 * mrSges(5,2) + t50 + (mrSges(4,1) * t148 + mrSges(4,2) * t150) * t149) * t149 + m(6) * (t108 * t149 + t28 * t94 + t29 * t97) + m(5) * (-t63 * t124 + t64 * t126 + t149 * t173) + m(4) * (t144 + (-t116 * t148 + t117 * t150) * t151) + m(3) * (qJ(2) * t147 + t144);
t183 = t9 * qJD(1);
t10 = -t123 * t114 - t125 * t115 + t167 * t73 - t96 * t74 + m(6) * (t167 * t29 - t28 * t96) + m(5) * (-t123 * t64 - t125 * t63) + (m(4) * (-t116 * t150 - t117 * t148) - t148 * t134 - t150 * t135) * t149;
t181 = qJD(1) * t10;
t158 = (t112 * t167 - t166 * t96) * t221 + (-t123 * t132 - t125 * t131) * t222;
t18 = (-m(6) / 0.2e1 - m(5) / 0.2e1 + (-t224 / 0.2e1 - t223 / 0.2e1 - 0.1e1 / 0.2e1) * m(4)) * t149 + t158;
t176 = t18 * qJD(1);
t19 = (t188 / 0.2e1 - t186 / 0.2e1 - t125 / 0.2e1) * t220 - t168 - t169;
t175 = t19 * qJD(1);
t34 = 0.2e1 * t170 + 0.2e1 * t236;
t174 = t34 * qJD(1);
t171 = m(4) * t212;
t156 = (t123 * t229 + t125 * t213) * mrSges(5,3) + ((t29 + t32) * t166 + (-t28 + t33) * t112) * t221 + t73 * t215 + t114 * t229 + t115 * t213 + t233;
t3 = t156 + t234;
t165 = t3 * qJD(1);
t136 = (mrSges(6,1) * t152 + mrSges(6,2) * t154) * pkin(4);
t23 = (t216 + t215) * mrSges(6,2) + (t214 + t112 / 0.2e1) * mrSges(6,1);
t157 = (t154 * t219 + t74 * t210 + (t154 * t218 + t210 * t96) * mrSges(6,3)) * pkin(4) - t164;
t6 = (-t28 / 0.2e1 + t33 / 0.2e1) * mrSges(6,2) + (-t29 / 0.2e1 - t32 / 0.2e1) * mrSges(6,1) + t157 + t164;
t161 = -qJD(1) * t6 - qJD(2) * t23 + qJD(4) * t136;
t133 = t136 * qJD(5);
t43 = (t125 - t186 + t188) * t172;
t35 = t170 + t236 - t88 / 0.2e1 - t200 / 0.2e1;
t17 = t171 + (-t223 - t224) * t171 + t158 + (m(6) + m(5)) * t212;
t8 = t159 + t197;
t5 = -t206 / 0.2e1 - t204 / 0.2e1 - t201 / 0.2e1 + t202 / 0.2e1 + t157 - t164;
t4 = t156 - t234;
t11 = [qJD(2) * t9 + qJD(3) * t10 - qJD(4) * t1 + qJD(5) * t2, t183 + 0.2e1 * ((t112 * t97 + t166 * t94) * t221 + (-t124 * t131 + t126 * t132) * t222) * qJD(2) + t17 * qJD(3) + t4 * qJD(4) + t8 * qJD(5), qJD(2) * t17 + qJD(4) * t43 + qJD(5) * t35 + t181, t4 * qJD(2) + t43 * qJD(3) + t5 * qJD(5) - t195 + (-t64 * mrSges(5,1) - t63 * mrSges(5,2) - Ifges(5,5) * t123 + t162 - t201 + t202 + (m(6) * (t152 * t33 + t154 * t32) + (-t152 * t96 - t154 * t167) * mrSges(6,3)) * pkin(4)) * qJD(4), t185 + t8 * qJD(2) + t35 * qJD(3) + t5 * qJD(4) + (t196 - t204 - t206) * qJD(5); qJD(3) * t18 + qJD(4) * t3 + qJD(5) * t7 - t183, 0, t176, (-t132 * mrSges(5,1) - t131 * mrSges(5,2) + (-t112 * t154 + t152 * t166) * t220 + t24) * qJD(4) + t235 + t165, t24 * qJD(4) + t184 + t235; -qJD(2) * t18 - qJD(4) * t19 + qJD(5) * t34 - t181, -t176, 0, -t175, t174; -qJD(2) * t3 + qJD(3) * t19 + qJD(5) * t6 + t195, qJD(5) * t23 - t165, t175, -t133, -t133 - t161; -qJD(2) * t7 - qJD(3) * t34 - qJD(4) * t6 - t185, -t23 * qJD(4) - t184, -t174, t161, 0;];
Cq = t11;
