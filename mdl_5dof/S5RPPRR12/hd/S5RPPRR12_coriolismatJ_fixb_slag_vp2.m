% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-31 18:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPPRR12_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR12_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR12_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR12_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR12_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR12_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR12_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:06:59
% EndTime: 2019-12-31 18:07:02
% DurationCPUTime: 1.42s
% Computational Cost: add. (4155->220), mult. (7815->307), div. (0->0), fcn. (8034->6), ass. (0->134)
t146 = sin(qJ(5));
t141 = t146 ^ 2;
t148 = cos(qJ(5));
t142 = t148 ^ 2;
t245 = t141 + t142;
t246 = -t245 + 0.1e1;
t179 = mrSges(6,3) * (t142 / 0.2e1 + t141 / 0.2e1);
t225 = Ifges(6,4) * t146;
t131 = Ifges(6,2) * t148 + t225;
t138 = Ifges(6,4) * t148;
t132 = Ifges(6,1) * t146 + t138;
t232 = -t148 / 0.2e1;
t233 = t146 / 0.2e1;
t244 = t131 * t233 + t132 * t232;
t137 = Ifges(6,5) * t148;
t243 = -t146 * Ifges(6,6) + t137;
t223 = Ifges(6,2) * t146;
t242 = t223 - t138;
t241 = mrSges(6,3) * t245;
t143 = sin(pkin(8));
t144 = cos(pkin(8));
t147 = sin(qJ(4));
t230 = cos(qJ(4));
t123 = t143 * t147 - t230 * t144;
t101 = t123 ^ 2;
t124 = t230 * t143 + t147 * t144;
t240 = t124 ^ 2;
t239 = -m(6) / 0.2e1;
t238 = m(6) / 0.2e1;
t205 = t123 * t148;
t89 = t124 * mrSges(6,1) + mrSges(6,3) * t205;
t237 = -t89 / 0.2e1;
t236 = -t123 / 0.2e1;
t235 = -t124 / 0.2e1;
t212 = t148 * mrSges(6,2);
t218 = t146 * mrSges(6,1);
t130 = t212 + t218;
t234 = t130 / 0.2e1;
t231 = t148 / 0.2e1;
t229 = pkin(4) * t123;
t228 = pkin(4) * t124;
t227 = pkin(7) * t124;
t145 = -pkin(1) - qJ(3);
t226 = -pkin(6) + t145;
t224 = Ifges(6,5) * t124;
t222 = Ifges(6,6) * t124;
t197 = t143 * pkin(3) + qJ(2);
t87 = pkin(7) * t123 + t197 + t228;
t128 = t226 * t143;
t184 = t226 * t144;
t93 = t230 * t128 + t147 * t184;
t37 = -t146 * t93 + t148 * t87;
t38 = t146 * t87 + t148 * t93;
t176 = t37 * t146 - t38 * t148;
t196 = t143 ^ 2 + t144 ^ 2;
t185 = m(4) * t196;
t202 = t146 * t123;
t194 = mrSges(6,3) * t202;
t88 = -t124 * mrSges(6,2) + t194;
t211 = t148 * t88;
t214 = t146 * t89;
t92 = t128 * t147 - t230 * t184;
t220 = t123 * t92;
t7 = t196 * mrSges(4,3) + (mrSges(5,3) * t124 - t211 + t214) * t124 + (mrSges(5,3) + t130) * t101 + m(6) * (t176 * t124 - t220) + m(5) * (-t124 * t93 - t220) - t145 * t185;
t221 = qJD(1) * t7;
t217 = t146 * mrSges(6,2);
t215 = t146 * t88;
t213 = t148 * mrSges(6,1);
t210 = t148 * t89;
t129 = -t213 + t217;
t2 = -t129 * t220 + t38 * t89 + ((-Ifges(6,4) * t202 - t224) * t146 + (-t38 * mrSges(6,3) - t222 + (t138 + (Ifges(6,1) - Ifges(6,2)) * t146) * t123) * t148) * t123 + (-t88 + t194) * t37;
t209 = t2 * qJD(1);
t163 = t101 * t129;
t166 = t213 / 0.2e1 - t217 / 0.2e1;
t192 = -t215 / 0.2e1;
t168 = t192 - t210 / 0.2e1;
t9 = -t163 / 0.2e1 + (-t123 * t179 - t168) * t124 + t166;
t208 = t9 * qJD(1);
t182 = -t124 * mrSges(5,1) + t123 * mrSges(5,2);
t15 = t143 * mrSges(4,1) + t144 * mrSges(4,2) + t215 + t210 + mrSges(3,3) + (m(3) + m(4)) * qJ(2) + m(6) * (t38 * t146 + t37 * t148) + m(5) * t197 - t182;
t207 = qJD(1) * t15;
t201 = t146 * t124;
t169 = t123 * mrSges(6,2) + mrSges(6,3) * t201;
t200 = t148 * t124;
t170 = -t123 * mrSges(6,1) + mrSges(6,3) * t200;
t94 = t227 - t229;
t151 = t238 * t245 * t94 + t169 * t233 + t170 * t231;
t154 = (-t227 * t245 + t229) * t238 + t129 * t236;
t183 = -t123 * mrSges(5,1) - t124 * mrSges(5,2);
t11 = t235 * t241 - t151 + t154 - t183;
t206 = t11 * qJD(1);
t204 = t124 * t123;
t203 = t124 * t129;
t165 = t212 / 0.2e1 + t218 / 0.2e1;
t158 = t165 * t124;
t167 = t214 / 0.2e1 - t211 / 0.2e1;
t16 = t158 + t167;
t199 = t16 * qJD(1);
t152 = (-t240 * t245 - t101) * t238 + m(5) * (-t101 - t240) / 0.2e1 - t185 / 0.2e1;
t161 = -m(4) / 0.2e1 - m(5) / 0.2e1 + t245 * t239;
t21 = t152 + t161;
t198 = t21 * qJD(1);
t191 = t92 * t234;
t189 = -t200 / 0.2e1;
t181 = t245 * t123;
t178 = Ifges(6,1) * t148 - t225;
t162 = t124 * t130;
t171 = -Ifges(5,4) + t243;
t174 = t146 * t92 + t148 * t94;
t175 = t146 * t94 - t148 * t92;
t1 = t92 * t162 - m(6) * (t37 * t174 + t38 * t175 + t92 * t93) - t175 * t88 - t38 * t169 - t174 * t89 - t37 * t170 - t197 * t183 + t171 * t240 + ((-t142 * Ifges(6,1) - Ifges(5,1) + Ifges(5,2) + Ifges(6,3) + (-t223 + 0.2e1 * t138) * t146) * t124 + t130 * t93 - t171 * t123) * t123;
t159 = t148 * t169;
t160 = t146 * t170;
t6 = t123 * t162 + ((t176 + t93) * t123 + t246 * t124 * t92) * t239 + t88 * t205 / 0.2e1 + t159 * t235 + t202 * t237 + t124 * t160 / 0.2e1;
t177 = -t1 * qJD(1) - t6 * qJD(2);
t24 = m(6) * t246 * t204;
t173 = -t6 * qJD(1) + t24 * qJD(2);
t172 = pkin(7) * t179;
t164 = t123 * (Ifges(6,5) * t146 + Ifges(6,6) * t148);
t153 = t175 * mrSges(6,2) / 0.2e1 - t174 * mrSges(6,1) / 0.2e1;
t155 = t138 / 0.4e1 + t132 / 0.4e1 - pkin(4) * mrSges(6,2) / 0.2e1 + (-Ifges(6,2) / 0.2e1 + Ifges(6,1) / 0.4e1) * t146;
t156 = t131 / 0.4e1 + pkin(4) * mrSges(6,1) / 0.2e1 + (-Ifges(6,1) / 0.2e1 + Ifges(6,2) / 0.4e1) * t148;
t3 = t191 + t168 * pkin(7) + t243 * t124 + (Ifges(6,3) / 0.2e1 + t172 + t156 * t148 + (0.5e1 / 0.4e1 * t138 + t155) * t146) * t123 + t153;
t52 = pkin(4) * t130 - t242 * t232 - t146 * t178 / 0.2e1 + t244;
t53 = (-t130 / 0.2e1 + t165) * t123;
t157 = t3 * qJD(1) - t53 * qJD(2) - t52 * qJD(4);
t54 = (t165 + t234) * t123;
t20 = t152 - t161;
t17 = t158 - t167;
t13 = -t124 * t179 + t151 + t154;
t10 = t163 / 0.2e1 + t124 * t192 + t89 * t189 + t166 + t204 * t241 / 0.2e1;
t5 = t6 * qJD(4);
t4 = t191 + t124 * t137 / 0.4e1 + Ifges(6,6) * t201 / 0.2e1 + Ifges(6,5) * t189 + Ifges(6,3) * t236 + t123 * t172 + (-pkin(7) * t88 / 0.2e1 - t222 / 0.2e1 + t155 * t123) * t146 + (pkin(7) * t237 + t224 / 0.4e1 + (0.5e1 / 0.4e1 * t225 + t156) * t123) * t148 - t153;
t8 = [qJD(2) * t15 + qJD(3) * t7 - qJD(4) * t1 - qJD(5) * t2, qJD(3) * t20 + qJD(5) * t10 + t207 - t5, qJD(2) * t20 + qJD(4) * t13 + qJD(5) * t17 + t221, t13 * qJD(3) + t4 * qJD(5) + t177 + (pkin(4) * t162 - t164 / 0.2e1 + t92 * mrSges(5,2) + (-Ifges(6,5) * t233 - Ifges(6,6) * t231 + Ifges(5,6)) * t123 + (-t178 * t233 + t242 * t231 - Ifges(5,5) + t244) * t124 + (-m(6) * pkin(4) - mrSges(5,1) + t129) * t93 + (-m(6) * t245 * t92 + t159 - t160) * pkin(7) + (-t174 * t146 + t175 * t148) * mrSges(6,3)) * qJD(4), -t209 + t10 * qJD(2) + t17 * qJD(3) + t4 * qJD(4) + (-mrSges(6,1) * t38 - mrSges(6,2) * t37 + t164) * qJD(5); qJD(3) * t21 - qJD(5) * t9 - t207 - t5, t24 * qJD(4), t198, (t203 + m(6) * (-pkin(7) * t181 - t228) - mrSges(6,3) * t181 + t182) * qJD(4) + t54 * qJD(5) + t173, t54 * qJD(4) + qJD(5) * t203 - t208; -qJD(2) * t21 - qJD(4) * t11 - qJD(5) * t16 - t221, -t198, 0, -t206, -qJD(5) * t130 - t199; qJD(3) * t11 + qJD(5) * t3 - t177, -qJD(5) * t53 - t173, t206, -t52 * qJD(5), (t129 * pkin(7) + t243) * qJD(5) + t157; qJD(2) * t9 + qJD(3) * t16 - qJD(4) * t3 + t209, qJD(4) * t53 + t208, t199, -t157, 0;];
Cq = t8;
