% Calculate matrix of centrifugal and coriolis load on the joints for
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
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RRRP6_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP6_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP6_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP6_coriolismatJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP6_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP6_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP6_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:18:11
% EndTime: 2019-12-31 17:18:16
% DurationCPUTime: 2.16s
% Computational Cost: add. (2142->276), mult. (4915->393), div. (0->0), fcn. (3804->4), ass. (0->134)
t163 = cos(qJ(3));
t241 = Ifges(4,6) + Ifges(5,6);
t246 = t241 * t163;
t242 = Ifges(4,5) + Ifges(5,5);
t162 = sin(qJ(2));
t164 = cos(qJ(2));
t156 = Ifges(5,4) * t163;
t161 = sin(qJ(3));
t175 = -Ifges(5,2) * t161 + t156;
t69 = -t164 * Ifges(5,6) + t162 * t175;
t157 = Ifges(4,4) * t163;
t176 = -Ifges(4,2) * t161 + t157;
t71 = -t164 * Ifges(4,6) + t162 * t176;
t245 = t69 + t71;
t244 = t242 * t161 + t246;
t240 = Ifges(4,3) + Ifges(5,3);
t152 = t161 * mrSges(5,1);
t153 = t163 * mrSges(5,2);
t239 = t153 + t152;
t154 = Ifges(5,5) * t163;
t155 = Ifges(4,5) * t163;
t238 = t154 + t155;
t133 = t161 * Ifges(5,1) + t156;
t134 = t161 * Ifges(4,1) + t157;
t137 = t162 * pkin(2) - pkin(6) * t164;
t201 = t161 * t162;
t61 = pkin(5) * t201 + t163 * t137;
t199 = t162 * t163;
t62 = -pkin(5) * t199 + t161 * t137;
t237 = -t161 * t61 + t163 * t62;
t236 = m(5) / 0.2e1;
t235 = -pkin(6) / 0.2e1;
t234 = -mrSges(4,1) / 0.2e1;
t233 = mrSges(4,2) / 0.2e1;
t126 = -pkin(2) * t164 - t162 * pkin(6) - pkin(1);
t198 = t163 * t164;
t59 = pkin(5) * t198 + t161 * t126;
t32 = -qJ(4) * t201 + t59;
t232 = m(5) * t32;
t29 = t162 * pkin(3) - qJ(4) * t198 + t61;
t231 = pkin(3) * t29;
t217 = -qJ(4) - pkin(6);
t127 = t217 * t161;
t230 = t127 / 0.2e1;
t129 = t217 * t163;
t229 = t129 / 0.2e1;
t228 = -t161 / 0.2e1;
t226 = t161 / 0.2e1;
t225 = t162 / 0.2e1;
t223 = t163 / 0.2e1;
t219 = pkin(3) * t163;
t148 = -pkin(2) - t219;
t221 = m(5) * t148;
t220 = pkin(3) * t161;
t216 = mrSges(4,1) * t163;
t215 = mrSges(5,3) * t163;
t214 = Ifges(4,4) * t161;
t213 = Ifges(5,4) * t161;
t113 = t163 * t126;
t181 = -qJ(4) * t199 + t113;
t28 = (-pkin(5) * t161 - pkin(3)) * t164 + t181;
t212 = t129 * t28;
t200 = t161 * t164;
t195 = pkin(5) * t200;
t31 = t181 - t195;
t211 = t129 * t31;
t206 = t164 * mrSges(5,1);
t205 = t164 * mrSges(5,2);
t130 = mrSges(4,1) * t161 + mrSges(4,2) * t163;
t100 = t164 * t130;
t193 = mrSges(5,3) * t201;
t116 = -t193 + t205;
t194 = mrSges(4,3) * t201;
t117 = mrSges(4,2) * t164 - t194;
t118 = -t162 * mrSges(5,2) - mrSges(5,3) * t200;
t119 = -t162 * mrSges(4,2) - mrSges(4,3) * t200;
t120 = -mrSges(5,3) * t199 - t206;
t121 = -mrSges(4,1) * t164 - mrSges(4,3) * t199;
t122 = t162 * mrSges(5,1) - mrSges(5,3) * t198;
t123 = t162 * mrSges(4,1) - mrSges(4,3) * t198;
t187 = pkin(5) + t220;
t124 = t187 * t162;
t125 = t187 * t164;
t188 = -Ifges(5,6) / 0.2e1 - Ifges(4,6) / 0.2e1;
t189 = -Ifges(5,5) / 0.2e1 - Ifges(4,5) / 0.2e1;
t34 = -qJ(4) * t200 + t62;
t58 = t113 - t195;
t70 = Ifges(5,6) * t162 + t164 * t175;
t72 = Ifges(4,6) * t162 + t164 * t176;
t177 = Ifges(5,1) * t163 - t213;
t73 = -t164 * Ifges(5,5) + t162 * t177;
t74 = Ifges(5,5) * t162 + t164 * t177;
t178 = Ifges(4,1) * t163 - t214;
t75 = -t164 * Ifges(4,5) + t162 * t178;
t76 = Ifges(4,5) * t162 + t164 * t178;
t98 = t239 * t162;
t99 = t239 * t164;
t3 = t34 * t116 + t62 * t117 + t32 * t118 + t59 * t119 + t29 * t120 + t61 * t121 + t28 * t122 + t58 * t123 + t124 * t99 + t125 * t98 + m(4) * (t58 * t61 + t59 * t62) + m(5) * (t124 * t125 + t28 * t29 + t32 * t34) + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t164 + (t73 / 0.2e1 + t75 / 0.2e1 + t189 * t164) * t163 + (-t69 / 0.2e1 - t71 / 0.2e1 - t188 * t164) * t161) * t164 + (-pkin(1) * mrSges(3,1) - Ifges(3,4) * t162 + pkin(5) * t100 + (t74 / 0.2e1 + t76 / 0.2e1 - t189 * t162) * t163 + (-t70 / 0.2e1 - t72 / 0.2e1 + t188 * t162) * t161 + (Ifges(3,1) - Ifges(3,2) + (m(4) * pkin(5) + t130) * pkin(5) - t240) * t164) * t162;
t204 = t3 * qJD(1);
t131 = t163 * Ifges(5,2) + t213;
t132 = t163 * Ifges(4,2) + t214;
t180 = -mrSges(4,2) * t161 + t216;
t146 = mrSges(5,1) * t199;
t182 = -mrSges(5,2) * t201 + t146;
t196 = pkin(3) * t199;
t4 = -t98 * t196 - t124 * t182 - t31 * t116 + t32 * t120 - t28 * t193 + t59 * t121 - (-t28 + t31) * t232 + (-t117 - t194) * t58 + (t32 * t215 - m(5) * t124 * t219 + t59 * t163 * mrSges(4,3) - t244 * t164 / 0.2e1 + (t73 + t75) * t226 + t245 * t223 + (-pkin(5) * t180 + (t131 + t132) * t228 - (-t134 - t133) * t163 / 0.2e1) * t162) * t162;
t203 = t4 * qJD(1);
t10 = (t161 * t116 - m(5) * (-t161 * t32 - t28 * t163) + t163 * t120) * t162;
t202 = qJD(1) * t10;
t192 = t125 * t236;
t191 = t220 / 0.2e1;
t190 = t117 * t235;
t184 = -t131 / 0.2e1 - t132 / 0.2e1;
t183 = t133 / 0.2e1 + t134 / 0.2e1;
t128 = -mrSges(5,1) * t163 + mrSges(5,2) * t161;
t159 = t161 ^ 2;
t160 = t163 ^ 2;
t169 = -Ifges(4,2) / 0.4e1 - Ifges(5,2) / 0.4e1 + Ifges(4,1) / 0.4e1 + Ifges(5,1) / 0.4e1;
t165 = pkin(5) * t130 / 0.2e1 + (-t160 / 0.2e1 - t159 / 0.2e1) * pkin(6) * mrSges(4,3) + (pkin(2) * t233 - t148 * mrSges(5,2) / 0.2e1 - t134 / 0.4e1 - t133 / 0.4e1 - t157 / 0.4e1 - t156 / 0.4e1 + mrSges(5,3) * t230 - t169 * t161) * t161 + (pkin(2) * t234 - t132 / 0.4e1 - t131 / 0.4e1 + mrSges(5,3) * t229 + t169 * t163 + (-0.3e1 / 0.4e1 * Ifges(4,4) - 0.3e1 / 0.4e1 * Ifges(5,4)) * t161 + (t221 / 0.2e1 + t128 / 0.2e1) * pkin(3)) * t163;
t166 = (t75 / 0.4e1 + t73 / 0.4e1 + t121 * t235 + (t31 / 0.2e1 - t28 / 0.2e1) * mrSges(5,3)) * t163 + t124 * t239 / 0.2e1 + t116 * t230 + t120 * t229 + t148 * t146 / 0.2e1;
t167 = -pkin(3) * t122 / 0.2e1 - t29 * mrSges(5,1) / 0.2e1 + t34 * mrSges(5,2) / 0.2e1 + t61 * t234 + t62 * t233;
t2 = (-t71 / 0.4e1 - t69 / 0.4e1 + t190 + pkin(3) * t98 / 0.2e1) * t161 + (-t231 / 0.2e1 + t124 * t191 + t212 / 0.2e1 - t211 / 0.2e1) * m(5) + (-t155 / 0.4e1 - t154 / 0.4e1 + t189 * t163 + (0.3e1 / 0.4e1 * Ifges(4,6) + 0.3e1 / 0.4e1 * Ifges(5,6)) * t161) * t164 + (-Ifges(4,3) / 0.2e1 - Ifges(5,3) / 0.2e1 + t165) * t162 + t166 + t167;
t5 = pkin(2) * t130 - t128 * t220 - t148 * t239 + (-t156 / 0.2e1 - t157 / 0.2e1 - t183) * t163 + (-pkin(3) * t221 + (Ifges(5,4) / 0.2e1 + Ifges(4,4) / 0.2e1) * t161 + (-Ifges(5,1) / 0.2e1 + Ifges(5,2) / 0.2e1 - Ifges(4,1) / 0.2e1 + Ifges(4,2) / 0.2e1) * t163 - t184) * t161;
t172 = t2 * qJD(1) - t5 * qJD(2);
t18 = m(5) * (-t127 * t161 - t129 * t163) + (t160 + t159) * mrSges(5,3);
t168 = m(5) * ((-t127 * t162 + t32) * t163 + (t129 * t162 - t28) * t161);
t7 = (t205 / 0.2e1 - t116 / 0.2e1) * t163 + (t206 / 0.2e1 + t120 / 0.2e1) * t161 + t192 - t168 / 0.2e1;
t171 = -qJD(1) * t7 + qJD(2) * t18;
t101 = -m(5) * t220 - t239;
t52 = -m(5) * t196 - t182;
t170 = qJD(1) * t52 + qJD(2) * t101;
t8 = t116 * t223 + t168 / 0.2e1 + t120 * t228 + t192 + (t153 / 0.2e1 + t152 / 0.2e1) * t164;
t1 = t161 * t190 + t165 * t162 + t98 * t191 + t166 - t167 + (t124 * t220 - t211 + t212 + t231) * t236 - t245 * t161 / 0.4e1 + t240 * t225 - (-t241 * t161 + t238) * t164 / 0.4e1 - t241 * t200 / 0.2e1 + t242 * t198 / 0.2e1;
t6 = [qJD(2) * t3 - qJD(3) * t4 - qJD(4) * t10, t1 * qJD(3) + t8 * qJD(4) + t204 + (t148 * t99 + t127 * t122 + t125 * t128 - t129 * t118 + t34 * t215 - t29 * t161 * mrSges(5,3) + m(5) * (t125 * t148 + t127 * t29 - t129 * t34) - pkin(2) * t100 + (m(4) * t237 + t163 * t119 - t161 * t123) * pkin(6) + (Ifges(3,5) + t183 * t163 + t184 * t161 + (-m(4) * pkin(2) - mrSges(3,1) - t180) * pkin(5)) * t164 + (t74 + t76) * t226 + t244 * t225 + (t70 + t72) * t223 + (pkin(5) * mrSges(3,2) - Ifges(3,6)) * t162 + t237 * mrSges(4,3)) * qJD(2), t1 * qJD(2) - t203 + (-mrSges(4,1) * t59 - mrSges(5,1) * t32 - mrSges(4,2) * t58 - mrSges(5,2) * t31 - pkin(3) * t232 + (-t246 + (mrSges(5,3) * pkin(3) - t242) * t161) * t162) * qJD(3), qJD(2) * t8 - t202; qJD(3) * t2 - qJD(4) * t7 - t204, -qJD(3) * t5 + qJD(4) * t18, t172 + (mrSges(5,1) * t129 - mrSges(5,2) * t127 - pkin(6) * t216 + (mrSges(4,2) * pkin(6) - t241) * t161 + (m(5) * t129 - t215) * pkin(3) + t238) * qJD(3), t171; -qJD(2) * t2 + qJD(4) * t52 + t203, qJD(4) * t101 - t172, 0, t170; qJD(2) * t7 - qJD(3) * t52 + t202, -qJD(3) * t101 - t171, -t170, 0;];
Cq = t6;
