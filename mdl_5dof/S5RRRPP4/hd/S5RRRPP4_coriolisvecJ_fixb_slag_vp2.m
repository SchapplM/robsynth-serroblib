% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
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
% Datum: 2019-12-31 20:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRPP4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP4_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP4_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP4_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP4_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:54:53
% EndTime: 2019-12-31 20:55:02
% DurationCPUTime: 3.85s
% Computational Cost: add. (3879->347), mult. (10577->458), div. (0->0), fcn. (7252->6), ass. (0->168)
t238 = Ifges(5,1) + Ifges(6,1);
t243 = Ifges(6,4) + Ifges(5,5);
t246 = -mrSges(6,1) - mrSges(5,1);
t171 = sin(qJ(3));
t172 = sin(qJ(2));
t173 = cos(qJ(3));
t174 = cos(qJ(2));
t145 = t171 * t174 + t173 * t172;
t135 = t145 * qJD(1);
t196 = t173 * t174;
t144 = -t171 * t172 + t196;
t134 = t144 * qJD(1);
t170 = sin(pkin(8));
t204 = cos(pkin(8));
t90 = -t204 * t134 + t135 * t170;
t231 = t90 / 0.2e1;
t245 = t90 * Ifges(5,4);
t244 = t90 * Ifges(6,5);
t242 = -Ifges(5,6) + Ifges(6,6);
t169 = qJD(2) + qJD(3);
t178 = t170 * t134 + t135 * t204;
t241 = t243 * t169 + t238 * t178 + t244 - t245;
t228 = t178 / 0.2e1;
t227 = -pkin(7) - pkin(6);
t158 = t227 * t174;
t150 = qJD(1) * t158;
t139 = t173 * t150;
t157 = t227 * t172;
t149 = qJD(1) * t157;
t142 = qJD(2) * pkin(2) + t149;
t101 = t142 * t171 - t139;
t201 = qJ(4) * t134;
t75 = t101 + t201;
t205 = t170 * t75;
t136 = t171 * t150;
t100 = t173 * t142 + t136;
t129 = t135 * qJ(4);
t74 = t100 - t129;
t68 = pkin(3) * t169 + t74;
t26 = t204 * t68 - t205;
t240 = t26 * mrSges(5,3);
t239 = mrSges(6,2) + mrSges(5,3);
t212 = -Ifges(5,4) + Ifges(6,5);
t78 = -mrSges(5,2) * t169 - mrSges(5,3) * t90;
t81 = -mrSges(6,2) * t90 + mrSges(6,3) * t169;
t211 = t78 + t81;
t237 = mrSges(6,2) * t178;
t210 = mrSges(5,3) * t178 + t169 * t246 + t237;
t236 = Ifges(5,4) * t178;
t235 = Ifges(6,5) * t178;
t105 = -t149 * t171 + t139;
t180 = t105 - t201;
t183 = t204 * t171;
t106 = t173 * t149 + t136;
t77 = -t129 + t106;
t234 = -t170 * t77 + t180 * t204 + (t170 * t173 + t183) * qJD(3) * pkin(2);
t146 = t171 * t157;
t111 = -t173 * t158 + t146;
t232 = -t90 / 0.2e1;
t226 = (pkin(1) * mrSges(3,1));
t225 = (pkin(1) * mrSges(3,2));
t110 = t173 * t157 + t158 * t171;
t179 = -qJ(4) * t145 + t110;
t87 = qJ(4) * t144 + t111;
t38 = t170 * t87 - t179 * t204;
t189 = qJD(2) * t227;
t191 = qJD(3) * t173;
t192 = qJD(3) * t171;
t59 = t135 * t189 + t142 * t191 + t150 * t192;
t108 = t169 * t145;
t97 = t108 * qJD(1);
t15 = -qJ(4) * t97 + qJD(4) * t134 + t59;
t60 = -t101 * qJD(3) + (t196 * t227 - t146) * qJD(2) * qJD(1);
t107 = t169 * t144;
t96 = t107 * qJD(1);
t176 = -t96 * qJ(4) - t135 * qJD(4) + t60;
t4 = t15 * t170 - t176 * t204;
t224 = t38 * t4;
t222 = t134 / 0.2e1;
t221 = t135 / 0.2e1;
t220 = -t169 / 0.2e1;
t165 = -pkin(2) * t174 - pkin(1);
t156 = qJD(1) * t165;
t218 = m(4) * t156;
t217 = pkin(3) * t135;
t216 = pkin(3) * t170;
t24 = -t169 * pkin(4) + qJD(5) - t26;
t215 = t24 * t90;
t69 = t204 * t75;
t27 = t170 * t68 + t69;
t214 = t27 * t178;
t57 = -t170 * t97 + t204 * t96;
t213 = t57 * mrSges(6,2);
t5 = t204 * t15 + t170 * t176;
t209 = mrSges(4,3) * t134;
t208 = Ifges(3,4) * t172;
t207 = t135 * mrSges(4,3);
t206 = t135 * Ifges(4,4);
t203 = Ifges(3,5) * qJD(2);
t202 = Ifges(3,6) * qJD(2);
t200 = qJD(2) * mrSges(3,1);
t199 = qJD(2) * mrSges(3,2);
t164 = pkin(2) * t173 + pkin(3);
t128 = pkin(2) * t183 + t170 * t164;
t195 = qJD(1) * t172;
t194 = qJD(1) * t174;
t193 = qJD(2) * t172;
t190 = t170 * t171 * pkin(2);
t167 = pkin(2) * t195;
t188 = t204 * pkin(3);
t187 = t203 / 0.2e1;
t186 = -t202 / 0.2e1;
t56 = t170 * t96 + t204 * t97;
t185 = t56 * mrSges(5,1) + t57 * mrSges(5,2);
t184 = t56 * mrSges(6,1) - t57 * mrSges(6,3);
t82 = pkin(3) * t97 + qJD(2) * t167;
t93 = pkin(2) * t193 + pkin(3) * t108;
t115 = -t144 * pkin(3) + t165;
t126 = t204 * pkin(2) * t191 - qJD(3) * t190;
t151 = t172 * t189;
t152 = t174 * t189;
t64 = t173 * t151 + t171 * t152 + t157 * t191 + t158 * t192;
t36 = pkin(4) * t178 + qJ(5) * t90 + t217;
t127 = t164 * t204 - t190;
t109 = -t134 * pkin(3) + qJD(4) + t156;
t65 = -qJD(3) * t111 - t151 * t171 + t173 * t152;
t177 = -qJ(4) * t107 - qJD(4) * t145 + t65;
t130 = Ifges(4,4) * t134;
t2 = qJD(5) * t169 + t5;
t25 = qJ(5) * t169 + t27;
t31 = t90 * pkin(4) - qJ(5) * t178 + t109;
t42 = t169 * Ifges(6,6) + t90 * Ifges(6,3) + t235;
t43 = -t90 * Ifges(5,2) + t169 * Ifges(5,6) + t236;
t85 = t134 * Ifges(4,2) + t169 * Ifges(4,6) + t206;
t86 = t135 * Ifges(4,1) + t169 * Ifges(4,5) + t130;
t175 = t246 * t4 + t2 * mrSges(6,3) - t5 * mrSges(5,2) - t135 * (Ifges(4,1) * t134 - t206) / 0.2e1 - t156 * (mrSges(4,1) * t135 + mrSges(4,2) * t134) + t243 * t57 + (Ifges(4,5) * t134 - Ifges(4,6) * t135 + t178 * t242 - t243 * t90) * t220 - (-t238 * t90 + t235 - t236 + t42) * t178 / 0.2e1 + (Ifges(6,3) * t178 - t244) * t232 + t25 * t237 + t242 * t56 + t43 * t228 + Ifges(4,5) * t96 - Ifges(4,6) * t97 - (-Ifges(4,2) * t135 + t130 + t86) * t134 / 0.2e1 - t31 * (mrSges(6,1) * t178 + mrSges(6,3) * t90) - t109 * (mrSges(5,1) * t178 - mrSges(5,2) * t90) - t90 * t240 - t59 * mrSges(4,2) + t60 * mrSges(4,1) + t101 * t207 + t100 * t209 + t85 * t221 + (-Ifges(5,2) * t178 + t241 - t245) * t231;
t166 = Ifges(3,4) * t194;
t163 = -t188 - pkin(4);
t162 = qJ(5) + t216;
t154 = mrSges(3,3) * t194 - t199;
t153 = -mrSges(3,3) * t195 + t200;
t133 = Ifges(3,1) * t195 + t166 + t203;
t132 = t202 + (t174 * Ifges(3,2) + t208) * qJD(1);
t121 = -pkin(4) - t127;
t120 = qJ(5) + t128;
t117 = qJD(5) + t126;
t114 = mrSges(4,1) * t169 - t207;
t113 = -mrSges(4,2) * t169 + t209;
t112 = t167 + t217;
t104 = t170 * t144 + t145 * t204;
t103 = -t144 * t204 + t145 * t170;
t99 = -mrSges(4,1) * t134 + mrSges(4,2) * t135;
t62 = t107 * t204 - t170 * t108;
t61 = t107 * t170 + t108 * t204;
t49 = mrSges(5,1) * t90 + mrSges(5,2) * t178;
t48 = mrSges(6,1) * t90 - mrSges(6,3) * t178;
t39 = t170 * t179 + t204 * t87;
t37 = t103 * pkin(4) - t104 * qJ(5) + t115;
t35 = t167 + t36;
t34 = t170 * t180 + t204 * t77;
t32 = -qJ(4) * t108 + qJD(4) * t144 + t64;
t29 = t204 * t74 - t205;
t28 = t170 * t74 + t69;
t9 = pkin(4) * t61 - qJ(5) * t62 - qJD(5) * t104 + t93;
t8 = t170 * t177 + t204 * t32;
t7 = t170 * t32 - t177 * t204;
t6 = pkin(4) * t56 - qJ(5) * t57 - qJD(5) * t178 + t82;
t1 = [(-t27 * mrSges(5,3) - t25 * mrSges(6,2) - t43 / 0.2e1 + t31 * mrSges(6,1) + t109 * mrSges(5,1) + t42 / 0.2e1 + Ifges(6,3) * t231 - Ifges(5,2) * t232) * t61 + (t107 * t221 + t96 * t145) * Ifges(4,1) + m(5) * (t109 * t93 + t115 * t82 - t26 * t7 + t27 * t8 + t39 * t5 + t224) + m(6) * (t2 * t39 + t24 * t7 + t25 * t8 + t31 * t9 + t37 * t6 + t224) + (-pkin(6) * t154 - t132 / 0.2e1 + t186 + (-(2 * t226) - 0.3e1 / 0.2e1 * t208 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t174) * qJD(1) + (t99 + qJD(1) * (-mrSges(4,1) * t144 + mrSges(4,2) * t145) + 0.2e1 * t218) * pkin(2)) * t193 + t210 * t7 + t211 * t8 + (mrSges(5,1) * t82 + mrSges(6,1) * t6 - mrSges(6,2) * t2 - mrSges(5,3) * t5 + t212 * t57 + (Ifges(5,2) + Ifges(6,3)) * t56) * t103 + (-pkin(6) * t153 + t133 / 0.2e1 + t187 + (-(2 * t225) + 0.3e1 / 0.2e1 * Ifges(3,4) * t174) * qJD(1)) * t174 * qJD(2) + (Ifges(4,5) * t107 - Ifges(4,6) * t108 + t242 * t61 + t243 * t62) * t169 / 0.2e1 + t165 * (mrSges(4,1) * t97 + mrSges(4,2) * t96) + (-t108 * t222 - t144 * t97) * Ifges(4,2) + (t107 * t222 - t108 * t221 + t144 * t96 - t97 * t145) * Ifges(4,4) + (-t100 * t107 - t101 * t108 - t110 * t96 - t111 * t97 + t144 * t59 - t145 * t60) * mrSges(4,3) + t37 * t184 + t115 * t185 + t241 * t62 / 0.2e1 + m(4) * (t100 * t65 + t101 * t64 + t110 * t60 + t111 * t59) + (t212 * t61 + t238 * t62) * t228 + t239 * (t38 * t57 - t39 * t56) + (mrSges(5,2) * t82 - mrSges(6,3) * t6 + t212 * t56 + t238 * t57 + t239 * t4) * t104 + (t109 * mrSges(5,2) + t24 * mrSges(6,2) - t31 * mrSges(6,3) + Ifges(5,4) * t232 + Ifges(6,5) * t231 - t240) * t62 + t93 * t49 + t107 * t86 / 0.2e1 - t108 * t85 / 0.2e1 + t64 * t113 + t65 * t114 + t9 * t48 + t156 * (mrSges(4,1) * t108 + mrSges(4,2) * t107); ((t113 * t173 - t114 * t171) * qJD(3) + (-t171 * t97 - t173 * t96) * mrSges(4,3)) * pkin(2) + t175 - t211 * t34 + (-t127 * t57 - t128 * t56 + t214) * mrSges(5,3) + (-t120 * t56 + t121 * t57 + t215) * mrSges(6,2) + ((t187 - t133 / 0.2e1 - t166 / 0.2e1 + qJD(1) * t225 + (t153 - t200) * pkin(6)) * t174 + (t186 + t132 / 0.2e1 + (t226 + t208 / 0.2e1 + (Ifges(3,2) / 0.2e1 - Ifges(3,1) / 0.2e1) * t174) * qJD(1) + (t154 + t199) * pkin(6) + (-t99 - t218) * pkin(2)) * t172) * qJD(1) + t126 * t78 - t112 * t49 - t106 * t113 - t105 * t114 + t117 * t81 - t35 * t48 + t234 * t210 + (t2 * t120 + t4 * t121 - t31 * t35 + (-t34 + t117) * t25 + t234 * t24) * m(6) + (-t109 * t112 - t127 * t4 + t128 * t5 + (t126 - t34) * t27 - t234 * t26) * m(5) + (-t100 * t105 - t101 * t106 + (t171 * t59 + t173 * t60 + (-t100 * t171 + t101 * t173) * qJD(3)) * pkin(2)) * m(4); qJD(5) * t81 - t100 * t113 + t101 * t114 + t163 * t213 - t49 * t217 - t36 * t48 + t175 - t211 * t29 - t210 * t28 + (-t162 * t56 + t215) * mrSges(6,2) + (t162 * t2 + t163 * t4 - t24 * t28 - t31 * t36 + (qJD(5) - t29) * t25) * m(6) + (-t109 * t217 + t26 * t28 - t27 * t29 + (t170 * t5 - t204 * t4) * pkin(3)) * m(5) + (-t188 * t57 - t216 * t56 + t214) * mrSges(5,3); -t210 * t178 + t211 * t90 + t184 + t185 + (-t178 * t24 + t25 * t90 + t6) * m(6) + (t178 * t26 + t27 * t90 + t82) * m(5); t213 - t169 * t81 + t178 * t48 + 0.2e1 * (t4 / 0.2e1 + t25 * t220 + t31 * t228) * m(6);];
tauc = t1(:);
