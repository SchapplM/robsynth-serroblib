% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRPRR4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR4_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR4_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR4_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR4_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR4_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR4_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:49:51
% EndTime: 2019-12-05 15:49:57
% DurationCPUTime: 1.63s
% Computational Cost: add. (3190->244), mult. (8478->385), div. (0->0), fcn. (8510->10), ass. (0->160)
t155 = sin(qJ(4));
t260 = m(6) * t155;
t157 = cos(qJ(5));
t146 = Ifges(6,5) * t157;
t154 = sin(qJ(5));
t233 = Ifges(6,6) * t154;
t259 = Ifges(5,4) - t146 / 0.2e1 + t233 / 0.2e1;
t147 = Ifges(6,4) * t157;
t184 = -Ifges(6,2) * t154 + t147;
t131 = Ifges(6,1) * t154 + t147;
t158 = cos(qJ(4));
t239 = t155 * pkin(4);
t135 = -pkin(8) * t158 + t239;
t153 = sin(pkin(10));
t143 = pkin(2) * t153 + pkin(7);
t210 = t154 * t155;
t78 = t135 * t157 + t143 * t210;
t209 = t155 * t157;
t79 = t135 * t154 - t143 * t209;
t258 = -t154 * t78 + t157 * t79;
t128 = -mrSges(6,1) * t157 + mrSges(6,2) * t154;
t257 = -m(6) * pkin(4) - mrSges(5,1) + t128;
t256 = -m(6) / 0.2e1;
t255 = m(6) / 0.2e1;
t219 = sin(pkin(5));
t220 = cos(pkin(10));
t187 = t220 * t219;
t156 = sin(qJ(2));
t193 = t156 * t219;
t240 = cos(qJ(2));
t101 = t153 * t193 - t240 * t187;
t188 = t240 * t219;
t102 = t153 * t188 + t156 * t187;
t221 = cos(pkin(5));
t72 = t102 * t158 + t221 * t155;
t43 = t101 * t154 + t157 * t72;
t254 = -t43 / 0.2e1;
t71 = t102 * t155 - t221 * t158;
t253 = t71 / 0.2e1;
t205 = t157 * t158;
t49 = -t101 * t205 + t102 * t154;
t222 = t49 * t157;
t204 = t158 * t154;
t48 = t101 * t204 + t102 * t157;
t223 = t48 * t154;
t182 = t222 - t223;
t215 = t101 * t158;
t252 = (t182 + t215) * t260;
t149 = t154 ^ 2;
t151 = t157 ^ 2;
t202 = t149 + t151;
t201 = -0.1e1 + t202;
t251 = t201 * t158 * t260;
t114 = t155 * t128;
t250 = t114 / 0.2e1;
t198 = mrSges(6,3) * t210;
t124 = mrSges(6,2) * t158 - t198;
t249 = -t124 / 0.2e1;
t248 = t143 / 0.2e1;
t247 = -t154 / 0.2e1;
t246 = t154 / 0.2e1;
t245 = t155 / 0.2e1;
t244 = -t157 / 0.2e1;
t243 = t157 / 0.2e1;
t242 = -t158 / 0.2e1;
t241 = -t158 / 0.4e1;
t228 = t157 * mrSges(6,2);
t231 = t154 * mrSges(6,1);
t129 = t228 + t231;
t115 = t129 * t155;
t150 = t155 ^ 2;
t152 = t158 ^ 2;
t116 = t158 * t129;
t207 = t157 * t124;
t197 = mrSges(6,3) * t209;
t126 = -mrSges(6,1) * t158 - t197;
t212 = t154 * t126;
t166 = -t116 / 0.2e1 - t212 / 0.2e1 + t207 / 0.2e1;
t174 = m(6) * t258;
t144 = -t220 * pkin(2) - pkin(3);
t117 = -t158 * pkin(4) - t155 * pkin(8) + t144;
t68 = t154 * t117 + t143 * t205;
t227 = t157 * t68;
t67 = t157 * t117 - t143 * t204;
t181 = -t154 * t67 + t227;
t125 = -t155 * mrSges(6,2) - mrSges(6,3) * t204;
t206 = t157 * t125;
t127 = t155 * mrSges(6,1) - mrSges(6,3) * t205;
t211 = t154 * t127;
t18 = m(6) * (t150 - t152) * t248 + (t181 * t255 + t166) * t158 + (t174 / 0.2e1 + t206 / 0.2e1 - t211 / 0.2e1 + t115 / 0.2e1) * t155;
t173 = t124 * t247 + t126 * t244;
t189 = mrSges(6,3) * (-t151 / 0.2e1 - t149 / 0.2e1);
t23 = t150 * t189 + t173 * t155 + t158 * t250;
t238 = t18 * qJD(4) + t23 * qJD(5);
t236 = Ifges(6,4) * t154;
t235 = Ifges(6,5) * t158;
t232 = Ifges(6,6) * t158;
t229 = t155 * t71;
t42 = t101 * t157 - t154 * t72;
t225 = t42 * t154;
t224 = t43 * t157;
t218 = t101 * t150;
t217 = t101 * t152;
t216 = t101 * t155;
t214 = t143 * t155;
t103 = t184 * t155 - t232;
t213 = t154 * t103;
t185 = Ifges(6,1) * t157 - t236;
t105 = t185 * t155 - t235;
t208 = t157 * t105;
t203 = t23 * qJD(2);
t200 = -t252 / 0.2e1;
t199 = t252 / 0.2e1;
t196 = pkin(4) * t250;
t195 = Ifges(6,2) / 0.4e1 - Ifges(6,1) / 0.4e1;
t130 = Ifges(6,2) * t157 + t236;
t194 = t130 * t247;
t191 = t146 - t233;
t190 = t202 * t158;
t186 = t155 * mrSges(5,1) + t158 * mrSges(5,2);
t183 = Ifges(6,5) * t154 + Ifges(6,6) * t157;
t177 = -t224 + t72 + t225;
t12 = m(6) * t177 * t71;
t16 = (-t177 * t158 - t201 * t229) * t255;
t180 = t12 * qJD(1) + t16 * qJD(3);
t179 = t48 * mrSges(6,1) / 0.2e1 - t49 * mrSges(6,2) / 0.2e1;
t178 = -t78 * mrSges(6,1) / 0.2e1 + t79 * mrSges(6,2) / 0.2e1;
t175 = t228 / 0.2e1 + t231 / 0.2e1;
t172 = t131 * t243 + t194;
t7 = m(5) * (-t158 * t72 + t102 - t229) * t101 + (-t229 * t101 + t42 * t48 + t43 * t49) * m(6);
t171 = -t7 * qJD(1) + qJD(3) * t200;
t160 = (t222 / 0.2e1 - t223 / 0.2e1) * mrSges(6,3) + (pkin(4) * t216 + t182 * pkin(8)) * t255 - t101 * t114 / 0.2e1;
t162 = -t42 * t127 / 0.2e1 + t125 * t254 - t72 * t115 / 0.2e1;
t165 = t72 * t214 + t78 * t42 + t79 * t43;
t169 = t143 * t158 - t181;
t1 = t165 * t256 + (t169 * t256 + t166) * t71 + t160 + t162;
t104 = Ifges(6,6) * t155 + t184 * t158;
t106 = Ifges(6,5) * t155 + t185 * t158;
t6 = m(6) * (t67 * t78 + t68 * t79) + t79 * t124 + t68 * t125 + t78 * t126 + t67 * t127 + (t208 / 0.2e1 - t213 / 0.2e1 + t143 * t115 + t144 * mrSges(5,2) + t259 * t158) * t158 + (t106 * t243 + t104 * t247 + t143 * t116 + t144 * mrSges(5,1) - t259 * t155 + (m(6) * t143 ^ 2 + Ifges(5,1) - Ifges(5,2) - Ifges(6,3)) * t158) * t155;
t170 = -t1 * qJD(1) + t6 * qJD(2) + t18 * qJD(3);
t168 = -t129 / 0.2e1 + t175;
t10 = t114 * t214 + t68 * t126 + (t105 * t246 + t103 * t243 + mrSges(6,3) * t227 + t183 * t242 + (-t131 * t244 + t194) * t155) * t155 + (-t124 - t198) * t67;
t163 = t42 * t249 + t43 * t126 / 0.2e1 + t71 * t250;
t3 = (t224 / 0.2e1 - t225 / 0.2e1) * t155 * mrSges(6,3) + t163 + t179;
t167 = -t3 * qJD(1) - t10 * qJD(2) + t23 * qJD(3);
t164 = t16 * qJD(1) + t18 * qJD(2) + qJD(3) * t251;
t13 = t168 * t71;
t37 = pkin(4) * t129 + t184 * t244 + t185 * t247 - t172;
t69 = t168 * t158;
t159 = t129 * t248 + pkin(8) * t189 + (-t131 / 0.4e1 - t147 / 0.4e1 + t195 * t154) * t154 + (-0.3e1 / 0.4e1 * t236 - t130 / 0.4e1 - t195 * t157) * t157;
t9 = t196 + t146 * t241 + (-pkin(8) * t126 / 0.2e1 + t105 / 0.4e1 - t235 / 0.2e1) * t157 + (pkin(8) * t249 - t103 / 0.4e1 + 0.3e1 / 0.4e1 * t232) * t154 + (-Ifges(6,3) / 0.2e1 + t159) * t155 + t178;
t161 = t13 * qJD(1) - t9 * qJD(2) - t69 * qJD(3) + t37 * qJD(4);
t73 = t143 * t218;
t70 = t129 * t242 - t175 * t158;
t14 = t129 * t253 + t175 * t71;
t8 = t208 / 0.4e1 - t213 / 0.4e1 + t196 + t191 * t241 + Ifges(6,5) * t205 / 0.2e1 - Ifges(6,6) * t204 / 0.2e1 + Ifges(6,3) * t245 + t173 * pkin(8) + t159 * t155 - t178;
t5 = qJD(2) * t199 + t16 * qJD(4);
t4 = t197 * t254 + t42 * t198 / 0.2e1 - t163 + t179;
t2 = -t71 * t207 / 0.2e1 + (t169 * t71 + t165) * t255 + t101 * t186 / 0.2e1 + mrSges(5,2) * t215 / 0.2e1 + mrSges(5,1) * t216 / 0.2e1 + t160 - t162 + (t212 + t116) * t253;
t11 = [t7 * qJD(2) + t12 * qJD(4), (-mrSges(3,2) * t188 - mrSges(3,1) * t193 + t49 * t124 + t48 * t126 + m(6) * (t67 * t48 + t68 * t49 - t73) - t115 * t216 + m(5) * (t102 * t144 - t143 * t217 - t73) + t102 * (-t158 * mrSges(5,1) + t155 * mrSges(5,2)) + m(4) * (-t101 * t153 - t220 * t102) * pkin(2) + t101 * mrSges(4,2) - t102 * mrSges(4,1) + (-t217 - t218) * mrSges(5,3)) * qJD(2) + t2 * qJD(4) + t4 * qJD(5) - t171, t5, t2 * qJD(2) + (t257 * t72 + (mrSges(5,2) - (m(6) * pkin(8) + mrSges(6,3)) * t202) * t71) * qJD(4) + t14 * qJD(5) + t180, t4 * qJD(2) + t14 * qJD(4) + (-mrSges(6,1) * t43 - mrSges(6,2) * t42) * qJD(5); -t1 * qJD(4) - t3 * qJD(5) + t171, qJD(4) * t6 - qJD(5) * t10, qJD(1) * t200 + t238, t8 * qJD(5) + t170 + (mrSges(5,2) * t214 - Ifges(5,6) * t155 - pkin(4) * t116 + t104 * t243 + t106 * t246 + t183 * t245 + (t174 + t206 - t211) * pkin(8) + (t257 * t143 + Ifges(5,5) + t172) * t158 + t258 * mrSges(6,3)) * qJD(4), t8 * qJD(4) + (-mrSges(6,1) * t68 - mrSges(6,2) * t67 - t155 * t183) * qJD(5) + t167; t5, qJD(1) * t199 + t238, qJD(4) * t251, (t114 + m(6) * (pkin(8) * t190 - t239) + mrSges(6,3) * t190 - t186) * qJD(4) + t70 * qJD(5) + t164, t70 * qJD(4) + t114 * qJD(5) + t203; qJD(2) * t1 - qJD(5) * t13 - t180, qJD(5) * t9 - t170, t69 * qJD(5) - t164, -t37 * qJD(5), (t128 * pkin(8) + t191) * qJD(5) - t161; t3 * qJD(2) + t13 * qJD(4), -qJD(4) * t9 - t167, -t69 * qJD(4) - t203, t161, 0;];
Cq = t11;
