% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S4RPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
% 
% Output:
% tauB_reg [6x(5*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S4RPPR5_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR5_invdynB_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR5_invdynB_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR5_invdynB_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR5_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR5_invdynB_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:39:50
% EndTime: 2019-12-31 16:39:52
% DurationCPUTime: 1.36s
% Computational Cost: add. (2804->215), mult. (5021->290), div. (0->0), fcn. (2300->6), ass. (0->135)
t190 = sin(pkin(6));
t191 = cos(pkin(6));
t198 = qJD(1) ^ 2;
t162 = -t190 * qJDD(1) + t191 * t198;
t163 = t191 * qJDD(1) + t190 * t198;
t194 = sin(qJ(1));
t196 = cos(qJ(1));
t123 = t196 * t162 + t194 * t163;
t189 = g(3) + qJDD(3);
t147 = qJ(3) * t163 + t190 * t189;
t203 = qJ(3) * t162 + t191 * t189;
t230 = -pkin(4) * t123 + t194 * t147 + t196 * t203;
t207 = -t194 * t162 + t196 * t163;
t228 = -pkin(4) * t207 + t196 * t147 - t194 * t203;
t186 = qJDD(1) * qJ(2);
t172 = t196 * g(1) + t194 * g(2);
t202 = (2 * qJD(2) * qJD(1)) - t172;
t201 = t186 + t202;
t224 = pkin(1) + pkin(2);
t143 = -t224 * t198 + t201;
t171 = t194 * g(1) - t196 * g(2);
t204 = -qJDD(2) + t171;
t200 = -t198 * qJ(2) - t204;
t199 = -t224 * qJDD(1) + t200;
t104 = t190 * t143 - t191 * t199;
t105 = t191 * t143 + t190 * t199;
t87 = t191 * t104 - t190 * t105;
t88 = t190 * t104 + t191 * t105;
t74 = t194 * t88 + t196 * t87;
t227 = t194 * t87 - t196 * t88;
t223 = qJDD(1) * pkin(1);
t193 = sin(qJ(4));
t187 = t193 ^ 2;
t222 = t187 * t198;
t102 = qJDD(1) * pkin(3) - t198 * pkin(5) + t104;
t221 = t193 * t102;
t195 = cos(qJ(4));
t177 = t195 * t198 * t193;
t169 = qJDD(4) + t177;
t220 = t193 * t169;
t170 = qJDD(4) - t177;
t219 = t193 * t170;
t216 = t195 * t102;
t215 = t195 * t169;
t214 = t195 * t170;
t103 = -t198 * pkin(3) - qJDD(1) * pkin(5) + t105;
t98 = t195 * t103 + t193 * t189;
t188 = t195 ^ 2;
t213 = -t187 - t188;
t212 = qJD(1) * qJD(4);
t211 = t193 * qJDD(1);
t210 = t195 * qJDD(1);
t209 = t193 * t212;
t208 = t195 * t212;
t97 = t193 * t103 - t195 * t189;
t150 = -t198 * pkin(1) + t201;
t151 = -t200 + t223;
t113 = t196 * t150 - t194 * t151;
t131 = -t194 * t171 - t196 * t172;
t206 = t190 * t177;
t205 = t191 * t177;
t165 = t194 * qJDD(1) + t196 * t198;
t153 = -pkin(4) * t165 + t196 * g(3);
t166 = t196 * qJDD(1) - t194 * t198;
t152 = pkin(4) * t166 + t194 * g(3);
t80 = t193 * t98 - t195 * t97;
t82 = t193 * t97 + t195 * t98;
t112 = t194 * t150 + t196 * t151;
t130 = t196 * t171 - t194 * t172;
t197 = qJD(4) ^ 2;
t183 = t188 * t198;
t176 = -t183 - t197;
t175 = t183 - t197;
t174 = -t197 - t222;
t173 = t197 - t222;
t168 = t183 - t222;
t167 = t183 + t222;
t164 = t213 * qJDD(1);
t161 = t209 - t210;
t160 = 0.2e1 * t209 - t210;
t159 = -t208 - t211;
t158 = 0.2e1 * t208 + t211;
t156 = t213 * t212;
t145 = t195 * t159 + t187 * t212;
t144 = -t193 * t161 + t188 * t212;
t141 = t190 * qJDD(4) + t191 * t156;
t140 = t191 * qJDD(4) - t190 * t156;
t137 = -t193 * t174 - t214;
t136 = -t193 * t173 + t215;
t135 = t195 * t176 - t220;
t134 = t195 * t175 - t219;
t133 = t195 * t174 - t219;
t132 = t193 * t176 + t215;
t128 = t191 * t164 - t190 * t167;
t127 = t190 * t164 + t191 * t167;
t122 = t193 * t158 + t195 * t160;
t121 = t191 * t145 - t206;
t120 = t191 * t144 + t206;
t119 = -t190 * t145 - t205;
t118 = -t190 * t144 + t205;
t117 = t191 * t136 - t190 * t211;
t116 = t191 * t134 - t190 * t210;
t115 = -t190 * t136 - t191 * t211;
t114 = -t190 * t134 - t191 * t210;
t111 = t191 * t137 - t190 * t158;
t110 = t191 * t135 - t190 * t160;
t109 = t190 * t137 + t191 * t158;
t108 = t190 * t135 + t191 * t160;
t107 = t191 * t122 - t190 * t168;
t106 = -t190 * t122 - t191 * t168;
t100 = t194 * t127 + t196 * t128;
t99 = -t196 * t127 + t194 * t128;
t96 = -pkin(5) * t133 + t216;
t95 = -pkin(5) * t132 + t221;
t94 = -pkin(3) * t133 + t98;
t93 = -pkin(3) * t132 + t97;
t92 = t194 * t109 + t196 * t111;
t91 = t194 * t108 + t196 * t110;
t90 = -t196 * t109 + t194 * t111;
t89 = -t196 * t108 + t194 * t110;
t84 = qJ(2) * t189 + qJ(3) * t87;
t83 = -qJ(3) * t88 + t224 * t189;
t79 = -qJ(3) * t127 - t191 * t80;
t78 = -qJ(3) * t128 + t190 * t80;
t77 = t190 * t102 + t191 * t82;
t76 = -t191 * t102 + t190 * t82;
t73 = qJ(2) * t133 - qJ(3) * t109 - t190 * t94 + t191 * t96;
t72 = qJ(2) * t132 - qJ(3) * t108 - t190 * t93 + t191 * t95;
t71 = -qJ(3) * t111 + t224 * t133 - t190 * t96 - t191 * t94;
t70 = -qJ(3) * t110 + t224 * t132 - t190 * t95 - t191 * t93;
t69 = t194 * t76 + t196 * t77;
t68 = t194 * t77 - t196 * t76;
t67 = -qJ(3) * t76 + (pkin(3) * t190 - pkin(5) * t191 + qJ(2)) * t80;
t66 = -qJ(3) * t77 + (pkin(3) * t191 + pkin(5) * t190 + t224) * t80;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t165, -t166, 0, t131, 0, 0, 0, 0, 0, 0, -t165, 0, t166, t113, 0, 0, 0, 0, 0, 0, -t123, t207, 0, -t227, 0, 0, 0, 0, 0, 0, t91, t92, t100, t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t166, -t165, 0, t130, 0, 0, 0, 0, 0, 0, t166, 0, t165, t112, 0, 0, 0, 0, 0, 0, t207, t123, 0, t74, 0, 0, 0, 0, 0, 0, t89, t90, t99, t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t189, 0, 0, 0, 0, 0, 0, -t132, -t133, 0, -t80; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t166, 0, -t165, 0, -t152, -t153, -t130, -pkin(4) * t130, 0, t166, 0, 0, t165, 0, -t152, -t112, t153, -pkin(4) * t112 + (-pkin(1) * t194 + qJ(2) * t196) * g(3), 0, 0, -t207, 0, -t123, 0, t228, t230, t74, -pkin(4) * t74 - t194 * t83 + t196 * t84, -t194 * t119 + t196 * t121, -t194 * t106 + t196 * t107, -t194 * t115 + t196 * t117, -t194 * t118 + t196 * t120, -t194 * t114 + t196 * t116, -t194 * t140 + t196 * t141, -pkin(4) * t89 - t194 * t70 + t196 * t72, -pkin(4) * t90 - t194 * t71 + t196 * t73, -pkin(4) * t99 - t194 * t78 + t196 * t79, -pkin(4) * t68 - t194 * t66 + t196 * t67; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t165, 0, t166, 0, t153, -t152, t131, pkin(4) * t131, 0, t165, 0, 0, -t166, 0, t153, t113, t152, pkin(4) * t113 + (pkin(1) * t196 + qJ(2) * t194) * g(3), 0, 0, -t123, 0, t207, 0, t230, -t228, t227, -pkin(4) * t227 + t194 * t84 + t196 * t83, t196 * t119 + t194 * t121, t196 * t106 + t194 * t107, t196 * t115 + t194 * t117, t196 * t118 + t194 * t120, t196 * t114 + t194 * t116, t196 * t140 + t194 * t141, pkin(4) * t91 + t194 * t72 + t196 * t70, pkin(4) * t92 + t194 * t73 + t196 * t71, pkin(4) * t100 + t194 * t79 + t196 * t78, pkin(4) * t69 + t194 * t67 + t196 * t66; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t171, t172, 0, 0, 0, 0, 0, qJDD(1), 0, 0, t204 + 0.2e1 * t223, 0, 0.2e1 * t186 + t202, pkin(1) * t151 + qJ(2) * t150, 0, 0, 0, 0, 0, qJDD(1), -qJ(2) * t162 + t224 * t163 + t104, qJ(2) * t163 + t224 * t162 + t105, 0, qJ(2) * t88 + t224 * t87, (-t159 + t208) * t193, t195 * t158 - t193 * t160, -t195 * t173 - t220, (-t161 - t209) * t195, -t193 * t175 - t214, 0, -pkin(3) * t160 - pkin(5) * t135 + qJ(2) * t110 - t224 * t108 + t216, -pkin(3) * t158 - pkin(5) * t137 + qJ(2) * t111 - t224 * t109 - t221, -pkin(3) * t167 - pkin(5) * t164 + qJ(2) * t128 - t224 * t127 - t82, pkin(3) * t102 - pkin(5) * t82 + qJ(2) * t77 - t224 * t76;];
tauB_reg = t1;
