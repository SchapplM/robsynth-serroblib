% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5PPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5PPRPR5_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR5_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR5_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRPR5_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR5_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRPR5_invdynB_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:33:31
% EndTime: 2019-12-31 17:33:35
% DurationCPUTime: 1.90s
% Computational Cost: add. (2714->242), mult. (4710->317), div. (0->0), fcn. (2996->6), ass. (0->147)
t196 = sin(qJ(3));
t200 = qJD(3) ^ 2;
t198 = cos(qJ(3));
t217 = t198 * qJDD(3);
t171 = t196 * t200 - t217;
t189 = g(3) - qJDD(1);
t149 = pkin(5) * t171 + t189 * t196;
t191 = sin(pkin(7));
t192 = cos(pkin(7));
t219 = t196 * qJDD(3);
t170 = t198 * t200 + t219;
t205 = pkin(5) * t170 + t189 * t198;
t210 = -t170 * t191 + t171 * t192;
t263 = -qJ(1) * t210 + t192 * t149 - t191 * t205;
t212 = t170 * t192 + t171 * t191;
t262 = qJ(1) * t212 - t191 * t149 - t192 * t205;
t174 = g(1) * t191 - g(2) * t192;
t167 = -qJDD(2) + t174;
t175 = g(1) * t192 + g(2) * t191;
t131 = -t196 * t167 - t198 * t175;
t201 = -(2 * qJD(4) * qJD(3)) - t131;
t221 = qJDD(3) * qJ(4);
t109 = -pkin(3) * t200 - t201 + t221;
t130 = t167 * t198 - t196 * t175;
t190 = qJDD(3) * pkin(3);
t110 = -qJ(4) * t200 + qJDD(4) + t130 - t190;
t88 = t109 * t196 - t110 * t198;
t90 = t109 * t198 + t110 * t196;
t66 = t191 * t90 - t192 * t88;
t67 = t191 * t88 + t192 * t90;
t95 = t130 * t198 - t131 * t196;
t96 = t130 * t196 + t131 * t198;
t72 = t191 * t96 + t192 * t95;
t251 = t191 * t95 - t192 * t96;
t243 = pkin(1) + pkin(2);
t246 = qJ(2) * t171 + t170 * t243;
t245 = -qJ(2) * t170 + t171 * t243 + t130;
t108 = -qJDD(3) * pkin(6) + t110;
t195 = sin(qJ(5));
t197 = cos(qJ(5));
t100 = t108 * t195 + t189 * t197;
t213 = t108 * t197 - t189 * t195;
t79 = t197 * t100 - t195 * t213;
t78 = t195 * t100 + t197 * t213;
t242 = pkin(3) + pkin(6);
t187 = t195 ^ 2;
t240 = t187 * t200;
t188 = t197 ^ 2;
t239 = t188 * t200;
t234 = t191 * t189;
t181 = t192 * t189;
t107 = -pkin(6) * t200 + t109;
t230 = t195 * t107;
t216 = t197 * t200 * t195;
t176 = qJDD(5) + t216;
t229 = t195 * t176;
t223 = t187 + t188;
t169 = t223 * qJDD(3);
t228 = t196 * t169;
t226 = t197 * t107;
t225 = t197 * t176;
t224 = t198 * t169;
t222 = qJD(3) * qJD(5);
t220 = t195 * qJDD(3);
t218 = t197 * qJDD(3);
t215 = t195 * t222;
t214 = t197 * t222;
t159 = t191 * t175;
t120 = t167 * t192 - t159;
t134 = t174 * t192 - t159;
t160 = t192 * t175;
t121 = -t167 * t191 - t160;
t135 = -t174 * t191 - t160;
t199 = qJD(5) ^ 2;
t208 = -t199 - t239;
t207 = t196 * t216;
t206 = t198 * t216;
t204 = qJDD(5) - t216;
t202 = t197 * t204;
t179 = t199 - t239;
t178 = -t199 - t240;
t177 = -t199 + t240;
t173 = (-t187 + t188) * t200;
t172 = t223 * t200;
t166 = -0.2e1 * t215 + t218;
t165 = -t215 + t218;
t164 = -t214 - t220;
t163 = 0.2e1 * t214 + t220;
t162 = t195 * t204;
t161 = t223 * t222;
t147 = qJDD(5) * t198 - t161 * t196;
t146 = -qJDD(5) * t196 - t161 * t198;
t145 = -t165 * t195 - t188 * t222;
t144 = -t164 * t197 - t187 * t222;
t141 = -t195 * t208 - t225;
t140 = t178 * t197 - t162;
t139 = t197 * t208 - t229;
t138 = -t179 * t197 - t162;
t137 = t195 * t178 + t202;
t136 = -t177 * t195 - t225;
t133 = -t172 * t198 - t228;
t132 = -t172 * t196 + t224;
t119 = t163 * t195 - t166 * t197;
t118 = -t144 * t196 - t206;
t117 = -t145 * t196 + t206;
t116 = -t145 * t198 - t207;
t115 = -t144 * t198 + t207;
t114 = -t138 * t196 + t197 * t217;
t113 = -t136 * t196 - t195 * t217;
t112 = -t138 * t198 - t196 * t218;
t111 = -t136 * t198 + t195 * t219;
t106 = t139 * t196 + t166 * t198;
t105 = t137 * t196 + t163 * t198;
t104 = -t139 * t198 + t166 * t196;
t103 = -t137 * t198 + t163 * t196;
t102 = -t119 * t196 + t173 * t198;
t101 = -t119 * t198 - t173 * t196;
t98 = t132 * t191 + t133 * t192;
t97 = -t132 * t192 + t133 * t191;
t92 = pkin(5) * t95 + qJ(2) * t189;
t91 = -pkin(5) * t96 + t189 * t243;
t86 = t104 * t191 + t106 * t192;
t85 = t103 * t191 + t105 * t192;
t84 = -t104 * t192 + t106 * t191;
t83 = -t103 * t192 + t105 * t191;
t82 = pkin(4) * t139 - qJ(4) * t141 - t100;
t81 = pkin(4) * t137 - qJ(4) * t140 + t213;
t80 = -pkin(5) * t88 + (pkin(3) * t196 - qJ(4) * t198 + qJ(2)) * t189;
t77 = -pkin(5) * t90 + (pkin(3) * t198 + qJ(4) * t196 + t243) * t189;
t76 = pkin(4) * t163 - t140 * t242 + t226;
t75 = pkin(4) * t166 - t141 * t242 - t230;
t74 = -pkin(4) * t172 - t79;
t71 = t107 * t198 + t196 * t78;
t70 = t107 * t196 - t198 * t78;
t69 = -pkin(4) * t224 - pkin(5) * t132 - t196 * t74;
t68 = pkin(4) * t228 - pkin(5) * t133 - t198 * t74;
t65 = pkin(4) * t78 - qJ(4) * t79;
t64 = pkin(4) * t107 - t242 * t79;
t63 = -pkin(5) * t104 + qJ(2) * t141 - t196 * t75 + t198 * t82;
t62 = -pkin(5) * t103 + qJ(2) * t140 - t196 * t76 + t198 * t81;
t61 = -pkin(5) * t106 + t141 * t243 - t196 * t82 - t198 * t75;
t60 = -pkin(5) * t105 + t140 * t243 - t196 * t81 - t198 * t76;
t59 = t191 * t70 + t192 * t71;
t58 = t191 * t71 - t192 * t70;
t57 = -pkin(5) * t70 + qJ(2) * t79 - t196 * t64 + t198 * t65;
t56 = -pkin(5) * t71 - t196 * t65 - t198 * t64 + t243 * t79;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t135, 0, 0, 0, 0, 0, 0, 0, 0, 0, t121, 0, 0, 0, 0, 0, 0, -t212, t210, 0, -t251, 0, 0, 0, 0, 0, 0, 0, t212, -t210, t67, 0, 0, 0, 0, 0, 0, t85, t86, t98, t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t134, 0, 0, 0, 0, 0, 0, 0, 0, 0, t120, 0, 0, 0, 0, 0, 0, t210, t212, 0, t72, 0, 0, 0, 0, 0, 0, 0, -t210, -t212, t66, 0, 0, 0, 0, 0, 0, t83, t84, t97, t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t189, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t189, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t189, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t189, 0, 0, 0, 0, 0, 0, -t140, -t141, 0, -t79; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t234, -t181, -t134, -qJ(1) * t134, 0, 0, 0, 0, 0, 0, -t234, -t120, t181, -qJ(1) * t120 + (-pkin(1) * t191 + qJ(2) * t192) * t189, 0, 0, -t210, 0, -t212, 0, t263, -t262, t72, -qJ(1) * t72 - t191 * t91 + t192 * t92, 0, t210, t212, 0, 0, 0, t66, -t263, t262, -qJ(1) * t66 - t191 * t77 + t192 * t80, -t116 * t191 + t117 * t192, -t101 * t191 + t102 * t192, -t112 * t191 + t114 * t192, -t115 * t191 + t118 * t192, -t111 * t191 + t113 * t192, -t146 * t191 + t147 * t192, -qJ(1) * t83 - t191 * t60 + t192 * t62, -qJ(1) * t84 - t191 * t61 + t192 * t63, -qJ(1) * t97 - t191 * t68 + t192 * t69, -qJ(1) * t58 - t191 * t56 + t192 * t57; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t181, -t234, t135, qJ(1) * t135, 0, 0, 0, 0, 0, 0, t181, t121, t234, qJ(1) * t121 + (pkin(1) * t192 + qJ(2) * t191) * t189, 0, 0, -t212, 0, t210, 0, -t262, -t263, t251, -qJ(1) * t251 + t191 * t92 + t192 * t91, 0, t212, -t210, 0, 0, 0, -t67, t262, t263, qJ(1) * t67 + t191 * t80 + t192 * t77, t116 * t192 + t117 * t191, t101 * t192 + t102 * t191, t112 * t192 + t114 * t191, t115 * t192 + t118 * t191, t111 * t192 + t113 * t191, t146 * t192 + t147 * t191, qJ(1) * t85 + t191 * t62 + t192 * t60, qJ(1) * t86 + t191 * t63 + t192 * t61, qJ(1) * t98 + t191 * t69 + t192 * t68, qJ(1) * t59 + t191 * t57 + t192 * t56; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t174, t175, 0, 0, 0, 0, 0, 0, 0, 0, t167, 0, -t175, pkin(1) * t167 - qJ(2) * t175, 0, 0, 0, 0, 0, -qJDD(3), t245, t131 + t246, 0, qJ(2) * t96 + t243 * t95, -qJDD(3), 0, 0, 0, 0, 0, 0, -qJDD(4) + 0.2e1 * t190 - t245, t201 - 0.2e1 * t221 - t246, pkin(3) * t110 + qJ(2) * t90 - qJ(4) * t109 - t243 * t88, (-t165 + t215) * t197, t163 * t197 + t166 * t195, t195 * t179 - t202, (t164 - t214) * t195, -t177 * t197 + t229, 0, qJ(2) * t105 - qJ(4) * t163 - t103 * t243 + t137 * t242 - t230, qJ(2) * t106 - qJ(4) * t166 - t104 * t243 + t139 * t242 - t226, qJ(2) * t133 + qJ(4) * t172 - t132 * t243 - t169 * t242 + t78, qJ(2) * t71 - qJ(4) * t107 + t242 * t78 - t243 * t70;];
tauB_reg = t1;
