% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRPRPR6
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S6PRPRPR6_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta5]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:39
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPRPR6_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR6_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR6_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRPR6_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR6_jacobiaD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:39:09
	% EndTime: 2019-10-09 21:39:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:39:09
	% EndTime: 2019-10-09 21:39:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:39:09
	% EndTime: 2019-10-09 21:39:09
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (46->7), mult. (159->21), div. (18->4), fcn. (175->5), ass. (0->15)
	t39 = cos(pkin(10));
	t41 = sin(qJ(2));
	t42 = cos(qJ(2));
	t45 = sin(pkin(10)) * cos(pkin(6));
	t37 = t39 * t42 - t41 * t45;
	t34 = 0.1e1 / t37 ^ 2;
	t49 = qJD(2) * t34;
	t36 = t39 * t41 + t42 * t45;
	t33 = t36 ^ 2;
	t30 = t33 * t34 + 0.1e1;
	t46 = t37 * t49;
	t47 = t36 / t37 * t49;
	t48 = (t33 * t47 + t36 * t46) / t30 ^ 2;
	t28 = 0.1e1 / t30;
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, -0.2e1 * t48 + 0.2e1 * (t28 * t46 + (t28 * t47 - t34 * t48) * t36) * t36, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:39:09
	% EndTime: 2019-10-09 21:39:09
	% DurationCPUTime: 0.47s
	% Computational Cost: add. (536->40), mult. (1574->99), div. (384->14), fcn. (2127->9), ass. (0->50)
	t90 = sin(pkin(6));
	t85 = 0.1e1 / t90 ^ 2;
	t89 = sin(pkin(10));
	t118 = 0.1e1 / t89 ^ 2 * t85;
	t92 = cos(qJ(2));
	t110 = t90 * t92;
	t108 = cos(pkin(10));
	t109 = cos(pkin(6));
	t104 = t109 * t108;
	t91 = sin(qJ(2));
	t78 = -t92 * t104 + t89 * t91;
	t67 = atan2(-t78, -t110);
	t65 = sin(t67);
	t66 = cos(t67);
	t63 = -t66 * t110 - t65 * t78;
	t60 = 0.1e1 / t63;
	t84 = 0.1e1 / t90;
	t86 = 0.1e1 / t92;
	t61 = 0.1e1 / t63 ^ 2;
	t87 = 0.1e1 / t92 ^ 2;
	t105 = t89 * t109;
	t101 = -t92 * t105 - t108 * t91;
	t115 = -0.2e1 * t101;
	t103 = t65 * t110 - t66 * t78;
	t107 = t66 * t90 * t91;
	t111 = t87 * t91;
	t106 = t78 * t111;
	t76 = t78 ^ 2;
	t71 = t76 * t85 * t87 + 0.1e1;
	t68 = 0.1e1 / t71;
	t112 = t68 * t84;
	t80 = t91 * t104 + t89 * t92;
	t73 = t80 * qJD(2);
	t55 = (qJD(2) * t106 + t73 * t86) * t112;
	t53 = qJD(2) * t107 + t103 * t55 - t65 * t73;
	t114 = t53 * t60 * t61;
	t113 = t61 * t101;
	t102 = t80 * t86 + t106;
	t82 = -t91 * t105 + t108 * t92;
	t88 = t86 * t87;
	t77 = t101 ^ 2;
	t75 = t82 * qJD(2);
	t74 = t101 * qJD(2);
	t72 = qJD(2) * t78;
	t70 = t82 ^ 2 * t118 + 0.1e1;
	t59 = t77 * t61 + 0.1e1;
	t56 = t102 * t112;
	t54 = t103 * t56 - t65 * t80 + t107;
	t52 = (-0.2e1 * t102 / t71 ^ 2 * (qJD(2) * t76 * t88 * t91 + t73 * t78 * t87) * t85 + (t73 * t111 - t72 * t86 + (t80 * t111 + (0.2e1 * t88 * t91 ^ 2 + t86) * t78) * qJD(2)) * t68) * t84;
	t1 = [0, t52, 0, 0, 0, 0; 0, 0.2e1 * (-t54 * t113 - t60 * t82) / t59 ^ 2 * (-t75 * t113 - t77 * t114) + (t54 * t114 * t115 + t74 * t60 + (-t82 * t53 - t54 * t75 - (-(qJD(2) * t110 - t52 * t78 - t56 * t73 + (t56 * t110 - t80) * t55) * t66 - (t55 * t56 * t78 + t72 + (t52 * t92 + (-qJD(2) * t56 - t55) * t91) * t90) * t65) * t101) * t61) / t59, 0, 0, 0, 0; 0, (-t75 / t70 + 0.1e1 / t70 ^ 2 * t82 * t74 * t115 * t118) * t84 / t89, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:39:09
	% EndTime: 2019-10-09 21:39:09
	% DurationCPUTime: 0.53s
	% Computational Cost: add. (682->54), mult. (2271->134), div. (423->14), fcn. (2956->11), ass. (0->64)
	t135 = sin(pkin(10));
	t137 = cos(pkin(10));
	t140 = sin(qJ(2));
	t138 = cos(pkin(6));
	t142 = cos(qJ(2));
	t155 = t138 * t142;
	t123 = t135 * t140 - t137 * t155;
	t156 = t138 * t140;
	t124 = t135 * t142 + t137 * t156;
	t136 = sin(pkin(6));
	t157 = t136 * t140;
	t114 = atan2(-t124, t157);
	t110 = sin(t114);
	t111 = cos(t114);
	t97 = -t110 * t124 + t111 * t157;
	t94 = 0.1e1 / t97;
	t126 = t135 * t155 + t137 * t140;
	t139 = sin(qJ(4));
	t141 = cos(qJ(4));
	t159 = t135 * t136;
	t109 = t126 * t139 + t141 * t159;
	t105 = 0.1e1 / t109;
	t132 = 0.1e1 / t140;
	t106 = 0.1e1 / t109 ^ 2;
	t133 = 0.1e1 / t140 ^ 2;
	t95 = 0.1e1 / t97 ^ 2;
	t108 = -t126 * t141 + t139 * t159;
	t104 = t108 ^ 2;
	t101 = t104 * t106 + 0.1e1;
	t127 = -t135 * t156 + t137 * t142;
	t120 = t127 * qJD(2);
	t103 = t109 * qJD(4) - t120 * t141;
	t162 = t106 * t108;
	t153 = qJD(4) * t108;
	t102 = t120 * t139 - t153;
	t163 = t102 * t105 * t106;
	t166 = 0.1e1 / t101 ^ 2 * (t103 * t162 - t104 * t163);
	t160 = t133 * t142;
	t152 = t124 * t160;
	t149 = t123 * t132 + t152;
	t121 = t124 ^ 2;
	t131 = 0.1e1 / t136 ^ 2;
	t115 = t121 * t131 * t133 + 0.1e1;
	t112 = 0.1e1 / t115;
	t130 = 0.1e1 / t136;
	t161 = t112 * t130;
	t90 = t149 * t161;
	t165 = t124 * t90;
	t164 = t127 * t95;
	t154 = qJD(2) * t142;
	t150 = t105 * t141 + t139 * t162;
	t134 = t132 * t133;
	t122 = t127 ^ 2;
	t119 = t126 * qJD(2);
	t118 = t124 * qJD(2);
	t117 = t123 * qJD(2);
	t99 = 0.1e1 / t101;
	t96 = t94 * t95;
	t93 = t122 * t95 + 0.1e1;
	t89 = (qJD(2) * t152 + t117 * t132) * t161;
	t87 = (t136 * t142 - t165) * t111 + (-t90 * t157 + t123) * t110;
	t86 = (-t124 * t89 + t136 * t154) * t111 + (-t89 * t157 + t117) * t110;
	t85 = (-0.2e1 * t149 * (-t117 * t124 * t133 - t121 * t134 * t154) * t131 / t115 ^ 2 + (-t117 * t160 + t118 * t132 + (-t123 * t160 + (-0.2e1 * t134 * t142 ^ 2 - t132) * t124) * qJD(2)) * t112) * t130;
	t1 = [0, t85, 0, 0, 0, 0; 0, 0.2e1 * (t126 * t94 + t87 * t164) / t93 ^ 2 * (-t122 * t96 * t86 - t119 * t164) + (t87 * t119 * t95 - t120 * t94 + (0.2e1 * t87 * t127 * t96 + t126 * t95) * t86 + (-(t117 * t90 + t123 * t89 - t124 * t85 + (-t89 * t90 - qJD(2)) * t157) * t111 - (t89 * t165 + t118 + (-t140 * t85 + (-qJD(2) * t90 - t89) * t142) * t136) * t110) * t164) / t93, 0, 0, 0, 0; 0, t150 * t99 * t119 + (0.2e1 * t150 * t166 + ((qJD(4) * t105 + 0.2e1 * t108 * t163) * t139 + (-t103 * t139 + (t102 - t153) * t141) * t106) * t99) * t127, 0, -0.2e1 * t166 + 0.2e1 * (t103 * t106 * t99 + (-t106 * t166 - t99 * t163) * t108) * t108, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:39:09
	% EndTime: 2019-10-09 21:39:10
	% DurationCPUTime: 1.18s
	% Computational Cost: add. (2688->100), mult. (8196->220), div. (535->12), fcn. (10572->13), ass. (0->99)
	t187 = sin(pkin(10));
	t190 = cos(pkin(10));
	t195 = cos(qJ(2));
	t191 = cos(pkin(6));
	t193 = sin(qJ(2));
	t218 = t191 * t193;
	t180 = t187 * t195 + t190 * t218;
	t174 = t180 * qJD(2);
	t188 = sin(pkin(6));
	t194 = cos(qJ(4));
	t192 = sin(qJ(4));
	t217 = t191 * t195;
	t205 = -t187 * t193 + t190 * t217;
	t202 = t205 * t192;
	t146 = qJD(4) * t202 + (qJD(4) * t188 * t190 + t174) * t194;
	t222 = t188 * t192;
	t165 = t190 * t222 - t205 * t194;
	t162 = t165 ^ 2;
	t219 = t188 * t195;
	t183 = t191 * t192 + t194 * t219;
	t178 = 0.1e1 / t183 ^ 2;
	t157 = t162 * t178 + 0.1e1;
	t155 = 0.1e1 / t157;
	t184 = t191 * t194 - t192 * t219;
	t221 = t188 * t193;
	t208 = qJD(2) * t221;
	t168 = t184 * qJD(4) - t194 * t208;
	t177 = 0.1e1 / t183;
	t227 = t165 * t178;
	t127 = (t146 * t177 - t168 * t227) * t155;
	t158 = atan2(t165, t183);
	t153 = sin(t158);
	t154 = cos(t158);
	t207 = -t153 * t183 + t154 * t165;
	t123 = t207 * t127 + t153 * t146 + t154 * t168;
	t137 = t153 * t165 + t154 * t183;
	t134 = 0.1e1 / t137;
	t135 = 0.1e1 / t137 ^ 2;
	t238 = t123 * t134 * t135;
	t181 = t187 * t217 + t190 * t193;
	t163 = -t181 * t194 + t187 * t222;
	t237 = 0.2e1 * t163 * t238;
	t225 = t168 * t177 * t178;
	t236 = (t146 * t227 - t162 * t225) / t157 ^ 2;
	t210 = t165 * t221;
	t203 = t177 * t180 + t178 * t210;
	t235 = t194 * t203;
	t220 = t188 * t194;
	t164 = t181 * t192 + t187 * t220;
	t209 = t187 * t218;
	t182 = t190 * t195 - t209;
	t186 = sin(pkin(11));
	t189 = cos(pkin(11));
	t152 = t164 * t189 + t182 * t186;
	t148 = 0.1e1 / t152;
	t149 = 0.1e1 / t152 ^ 2;
	t234 = t135 * t163;
	t216 = qJD(2) * t195;
	t176 = -qJD(2) * t209 + t190 * t216;
	t143 = -t163 * qJD(4) + t176 * t192;
	t175 = t181 * qJD(2);
	t142 = t143 * t189 - t175 * t186;
	t233 = t142 * t148 * t149;
	t232 = t148 * t186;
	t151 = t164 * t186 - t182 * t189;
	t231 = t149 * t151;
	t230 = t151 * t189;
	t229 = t153 * t163;
	t228 = t154 * t163;
	t226 = t165 * t184;
	t224 = t182 * t192;
	t223 = t182 * t194;
	t215 = qJD(4) * t192;
	t161 = t163 ^ 2;
	t133 = t135 * t161 + 0.1e1;
	t144 = t164 * qJD(4) - t176 * t194;
	t214 = 0.2e1 * (t144 * t234 - t161 * t238) / t133 ^ 2;
	t147 = t151 ^ 2;
	t140 = t147 * t149 + 0.1e1;
	t141 = t143 * t186 + t175 * t189;
	t213 = 0.2e1 * (t141 * t231 - t147 * t233) / t140 ^ 2;
	t211 = t151 * t233;
	t166 = t190 * t220 + t202;
	t206 = -t166 * t177 + t178 * t226;
	t204 = qJD(4) * t223 - t175 * t192;
	t173 = t205 * qJD(2);
	t167 = -t183 * qJD(4) + t192 * t208;
	t160 = -t181 * t186 + t189 * t224;
	t159 = t181 * t189 + t186 * t224;
	t145 = t165 * qJD(4) + t174 * t192;
	t138 = 0.1e1 / t140;
	t130 = 0.1e1 / t133;
	t129 = t155 * t235;
	t128 = t206 * t155;
	t125 = (t153 * t180 - t154 * t221) * t194 + t207 * t129;
	t124 = -t207 * t128 + t153 * t166 + t154 * t184;
	t122 = 0.2e1 * t206 * t236 + (0.2e1 * t225 * t226 - t145 * t177 + (-t146 * t184 - t165 * t167 - t166 * t168) * t178) * t155;
	t120 = -0.2e1 * t235 * t236 + (-t203 * t215 + (-0.2e1 * t210 * t225 + t173 * t177 + (-t168 * t180 + (t146 * t193 + t165 * t216) * t188) * t178) * t194) * t155;
	t1 = [0, t120, 0, t122, 0, 0; 0, (t125 * t234 + t134 * t223) * t214 + ((t175 * t194 + t182 * t215) * t134 + t125 * t237 + (-t125 * t144 + t223 * t123 - (t120 * t165 + t129 * t146 + (t193 * t215 - t194 * t216) * t188 + (-t129 * t183 + t180 * t194) * t127) * t228 - (-t180 * t215 - t120 * t183 - t129 * t168 + t173 * t194 + (-t129 * t165 + t193 * t220) * t127) * t229) * t135) * t130, 0, (t124 * t234 - t134 * t164) * t214 + (t124 * t237 + t143 * t134 + (-t164 * t123 - t124 * t144 - (t122 * t165 - t128 * t146 + t167 + (t128 * t183 + t166) * t127) * t228 - (-t122 * t183 + t128 * t168 - t145 + (t128 * t165 - t184) * t127) * t229) * t135) * t130, 0, 0; 0, (-t148 * t159 + t160 * t231) * t213 + ((t176 * t189 + t204 * t186) * t148 + 0.2e1 * t160 * t211 + (-t159 * t142 - (-t176 * t186 + t204 * t189) * t151 - t160 * t141) * t149) * t138, 0, (-t149 * t230 + t232) * t163 * t213 + (-0.2e1 * t163 * t189 * t211 - t144 * t232 + (t144 * t230 + (t141 * t189 + t142 * t186) * t163) * t149) * t138, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:39:09
	% EndTime: 2019-10-09 21:39:10
	% DurationCPUTime: 1.33s
	% Computational Cost: add. (3313->109), mult. (9085->231), div. (559->12), fcn. (11668->13), ass. (0->105)
	t217 = sin(pkin(10));
	t219 = cos(pkin(10));
	t224 = cos(qJ(2));
	t220 = cos(pkin(6));
	t222 = sin(qJ(2));
	t250 = t220 * t222;
	t208 = t217 * t224 + t219 * t250;
	t202 = t208 * qJD(2);
	t218 = sin(pkin(6));
	t223 = cos(qJ(4));
	t221 = sin(qJ(4));
	t249 = t220 * t224;
	t234 = -t217 * t222 + t219 * t249;
	t232 = t234 * t221;
	t179 = qJD(4) * t232 + (t219 * t218 * qJD(4) + t202) * t223;
	t254 = t218 * t221;
	t192 = t219 * t254 - t234 * t223;
	t189 = t192 ^ 2;
	t251 = t218 * t224;
	t211 = t220 * t221 + t223 * t251;
	t206 = 0.1e1 / t211 ^ 2;
	t184 = t189 * t206 + 0.1e1;
	t182 = 0.1e1 / t184;
	t212 = t220 * t223 - t221 * t251;
	t253 = t218 * t222;
	t239 = qJD(2) * t253;
	t195 = t212 * qJD(4) - t223 * t239;
	t205 = 0.1e1 / t211;
	t259 = t192 * t206;
	t154 = (t179 * t205 - t195 * t259) * t182;
	t185 = atan2(t192, t211);
	t180 = sin(t185);
	t181 = cos(t185);
	t237 = -t180 * t211 + t181 * t192;
	t150 = t237 * t154 + t180 * t179 + t181 * t195;
	t166 = t180 * t192 + t181 * t211;
	t163 = 0.1e1 / t166;
	t164 = 0.1e1 / t166 ^ 2;
	t272 = t150 * t163 * t164;
	t209 = t217 * t249 + t219 * t222;
	t190 = -t209 * t223 + t217 * t254;
	t271 = 0.2e1 * t190 * t272;
	t257 = t195 * t205 * t206;
	t270 = (t179 * t259 - t189 * t257) / t184 ^ 2;
	t241 = t192 * t253;
	t233 = t205 * t208 + t206 * t241;
	t269 = t223 * t233;
	t252 = t218 * t223;
	t191 = t209 * t221 + t217 * t252;
	t240 = t217 * t250;
	t210 = t219 * t224 - t240;
	t216 = pkin(11) + qJ(6);
	t214 = sin(t216);
	t215 = cos(t216);
	t175 = t191 * t215 + t210 * t214;
	t171 = 0.1e1 / t175;
	t172 = 0.1e1 / t175 ^ 2;
	t248 = qJD(2) * t224;
	t204 = -qJD(2) * t240 + t219 * t248;
	t176 = -t190 * qJD(4) + t204 * t221;
	t203 = t209 * qJD(2);
	t161 = t175 * qJD(6) + t176 * t214 + t203 * t215;
	t174 = t191 * t214 - t210 * t215;
	t170 = t174 ^ 2;
	t169 = t170 * t172 + 0.1e1;
	t264 = t172 * t174;
	t246 = qJD(6) * t174;
	t162 = t176 * t215 - t203 * t214 - t246;
	t267 = t162 * t171 * t172;
	t268 = (t161 * t264 - t170 * t267) / t169 ^ 2;
	t266 = t164 * t190;
	t265 = t171 * t214;
	t263 = t174 * t215;
	t177 = t191 * qJD(4) - t204 * t223;
	t262 = t177 * t164;
	t261 = t180 * t190;
	t260 = t181 * t190;
	t258 = t192 * t212;
	t256 = t210 * t221;
	t255 = t210 * t223;
	t247 = qJD(4) * t221;
	t188 = t190 ^ 2;
	t160 = t188 * t164 + 0.1e1;
	t245 = 0.2e1 * (-t188 * t272 + t190 * t262) / t160 ^ 2;
	t244 = -0.2e1 * t268;
	t242 = t174 * t267;
	t238 = qJD(6) * t256 + t204;
	t236 = t172 * t263 - t265;
	t193 = t219 * t252 + t232;
	t235 = -t193 * t205 + t206 * t258;
	t231 = qJD(4) * t255 - qJD(6) * t209 - t203 * t221;
	t201 = t234 * qJD(2);
	t194 = -t211 * qJD(4) + t221 * t239;
	t187 = -t209 * t214 + t215 * t256;
	t186 = t209 * t215 + t214 * t256;
	t178 = t192 * qJD(4) + t202 * t221;
	t167 = 0.1e1 / t169;
	t157 = 0.1e1 / t160;
	t156 = t182 * t269;
	t155 = t235 * t182;
	t152 = (t180 * t208 - t181 * t253) * t223 + t237 * t156;
	t151 = -t237 * t155 + t180 * t193 + t181 * t212;
	t149 = 0.2e1 * t235 * t270 + (0.2e1 * t257 * t258 - t178 * t205 + (-t179 * t212 - t192 * t194 - t193 * t195) * t206) * t182;
	t147 = -0.2e1 * t269 * t270 + (-t233 * t247 + (-0.2e1 * t241 * t257 + t201 * t205 + (-t195 * t208 + (t179 * t222 + t192 * t248) * t218) * t206) * t223) * t182;
	t1 = [0, t147, 0, t149, 0, 0; 0, (t152 * t266 + t163 * t255) * t245 + ((t203 * t223 + t210 * t247) * t163 + (-t262 + t271) * t152 + (t255 * t150 - (t147 * t192 + t156 * t179 + (t222 * t247 - t223 * t248) * t218 + (-t156 * t211 + t208 * t223) * t154) * t260 - (-t208 * t247 - t147 * t211 - t156 * t195 + t201 * t223 + (-t156 * t192 + t222 * t252) * t154) * t261) * t164) * t157, 0, (t151 * t266 - t163 * t191) * t245 + (t151 * t271 + t176 * t163 + (-t191 * t150 - t151 * t177 - (t149 * t192 - t155 * t179 + t194 + (t155 * t211 + t193) * t154) * t260 - (-t149 * t211 + t155 * t195 - t178 + (t155 * t192 - t212) * t154) * t261) * t164) * t157, 0, 0; 0, 0.2e1 * (-t171 * t186 + t187 * t264) * t268 + (0.2e1 * t187 * t242 + t238 * t171 * t215 + t231 * t265 + (t238 * t174 * t214 - t187 * t161 - t186 * t162 - t231 * t263) * t172) * t167, 0, t236 * t190 * t244 + (t236 * t177 + ((-qJD(6) * t171 - 0.2e1 * t242) * t215 + (t161 * t215 + (t162 - t246) * t214) * t172) * t190) * t167, 0, t244 + 0.2e1 * (t161 * t172 * t167 + (-t167 * t267 - t172 * t268) * t174) * t174;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end