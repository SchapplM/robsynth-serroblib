% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRRPR10
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S5RRRPR10_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RRRPR10_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR10_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR10_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPR10_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR10_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:31:28
	% EndTime: 2019-12-31 21:31:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:31:28
	% EndTime: 2019-12-31 21:31:28
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:31:28
	% EndTime: 2019-12-31 21:31:29
	% DurationCPUTime: 0.29s
	% Computational Cost: add. (215->39), mult. (853->106), div. (126->12), fcn. (1047->9), ass. (0->54)
	t99 = sin(pkin(5));
	t93 = t99 ^ 2;
	t100 = cos(pkin(5));
	t95 = 0.1e1 / t100 ^ 2;
	t104 = cos(qJ(1));
	t98 = t104 ^ 2;
	t89 = t98 * t93 * t95 + 0.1e1;
	t102 = sin(qJ(1));
	t97 = t102 ^ 2;
	t126 = 0.1e1 / t89 ^ 2 * t97;
	t131 = t126 * t95;
	t122 = t104 * t99;
	t88 = atan2(t122, t100);
	t84 = sin(t88);
	t85 = cos(t88);
	t72 = t85 * t100 + t84 * t122;
	t67 = 0.1e1 / t72;
	t103 = cos(qJ(2));
	t118 = t104 * t103;
	t101 = sin(qJ(2));
	t121 = t102 * t101;
	t113 = t100 * t121 - t118;
	t77 = 0.1e1 / t113;
	t94 = 0.1e1 / t100;
	t68 = 0.1e1 / t72 ^ 2;
	t78 = 0.1e1 / t113 ^ 2;
	t119 = t104 * t101;
	t120 = t102 * t103;
	t81 = -t100 * t119 - t120;
	t82 = t100 * t120 + t119;
	t71 = t81 * qJD(1) - t82 * qJD(2);
	t128 = t71 * t77 * t78;
	t115 = t100 * t118;
	t70 = -qJD(1) * t115 - qJD(2) * t118 + (qJD(2) * t100 + qJD(1)) * t121;
	t129 = t70 * t78;
	t76 = t82 ^ 2;
	t75 = t76 * t78 + 0.1e1;
	t130 = (t76 * t128 - t82 * t129) / t75 ^ 2;
	t127 = t81 * t82;
	t125 = t93 * t94;
	t124 = t102 * t68;
	t123 = t104 * t68;
	t117 = qJD(1) * t104;
	t86 = 0.1e1 / t89;
	t116 = (t86 - 0.1e1) * t99;
	t114 = -0.2e1 * t94 * t131;
	t80 = t115 - t121;
	t63 = (-t104 * t85 * t86 * t125 + t84 * t116) * t102;
	t92 = t99 * t93;
	t73 = 0.1e1 / t75;
	t69 = t67 * t68;
	t66 = t97 * t93 * t68 + 0.1e1;
	t62 = qJD(1) * t63;
	t1 = [(-t86 * t94 * t99 + t92 * t114) * t117, 0, 0, 0, 0; (0.2e1 * (-t104 * t67 + t63 * t124) / t66 ^ 2 * (-t62 * t69 * t97 + t117 * t124) * t93 + ((0.2e1 * t102 * t63 * t69 - t123) * t62 + (-t63 * t123 + (-t67 + (-t92 * t131 - t116) * t84 * t123 - (t93 ^ 2 * t98 * t114 + (-t126 + (0.2e1 * t97 - t98) * t86) * t125) * t68 * t85) * t102) * qJD(1)) / t66) * t99, 0, 0, 0, 0; 0.2e1 * (t78 * t127 + t77 * t80) * t130 + (-(-t82 * qJD(1) + t81 * qJD(2)) * t77 - 0.2e1 * t127 * t128 + (-t80 * t71 - (t113 * qJD(1) - t80 * qJD(2)) * t82 + t81 * t70) * t78) * t73, -0.2e1 * t130 + 0.2e1 * (-t73 * t129 + (t73 * t128 - t78 * t130) * t82) * t82, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:31:29
	% EndTime: 2019-12-31 21:31:29
	% DurationCPUTime: 0.76s
	% Computational Cost: add. (1479->91), mult. (4303->201), div. (668->14), fcn. (5516->11), ass. (0->91)
	t171 = sin(qJ(2));
	t172 = sin(qJ(1));
	t225 = cos(pkin(5));
	t195 = t172 * t225;
	t193 = t171 * t195;
	t174 = cos(qJ(2));
	t175 = cos(qJ(1));
	t209 = t175 * t174;
	t157 = -t193 + t209;
	t170 = sin(qJ(3));
	t173 = cos(qJ(3));
	t169 = sin(pkin(5));
	t213 = t169 * t172;
	t185 = -t157 * t170 + t173 * t213;
	t227 = t185 * qJD(3);
	t194 = t175 * t225;
	t192 = t174 * t194;
	t210 = t172 * t171;
	t153 = -t192 + t210;
	t212 = t169 * t174;
	t147 = atan2(-t153, -t212);
	t145 = sin(t147);
	t146 = cos(t147);
	t151 = t153 ^ 2;
	t165 = 0.1e1 / t169 ^ 2;
	t167 = 0.1e1 / t174 ^ 2;
	t150 = t151 * t165 * t167 + 0.1e1;
	t148 = 0.1e1 / t150;
	t164 = 0.1e1 / t169;
	t166 = 0.1e1 / t174;
	t199 = t153 * t164 * t166;
	t226 = (t146 * t199 - t145) * t148 + t145;
	t129 = -t145 * t153 - t146 * t212;
	t126 = 0.1e1 / t129;
	t144 = t157 * t173 + t170 * t213;
	t138 = 0.1e1 / t144;
	t127 = 0.1e1 / t129 ^ 2;
	t139 = 0.1e1 / t144 ^ 2;
	t182 = -t171 * t194 - t172 * t174;
	t183 = -t175 * t171 - t174 * t195;
	t135 = -t183 * qJD(1) - t182 * qJD(2);
	t207 = qJD(2) * t171;
	t196 = t167 * t207;
	t184 = t135 * t166 + t153 * t196;
	t215 = t148 * t164;
	t118 = t184 * t215;
	t188 = t145 * t212 - t146 * t153;
	t200 = t146 * t169 * t171;
	t114 = qJD(2) * t200 + t188 * t118 - t145 * t135;
	t224 = t114 * t126 * t127;
	t214 = t167 * t171;
	t187 = t153 * t214 - t166 * t182;
	t119 = t187 * t215;
	t115 = t188 * t119 + t145 * t182 + t200;
	t223 = t115 * t183;
	t134 = t182 * qJD(1) + t183 * qJD(2);
	t208 = qJD(1) * t169;
	t197 = t175 * t208;
	t124 = t144 * qJD(3) + t134 * t170 - t173 * t197;
	t137 = t185 ^ 2;
	t132 = t137 * t139 + 0.1e1;
	t218 = t139 * t185;
	t125 = t134 * t173 + t170 * t197 + t227;
	t220 = t125 * t138 * t139;
	t222 = (-t124 * t218 - t137 * t220) / t132 ^ 2;
	t168 = t166 * t167;
	t221 = (t135 * t153 * t167 + t151 * t168 * t207) * t165 / t150 ^ 2;
	t191 = qJD(2) * t225 + qJD(1);
	t206 = qJD(2) * t174;
	t133 = -qJD(1) * t192 - t175 * t206 + t191 * t210;
	t219 = t133 * t127;
	t217 = t145 * t183;
	t216 = t146 * t183;
	t211 = t169 * t175;
	t152 = t183 ^ 2;
	t122 = t127 * t152 + 0.1e1;
	t205 = 0.2e1 * (-t152 * t224 + t183 * t219) / t122 ^ 2;
	t204 = 0.2e1 * t224;
	t203 = 0.2e1 * t222;
	t202 = -0.2e1 * t221;
	t201 = t185 * t220;
	t198 = t172 * t208;
	t189 = t170 * t138 + t173 * t218;
	t186 = -t170 * t182 + t173 * t211;
	t142 = t170 * t211 + t173 * t182;
	t136 = -qJD(1) * t193 - t172 * t207 + t191 * t209;
	t130 = 0.1e1 / t132;
	t120 = 0.1e1 / t122;
	t117 = t226 * t183;
	t113 = (t187 * t202 + (t135 * t214 + t136 * t166 + (-t182 * t214 + (0.2e1 * t168 * t171 ^ 2 + t166) * t153) * qJD(2)) * t148) * t164;
	t1 = [(-t183 * t166 * t202 + (-t133 * t166 - t183 * t196) * t148) * t164, t113, 0, 0, 0; t153 * t126 * t205 + (-t135 * t126 + (t114 * t153 + t117 * t133) * t127) * t120 - ((t117 * t204 - t226 * t219) * t120 + (t117 * t205 + ((t118 * t148 * t199 + t202) * t217 + (0.2e1 * t199 * t221 - t118 + (-t184 * t164 + t118) * t148) * t216) * t120) * t127) * t183, (-t126 * t157 - t127 * t223) * t205 + (-t204 * t223 + t134 * t126 + (-t157 * t114 + t115 * t133 + (t169 * t206 - t113 * t153 - t119 * t135 + (t119 * t212 + t182) * t118) * t216 + (t118 * t119 * t153 - t136 + (t113 * t174 + (-qJD(2) * t119 - t118) * t171) * t169) * t217) * t127) * t120, 0, 0, 0; (t138 * t186 - t142 * t218) * t203 + ((t142 * qJD(3) - t136 * t170 + t173 * t198) * t138 - 0.2e1 * t142 * t201 + (t186 * t125 + (t186 * qJD(3) - t136 * t173 - t170 * t198) * t185 - t142 * t124) * t139) * t130, -t189 * t183 * t203 + (t189 * t133 - ((-qJD(3) * t138 + 0.2e1 * t201) * t173 + (t124 * t173 + (t125 + t227) * t170) * t139) * t183) * t130, -0.2e1 * t222 - 0.2e1 * (t124 * t139 * t130 - (-t130 * t220 - t139 * t222) * t185) * t185, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:31:29
	% EndTime: 2019-12-31 21:31:29
	% DurationCPUTime: 0.75s
	% Computational Cost: add. (1788->92), mult. (4303->201), div. (668->14), fcn. (5516->11), ass. (0->91)
	t179 = sin(qJ(2));
	t180 = sin(qJ(1));
	t231 = cos(pkin(5));
	t202 = t180 * t231;
	t200 = t179 * t202;
	t181 = cos(qJ(2));
	t182 = cos(qJ(1));
	t216 = t182 * t181;
	t163 = -t200 + t216;
	t174 = qJ(3) + pkin(10);
	t170 = sin(t174);
	t171 = cos(t174);
	t178 = sin(pkin(5));
	t220 = t178 * t180;
	t192 = -t163 * t170 + t171 * t220;
	t233 = t192 * qJD(3);
	t201 = t182 * t231;
	t199 = t181 * t201;
	t217 = t180 * t179;
	t159 = -t199 + t217;
	t219 = t178 * t181;
	t153 = atan2(-t159, -t219);
	t151 = sin(t153);
	t152 = cos(t153);
	t157 = t159 ^ 2;
	t173 = 0.1e1 / t178 ^ 2;
	t176 = 0.1e1 / t181 ^ 2;
	t156 = t157 * t173 * t176 + 0.1e1;
	t154 = 0.1e1 / t156;
	t172 = 0.1e1 / t178;
	t175 = 0.1e1 / t181;
	t206 = t159 * t172 * t175;
	t232 = (t152 * t206 - t151) * t154 + t151;
	t135 = -t151 * t159 - t152 * t219;
	t132 = 0.1e1 / t135;
	t150 = t163 * t171 + t170 * t220;
	t144 = 0.1e1 / t150;
	t133 = 0.1e1 / t135 ^ 2;
	t145 = 0.1e1 / t150 ^ 2;
	t189 = -t179 * t201 - t180 * t181;
	t190 = -t182 * t179 - t181 * t202;
	t141 = -t190 * qJD(1) - t189 * qJD(2);
	t214 = qJD(2) * t179;
	t203 = t176 * t214;
	t191 = t141 * t175 + t159 * t203;
	t222 = t154 * t172;
	t124 = t191 * t222;
	t195 = t151 * t219 - t152 * t159;
	t207 = t152 * t178 * t179;
	t120 = qJD(2) * t207 + t195 * t124 - t151 * t141;
	t230 = t120 * t132 * t133;
	t140 = t189 * qJD(1) + t190 * qJD(2);
	t215 = qJD(1) * t178;
	t204 = t182 * t215;
	t129 = t150 * qJD(3) + t140 * t170 - t171 * t204;
	t143 = t192 ^ 2;
	t138 = t143 * t145 + 0.1e1;
	t225 = t145 * t192;
	t130 = t140 * t171 + t170 * t204 + t233;
	t228 = t130 * t144 * t145;
	t229 = (-t129 * t225 - t143 * t228) / t138 ^ 2;
	t177 = t175 * t176;
	t227 = (t141 * t159 * t176 + t157 * t177 * t214) * t173 / t156 ^ 2;
	t226 = t133 * t190;
	t224 = t151 * t190;
	t223 = t152 * t190;
	t221 = t176 * t179;
	t218 = t178 * t182;
	t213 = qJD(2) * t181;
	t158 = t190 ^ 2;
	t128 = t158 * t133 + 0.1e1;
	t198 = qJD(2) * t231 + qJD(1);
	t139 = -qJD(1) * t199 - t182 * t213 + t198 * t217;
	t212 = 0.2e1 * (t139 * t226 - t158 * t230) / t128 ^ 2;
	t211 = 0.2e1 * t230;
	t210 = 0.2e1 * t229;
	t209 = -0.2e1 * t227;
	t208 = t192 * t228;
	t205 = t180 * t215;
	t196 = t170 * t144 + t171 * t225;
	t194 = t159 * t221 - t175 * t189;
	t193 = -t170 * t189 + t171 * t218;
	t148 = t170 * t218 + t171 * t189;
	t142 = -qJD(1) * t200 - t180 * t214 + t198 * t216;
	t136 = 0.1e1 / t138;
	t126 = 0.1e1 / t128;
	t125 = t194 * t222;
	t123 = t232 * t190;
	t121 = t195 * t125 + t151 * t189 + t207;
	t119 = (t194 * t209 + (t141 * t221 + t142 * t175 + (-t189 * t221 + (0.2e1 * t177 * t179 ^ 2 + t175) * t159) * qJD(2)) * t154) * t172;
	t1 = [(-t190 * t175 * t209 + (-t139 * t175 - t190 * t203) * t154) * t172, t119, 0, 0, 0; t159 * t132 * t212 + (-t141 * t132 + (t120 * t159 + t123 * t139) * t133) * t126 - (t123 * t211 * t126 + (t123 * t212 + ((t124 * t154 * t206 + t209) * t224 + (0.2e1 * t206 * t227 - t124 + (-t191 * t172 + t124) * t154) * t223 - t232 * t139) * t126) * t133) * t190, (-t121 * t226 - t132 * t163) * t212 + (-t121 * t190 * t211 + t140 * t132 + (-t163 * t120 + t121 * t139 + (t178 * t213 - t119 * t159 - t125 * t141 + (t125 * t219 + t189) * t124) * t223 + (t124 * t125 * t159 - t142 + (t119 * t181 + (-qJD(2) * t125 - t124) * t179) * t178) * t224) * t133) * t126, 0, 0, 0; (t144 * t193 - t148 * t225) * t210 + ((t148 * qJD(3) - t142 * t170 + t171 * t205) * t144 - 0.2e1 * t148 * t208 + (t193 * t130 + (t193 * qJD(3) - t142 * t171 - t170 * t205) * t192 - t148 * t129) * t145) * t136, -t196 * t190 * t210 + (t196 * t139 - ((-qJD(3) * t144 + 0.2e1 * t208) * t171 + (t129 * t171 + (t130 + t233) * t170) * t145) * t190) * t136, -0.2e1 * t229 - 0.2e1 * (t129 * t145 * t136 - (-t136 * t228 - t145 * t229) * t192) * t192, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:31:29
	% EndTime: 2019-12-31 21:31:30
	% DurationCPUTime: 1.47s
	% Computational Cost: add. (8428->149), mult. (13478->297), div. (726->12), fcn. (17045->13), ass. (0->128)
	t271 = cos(pkin(5));
	t273 = sin(qJ(2));
	t345 = sin(qJ(1));
	t306 = t345 * t273;
	t296 = t271 * t306;
	t300 = t345 * qJD(2);
	t275 = cos(qJ(2));
	t276 = cos(qJ(1));
	t322 = t276 * t275;
	t270 = sin(pkin(5));
	t326 = t270 * t276;
	t350 = -qJD(1) * t296 - t273 * t300 + (qJD(2) * t271 + qJD(1)) * t322 - qJD(3) * t326;
	t269 = qJ(3) + pkin(10);
	t267 = sin(t269);
	t268 = cos(t269);
	t305 = t345 * t275;
	t323 = t276 * t273;
	t288 = -t271 * t323 - t305;
	t240 = -t267 * t288 + t268 * t326;
	t328 = t270 * t273;
	t250 = t267 * t328 - t271 * t268;
	t229 = atan2(-t240, t250);
	t216 = sin(t229);
	t217 = cos(t229);
	t207 = -t216 * t240 + t217 * t250;
	t205 = 0.1e1 / t207 ^ 2;
	t256 = -t296 + t322;
	t307 = t270 * t345;
	t287 = -t256 * t267 + t268 * t307;
	t239 = t287 ^ 2;
	t201 = t239 * t205 + 0.1e1;
	t286 = -t271 * t305 - t323;
	t233 = t288 * qJD(1) + t286 * qJD(2);
	t246 = t256 * t268 + t267 * t307;
	t304 = qJD(1) * t326;
	t211 = t246 * qJD(3) + t233 * t267 - t268 * t304;
	t338 = t211 * t205;
	t238 = t240 ^ 2;
	t248 = 0.1e1 / t250 ^ 2;
	t228 = t238 * t248 + 0.1e1;
	t222 = 0.1e1 / t228;
	t301 = t345 * qJD(1);
	t295 = t270 * t301;
	t319 = qJD(3) * t268;
	t213 = t350 * t267 - t268 * t295 - t288 * t319;
	t251 = t271 * t267 + t268 * t328;
	t320 = qJD(2) * t275;
	t303 = t270 * t320;
	t236 = t251 * qJD(3) + t267 * t303;
	t247 = 0.1e1 / t250;
	t332 = t240 * t248;
	t292 = -t213 * t247 + t236 * t332;
	t195 = t292 * t222;
	t293 = -t216 * t250 - t217 * t240;
	t190 = t293 * t195 - t216 * t213 + t217 * t236;
	t204 = 0.1e1 / t207;
	t206 = t204 * t205;
	t343 = t190 * t206;
	t317 = 0.2e1 * (-t239 * t343 - t287 * t338) / t201 ^ 2;
	t349 = t236 * t248;
	t308 = t271 * t322;
	t253 = -t306 + t308;
	t327 = t270 * t275;
	t289 = -t247 * t253 + t327 * t332;
	t348 = t267 * t289;
	t214 = (qJD(3) * t288 + t295) * t267 + t350 * t268;
	t274 = cos(qJ(5));
	t272 = sin(qJ(5));
	t330 = t286 * t272;
	t227 = t246 * t274 - t330;
	t219 = 0.1e1 / t227;
	t220 = 0.1e1 / t227 ^ 2;
	t347 = -0.2e1 * t240;
	t346 = -0.2e1 * t287;
	t212 = t287 * qJD(3) + t233 * t268 + t267 * t304;
	t232 = -qJD(1) * t308 - t276 * t320 + (t271 * t300 + t301) * t273;
	t202 = t227 * qJD(5) + t212 * t272 + t232 * t274;
	t329 = t286 * t274;
	t226 = t246 * t272 + t329;
	t218 = t226 ^ 2;
	t210 = t218 * t220 + 0.1e1;
	t335 = t220 * t226;
	t318 = qJD(5) * t226;
	t203 = t212 * t274 - t232 * t272 - t318;
	t340 = t203 * t219 * t220;
	t342 = (t202 * t335 - t218 * t340) / t210 ^ 2;
	t334 = t247 * t349;
	t341 = (t213 * t332 - t238 * t334) / t228 ^ 2;
	t339 = t205 * t287;
	t337 = t216 * t287;
	t336 = t217 * t287;
	t333 = t240 * t247;
	t331 = t286 * t267;
	t325 = t272 * t219;
	t324 = t274 * t226;
	t321 = qJD(2) * t273;
	t316 = -0.2e1 * t342;
	t315 = 0.2e1 * t342;
	t314 = -0.2e1 * t341;
	t313 = t206 * t346;
	t312 = t247 * t341;
	t311 = t205 * t337;
	t310 = t205 * t336;
	t309 = t226 * t340;
	t299 = 0.2e1 * t309;
	t298 = t334 * t347;
	t242 = -t267 * t326 - t268 * t288;
	t294 = -qJD(5) * t268 * t286 + t233;
	t225 = -t242 * t274 + t253 * t272;
	t224 = -t242 * t272 - t253 * t274;
	t291 = t220 * t324 - t325;
	t290 = -t242 * t247 + t251 * t332;
	t284 = -t216 + (t217 * t333 + t216) * t222;
	t283 = -qJD(3) * t331 + qJD(5) * t256 + t232 * t268;
	t237 = -t250 * qJD(3) + t268 * t303;
	t234 = t286 * qJD(1) + t288 * qJD(2);
	t231 = t256 * t272 + t268 * t329;
	t230 = -t256 * t274 + t268 * t330;
	t208 = 0.1e1 / t210;
	t199 = 0.1e1 / t201;
	t198 = t222 * t348;
	t196 = t290 * t222;
	t194 = t284 * t287;
	t192 = (-t216 * t253 + t217 * t327) * t267 + t293 * t198;
	t191 = t293 * t196 - t216 * t242 + t217 * t251;
	t189 = t290 * t314 + (t251 * t298 - t214 * t247 + (t213 * t251 + t236 * t242 + t237 * t240) * t248) * t222;
	t187 = t314 * t348 + (t289 * t319 + (t298 * t327 - t234 * t247 + (t236 * t253 + (t213 * t275 - t240 * t321) * t270) * t248) * t267) * t222;
	t1 = [t312 * t346 + (-t211 * t247 - t287 * t349) * t222, t187, t189, 0, 0; t240 * t204 * t317 + (-t213 * t204 + (t190 * t240 + t194 * t211) * t205) * t199 - (-t194 * t205 * t317 + (-0.2e1 * t194 * t343 + (-t195 * t222 * t333 + t314) * t311 + (t312 * t347 - t195 + (t195 - t292) * t222) * t310 - t284 * t338) * t199) * t287, (-t192 * t339 - t204 * t331) * t317 + (-t192 * t338 + (t232 * t267 + t286 * t319) * t204 + (t192 * t313 - t205 * t331) * t190 + (-t187 * t240 - t198 * t213 + (-t267 * t321 + t275 * t319) * t270 + (-t198 * t250 - t253 * t267) * t195) * t310 + (-t253 * t319 - t187 * t250 - t198 * t236 - t234 * t267 + (t198 * t240 - t267 * t327) * t195) * t311) * t199, (-t191 * t339 - t204 * t246) * t317 + (t191 * t190 * t313 + t212 * t204 + (-t246 * t190 - t191 * t211 + (-t189 * t240 - t196 * t213 + t237 + (-t196 * t250 - t242) * t195) * t336 + (-t189 * t250 - t196 * t236 - t214 + (t196 * t240 - t251) * t195) * t337) * t205) * t199, 0, 0; (-t219 * t224 + t225 * t335) * t315 + ((t225 * qJD(5) - t214 * t272 - t234 * t274) * t219 + t225 * t299 + (-t224 * t203 - (-t224 * qJD(5) - t214 * t274 + t234 * t272) * t226 - t225 * t202) * t220) * t208, (-t219 * t230 + t231 * t335) * t315 + (t231 * t299 - t294 * t219 * t274 + t283 * t325 + (-t294 * t226 * t272 - t231 * t202 - t230 * t203 - t283 * t324) * t220) * t208, -t291 * t287 * t316 + (t291 * t211 - ((-qJD(5) * t219 - 0.2e1 * t309) * t274 + (t202 * t274 + (t203 - t318) * t272) * t220) * t287) * t208, 0, t316 + 0.2e1 * (t202 * t220 * t208 + (-t208 * t340 - t220 * t342) * t226) * t226;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end