% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRRPR8
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
%   Wie in S5RRRPR8_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 20:08
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RRRPR8_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR8_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR8_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPR8_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR8_jacobiaD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:07:57
	% EndTime: 2019-12-29 20:07:57
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:07:57
	% EndTime: 2019-12-29 20:07:57
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:07:57
	% EndTime: 2019-12-29 20:07:57
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:07:57
	% EndTime: 2019-12-29 20:07:57
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:07:58
	% EndTime: 2019-12-29 20:07:59
	% DurationCPUTime: 1.28s
	% Computational Cost: add. (3078->72), mult. (2858->158), div. (686->14), fcn. (3330->7), ass. (0->75)
	t122 = sin(qJ(1));
	t115 = t122 ^ 2;
	t121 = qJ(2) + qJ(3);
	t112 = sin(t121);
	t107 = t112 ^ 2;
	t113 = cos(t121);
	t110 = 0.1e1 / t113 ^ 2;
	t156 = t107 * t110;
	t102 = t115 * t156 + 0.1e1;
	t106 = t112 * t107;
	t108 = t113 ^ 2;
	t109 = 0.1e1 / t113;
	t114 = qJD(2) + qJD(3);
	t155 = t109 * t112;
	t131 = t114 * (t106 * t109 / t108 + t155);
	t123 = cos(qJ(1));
	t147 = qJD(1) * t123;
	t139 = t122 * t147;
	t164 = 0.1e1 / t102 ^ 2 * (t115 * t131 + t139 * t156);
	t172 = -0.2e1 * t164;
	t100 = 0.1e1 / t102;
	t137 = 0.1e1 + t156;
	t170 = t122 * t137;
	t95 = t100 * t170;
	t171 = t122 * t95 - 0.1e1;
	t151 = 0.1e1 / t122 * t123;
	t120 = t123 ^ 2;
	t169 = qJD(1) * (0.1e1 / t115 * t120 + 0.1e1) * t151;
	t149 = t122 * t112;
	t99 = atan2(-t149, -t113);
	t97 = sin(t99);
	t143 = t97 * t149;
	t98 = cos(t99);
	t94 = -t113 * t98 - t143;
	t91 = 0.1e1 / t94;
	t92 = 0.1e1 / t94 ^ 2;
	t153 = t113 * t114;
	t135 = t112 * t120 * t153;
	t152 = t114 * t122;
	t161 = t113 * t97;
	t140 = t110 * t152;
	t86 = (-(-t112 * t147 - t113 * t152) * t109 + t107 * t140) * t100;
	t81 = (t86 - t152) * t161 + (-t97 * t147 + (-t122 * t86 + t114) * t98) * t112;
	t167 = t81 * t91 * t92;
	t89 = t107 * t120 * t92 + 0.1e1;
	t168 = (t92 * t135 + (-t120 * t167 - t139 * t92) * t107) / t89 ^ 2;
	t87 = 0.1e1 / t89;
	t165 = t87 * t92;
	t163 = t112 * t97;
	t162 = t112 * t98;
	t160 = t114 * t95;
	t158 = t123 * t92;
	t157 = t107 * t109;
	t154 = t112 * t123;
	t117 = 0.1e1 / t122 ^ 2;
	t150 = t117 * t120;
	t148 = qJD(1) * t122;
	t146 = 0.2e1 * t167;
	t105 = t108 * t150 + 0.1e1;
	t145 = 0.2e1 / t105 ^ 2 * (-t108 * t169 - t117 * t135);
	t144 = t91 * t168;
	t142 = t87 * t153;
	t141 = t122 * t157;
	t138 = 0.2e1 * t92 * t168;
	t136 = 0.1e1 + t150;
	t134 = t137 * t123;
	t133 = t136 * t112;
	t130 = -t141 * t98 + t163;
	t103 = 0.1e1 / t105;
	t85 = (t100 * t130 - t163) * t123;
	t84 = t112 * t145 * t151 + (qJD(1) * t133 - t151 * t153) * t103;
	t83 = (-t122 + t95) * t161 - t171 * t162;
	t82 = t170 * t172 + (qJD(1) * t134 + 0.2e1 * t122 * t131) * t100;
	t79 = (-t91 * t87 * t148 + (-0.2e1 * t144 + (-t114 * t83 - t81) * t165) * t123) * t113 + (t83 * t123 * t138 + (-t123 * t114 * t91 - ((-t122 * t82 - t147 * t95) * t98 + (t171 * t86 + t152 - t160) * t97) * t92 * t154 + (t123 * t146 + t148 * t92) * t83 - ((t82 - t147) * t97 + (t86 * t95 + t114 + (-t86 - t160) * t122) * t98) * t113 * t158) * t87) * t112;
	t1 = [t109 * t154 * t172 + (t114 * t134 - t148 * t155) * t100, t82, t82, 0, 0; (-t91 * t142 + (0.2e1 * t144 + (qJD(1) * t85 + t81) * t165) * t112) * t122 + (-t85 * t92 * t142 + (t85 * t138 + (t85 * t146 + (t86 * t162 + t97 * t153 + 0.2e1 * t130 * t164 + ((-t141 * t86 - t153) * t97 + (t106 * t140 - (t86 - 0.2e1 * t152) * t112) * t98) * t100) * t158) * t87) * t112 + (-t91 + (-t143 + (t143 - (t115 - t120) * t98 * t157) * t100) * t92) * t112 * t87 * qJD(1)) * t123, t79, t79, 0, 0; t136 * t113 * t145 + (0.2e1 * t113 * t169 + t114 * t133) * t103, t84, t84, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:08:03
	% EndTime: 2019-12-29 20:08:04
	% DurationCPUTime: 1.62s
	% Computational Cost: add. (3360->95), mult. (3810->207), div. (753->12), fcn. (4455->9), ass. (0->95)
	t174 = sin(qJ(1));
	t172 = qJ(2) + qJ(3);
	t167 = sin(t172);
	t163 = 0.1e1 / t167 ^ 2;
	t168 = cos(t172);
	t166 = t168 ^ 2;
	t220 = t163 * t166;
	t195 = 0.1e1 + t220;
	t235 = t174 * t195;
	t170 = t174 ^ 2;
	t160 = t170 * t220 + 0.1e1;
	t158 = 0.1e1 / t160;
	t162 = 0.1e1 / t167;
	t176 = cos(qJ(1));
	t207 = qJD(1) * t176;
	t196 = t168 * t207;
	t169 = qJD(2) + qJD(3);
	t215 = t169 * t174;
	t198 = t163 * t215;
	t132 = ((t167 * t215 - t196) * t162 + t166 * t198) * t158;
	t234 = -t132 + t215;
	t191 = qJD(1) * t167 + qJD(5);
	t214 = t169 * t176;
	t233 = -t168 * t214 + t191 * t174;
	t213 = t174 * t168;
	t157 = atan2(-t213, t167);
	t156 = cos(t157);
	t155 = sin(t157);
	t200 = t155 * t213;
	t142 = t156 * t167 - t200;
	t139 = 0.1e1 / t142;
	t173 = sin(qJ(5));
	t210 = t176 * t173;
	t175 = cos(qJ(5));
	t211 = t174 * t175;
	t152 = t167 * t210 + t211;
	t148 = 0.1e1 / t152;
	t140 = 0.1e1 / t142 ^ 2;
	t149 = 0.1e1 / t152 ^ 2;
	t232 = t158 - 0.1e1;
	t222 = t156 * t168;
	t127 = (-t132 * t174 + t169) * t222 + (t234 * t167 - t196) * t155;
	t231 = t127 * t139 * t140;
	t192 = qJD(5) * t167 + qJD(1);
	t187 = t192 * t176;
	t136 = t173 * t187 + t233 * t175;
	t209 = t176 * t175;
	t212 = t174 * t173;
	t151 = -t167 * t209 + t212;
	t147 = t151 ^ 2;
	t146 = t147 * t149 + 0.1e1;
	t225 = t149 * t151;
	t137 = -t233 * t173 + t175 * t187;
	t229 = t137 * t148 * t149;
	t230 = (t136 * t225 - t147 * t229) / t146 ^ 2;
	t165 = t168 * t166;
	t221 = t162 * t168;
	t185 = t169 * (-t162 * t163 * t165 - t221);
	t218 = t166 * t174;
	t189 = t207 * t218;
	t228 = (t163 * t189 + t170 * t185) / t160 ^ 2;
	t227 = t140 * t168;
	t226 = t140 * t176;
	t224 = t151 * t173;
	t223 = t155 * t174;
	t171 = t176 ^ 2;
	t219 = t166 * t171;
	t217 = t167 * t169;
	t216 = t168 * t169;
	t208 = qJD(1) * t174;
	t135 = t140 * t219 + 0.1e1;
	t206 = 0.2e1 * (-t219 * t231 + (-t167 * t171 * t216 - t189) * t140) / t135 ^ 2;
	t205 = 0.2e1 * t231;
	t204 = 0.2e1 * t230;
	t203 = -0.2e1 * t228;
	t202 = t168 * t228;
	t201 = t168 * t226;
	t199 = t162 * t218;
	t194 = t168 * t206;
	t193 = 0.2e1 * t151 * t229;
	t190 = t156 * t158 * t162 * t166;
	t188 = t195 * t176;
	t186 = t148 * t175 + t149 * t224;
	t184 = t186 * t176;
	t154 = -t167 * t212 + t209;
	t153 = t167 * t211 + t210;
	t144 = 0.1e1 / t146;
	t143 = t158 * t235;
	t133 = 0.1e1 / t135;
	t131 = (t232 * t168 * t155 + t174 * t190) * t176;
	t129 = t167 * t223 + t222 + (-t155 * t167 - t156 * t213) * t143;
	t128 = t203 * t235 + (qJD(1) * t188 + 0.2e1 * t174 * t185) * t158;
	t125 = t168 * t184 * t204 + (t184 * t217 + (t186 * t208 + ((qJD(5) * t148 + t193) * t173 + (-t136 * t173 + (-qJD(5) * t151 + t137) * t175) * t149) * t176) * t168) * t144;
	t124 = (t129 * t227 + t139 * t167) * t176 * t206 + ((t139 * t208 + (t129 * t169 + t127) * t226) * t167 + (-t139 * t214 - (-t128 * t156 * t174 + t234 * t155 + (t132 * t223 - t155 * t169 - t156 * t207) * t143) * t201 + (t140 * t208 + t176 * t205) * t129 - ((-t128 + t207) * t155 + ((t143 * t174 - 0.1e1) * t169 + (-t143 + t174) * t132) * t156) * t167 * t226) * t168) * t133;
	t1 = [0.2e1 * t176 * t162 * t202 + (t169 * t188 + t208 * t221) * t158, t128, t128, 0, 0; (t139 * t194 + (t139 * t217 + (qJD(1) * t131 + t127) * t227) * t133) * t174 + (t140 * t194 * t131 + (-((-0.2e1 * t202 + t217 + (-t132 * t199 - t217) * t158) * t155 + (t199 * t203 - t132 * t168 + (-t165 * t198 + (t132 - 0.2e1 * t215) * t168) * t158) * t156) * t201 + (t140 * t217 + t168 * t205) * t131 + (-t139 + ((t170 - t171) * t190 + t232 * t200) * t140) * t168 * qJD(1)) * t133) * t176, t124, t124, 0, 0; (-t148 * t153 + t154 * t225) * t204 + (t154 * t193 + (-t154 * t136 - t153 * t137 + t192 * t151 * t211 - (-t169 * t213 - t191 * t176) * t224) * t149 + (t191 * t209 + (-t192 * t173 + t175 * t216) * t174) * t148) * t144, t125, t125, 0, -0.2e1 * t230 + 0.2e1 * (t136 * t149 * t144 + (-t144 * t229 - t149 * t230) * t151) * t151;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end