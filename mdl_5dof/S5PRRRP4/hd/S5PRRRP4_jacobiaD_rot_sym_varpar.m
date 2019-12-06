% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PRRRP4
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
%   Wie in S5PRRRP4_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5PRRRP4_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP4_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP4_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRP4_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP4_jacobiaD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:46:43
	% EndTime: 2019-12-05 16:46:43
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:46:43
	% EndTime: 2019-12-05 16:46:43
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:46:43
	% EndTime: 2019-12-05 16:46:43
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:46:43
	% EndTime: 2019-12-05 16:46:43
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:46:43
	% EndTime: 2019-12-05 16:46:44
	% DurationCPUTime: 0.39s
	% Computational Cost: add. (1992->45), mult. (2152->102), div. (440->12), fcn. (2434->9), ass. (0->60)
	t139 = qJ(2) + qJ(3);
	t134 = sin(t139);
	t130 = t134 ^ 2;
	t135 = cos(t139);
	t165 = t130 / t135 ^ 2;
	t155 = 0.1e1 + t165;
	t138 = qJD(2) + qJD(3);
	t140 = sin(pkin(8));
	t175 = t138 * t140;
	t142 = sin(qJ(4));
	t143 = cos(qJ(4));
	t141 = cos(pkin(8));
	t163 = t135 * t141;
	t119 = t140 * t142 + t143 * t163;
	t136 = t140 ^ 2;
	t125 = t136 * t165 + 0.1e1;
	t123 = 0.1e1 / t125;
	t107 = t155 * t140 * t123;
	t162 = t140 * t134;
	t122 = atan2(-t162, -t135);
	t121 = cos(t122);
	t120 = sin(t122);
	t166 = t120 * t135;
	t152 = t121 * t134 - t140 * t166;
	t158 = t121 * t162;
	t153 = -t158 + t166;
	t96 = t153 * t107 + t152;
	t174 = 0.2e1 * t96;
	t106 = -t120 * t162 - t121 * t135;
	t103 = 0.1e1 / t106;
	t115 = 0.1e1 / t119;
	t104 = 0.1e1 / t106 ^ 2;
	t116 = 0.1e1 / t119 ^ 2;
	t137 = t141 ^ 2;
	t101 = t137 * t130 * t104 + 0.1e1;
	t164 = t135 * t138;
	t169 = t104 * t134;
	t102 = t107 * t138;
	t95 = t153 * t102 + t152 * t138;
	t171 = t103 * t104 * t95;
	t173 = 0.1e1 / t101 ^ 2 * (-t130 * t171 + t164 * t169) * t137;
	t99 = 0.1e1 / t101;
	t172 = t104 * t99;
	t118 = -t140 * t143 + t142 * t163;
	t114 = t118 ^ 2;
	t110 = t114 * t116 + 0.1e1;
	t157 = t134 * t138 * t141;
	t111 = -t119 * qJD(4) + t142 * t157;
	t167 = t116 * t118;
	t160 = qJD(4) * t118;
	t112 = -t143 * t157 - t160;
	t168 = t112 * t115 * t116;
	t170 = 0.1e1 / t110 ^ 2 * (-t111 * t167 - t114 * t168);
	t159 = -0.2e1 * t170;
	t154 = -t115 * t142 + t143 * t167;
	t108 = 0.1e1 / t110;
	t97 = 0.2e1 * (t123 - t155 / t125 ^ 2 * t136) * t155 / t135 * t134 * t175;
	t93 = (t154 * t134 * t159 + (t154 * t164 + ((-qJD(4) * t115 - 0.2e1 * t118 * t168) * t143 + (-t111 * t143 + (t112 - t160) * t142) * t116) * t134) * t108) * t141;
	t92 = ((-0.2e1 * t103 * t173 + (-t138 * t96 - t95) * t172) * t135 + (t104 * t173 * t174 + (t97 * t104 * t158 + t171 * t174 - t138 * t103 - (t175 - t102 + (t102 * t140 - t138) * t107) * t120 * t169) * t99 - (t120 * t97 + ((-t107 * t140 + 0.1e1) * t138 + (t107 - t140) * t102) * t121) * t135 * t172) * t134) * t141;
	t1 = [0, t97, t97, 0, 0; 0, t92, t92, 0, 0; 0, t93, t93, t159 + 0.2e1 * (-t108 * t111 * t116 + (-t108 * t168 - t116 * t170) * t118) * t118, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:46:43
	% EndTime: 2019-12-05 16:46:44
	% DurationCPUTime: 0.73s
	% Computational Cost: add. (3924->85), mult. (5873->199), div. (1202->15), fcn. (7475->9), ass. (0->98)
	t175 = qJ(2) + qJ(3);
	t169 = cos(t175);
	t177 = cos(pkin(8));
	t179 = cos(qJ(4));
	t208 = t177 * t179;
	t176 = sin(pkin(8));
	t178 = sin(qJ(4));
	t212 = t176 * t178;
	t158 = t169 * t208 + t212;
	t152 = 0.1e1 / t158 ^ 2;
	t168 = sin(t175);
	t164 = t168 ^ 2;
	t170 = t177 ^ 2;
	t220 = t164 * t170;
	t194 = t152 * t220;
	t147 = 0.1e1 + t194;
	t209 = t177 * t178;
	t211 = t176 * t179;
	t157 = t169 * t209 - t211;
	t171 = qJD(2) + qJD(3);
	t214 = t171 * t179;
	t195 = t168 * t214;
	t139 = -t157 * qJD(4) - t177 * t195;
	t151 = 0.1e1 / t158;
	t227 = t139 * t151 * t152;
	t200 = t164 * t227;
	t215 = t169 * t171;
	t233 = (t152 * t168 * t215 - t200) * t170 / t147 ^ 2;
	t165 = 0.1e1 / t168;
	t154 = t169 * t212 + t208;
	t156 = t169 * t211 - t209;
	t172 = 0.1e1 / t178;
	t173 = 0.1e1 / t178 ^ 2;
	t213 = t173 * t179;
	t189 = t154 * t213 - t156 * t172;
	t232 = t165 * t189;
	t205 = qJD(4) * t179;
	t191 = t169 * t205;
	t217 = t168 * t176;
	t136 = -t176 * t191 + (qJD(4) * t177 + t171 * t217) * t178;
	t149 = t154 ^ 2;
	t166 = 0.1e1 / t168 ^ 2;
	t218 = t166 * t173;
	t148 = t149 * t218 + 0.1e1;
	t174 = t172 * t173;
	t192 = t166 * t205;
	t197 = t165 / t164 * t215;
	t199 = t154 * t218;
	t231 = -0.2e1 * (-t136 * t199 + (-t173 * t197 - t174 * t192) * t149) / t148 ^ 2;
	t216 = t168 * t178;
	t146 = atan2(-t154, t216);
	t141 = cos(t146);
	t140 = sin(t146);
	t226 = t140 * t154;
	t135 = t141 * t216 - t226;
	t132 = 0.1e1 / t135;
	t133 = 0.1e1 / t135 ^ 2;
	t230 = 0.2e1 * t157;
	t229 = t133 * t157;
	t196 = t168 * t171 * t177;
	t206 = qJD(4) * t178;
	t138 = -t176 * t206 - t177 * t191 + t178 * t196;
	t228 = t138 * t133;
	t225 = t140 * t157;
	t224 = t140 * t168;
	t223 = t141 * t154;
	t222 = t141 * t157;
	t221 = t141 * t169;
	t219 = t165 * t172;
	t210 = t177 * t132;
	t144 = 0.1e1 / t148;
	t198 = t166 * t169 * t172;
	t190 = t154 * t198 + t176;
	t131 = t190 * t144;
	t207 = -t131 + t176;
	t188 = t168 * t205 + t178 * t215;
	t124 = (t136 * t219 + t188 * t199) * t144;
	t187 = -t124 * t154 + t188;
	t121 = (-t124 * t216 + t136) * t140 + t187 * t141;
	t150 = t157 ^ 2;
	t129 = t150 * t133 + 0.1e1;
	t134 = t132 * t133;
	t204 = 0.2e1 * (-t150 * t134 * t121 - t157 * t228) / t129 ^ 2;
	t203 = t134 * t230;
	t202 = t168 * t230;
	t201 = t133 * t225;
	t193 = t168 * t210;
	t142 = 0.1e1 / t147;
	t137 = qJD(4) * t154 + t176 * t195;
	t127 = 0.1e1 / t129;
	t126 = t144 * t232;
	t123 = -t131 * t223 + (t207 * t224 + t221) * t178;
	t122 = t141 * t168 * t179 - t140 * t156 + (-t140 * t216 - t223) * t126;
	t120 = t190 * t231 + (-t136 * t198 + (-t171 * t219 + (-0.2e1 * t172 * t197 - t173 * t192) * t169) * t154) * t144;
	t119 = 0.2e1 * (t151 * t169 * t177 + t179 * t194) * t233 + (0.2e1 * t170 * t179 * t200 + t151 * t196 + (t206 * t220 + (t139 * t177 - 0.2e1 * t170 * t195) * t169) * t152) * t142;
	t117 = t231 * t232 + (-t189 * t166 * t215 + (-t136 * t213 + t137 * t172 + (t156 * t213 + (-0.2e1 * t174 * t179 ^ 2 - t172) * t154) * qJD(4)) * t165) * t144;
	t116 = t123 * t204 * t229 + (-(-t120 * t223 + (t124 * t226 + t136 * t141) * t131) * t229 + (t121 * t203 + t228) * t123 + (-t193 - (-t131 * t224 + t140 * t217 + t221) * t229) * t205) * t127 + (t193 * t204 + ((-t171 * t210 - (t207 * t171 - t124) * t201) * t169 + (t120 * t225 + t177 * t121 - (t207 * t124 - t171) * t222) * t133 * t168) * t127) * t178;
	t1 = [0, t120, t120, t117, 0; 0, t116, t116, (t122 * t229 - t132 * t158) * t204 + (t122 * t228 + t139 * t132 + (t122 * t203 - t158 * t133) * t121 - (-t168 * t206 + t169 * t214 - t117 * t154 + t126 * t136 + (-t126 * t216 - t156) * t124) * t133 * t222 - (t137 + (-t117 * t178 - t124 * t179) * t168 - t187 * t126) * t201) * t127, 0; 0, t119, t119, (t152 * t202 * t233 + (t202 * t227 + (t138 * t168 - t157 * t215) * t152) * t142) * t177, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end