% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRRP4
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
%   Wie in S6RPRRRP4_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:50
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRRP4_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP4_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP4_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRP4_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP4_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:50:39
	% EndTime: 2019-10-10 01:50:39
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:50:39
	% EndTime: 2019-10-10 01:50:39
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:50:39
	% EndTime: 2019-10-10 01:50:39
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:50:39
	% EndTime: 2019-10-10 01:50:39
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:50:39
	% EndTime: 2019-10-10 01:50:39
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:50:39
	% EndTime: 2019-10-10 01:50:40
	% DurationCPUTime: 1.13s
	% Computational Cost: add. (5587->97), mult. (3810->208), div. (753->12), fcn. (4455->9), ass. (0->95)
	t174 = sin(qJ(1));
	t236 = 0.2e1 * t174;
	t171 = t174 ^ 2;
	t169 = pkin(10) + qJ(3) + qJ(4);
	t166 = sin(t169);
	t162 = t166 ^ 2;
	t167 = cos(t169);
	t164 = 0.1e1 / t167 ^ 2;
	t221 = t162 * t164;
	t158 = t171 * t221 + 0.1e1;
	t156 = 0.1e1 / t158;
	t163 = 0.1e1 / t167;
	t176 = cos(qJ(1));
	t207 = qJD(1) * t176;
	t197 = t166 * t207;
	t170 = qJD(3) + qJD(4);
	t215 = t170 * t174;
	t200 = t164 * t215;
	t130 = (-(-t167 * t215 - t197) * t163 + t162 * t200) * t156;
	t235 = t130 - t215;
	t175 = cos(qJ(5));
	t209 = t175 * t176;
	t173 = sin(qJ(5));
	t211 = t174 * t173;
	t155 = t167 * t209 + t211;
	t212 = t174 * t166;
	t151 = atan2(-t212, -t167);
	t146 = cos(t151);
	t145 = sin(t151);
	t202 = t145 * t212;
	t140 = -t146 * t167 - t202;
	t137 = 0.1e1 / t140;
	t148 = 0.1e1 / t155;
	t138 = 0.1e1 / t140 ^ 2;
	t149 = 0.1e1 / t155 ^ 2;
	t234 = t156 - 0.1e1;
	t226 = t146 * t166;
	t125 = (-t130 * t174 + t170) * t226 + (t235 * t167 - t197) * t145;
	t233 = t125 * t137 * t138;
	t185 = t167 * t211 + t209;
	t214 = t170 * t176;
	t198 = t166 * t214;
	t135 = t185 * qJD(1) - t155 * qJD(5) + t173 * t198;
	t210 = t174 * t175;
	t213 = t173 * t176;
	t154 = t167 * t213 - t210;
	t147 = t154 ^ 2;
	t144 = t147 * t149 + 0.1e1;
	t224 = t149 * t154;
	t191 = -qJD(1) * t167 + qJD(5);
	t192 = qJD(5) * t167 - qJD(1);
	t136 = -t192 * t213 + (t191 * t174 - t198) * t175;
	t230 = t136 * t148 * t149;
	t232 = (-t135 * t224 - t147 * t230) / t144 ^ 2;
	t161 = t166 * t162;
	t218 = t163 * t166;
	t184 = t170 * (t161 * t163 * t164 + t218);
	t219 = t162 * t174;
	t189 = t207 * t219;
	t231 = (t164 * t189 + t171 * t184) / t158 ^ 2;
	t229 = t138 * t166;
	t228 = t138 * t176;
	t227 = t145 * t174;
	t225 = t148 * t173;
	t223 = t154 * t175;
	t222 = t162 * t163;
	t172 = t176 ^ 2;
	t220 = t162 * t172;
	t217 = t166 * t176;
	t216 = t167 * t170;
	t208 = qJD(1) * t174;
	t133 = t138 * t220 + 0.1e1;
	t206 = 0.2e1 * (-t220 * t233 + (t166 * t172 * t216 - t189) * t138) / t133 ^ 2;
	t205 = 0.2e1 * t233;
	t204 = -0.2e1 * t232;
	t203 = t138 * t217;
	t201 = t154 * t230;
	t196 = 0.1e1 + t221;
	t195 = t166 * t206;
	t194 = -0.2e1 * t166 * t231;
	t193 = t231 * t236;
	t190 = t146 * t156 * t222;
	t188 = t196 * t176;
	t187 = t191 * t176;
	t186 = t149 * t223 - t225;
	t153 = -t167 * t210 + t213;
	t142 = 0.1e1 / t144;
	t141 = t196 * t174 * t156;
	t131 = 0.1e1 / t133;
	t129 = (t234 * t166 * t145 - t174 * t190) * t176;
	t127 = -t167 * t227 + t226 + (t145 * t167 - t146 * t212) * t141;
	t126 = -t196 * t193 + (qJD(1) * t188 + t184 * t236) * t156;
	t123 = t186 * t204 * t217 + (t186 * t167 * t214 + (-t186 * t208 + ((-qJD(5) * t148 - 0.2e1 * t201) * t175 + (-t135 * t175 + (-qJD(5) * t154 + t136) * t173) * t149) * t176) * t166) * t142;
	t122 = (t127 * t229 - t137 * t167) * t176 * t206 + ((-t137 * t208 + (-t127 * t170 - t125) * t228) * t167 + (-t137 * t214 - (-t126 * t146 * t174 - t235 * t145 + (t130 * t227 - t145 * t170 - t146 * t207) * t141) * t203 + (t138 * t208 + t176 * t205) * t127 - ((t126 - t207) * t145 + ((-t141 * t174 + 0.1e1) * t170 + (t141 - t174) * t130) * t146) * t167 * t228) * t166) * t131;
	t1 = [t163 * t176 * t194 + (t170 * t188 - t208 * t218) * t156, 0, t126, t126, 0, 0; (t137 * t195 + (-t137 * t216 + (qJD(1) * t129 + t125) * t229) * t131) * t174 + (t138 * t195 * t129 + (-((t194 - t216 + (t130 * t163 * t219 + t216) * t156) * t145 + (t193 * t222 - t130 * t166 + (-t161 * t200 + (t130 - 0.2e1 * t215) * t166) * t156) * t146) * t203 + (-t138 * t216 + t166 * t205) * t129 + (-t137 + ((-t171 + t172) * t190 + t234 * t202) * t138) * t166 * qJD(1)) * t131) * t176, 0, t122, t122, 0, 0; 0.2e1 * (t148 * t185 + t153 * t224) * t232 + (0.2e1 * t153 * t201 - t192 * t148 * t210 + (t170 * t212 + t187) * t225 + (t153 * t135 + t185 * t136 - t187 * t223 - (t166 * t170 * t175 + t192 * t173) * t154 * t174) * t149) * t142, 0, t123, t123, t204 + 0.2e1 * (-t135 * t142 * t149 + (-t142 * t230 - t149 * t232) * t154) * t154, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:50:39
	% EndTime: 2019-10-10 01:50:40
	% DurationCPUTime: 1.13s
	% Computational Cost: add. (5587->97), mult. (3810->208), div. (753->12), fcn. (4455->9), ass. (0->95)
	t172 = sin(qJ(1));
	t234 = 0.2e1 * t172;
	t169 = t172 ^ 2;
	t167 = pkin(10) + qJ(3) + qJ(4);
	t164 = sin(t167);
	t160 = t164 ^ 2;
	t165 = cos(t167);
	t162 = 0.1e1 / t165 ^ 2;
	t219 = t160 * t162;
	t156 = t169 * t219 + 0.1e1;
	t154 = 0.1e1 / t156;
	t161 = 0.1e1 / t165;
	t174 = cos(qJ(1));
	t205 = qJD(1) * t174;
	t195 = t164 * t205;
	t168 = qJD(3) + qJD(4);
	t213 = t168 * t172;
	t198 = t162 * t213;
	t128 = (-(-t165 * t213 - t195) * t161 + t160 * t198) * t154;
	t233 = t128 - t213;
	t173 = cos(qJ(5));
	t207 = t173 * t174;
	t171 = sin(qJ(5));
	t209 = t172 * t171;
	t153 = t165 * t207 + t209;
	t210 = t172 * t164;
	t149 = atan2(-t210, -t165);
	t144 = cos(t149);
	t143 = sin(t149);
	t200 = t143 * t210;
	t138 = -t144 * t165 - t200;
	t135 = 0.1e1 / t138;
	t146 = 0.1e1 / t153;
	t136 = 0.1e1 / t138 ^ 2;
	t147 = 0.1e1 / t153 ^ 2;
	t232 = t154 - 0.1e1;
	t224 = t144 * t164;
	t123 = (-t128 * t172 + t168) * t224 + (t233 * t165 - t195) * t143;
	t231 = t123 * t135 * t136;
	t183 = t165 * t209 + t207;
	t212 = t168 * t174;
	t196 = t164 * t212;
	t133 = t183 * qJD(1) - t153 * qJD(5) + t171 * t196;
	t208 = t172 * t173;
	t211 = t171 * t174;
	t152 = t165 * t211 - t208;
	t145 = t152 ^ 2;
	t142 = t145 * t147 + 0.1e1;
	t222 = t147 * t152;
	t189 = -qJD(1) * t165 + qJD(5);
	t190 = qJD(5) * t165 - qJD(1);
	t134 = -t190 * t211 + (t189 * t172 - t196) * t173;
	t228 = t134 * t146 * t147;
	t230 = (-t133 * t222 - t145 * t228) / t142 ^ 2;
	t159 = t164 * t160;
	t216 = t161 * t164;
	t182 = t168 * (t159 * t161 * t162 + t216);
	t217 = t160 * t172;
	t187 = t205 * t217;
	t229 = (t162 * t187 + t169 * t182) / t156 ^ 2;
	t227 = t136 * t164;
	t226 = t136 * t174;
	t225 = t143 * t172;
	t223 = t146 * t171;
	t221 = t152 * t173;
	t220 = t160 * t161;
	t170 = t174 ^ 2;
	t218 = t160 * t170;
	t215 = t164 * t174;
	t214 = t165 * t168;
	t206 = qJD(1) * t172;
	t131 = t136 * t218 + 0.1e1;
	t204 = 0.2e1 * (-t218 * t231 + (t164 * t170 * t214 - t187) * t136) / t131 ^ 2;
	t203 = 0.2e1 * t231;
	t202 = -0.2e1 * t230;
	t201 = t136 * t215;
	t199 = t152 * t228;
	t194 = 0.1e1 + t219;
	t193 = t164 * t204;
	t192 = -0.2e1 * t164 * t229;
	t191 = t229 * t234;
	t188 = t144 * t154 * t220;
	t186 = t194 * t174;
	t185 = t189 * t174;
	t184 = t147 * t221 - t223;
	t151 = -t165 * t208 + t211;
	t140 = 0.1e1 / t142;
	t139 = t194 * t172 * t154;
	t129 = 0.1e1 / t131;
	t127 = (t232 * t164 * t143 - t172 * t188) * t174;
	t125 = -t165 * t225 + t224 + (t143 * t165 - t144 * t210) * t139;
	t124 = -t194 * t191 + (qJD(1) * t186 + t182 * t234) * t154;
	t121 = t184 * t202 * t215 + (t184 * t165 * t212 + (-t184 * t206 + ((-qJD(5) * t146 - 0.2e1 * t199) * t173 + (-t133 * t173 + (-qJD(5) * t152 + t134) * t171) * t147) * t174) * t164) * t140;
	t120 = (t125 * t227 - t135 * t165) * t174 * t204 + ((-t135 * t206 + (-t125 * t168 - t123) * t226) * t165 + (-t135 * t212 - (-t124 * t144 * t172 - t233 * t143 + (t128 * t225 - t143 * t168 - t144 * t205) * t139) * t201 + (t136 * t206 + t174 * t203) * t125 - ((t124 - t205) * t143 + ((-t139 * t172 + 0.1e1) * t168 + (t139 - t172) * t128) * t144) * t165 * t226) * t164) * t129;
	t1 = [t161 * t174 * t192 + (t168 * t186 - t206 * t216) * t154, 0, t124, t124, 0, 0; (t135 * t193 + (-t135 * t214 + (qJD(1) * t127 + t123) * t227) * t129) * t172 + (t136 * t193 * t127 + (-((t192 - t214 + (t128 * t161 * t217 + t214) * t154) * t143 + (t191 * t220 - t128 * t164 + (-t159 * t198 + (t128 - 0.2e1 * t213) * t164) * t154) * t144) * t201 + (-t136 * t214 + t164 * t203) * t127 + (-t135 + ((-t169 + t170) * t188 + t232 * t200) * t136) * t164 * qJD(1)) * t129) * t174, 0, t120, t120, 0, 0; 0.2e1 * (t146 * t183 + t151 * t222) * t230 + (0.2e1 * t151 * t199 - t190 * t146 * t208 + (t168 * t210 + t185) * t223 + (t151 * t133 + t183 * t134 - t185 * t221 - (t164 * t168 * t173 + t190 * t171) * t152 * t172) * t147) * t140, 0, t121, t121, t202 + 0.2e1 * (-t133 * t140 * t147 + (-t140 * t228 - t147 * t230) * t152) * t152, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end