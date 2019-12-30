% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRRRP6
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
%   Wie in S5RRRRP6_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 20:33
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RRRRP6_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP6_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP6_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRP6_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP6_jacobiaD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:33:16
	% EndTime: 2019-12-29 20:33:16
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:33:16
	% EndTime: 2019-12-29 20:33:16
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:33:09
	% EndTime: 2019-12-29 20:33:09
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:33:14
	% EndTime: 2019-12-29 20:33:14
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:33:10
	% EndTime: 2019-12-29 20:33:12
	% DurationCPUTime: 1.69s
	% Computational Cost: add. (3645->97), mult. (3810->208), div. (753->12), fcn. (4455->9), ass. (0->95)
	t173 = sin(qJ(1));
	t235 = 0.2e1 * t173;
	t169 = t173 ^ 2;
	t171 = qJ(2) + qJ(3);
	t166 = sin(t171);
	t162 = t166 ^ 2;
	t167 = cos(t171);
	t164 = 0.1e1 / t167 ^ 2;
	t220 = t162 * t164;
	t157 = t169 * t220 + 0.1e1;
	t155 = 0.1e1 / t157;
	t163 = 0.1e1 / t167;
	t175 = cos(qJ(1));
	t206 = qJD(1) * t175;
	t196 = t166 * t206;
	t168 = qJD(2) + qJD(3);
	t214 = t168 * t173;
	t199 = t164 * t214;
	t129 = (-(-t167 * t214 - t196) * t163 + t162 * t199) * t155;
	t234 = t129 - t214;
	t174 = cos(qJ(4));
	t208 = t174 * t175;
	t172 = sin(qJ(4));
	t210 = t173 * t172;
	t151 = t167 * t208 + t210;
	t211 = t173 * t166;
	t154 = atan2(-t211, -t167);
	t153 = cos(t154);
	t152 = sin(t154);
	t200 = t152 * t211;
	t139 = -t153 * t167 - t200;
	t136 = 0.1e1 / t139;
	t145 = 0.1e1 / t151;
	t137 = 0.1e1 / t139 ^ 2;
	t146 = 0.1e1 / t151 ^ 2;
	t233 = t155 - 0.1e1;
	t222 = t153 * t166;
	t124 = (-t129 * t173 + t168) * t222 + (t167 * t234 - t196) * t152;
	t232 = t124 * t136 * t137;
	t184 = t167 * t210 + t208;
	t213 = t168 * t175;
	t197 = t166 * t213;
	t133 = t184 * qJD(1) - qJD(4) * t151 + t172 * t197;
	t209 = t173 * t174;
	t212 = t172 * t175;
	t150 = t167 * t212 - t209;
	t144 = t150 ^ 2;
	t143 = t144 * t146 + 0.1e1;
	t225 = t146 * t150;
	t190 = -qJD(1) * t167 + qJD(4);
	t191 = qJD(4) * t167 - qJD(1);
	t134 = -t191 * t212 + (t190 * t173 - t197) * t174;
	t230 = t134 * t145 * t146;
	t231 = (-t133 * t225 - t144 * t230) / t143 ^ 2;
	t161 = t166 * t162;
	t217 = t163 * t166;
	t183 = t168 * (t161 * t163 * t164 + t217);
	t218 = t162 * t173;
	t188 = t206 * t218;
	t229 = (t164 * t188 + t169 * t183) / t157 ^ 2;
	t228 = t137 * t166;
	t227 = t137 * t175;
	t226 = t145 * t172;
	t224 = t150 * t174;
	t223 = t152 * t173;
	t221 = t162 * t163;
	t170 = t175 ^ 2;
	t219 = t162 * t170;
	t216 = t166 * t175;
	t215 = t167 * t168;
	t207 = qJD(1) * t173;
	t132 = t137 * t219 + 0.1e1;
	t205 = 0.2e1 * (-t219 * t232 + (t166 * t170 * t215 - t188) * t137) / t132 ^ 2;
	t204 = 0.2e1 * t232;
	t203 = -0.2e1 * t231;
	t202 = t150 * t230;
	t201 = t137 * t216;
	t195 = 0.1e1 + t220;
	t194 = t166 * t205;
	t193 = -0.2e1 * t166 * t229;
	t192 = t229 * t235;
	t189 = t153 * t155 * t221;
	t187 = t195 * t175;
	t186 = t190 * t175;
	t185 = t146 * t224 - t226;
	t149 = -t167 * t209 + t212;
	t141 = 0.1e1 / t143;
	t140 = t195 * t173 * t155;
	t130 = 0.1e1 / t132;
	t128 = (t233 * t166 * t152 - t173 * t189) * t175;
	t126 = -t167 * t223 + t222 + (t152 * t167 - t153 * t211) * t140;
	t125 = -t195 * t192 + (qJD(1) * t187 + t183 * t235) * t155;
	t122 = t185 * t203 * t216 + (t185 * t167 * t213 + (-t185 * t207 + ((-qJD(4) * t145 - 0.2e1 * t202) * t174 + (-t133 * t174 + (-qJD(4) * t150 + t134) * t172) * t146) * t175) * t166) * t141;
	t121 = (t126 * t228 - t136 * t167) * t175 * t205 + ((-t136 * t207 + (-t126 * t168 - t124) * t227) * t167 + (-t136 * t213 - (-t125 * t153 * t173 - t234 * t152 + (t129 * t223 - t152 * t168 - t153 * t206) * t140) * t201 + (t137 * t207 + t175 * t204) * t126 - ((t125 - t206) * t152 + ((-t140 * t173 + 0.1e1) * t168 + (t140 - t173) * t129) * t153) * t167 * t227) * t166) * t130;
	t1 = [t163 * t175 * t193 + (t168 * t187 - t207 * t217) * t155, t125, t125, 0, 0; (t136 * t194 + (-t136 * t215 + (qJD(1) * t128 + t124) * t228) * t130) * t173 + (t137 * t194 * t128 + (-((t193 - t215 + (t129 * t163 * t218 + t215) * t155) * t152 + (t192 * t221 - t129 * t166 + (-t161 * t199 + (t129 - 0.2e1 * t214) * t166) * t155) * t153) * t201 + (-t137 * t215 + t166 * t204) * t128 + (-t136 + ((-t169 + t170) * t189 + t233 * t200) * t137) * t166 * qJD(1)) * t130) * t175, t121, t121, 0, 0; 0.2e1 * (t145 * t184 + t149 * t225) * t231 + (0.2e1 * t149 * t202 - t191 * t145 * t209 + (t168 * t211 + t186) * t226 + (t149 * t133 + t184 * t134 - t186 * t224 - (t166 * t168 * t174 + t191 * t172) * t150 * t173) * t146) * t141, t122, t122, t203 + 0.2e1 * (-t133 * t141 * t146 + (-t141 * t230 - t146 * t231) * t150) * t150, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:33:14
	% EndTime: 2019-12-29 20:33:16
	% DurationCPUTime: 1.68s
	% Computational Cost: add. (3645->97), mult. (3810->208), div. (753->12), fcn. (4455->9), ass. (0->95)
	t171 = sin(qJ(1));
	t233 = 0.2e1 * t171;
	t167 = t171 ^ 2;
	t169 = qJ(2) + qJ(3);
	t164 = sin(t169);
	t160 = t164 ^ 2;
	t165 = cos(t169);
	t162 = 0.1e1 / t165 ^ 2;
	t218 = t160 * t162;
	t155 = t167 * t218 + 0.1e1;
	t153 = 0.1e1 / t155;
	t161 = 0.1e1 / t165;
	t173 = cos(qJ(1));
	t204 = qJD(1) * t173;
	t194 = t164 * t204;
	t166 = qJD(2) + qJD(3);
	t212 = t166 * t171;
	t197 = t162 * t212;
	t127 = (-(-t165 * t212 - t194) * t161 + t160 * t197) * t153;
	t232 = t127 - t212;
	t172 = cos(qJ(4));
	t206 = t172 * t173;
	t170 = sin(qJ(4));
	t208 = t171 * t170;
	t149 = t165 * t206 + t208;
	t209 = t171 * t164;
	t152 = atan2(-t209, -t165);
	t151 = cos(t152);
	t150 = sin(t152);
	t198 = t150 * t209;
	t137 = -t151 * t165 - t198;
	t134 = 0.1e1 / t137;
	t143 = 0.1e1 / t149;
	t135 = 0.1e1 / t137 ^ 2;
	t144 = 0.1e1 / t149 ^ 2;
	t231 = t153 - 0.1e1;
	t220 = t151 * t164;
	t122 = (-t127 * t171 + t166) * t220 + (t165 * t232 - t194) * t150;
	t230 = t122 * t134 * t135;
	t182 = t165 * t208 + t206;
	t211 = t166 * t173;
	t195 = t164 * t211;
	t131 = t182 * qJD(1) - qJD(4) * t149 + t170 * t195;
	t207 = t171 * t172;
	t210 = t170 * t173;
	t148 = t165 * t210 - t207;
	t142 = t148 ^ 2;
	t141 = t142 * t144 + 0.1e1;
	t223 = t144 * t148;
	t188 = -qJD(1) * t165 + qJD(4);
	t189 = qJD(4) * t165 - qJD(1);
	t132 = -t189 * t210 + (t188 * t171 - t195) * t172;
	t228 = t132 * t143 * t144;
	t229 = (-t131 * t223 - t142 * t228) / t141 ^ 2;
	t159 = t164 * t160;
	t215 = t161 * t164;
	t181 = t166 * (t159 * t161 * t162 + t215);
	t216 = t160 * t171;
	t186 = t204 * t216;
	t227 = (t162 * t186 + t167 * t181) / t155 ^ 2;
	t226 = t135 * t164;
	t225 = t135 * t173;
	t224 = t143 * t170;
	t222 = t148 * t172;
	t221 = t150 * t171;
	t219 = t160 * t161;
	t168 = t173 ^ 2;
	t217 = t160 * t168;
	t214 = t164 * t173;
	t213 = t165 * t166;
	t205 = qJD(1) * t171;
	t130 = t135 * t217 + 0.1e1;
	t203 = 0.2e1 * (-t217 * t230 + (t164 * t168 * t213 - t186) * t135) / t130 ^ 2;
	t202 = 0.2e1 * t230;
	t201 = -0.2e1 * t229;
	t200 = t148 * t228;
	t199 = t135 * t214;
	t193 = 0.1e1 + t218;
	t192 = t164 * t203;
	t191 = -0.2e1 * t164 * t227;
	t190 = t227 * t233;
	t187 = t151 * t153 * t219;
	t185 = t193 * t173;
	t184 = t188 * t173;
	t183 = t144 * t222 - t224;
	t147 = -t165 * t207 + t210;
	t139 = 0.1e1 / t141;
	t138 = t193 * t171 * t153;
	t128 = 0.1e1 / t130;
	t126 = (t231 * t164 * t150 - t171 * t187) * t173;
	t124 = -t165 * t221 + t220 + (t150 * t165 - t151 * t209) * t138;
	t123 = -t193 * t190 + (qJD(1) * t185 + t181 * t233) * t153;
	t120 = t183 * t201 * t214 + (t183 * t165 * t211 + (-t183 * t205 + ((-qJD(4) * t143 - 0.2e1 * t200) * t172 + (-t131 * t172 + (-qJD(4) * t148 + t132) * t170) * t144) * t173) * t164) * t139;
	t119 = (t124 * t226 - t134 * t165) * t173 * t203 + ((-t134 * t205 + (-t124 * t166 - t122) * t225) * t165 + (-t134 * t211 - (-t123 * t151 * t171 - t232 * t150 + (t127 * t221 - t150 * t166 - t151 * t204) * t138) * t199 + (t135 * t205 + t173 * t202) * t124 - ((t123 - t204) * t150 + ((-t138 * t171 + 0.1e1) * t166 + (t138 - t171) * t127) * t151) * t165 * t225) * t164) * t128;
	t1 = [t161 * t173 * t191 + (t166 * t185 - t205 * t215) * t153, t123, t123, 0, 0; (t134 * t192 + (-t134 * t213 + (qJD(1) * t126 + t122) * t226) * t128) * t171 + (t135 * t192 * t126 + (-((t191 - t213 + (t127 * t161 * t216 + t213) * t153) * t150 + (t190 * t219 - t127 * t164 + (-t159 * t197 + (t127 - 0.2e1 * t212) * t164) * t153) * t151) * t199 + (-t135 * t213 + t164 * t202) * t126 + (-t134 + ((-t167 + t168) * t187 + t231 * t198) * t135) * t164 * qJD(1)) * t128) * t173, t119, t119, 0, 0; 0.2e1 * (t143 * t182 + t147 * t223) * t229 + (0.2e1 * t147 * t200 - t189 * t143 * t207 + (t166 * t209 + t184) * t224 + (t147 * t131 + t182 * t132 - t184 * t222 - (t164 * t166 * t172 + t189 * t170) * t148 * t171) * t144) * t139, t120, t120, t201 + 0.2e1 * (-t131 * t139 * t144 + (-t139 * t228 - t144 * t229) * t148) * t148, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end