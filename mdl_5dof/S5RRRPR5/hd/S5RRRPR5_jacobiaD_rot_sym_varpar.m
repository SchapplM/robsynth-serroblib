% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRRPR5
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
%   Wie in S5RRRPR5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 19:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RRRPR5_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR5_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR5_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPR5_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR5_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:59:38
	% EndTime: 2019-12-29 19:59:38
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:59:38
	% EndTime: 2019-12-29 19:59:38
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:59:38
	% EndTime: 2019-12-29 19:59:38
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:59:38
	% EndTime: 2019-12-29 19:59:38
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:59:38
	% EndTime: 2019-12-29 19:59:38
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:59:39
	% EndTime: 2019-12-29 19:59:40
	% DurationCPUTime: 1.63s
	% Computational Cost: add. (5587->97), mult. (3810->208), div. (753->12), fcn. (4455->9), ass. (0->95)
	t179 = sin(qJ(1));
	t241 = 0.2e1 * t179;
	t176 = t179 ^ 2;
	t174 = qJ(2) + qJ(3) + pkin(9);
	t171 = sin(t174);
	t167 = t171 ^ 2;
	t172 = cos(t174);
	t169 = 0.1e1 / t172 ^ 2;
	t226 = t167 * t169;
	t163 = t176 * t226 + 0.1e1;
	t161 = 0.1e1 / t163;
	t168 = 0.1e1 / t172;
	t181 = cos(qJ(1));
	t212 = qJD(1) * t181;
	t202 = t171 * t212;
	t175 = qJD(2) + qJD(3);
	t220 = t175 * t179;
	t205 = t169 * t220;
	t135 = (-(-t172 * t220 - t202) * t168 + t167 * t205) * t161;
	t240 = t135 - t220;
	t180 = cos(qJ(5));
	t214 = t181 * t180;
	t178 = sin(qJ(5));
	t217 = t179 * t178;
	t160 = t172 * t214 + t217;
	t218 = t179 * t171;
	t156 = atan2(-t218, -t172);
	t151 = cos(t156);
	t150 = sin(t156);
	t207 = t150 * t218;
	t145 = -t151 * t172 - t207;
	t142 = 0.1e1 / t145;
	t153 = 0.1e1 / t160;
	t143 = 0.1e1 / t145 ^ 2;
	t154 = 0.1e1 / t160 ^ 2;
	t239 = t161 - 0.1e1;
	t231 = t151 * t171;
	t130 = (-t135 * t179 + t175) * t231 + (t240 * t172 - t202) * t150;
	t238 = t130 * t142 * t143;
	t190 = t172 * t217 + t214;
	t219 = t175 * t181;
	t203 = t171 * t219;
	t140 = qJD(1) * t190 - t160 * qJD(5) + t178 * t203;
	t215 = t181 * t178;
	t216 = t179 * t180;
	t159 = t172 * t215 - t216;
	t152 = t159 ^ 2;
	t149 = t152 * t154 + 0.1e1;
	t229 = t154 * t159;
	t196 = -qJD(1) * t172 + qJD(5);
	t197 = qJD(5) * t172 - qJD(1);
	t141 = -t197 * t215 + (t179 * t196 - t203) * t180;
	t235 = t141 * t153 * t154;
	t237 = (-t140 * t229 - t152 * t235) / t149 ^ 2;
	t166 = t171 * t167;
	t223 = t168 * t171;
	t189 = t175 * (t166 * t168 * t169 + t223);
	t224 = t167 * t179;
	t194 = t212 * t224;
	t236 = (t169 * t194 + t176 * t189) / t163 ^ 2;
	t234 = t143 * t171;
	t233 = t143 * t181;
	t232 = t150 * t179;
	t230 = t153 * t178;
	t228 = t159 * t180;
	t227 = t167 * t168;
	t177 = t181 ^ 2;
	t225 = t167 * t177;
	t222 = t171 * t181;
	t221 = t172 * t175;
	t213 = qJD(1) * t179;
	t138 = t143 * t225 + 0.1e1;
	t211 = 0.2e1 * (-t225 * t238 + (t171 * t177 * t221 - t194) * t143) / t138 ^ 2;
	t210 = 0.2e1 * t238;
	t209 = -0.2e1 * t237;
	t208 = t143 * t222;
	t206 = t159 * t235;
	t201 = 0.1e1 + t226;
	t200 = t171 * t211;
	t199 = -0.2e1 * t171 * t236;
	t198 = t236 * t241;
	t195 = t151 * t161 * t227;
	t193 = t201 * t181;
	t192 = t196 * t181;
	t191 = t154 * t228 - t230;
	t158 = -t172 * t216 + t215;
	t147 = 0.1e1 / t149;
	t146 = t201 * t179 * t161;
	t136 = 0.1e1 / t138;
	t134 = (t239 * t171 * t150 - t179 * t195) * t181;
	t132 = -t172 * t232 + t231 + (t150 * t172 - t151 * t218) * t146;
	t131 = -t201 * t198 + (qJD(1) * t193 + t189 * t241) * t161;
	t128 = t191 * t209 * t222 + (t191 * t172 * t219 + (-t191 * t213 + ((-qJD(5) * t153 - 0.2e1 * t206) * t180 + (-t140 * t180 + (-qJD(5) * t159 + t141) * t178) * t154) * t181) * t171) * t147;
	t127 = (t132 * t234 - t142 * t172) * t181 * t211 + ((-t142 * t213 + (-t132 * t175 - t130) * t233) * t172 + (-t142 * t219 - (-t131 * t151 * t179 - t240 * t150 + (t135 * t232 - t150 * t175 - t151 * t212) * t146) * t208 + (t143 * t213 + t181 * t210) * t132 - ((t131 - t212) * t150 + ((-t146 * t179 + 0.1e1) * t175 + (t146 - t179) * t135) * t151) * t172 * t233) * t171) * t136;
	t1 = [t181 * t168 * t199 + (t175 * t193 - t213 * t223) * t161, t131, t131, 0, 0; (t142 * t200 + (-t142 * t221 + (qJD(1) * t134 + t130) * t234) * t136) * t179 + (t143 * t200 * t134 + (-((t199 - t221 + (t135 * t168 * t224 + t221) * t161) * t150 + (t198 * t227 - t135 * t171 + (-t166 * t205 + (t135 - 0.2e1 * t220) * t171) * t161) * t151) * t208 + (-t143 * t221 + t171 * t210) * t134 + (-t142 + ((-t176 + t177) * t195 + t239 * t207) * t143) * t171 * qJD(1)) * t136) * t181, t127, t127, 0, 0; 0.2e1 * (t153 * t190 + t158 * t229) * t237 + (0.2e1 * t158 * t206 - t197 * t153 * t216 + (t175 * t218 + t192) * t230 + (t158 * t140 + t190 * t141 - t192 * t228 - (t171 * t175 * t180 + t178 * t197) * t159 * t179) * t154) * t147, t128, t128, 0, t209 + 0.2e1 * (-t140 * t154 * t147 + (-t147 * t235 - t154 * t237) * t159) * t159;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end