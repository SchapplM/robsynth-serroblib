% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPPRRR7
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
%   Wie in S6RPPRRR7_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:11
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPPRRR7_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR7_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR7_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRR7_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR7_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:11:46
	% EndTime: 2019-10-10 00:11:46
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:11:46
	% EndTime: 2019-10-10 00:11:46
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:11:46
	% EndTime: 2019-10-10 00:11:46
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:11:46
	% EndTime: 2019-10-10 00:11:46
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:11:46
	% EndTime: 2019-10-10 00:11:46
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:11:46
	% EndTime: 2019-10-10 00:11:46
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:11:47
	% EndTime: 2019-10-10 00:11:48
	% DurationCPUTime: 1.05s
	% Computational Cost: add. (5302->95), mult. (3810->207), div. (753->12), fcn. (4455->9), ass. (0->97)
	t170 = cos(qJ(1));
	t230 = 0.2e1 * t170;
	t166 = t170 ^ 2;
	t163 = pkin(10) + qJ(4) + qJ(5);
	t160 = sin(t163);
	t156 = 0.1e1 / t160 ^ 2;
	t161 = cos(t163);
	t159 = t161 ^ 2;
	t213 = t156 * t159;
	t154 = t166 * t213 + 0.1e1;
	t152 = 0.1e1 / t154;
	t155 = 0.1e1 / t160;
	t168 = sin(qJ(1));
	t202 = qJD(1) * t168;
	t194 = t161 * t202;
	t164 = qJD(4) + qJD(5);
	t208 = t164 * t170;
	t195 = t156 * t208;
	t126 = ((t160 * t208 + t194) * t155 + t159 * t195) * t152;
	t229 = -t126 + t208;
	t204 = t170 * t161;
	t151 = atan2(-t204, t160);
	t141 = sin(t151);
	t142 = cos(t151);
	t136 = -t141 * t204 + t142 * t160;
	t133 = 0.1e1 / t136;
	t169 = cos(qJ(6));
	t205 = t168 * t169;
	t167 = sin(qJ(6));
	t207 = t167 * t170;
	t148 = t160 * t205 + t207;
	t144 = 0.1e1 / t148;
	t134 = 0.1e1 / t136 ^ 2;
	t145 = 0.1e1 / t148 ^ 2;
	t165 = t168 ^ 2;
	t212 = t159 * t165;
	t129 = t134 * t212 + 0.1e1;
	t201 = qJD(1) * t170;
	t186 = t159 * t168 * t201;
	t210 = t161 * t164;
	t219 = t142 * t161;
	t225 = t126 * t170;
	t121 = (t164 - t225) * t219 + (t160 * t229 + t194) * t141;
	t227 = t121 * t133 * t134;
	t228 = (-t212 * t227 + (-t160 * t165 * t210 + t186) * t134) / t129 ^ 2;
	t189 = qJD(1) * t160 + qJD(6);
	t182 = t189 * t170;
	t190 = qJD(6) * t160 + qJD(1);
	t183 = t190 * t169;
	t209 = t164 * t168;
	t131 = t168 * t183 + (t161 * t209 + t182) * t167;
	t203 = t170 * t169;
	t206 = t168 * t167;
	t147 = t160 * t206 - t203;
	t143 = t147 ^ 2;
	t140 = t143 * t145 + 0.1e1;
	t217 = t145 * t147;
	t184 = t190 * t167;
	t132 = t169 * t182 + (t169 * t210 - t184) * t168;
	t223 = t132 * t144 * t145;
	t226 = (t131 * t217 - t143 * t223) / t140 ^ 2;
	t158 = t161 * t159;
	t214 = t155 * t161;
	t180 = t164 * (-t155 * t156 * t158 - t214);
	t224 = (-t156 * t186 + t166 * t180) / t154 ^ 2;
	t222 = t134 * t161;
	t221 = t134 * t168;
	t220 = t141 * t170;
	t218 = t144 * t167;
	t216 = t147 * t169;
	t215 = t155 * t159;
	t211 = t160 * t164;
	t200 = -0.2e1 * t227;
	t199 = 0.2e1 * t226;
	t198 = t161 * t228;
	t197 = t161 * t224;
	t196 = t161 * t221;
	t193 = 0.1e1 + t213;
	t192 = t224 * t230;
	t191 = 0.2e1 * t147 * t223;
	t188 = t142 * t152 * t215;
	t187 = (-t152 + 0.1e1) * t161 * t141;
	t185 = t193 * t168;
	t181 = t145 * t216 - t218;
	t179 = t181 * t168;
	t178 = t164 * t204 - t189 * t168;
	t150 = t160 * t203 - t206;
	t149 = t160 * t207 + t205;
	t138 = 0.1e1 / t140;
	t137 = t193 * t170 * t152;
	t127 = 0.1e1 / t129;
	t125 = (-t170 * t188 + t187) * t168;
	t123 = t160 * t220 + t219 + (-t141 * t160 - t142 * t204) * t137;
	t122 = -t193 * t192 + (-qJD(1) * t185 + t180 * t230) * t152;
	t119 = t161 * t179 * t199 + (t179 * t211 + (-t181 * t201 + ((qJD(6) * t144 + t191) * t169 + (-t131 * t169 + (qJD(6) * t147 - t132) * t167) * t145) * t168) * t161) * t138;
	t118 = 0.2e1 * (-t123 * t222 - t133 * t160) * t168 * t228 + ((t133 * t201 + (-t123 * t164 - t121) * t221) * t160 + (t133 * t209 + (-t122 * t142 * t170 + t229 * t141 + (t126 * t220 - t141 * t164 + t142 * t202) * t137) * t196 + (t134 * t201 + t168 * t200) * t123 + ((-t122 - t202) * t141 + ((t137 * t170 - 0.1e1) * t164 + (-t137 + t170) * t126) * t142) * t160 * t221) * t161) * t127;
	t1 = [-0.2e1 * t168 * t155 * t197 + (-t164 * t185 + t201 * t214) * t152, 0, 0, t122, t122, 0; (0.2e1 * t133 * t198 + (t133 * t211 + (qJD(1) * t125 + t121) * t222) * t127) * t170 + (-0.2e1 * t134 * t198 * t125 + (((0.2e1 * t197 - t211 + (t215 * t225 + t211) * t152) * t141 + (t192 * t215 + t126 * t161 + (t158 * t195 + (-t126 + 0.2e1 * t208) * t161) * t152) * t142) * t196 + (-t134 * t211 + t161 * t200) * t125 + (t133 + ((t165 - t166) * t188 + t170 * t187) * t134) * t161 * qJD(1)) * t127) * t168, 0, 0, t118, t118, 0; (-t144 * t149 + t150 * t217) * t199 + (t150 * t191 + t170 * t144 * t183 + t178 * t218 + (t170 * t147 * t184 - t150 * t131 - t149 * t132 - t178 * t216) * t145) * t138, 0, 0, t119, t119, -0.2e1 * t226 + 0.2e1 * (t131 * t145 * t138 + (-t138 * t223 - t145 * t226) * t147) * t147;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end