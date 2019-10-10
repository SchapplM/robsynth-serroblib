% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRRRR2
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
%   Wie in S5RRRRR2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:02
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RRRRR2_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR2_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR2_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR2_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5RRRRR2_jacobiaD_rot_sym_varpar: pkin has to be [2x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:02:47
	% EndTime: 2019-10-09 21:02:47
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:02:47
	% EndTime: 2019-10-09 21:02:47
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:02:47
	% EndTime: 2019-10-09 21:02:47
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (74->0), mult. (74->0), div. (30->0), fcn. (44->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:02:47
	% EndTime: 2019-10-09 21:02:47
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:02:47
	% EndTime: 2019-10-09 21:02:47
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:02:47
	% EndTime: 2019-10-09 21:02:48
	% DurationCPUTime: 1.26s
	% Computational Cost: add. (6736->101), mult. (4823->213), div. (942->12), fcn. (5705->9), ass. (0->100)
	t184 = qJ(1) + qJ(2);
	t178 = sin(t184);
	t183 = qJ(3) + qJ(4);
	t177 = sin(t183);
	t171 = t177 ^ 2;
	t179 = cos(t183);
	t174 = 0.1e1 / t179 ^ 2;
	t232 = t171 * t174;
	t203 = 0.1e1 + t232;
	t246 = t178 * t203;
	t180 = cos(t184);
	t186 = cos(qJ(5));
	t220 = t180 * t186;
	t185 = sin(qJ(5));
	t225 = t178 * t185;
	t160 = t179 * t220 + t225;
	t182 = qJD(1) + qJD(2);
	t200 = qJD(5) * t179 - t182;
	t181 = qJD(3) + qJD(4);
	t219 = t181 * t177;
	t245 = t200 * t185 + t186 * t219;
	t228 = t178 * t177;
	t163 = atan2(-t228, -t179);
	t161 = sin(t163);
	t162 = cos(t163);
	t148 = -t161 * t228 - t162 * t179;
	t145 = 0.1e1 / t148;
	t154 = 0.1e1 / t160;
	t173 = 0.1e1 / t179;
	t146 = 0.1e1 / t148 ^ 2;
	t155 = 0.1e1 / t160 ^ 2;
	t172 = t178 ^ 2;
	t166 = t172 * t232 + 0.1e1;
	t164 = 0.1e1 / t166;
	t222 = t180 * t182;
	t204 = t177 * t222;
	t227 = t178 * t181;
	t210 = t174 * t227;
	t223 = t179 * t181;
	t138 = (-(-t178 * t223 - t204) * t173 + t171 * t210) * t164;
	t233 = t162 * t177;
	t132 = (-t138 * t178 + t181) * t233 + (-t204 + (t138 - t227) * t179) * t161;
	t244 = t132 * t145 * t146;
	t206 = t179 * t225;
	t208 = t185 * t219;
	t139 = t182 * t206 + (t182 * t186 + t208) * t180 - t160 * qJD(5);
	t221 = t180 * t185;
	t224 = t178 * t186;
	t159 = t179 * t221 - t224;
	t153 = t159 ^ 2;
	t152 = t153 * t155 + 0.1e1;
	t235 = t155 * t159;
	t201 = -t179 * t182 + qJD(5);
	t196 = t201 * t186;
	t140 = t178 * t196 - t245 * t180;
	t240 = t140 * t154 * t155;
	t243 = (-t139 * t235 - t153 * t240) / t152 ^ 2;
	t230 = t171 * t178;
	t211 = t173 * t230;
	t199 = t162 * t211;
	t137 = (-t164 * t199 + (t164 - 0.1e1) * t177 * t161) * t180;
	t242 = t137 * t146;
	t241 = t138 * t161;
	t170 = t177 * t171;
	t194 = (t170 * t174 + t177) * t173 * t181;
	t198 = t222 * t230;
	t239 = (t172 * t194 + t174 * t198) / t166 ^ 2;
	t238 = t146 * t177;
	t237 = t146 * t180;
	t149 = t164 * t246;
	t236 = t149 * t178;
	t234 = t161 * t179;
	t176 = t180 ^ 2;
	t231 = t171 * t176;
	t229 = t177 * t180;
	t226 = t178 * t182;
	t218 = t181 * t180;
	t217 = t182 * t145;
	t143 = t146 * t231 + 0.1e1;
	t216 = 0.2e1 * (-t231 * t244 + (t176 * t179 * t219 - t198) * t146) / t143 ^ 2;
	t215 = 0.2e1 * t244;
	t214 = -0.2e1 * t243;
	t213 = -0.2e1 * t239;
	t212 = t159 * t240;
	t209 = t177 * t226;
	t202 = t177 * t216;
	t197 = t180 * t203;
	t195 = -t154 * t185 + t186 * t235;
	t158 = -t179 * t224 + t221;
	t157 = -t206 - t220;
	t150 = 0.1e1 / t152;
	t141 = 0.1e1 / t143;
	t136 = t173 * t213 * t229 + (-t173 * t209 + t181 * t197) * t164;
	t135 = -t178 * t234 + t233 + (-t162 * t228 + t234) * t149;
	t133 = t213 * t246 + (0.2e1 * t178 * t194 + t182 * t197) * t164;
	t130 = 0.2e1 * (-t154 * t157 + t158 * t235) * t243 + (0.2e1 * t158 * t212 + (t158 * t139 - t157 * t140 + (-t245 * t178 - t180 * t196) * t159) * t155 + (t201 * t221 + (-t200 * t186 + t208) * t178) * t154) * t150;
	t129 = t195 * t214 * t229 + (t195 * t179 * t218 + (-t195 * t226 + ((-qJD(5) * t154 - 0.2e1 * t212) * t186 + (-t139 * t186 + (-qJD(5) * t159 + t140) * t185) * t155) * t180) * t177) * t150;
	t128 = (t145 * t202 + (-t145 * t223 + (t137 * t182 + t132) * t238) * t141) * t178 + (t202 * t242 + (-t223 * t242 + (t137 * t215 - t217 + (-(-t161 * t223 + 0.2e1 * t199 * t239) * t180 + (-t161 * t226 - (-t138 * t162 + t161 * t213) * t180) * t177 + ((-(t138 * t211 + t223) * t180 + t209) * t161 + ((t170 * t210 - (t138 - 0.2e1 * t227) * t177) * t180 + (-t172 + t176) * t182 * t173 * t171) * t162) * t164) * t146) * t177) * t141) * t180;
	t127 = (t135 * t238 - t145 * t179) * t180 * t216 + ((-t178 * t217 + (-t135 * t181 - t132) * t237) * t179 + (-t145 * t218 - (-t133 * t162 * t178 + t161 * t227 + t236 * t241 - t241 + (-t161 * t181 - t162 * t222) * t149) * t146 * t229 + (t146 * t226 + t180 * t215) * t135 - ((t133 - t222) * t161 + ((0.1e1 - t236) * t181 + (t149 - t178) * t138) * t162) * t179 * t237) * t177) * t141;
	t1 = [t136, t136, t133, t133, 0; t128, t128, t127, t127, 0; t130, t130, t129, t129, t214 + 0.2e1 * (-t139 * t155 * t150 + (-t150 * t240 - t155 * t243) * t159) * t159;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end