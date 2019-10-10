% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRPRR1
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
%   Wie in S6RPRPRR1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:46
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPRR1_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR1_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR1_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR1_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR1_jacobiaD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:46:14
	% EndTime: 2019-10-10 00:46:14
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:46:14
	% EndTime: 2019-10-10 00:46:14
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:46:14
	% EndTime: 2019-10-10 00:46:14
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (31->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:46:14
	% EndTime: 2019-10-10 00:46:14
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:46:14
	% EndTime: 2019-10-10 00:46:14
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:46:14
	% EndTime: 2019-10-10 00:46:14
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:46:14
	% EndTime: 2019-10-10 00:46:15
	% DurationCPUTime: 1.15s
	% Computational Cost: add. (7044->97), mult. (3810->206), div. (753->12), fcn. (4455->9), ass. (0->96)
	t182 = qJ(1) + pkin(10);
	t178 = sin(t182);
	t245 = 0.2e1 * t178;
	t176 = t178 ^ 2;
	t180 = qJ(3) + pkin(11) + qJ(5);
	t174 = sin(t180);
	t170 = t174 ^ 2;
	t175 = cos(t180);
	t172 = 0.1e1 / t175 ^ 2;
	t231 = t170 * t172;
	t165 = t176 * t231 + 0.1e1;
	t163 = 0.1e1 / t165;
	t171 = 0.1e1 / t175;
	t179 = cos(t182);
	t216 = qJD(1) * t179;
	t205 = t174 * t216;
	t181 = qJD(3) + qJD(5);
	t223 = t178 * t181;
	t209 = t172 * t223;
	t137 = (-(-t175 * t223 - t205) * t171 + t170 * t209) * t163;
	t244 = t137 - t223;
	t184 = cos(qJ(6));
	t219 = t179 * t184;
	t183 = sin(qJ(6));
	t222 = t178 * t183;
	t159 = t175 * t219 + t222;
	t200 = qJD(6) * t175 - qJD(1);
	t226 = t174 * t181;
	t243 = t200 * t183 + t184 * t226;
	t224 = t178 * t174;
	t162 = atan2(-t224, -t175);
	t161 = cos(t162);
	t160 = sin(t162);
	t210 = t160 * t224;
	t147 = -t161 * t175 - t210;
	t144 = 0.1e1 / t147;
	t153 = 0.1e1 / t159;
	t145 = 0.1e1 / t147 ^ 2;
	t154 = 0.1e1 / t159 ^ 2;
	t242 = t163 - 0.1e1;
	t233 = t161 * t174;
	t132 = (-t137 * t178 + t181) * t233 + (t244 * t175 - t205) * t160;
	t241 = t132 * t144 * t145;
	t193 = t175 * t222 + t219;
	t208 = t183 * t226;
	t141 = t193 * qJD(1) - t159 * qJD(6) + t179 * t208;
	t220 = t179 * t183;
	t221 = t178 * t184;
	t158 = t175 * t220 - t221;
	t152 = t158 ^ 2;
	t151 = t152 * t154 + 0.1e1;
	t235 = t154 * t158;
	t199 = -qJD(1) * t175 + qJD(6);
	t195 = t199 * t184;
	t142 = t178 * t195 - t243 * t179;
	t239 = t142 * t153 * t154;
	t240 = (-t141 * t235 - t152 * t239) / t151 ^ 2;
	t169 = t174 * t170;
	t228 = t171 * t174;
	t192 = (t169 * t171 * t172 + t228) * t181;
	t229 = t170 * t178;
	t197 = t216 * t229;
	t238 = (t172 * t197 + t176 * t192) / t165 ^ 2;
	t237 = t145 * t174;
	t236 = t145 * t179;
	t234 = t160 * t178;
	t232 = t170 * t171;
	t177 = t179 ^ 2;
	t230 = t170 * t177;
	t227 = t174 * t179;
	t225 = t175 * t181;
	t218 = t181 * t144;
	t217 = qJD(1) * t178;
	t140 = t145 * t230 + 0.1e1;
	t215 = 0.2e1 * (-t230 * t241 + (t174 * t177 * t225 - t197) * t145) / t140 ^ 2;
	t214 = 0.2e1 * t241;
	t213 = -0.2e1 * t240;
	t212 = t145 * t227;
	t211 = t158 * t239;
	t204 = 0.1e1 + t231;
	t203 = t174 * t215;
	t202 = -0.2e1 * t174 * t238;
	t201 = t238 * t245;
	t198 = t161 * t163 * t232;
	t196 = t204 * t179;
	t194 = -t153 * t183 + t184 * t235;
	t157 = -t175 * t221 + t220;
	t149 = 0.1e1 / t151;
	t148 = t204 * t178 * t163;
	t138 = 0.1e1 / t140;
	t136 = (t242 * t174 * t160 - t178 * t198) * t179;
	t134 = -t175 * t234 + t233 + (t160 * t175 - t161 * t224) * t148;
	t133 = -t204 * t201 + (qJD(1) * t196 + t192 * t245) * t163;
	t130 = t194 * t213 * t227 + (t194 * t179 * t225 + (-t194 * t217 + ((-qJD(6) * t153 - 0.2e1 * t211) * t184 + (-t141 * t184 + (-qJD(6) * t158 + t142) * t183) * t154) * t179) * t174) * t149;
	t129 = (t134 * t237 - t144 * t175) * t179 * t215 + ((-t144 * t217 + (-t134 * t181 - t132) * t236) * t175 + (-t179 * t218 - (-t133 * t161 * t178 - t244 * t160 + (t137 * t234 - t160 * t181 - t161 * t216) * t148) * t212 + (t145 * t217 + t179 * t214) * t134 - ((t133 - t216) * t160 + ((-t148 * t178 + 0.1e1) * t181 + (t148 - t178) * t137) * t161) * t175 * t236) * t174) * t138;
	t1 = [t179 * t171 * t202 + (t181 * t196 - t217 * t228) * t163, 0, t133, 0, t133, 0; (t144 * t203 + (-t175 * t218 + (qJD(1) * t136 + t132) * t237) * t138) * t178 + (t145 * t203 * t136 + (-((t202 - t225 + (t137 * t171 * t229 + t225) * t163) * t160 + (t201 * t232 - t137 * t174 + (-t169 * t209 + (t137 - 0.2e1 * t223) * t174) * t163) * t161) * t212 + (-t145 * t225 + t174 * t214) * t136 + (-t144 + ((-t176 + t177) * t198 + t242 * t210) * t145) * t174 * qJD(1)) * t138) * t179, 0, t129, 0, t129, 0; 0.2e1 * (t153 * t193 + t157 * t235) * t240 + (0.2e1 * t157 * t211 + (t157 * t141 + t193 * t142 + (-t243 * t178 - t179 * t195) * t158) * t154 + (t199 * t220 + (-t200 * t184 + t208) * t178) * t153) * t149, 0, t130, 0, t130, t213 + 0.2e1 * (-t141 * t154 * t149 + (-t149 * t239 - t154 * t240) * t158) * t158;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end