% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRRPRR1
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
%   Wie in S6RRRPRR1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:55
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPRR1_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR1_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR1_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR1_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR1_jacobiaD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:55:18
	% EndTime: 2019-10-10 11:55:18
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:55:18
	% EndTime: 2019-10-10 11:55:18
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:55:18
	% EndTime: 2019-10-10 11:55:18
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:55:18
	% EndTime: 2019-10-10 11:55:18
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:55:18
	% EndTime: 2019-10-10 11:55:18
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:55:18
	% EndTime: 2019-10-10 11:55:18
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:55:18
	% EndTime: 2019-10-10 11:55:19
	% DurationCPUTime: 1.21s
	% Computational Cost: add. (10386->97), mult. (5101->208), div. (1026->12), fcn. (5942->9), ass. (0->95)
	t183 = sin(qJ(1));
	t245 = 0.2e1 * t183;
	t180 = t183 ^ 2;
	t177 = qJ(2) + qJ(3) + pkin(11) + qJ(5);
	t174 = sin(t177);
	t170 = t174 ^ 2;
	t175 = cos(t177);
	t172 = 0.1e1 / t175 ^ 2;
	t230 = t170 * t172;
	t167 = t180 * t230 + 0.1e1;
	t161 = 0.1e1 / t167;
	t171 = 0.1e1 / t175;
	t185 = cos(qJ(1));
	t216 = qJD(1) * t185;
	t206 = t174 * t216;
	t179 = qJD(2) + qJD(3) + qJD(5);
	t224 = t179 * t183;
	t209 = t172 * t224;
	t139 = (-(-t175 * t224 - t206) * t171 + t170 * t209) * t161;
	t244 = t139 - t224;
	t184 = cos(qJ(6));
	t218 = t184 * t185;
	t182 = sin(qJ(6));
	t220 = t183 * t182;
	t166 = t175 * t218 + t220;
	t221 = t183 * t174;
	t156 = atan2(-t221, -t175);
	t155 = cos(t156);
	t154 = sin(t156);
	t210 = t154 * t221;
	t147 = -t155 * t175 - t210;
	t144 = 0.1e1 / t147;
	t158 = 0.1e1 / t166;
	t145 = 0.1e1 / t147 ^ 2;
	t159 = 0.1e1 / t166 ^ 2;
	t243 = t161 - 0.1e1;
	t235 = t155 * t174;
	t134 = (-t139 * t183 + t179) * t235 + (t244 * t175 - t206) * t154;
	t242 = t134 * t144 * t145;
	t194 = t175 * t220 + t218;
	t223 = t179 * t185;
	t208 = t174 * t223;
	t148 = t194 * qJD(1) - t166 * qJD(6) + t182 * t208;
	t219 = t183 * t184;
	t222 = t182 * t185;
	t165 = t175 * t222 - t219;
	t157 = t165 ^ 2;
	t153 = t157 * t159 + 0.1e1;
	t233 = t159 * t165;
	t200 = -qJD(1) * t175 + qJD(6);
	t201 = qJD(6) * t175 - qJD(1);
	t149 = -t201 * t222 + (t200 * t183 - t208) * t184;
	t237 = t149 * t158 * t159;
	t241 = (-t148 * t233 - t157 * t237) / t153 ^ 2;
	t169 = t174 * t170;
	t227 = t171 * t174;
	t193 = t179 * (t169 * t171 * t172 + t227);
	t228 = t170 * t183;
	t198 = t216 * t228;
	t240 = (t172 * t198 + t180 * t193) / t167 ^ 2;
	t239 = t145 * t174;
	t238 = t145 * t185;
	t236 = t154 * t183;
	t234 = t158 * t182;
	t232 = t165 * t184;
	t231 = t170 * t171;
	t181 = t185 ^ 2;
	t229 = t170 * t181;
	t226 = t174 * t185;
	t225 = t175 * t179;
	t217 = qJD(1) * t183;
	t142 = t145 * t229 + 0.1e1;
	t215 = 0.2e1 * (-t229 * t242 + (t174 * t181 * t225 - t198) * t145) / t142 ^ 2;
	t214 = 0.2e1 * t242;
	t213 = -0.2e1 * t241;
	t212 = t145 * t226;
	t211 = t165 * t237;
	t205 = 0.1e1 + t230;
	t204 = t174 * t215;
	t203 = -0.2e1 * t174 * t240;
	t202 = t240 * t245;
	t199 = t155 * t161 * t231;
	t197 = t205 * t185;
	t196 = t200 * t185;
	t195 = t159 * t232 - t234;
	t164 = -t175 * t219 + t222;
	t151 = 0.1e1 / t153;
	t150 = t205 * t183 * t161;
	t140 = 0.1e1 / t142;
	t138 = (t243 * t174 * t154 - t183 * t199) * t185;
	t136 = -t175 * t236 + t235 + (t154 * t175 - t155 * t221) * t150;
	t135 = -t205 * t202 + (qJD(1) * t197 + t193 * t245) * t161;
	t132 = t195 * t213 * t226 + (t195 * t175 * t223 + (-t195 * t217 + ((-qJD(6) * t158 - 0.2e1 * t211) * t184 + (-t148 * t184 + (-qJD(6) * t165 + t149) * t182) * t159) * t185) * t174) * t151;
	t131 = (t136 * t239 - t144 * t175) * t185 * t215 + ((-t144 * t217 + (-t136 * t179 - t134) * t238) * t175 + (-t144 * t223 - (-t135 * t155 * t183 - t244 * t154 + (t139 * t236 - t154 * t179 - t155 * t216) * t150) * t212 + (t145 * t217 + t185 * t214) * t136 - ((t135 - t216) * t154 + ((-t150 * t183 + 0.1e1) * t179 + (t150 - t183) * t139) * t155) * t175 * t238) * t174) * t140;
	t1 = [t171 * t185 * t203 + (t179 * t197 - t217 * t227) * t161, t135, t135, 0, t135, 0; (t144 * t204 + (-t144 * t225 + (qJD(1) * t138 + t134) * t239) * t140) * t183 + (t145 * t204 * t138 + (-((t203 - t225 + (t139 * t171 * t228 + t225) * t161) * t154 + (t202 * t231 - t139 * t174 + (-t169 * t209 + (t139 - 0.2e1 * t224) * t174) * t161) * t155) * t212 + (-t145 * t225 + t174 * t214) * t138 + (-t144 + ((-t180 + t181) * t199 + t243 * t210) * t145) * t174 * qJD(1)) * t140) * t185, t131, t131, 0, t131, 0; 0.2e1 * (t158 * t194 + t164 * t233) * t241 + (0.2e1 * t164 * t211 - t201 * t158 * t219 + (t179 * t221 + t196) * t234 + (t164 * t148 + t194 * t149 - t196 * t232 - (t174 * t179 * t184 + t201 * t182) * t165 * t183) * t159) * t151, t132, t132, 0, t132, t213 + 0.2e1 * (-t148 * t151 * t159 + (-t151 * t237 - t159 * t241) * t165) * t165;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end