% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRPRR7
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
%   Wie in S6RPRPRR7_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:56
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPRR7_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR7_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR7_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR7_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR7_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:56:39
	% EndTime: 2019-10-10 00:56:39
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:56:39
	% EndTime: 2019-10-10 00:56:39
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:56:39
	% EndTime: 2019-10-10 00:56:39
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:56:39
	% EndTime: 2019-10-10 00:56:39
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:56:39
	% EndTime: 2019-10-10 00:56:39
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:56:39
	% EndTime: 2019-10-10 00:56:39
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:56:39
	% EndTime: 2019-10-10 00:56:40
	% DurationCPUTime: 1.05s
	% Computational Cost: add. (5302->95), mult. (3810->207), div. (753->12), fcn. (4455->9), ass. (0->97)
	t169 = cos(qJ(1));
	t229 = 0.2e1 * t169;
	t165 = t169 ^ 2;
	t162 = qJ(3) + pkin(10) + qJ(5);
	t159 = sin(t162);
	t155 = 0.1e1 / t159 ^ 2;
	t160 = cos(t162);
	t158 = t160 ^ 2;
	t212 = t155 * t158;
	t153 = t165 * t212 + 0.1e1;
	t151 = 0.1e1 / t153;
	t154 = 0.1e1 / t159;
	t167 = sin(qJ(1));
	t201 = qJD(1) * t167;
	t193 = t160 * t201;
	t163 = qJD(3) + qJD(5);
	t207 = t163 * t169;
	t194 = t155 * t207;
	t125 = ((t159 * t207 + t193) * t154 + t158 * t194) * t151;
	t228 = -t125 + t207;
	t203 = t169 * t160;
	t150 = atan2(-t203, t159);
	t140 = sin(t150);
	t141 = cos(t150);
	t135 = -t140 * t203 + t141 * t159;
	t132 = 0.1e1 / t135;
	t168 = cos(qJ(6));
	t204 = t167 * t168;
	t166 = sin(qJ(6));
	t206 = t166 * t169;
	t147 = t159 * t204 + t206;
	t143 = 0.1e1 / t147;
	t133 = 0.1e1 / t135 ^ 2;
	t144 = 0.1e1 / t147 ^ 2;
	t164 = t167 ^ 2;
	t211 = t158 * t164;
	t128 = t133 * t211 + 0.1e1;
	t200 = qJD(1) * t169;
	t185 = t158 * t167 * t200;
	t209 = t160 * t163;
	t218 = t141 * t160;
	t224 = t125 * t169;
	t120 = (t163 - t224) * t218 + (t159 * t228 + t193) * t140;
	t226 = t120 * t132 * t133;
	t227 = (-t211 * t226 + (-t159 * t164 * t209 + t185) * t133) / t128 ^ 2;
	t188 = qJD(1) * t159 + qJD(6);
	t181 = t188 * t169;
	t189 = qJD(6) * t159 + qJD(1);
	t182 = t189 * t168;
	t208 = t163 * t167;
	t130 = t167 * t182 + (t160 * t208 + t181) * t166;
	t202 = t169 * t168;
	t205 = t167 * t166;
	t146 = t159 * t205 - t202;
	t142 = t146 ^ 2;
	t139 = t142 * t144 + 0.1e1;
	t216 = t144 * t146;
	t183 = t189 * t166;
	t131 = t168 * t181 + (t168 * t209 - t183) * t167;
	t222 = t131 * t143 * t144;
	t225 = (t130 * t216 - t142 * t222) / t139 ^ 2;
	t157 = t160 * t158;
	t213 = t154 * t160;
	t179 = t163 * (-t154 * t155 * t157 - t213);
	t223 = (-t155 * t185 + t165 * t179) / t153 ^ 2;
	t221 = t133 * t160;
	t220 = t133 * t167;
	t219 = t140 * t169;
	t217 = t143 * t166;
	t215 = t146 * t168;
	t214 = t154 * t158;
	t210 = t159 * t163;
	t199 = -0.2e1 * t226;
	t198 = 0.2e1 * t225;
	t197 = t160 * t227;
	t196 = t160 * t223;
	t195 = t160 * t220;
	t192 = 0.1e1 + t212;
	t191 = t223 * t229;
	t190 = 0.2e1 * t146 * t222;
	t187 = t141 * t151 * t214;
	t186 = (-t151 + 0.1e1) * t160 * t140;
	t184 = t192 * t167;
	t180 = t144 * t215 - t217;
	t178 = t180 * t167;
	t177 = t163 * t203 - t188 * t167;
	t149 = t159 * t202 - t205;
	t148 = t159 * t206 + t204;
	t137 = 0.1e1 / t139;
	t136 = t192 * t169 * t151;
	t126 = 0.1e1 / t128;
	t124 = (-t169 * t187 + t186) * t167;
	t122 = t159 * t219 + t218 + (-t140 * t159 - t141 * t203) * t136;
	t121 = -t192 * t191 + (-qJD(1) * t184 + t179 * t229) * t151;
	t118 = t160 * t178 * t198 + (t178 * t210 + (-t180 * t200 + ((qJD(6) * t143 + t190) * t168 + (-t130 * t168 + (qJD(6) * t146 - t131) * t166) * t144) * t167) * t160) * t137;
	t117 = 0.2e1 * (-t122 * t221 - t132 * t159) * t167 * t227 + ((t132 * t200 + (-t122 * t163 - t120) * t220) * t159 + (t132 * t208 + (-t121 * t141 * t169 + t228 * t140 + (t125 * t219 - t140 * t163 + t141 * t201) * t136) * t195 + (t133 * t200 + t167 * t199) * t122 + ((-t121 - t201) * t140 + ((t136 * t169 - 0.1e1) * t163 + (-t136 + t169) * t125) * t141) * t159 * t220) * t160) * t126;
	t1 = [-0.2e1 * t167 * t154 * t196 + (-t163 * t184 + t200 * t213) * t151, 0, t121, 0, t121, 0; (0.2e1 * t132 * t197 + (t132 * t210 + (qJD(1) * t124 + t120) * t221) * t126) * t169 + (-0.2e1 * t133 * t197 * t124 + (((0.2e1 * t196 - t210 + (t214 * t224 + t210) * t151) * t140 + (t191 * t214 + t125 * t160 + (t157 * t194 + (-t125 + 0.2e1 * t207) * t160) * t151) * t141) * t195 + (-t133 * t210 + t160 * t199) * t124 + (t132 + ((t164 - t165) * t187 + t169 * t186) * t133) * t160 * qJD(1)) * t126) * t167, 0, t117, 0, t117, 0; (-t143 * t148 + t149 * t216) * t198 + (t149 * t190 + t169 * t143 * t182 + t177 * t217 + (t169 * t146 * t183 - t149 * t130 - t148 * t131 - t177 * t215) * t144) * t137, 0, t118, 0, t118, -0.2e1 * t225 + 0.2e1 * (t130 * t137 * t144 + (-t137 * t222 - t144 * t225) * t146) * t146;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end