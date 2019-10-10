% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRPR7
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
%   Wie in S6RPRRPR7_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:33
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRPR7_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR7_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR7_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR7_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR7_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:33:56
	% EndTime: 2019-10-10 01:33:56
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:33:56
	% EndTime: 2019-10-10 01:33:56
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:33:56
	% EndTime: 2019-10-10 01:33:56
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:33:56
	% EndTime: 2019-10-10 01:33:56
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:33:56
	% EndTime: 2019-10-10 01:33:56
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:33:56
	% EndTime: 2019-10-10 01:33:56
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:33:56
	% EndTime: 2019-10-10 01:33:57
	% DurationCPUTime: 1.05s
	% Computational Cost: add. (5302->95), mult. (3810->207), div. (753->12), fcn. (4455->9), ass. (0->97)
	t173 = cos(qJ(1));
	t233 = 0.2e1 * t173;
	t169 = t173 ^ 2;
	t166 = qJ(3) + qJ(4) + pkin(10);
	t163 = sin(t166);
	t159 = 0.1e1 / t163 ^ 2;
	t164 = cos(t166);
	t162 = t164 ^ 2;
	t216 = t159 * t162;
	t157 = t169 * t216 + 0.1e1;
	t155 = 0.1e1 / t157;
	t158 = 0.1e1 / t163;
	t171 = sin(qJ(1));
	t205 = qJD(1) * t171;
	t197 = t164 * t205;
	t167 = qJD(3) + qJD(4);
	t211 = t167 * t173;
	t198 = t159 * t211;
	t129 = ((t163 * t211 + t197) * t158 + t162 * t198) * t155;
	t232 = -t129 + t211;
	t207 = t173 * t164;
	t154 = atan2(-t207, t163);
	t144 = sin(t154);
	t145 = cos(t154);
	t139 = -t144 * t207 + t145 * t163;
	t136 = 0.1e1 / t139;
	t172 = cos(qJ(6));
	t208 = t171 * t172;
	t170 = sin(qJ(6));
	t210 = t170 * t173;
	t151 = t163 * t208 + t210;
	t147 = 0.1e1 / t151;
	t137 = 0.1e1 / t139 ^ 2;
	t148 = 0.1e1 / t151 ^ 2;
	t168 = t171 ^ 2;
	t215 = t162 * t168;
	t132 = t137 * t215 + 0.1e1;
	t204 = qJD(1) * t173;
	t189 = t162 * t171 * t204;
	t213 = t164 * t167;
	t222 = t145 * t164;
	t228 = t129 * t173;
	t124 = (t167 - t228) * t222 + (t232 * t163 + t197) * t144;
	t230 = t124 * t136 * t137;
	t231 = (-t215 * t230 + (-t163 * t168 * t213 + t189) * t137) / t132 ^ 2;
	t192 = qJD(1) * t163 + qJD(6);
	t185 = t192 * t173;
	t193 = qJD(6) * t163 + qJD(1);
	t186 = t193 * t172;
	t212 = t167 * t171;
	t134 = t171 * t186 + (t164 * t212 + t185) * t170;
	t206 = t173 * t172;
	t209 = t171 * t170;
	t150 = t163 * t209 - t206;
	t146 = t150 ^ 2;
	t143 = t146 * t148 + 0.1e1;
	t220 = t148 * t150;
	t187 = t193 * t170;
	t135 = t172 * t185 + (t172 * t213 - t187) * t171;
	t226 = t135 * t147 * t148;
	t229 = (t134 * t220 - t146 * t226) / t143 ^ 2;
	t161 = t164 * t162;
	t217 = t158 * t164;
	t183 = t167 * (-t158 * t159 * t161 - t217);
	t227 = (-t159 * t189 + t169 * t183) / t157 ^ 2;
	t225 = t137 * t164;
	t224 = t137 * t171;
	t223 = t144 * t173;
	t221 = t147 * t170;
	t219 = t150 * t172;
	t218 = t158 * t162;
	t214 = t163 * t167;
	t203 = -0.2e1 * t230;
	t202 = 0.2e1 * t229;
	t201 = t164 * t231;
	t200 = t164 * t227;
	t199 = t164 * t224;
	t196 = 0.1e1 + t216;
	t195 = t227 * t233;
	t194 = 0.2e1 * t150 * t226;
	t191 = t145 * t155 * t218;
	t190 = (-t155 + 0.1e1) * t164 * t144;
	t188 = t196 * t171;
	t184 = t148 * t219 - t221;
	t182 = t184 * t171;
	t181 = t167 * t207 - t192 * t171;
	t153 = t163 * t206 - t209;
	t152 = t163 * t210 + t208;
	t141 = 0.1e1 / t143;
	t140 = t196 * t173 * t155;
	t130 = 0.1e1 / t132;
	t128 = (-t173 * t191 + t190) * t171;
	t126 = t163 * t223 + t222 + (-t144 * t163 - t145 * t207) * t140;
	t125 = -t196 * t195 + (-qJD(1) * t188 + t183 * t233) * t155;
	t122 = t164 * t182 * t202 + (t182 * t214 + (-t184 * t204 + ((qJD(6) * t147 + t194) * t172 + (-t134 * t172 + (qJD(6) * t150 - t135) * t170) * t148) * t171) * t164) * t141;
	t121 = 0.2e1 * (-t126 * t225 - t136 * t163) * t171 * t231 + ((t136 * t204 + (-t126 * t167 - t124) * t224) * t163 + (t136 * t212 + (-t125 * t145 * t173 + t232 * t144 + (t129 * t223 - t144 * t167 + t145 * t205) * t140) * t199 + (t137 * t204 + t171 * t203) * t126 + ((-t125 - t205) * t144 + ((t140 * t173 - 0.1e1) * t167 + (-t140 + t173) * t129) * t145) * t163 * t224) * t164) * t130;
	t1 = [-0.2e1 * t171 * t158 * t200 + (-t167 * t188 + t204 * t217) * t155, 0, t125, t125, 0, 0; (0.2e1 * t136 * t201 + (t136 * t214 + (qJD(1) * t128 + t124) * t225) * t130) * t173 + (-0.2e1 * t137 * t201 * t128 + (((0.2e1 * t200 - t214 + (t218 * t228 + t214) * t155) * t144 + (t195 * t218 + t129 * t164 + (t161 * t198 + (-t129 + 0.2e1 * t211) * t164) * t155) * t145) * t199 + (-t137 * t214 + t164 * t203) * t128 + (t136 + ((t168 - t169) * t191 + t173 * t190) * t137) * t164 * qJD(1)) * t130) * t171, 0, t121, t121, 0, 0; (-t147 * t152 + t153 * t220) * t202 + (t153 * t194 + t173 * t147 * t186 + t181 * t221 + (t173 * t150 * t187 - t153 * t134 - t152 * t135 - t181 * t219) * t148) * t141, 0, t122, t122, 0, -0.2e1 * t229 + 0.2e1 * (t134 * t141 * t148 + (-t141 * t226 - t148 * t229) * t150) * t150;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end