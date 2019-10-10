% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PRRRR1
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
%   Wie in S5PRRRR1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 20:54
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5PRRRR1_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR1_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRR1_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_jacobiaD_rot_sym_varpar: pkin has to be [2x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:54:42
	% EndTime: 2019-10-09 20:54:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:54:42
	% EndTime: 2019-10-09 20:54:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:54:42
	% EndTime: 2019-10-09 20:54:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:54:42
	% EndTime: 2019-10-09 20:54:42
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (30->4), mult. (82->10), div. (42->6), fcn. (80->4), ass. (0->8)
	t50 = sin(qJ(2));
	t52 = cos(qJ(2));
	t61 = t50 ^ 2 / t52 ^ 2;
	t64 = 0.1e1 + t61;
	t38 = cos(atan2(0, -t52));
	t36 = 0.1e1 / t38 ^ 2;
	t32 = t36 * t61 + 0.1e1;
	t1 = [0, 0, 0, 0, 0; 0, 0.2e1 * (-0.1e1 / t32 + t64 / t32 ^ 2 * t36) * t64 / t52 * t50 * qJD(2) / t38, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:54:42
	% EndTime: 2019-10-09 20:54:42
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (95->4), mult. (119->10), div. (57->6), fcn. (102->4), ass. (0->8)
	t69 = sin(qJ(2));
	t70 = cos(qJ(2));
	t79 = t69 ^ 2 / t70 ^ 2;
	t82 = 0.1e1 + t79;
	t54 = cos(atan2(0, -t70));
	t53 = 0.1e1 / t54 ^ 2;
	t48 = t53 * t79 + 0.1e1;
	t1 = [0, 0, 0, 0, 0; 0, 0.2e1 * (-0.1e1 / t48 + t82 / t48 ^ 2 * t53) * t82 / t70 * t69 * qJD(2) / t54, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:54:42
	% EndTime: 2019-10-09 20:54:43
	% DurationCPUTime: 1.13s
	% Computational Cost: add. (3339->97), mult. (3536->225), div. (1131->15), fcn. (4495->9), ass. (0->97)
	t167 = qJ(3) + qJ(4);
	t160 = sin(t167);
	t155 = 0.1e1 / t160;
	t169 = sin(qJ(2));
	t163 = 0.1e1 / t169;
	t164 = 0.1e1 / t169 ^ 2;
	t165 = t163 * t164;
	t171 = cos(qJ(2));
	t166 = t171 ^ 2;
	t232 = t155 * (0.2e1 * t165 * t166 + t163);
	t156 = 0.1e1 / t160 ^ 2;
	t161 = cos(t167);
	t159 = t161 ^ 2;
	t216 = t159 * t164;
	t150 = t156 * t216 + 0.1e1;
	t148 = 0.1e1 / t150;
	t162 = qJD(3) + qJD(4);
	t201 = qJD(2) * t171;
	t192 = t160 * t201;
	t210 = t162 * t169;
	t122 = (t162 * t163 + (t161 * t210 + t192) * t161 * t164 * t156) * t148;
	t231 = t122 - t210;
	t170 = cos(qJ(5));
	t204 = t170 * t171;
	t168 = sin(qJ(5));
	t206 = t169 * t168;
	t144 = t161 * t204 + t206;
	t207 = t169 * t160;
	t147 = atan2(-t161, t207);
	t146 = cos(t147);
	t195 = t146 * t207;
	t145 = sin(t147);
	t221 = t145 * t161;
	t132 = t195 - t221;
	t129 = 0.1e1 / t132;
	t138 = 0.1e1 / t144;
	t130 = 0.1e1 / t132 ^ 2;
	t139 = 0.1e1 / t144 ^ 2;
	t222 = t145 * t160;
	t117 = (-t122 * t169 + t162) * t222 + (-t231 * t161 + t192) * t146;
	t230 = t117 * t129 * t130;
	t203 = qJD(2) * t161;
	t187 = qJD(5) - t203;
	t188 = qJD(5) * t161 - qJD(2);
	t209 = t162 * t171;
	t193 = t160 * t209;
	t208 = t168 * t171;
	t127 = -t188 * t208 + (t187 * t169 - t193) * t170;
	t229 = t127 * t138 * t139;
	t154 = t160 ^ 2;
	t158 = t161 * t159;
	t180 = t162 * (-t161 - 0.1e1 / t154 * t158) * t155;
	t217 = t156 * t159;
	t228 = (-t165 * t201 * t217 + t164 * t180) / t150 ^ 2;
	t227 = t130 * t160;
	t181 = t161 * t206 + t204;
	t126 = t181 * qJD(2) - t144 * qJD(5) + t168 * t193;
	t205 = t169 * t170;
	t143 = t161 * t208 - t205;
	t137 = t143 ^ 2;
	t136 = t137 * t139 + 0.1e1;
	t224 = t139 * t143;
	t226 = 0.1e1 / t136 ^ 2 * (-t126 * t224 - t137 * t229);
	t225 = t138 * t168;
	t223 = t143 * t170;
	t220 = t146 * t169;
	t219 = t148 * t163;
	t218 = t154 * t166;
	t215 = t160 * t162;
	t214 = t160 * t171;
	t213 = t161 * t162;
	t212 = t161 * t171;
	t211 = t162 * t164;
	t202 = qJD(2) * t169;
	t125 = t130 * t218 + 0.1e1;
	t200 = 0.2e1 * (-t218 * t230 + (-t154 * t169 * t201 + t160 * t166 * t213) * t130) / t125 ^ 2;
	t199 = 0.2e1 * t230;
	t198 = -0.2e1 * t226;
	t197 = t143 * t229;
	t196 = t130 * t214;
	t191 = 0.1e1 + t217;
	t190 = t160 * t200;
	t186 = t155 * t164 * t228;
	t185 = t187 * t171;
	t184 = t148 * t155 * t216 - t160;
	t183 = t164 * t171 * t191;
	t182 = t139 * t223 - t225;
	t142 = -t161 * t205 + t208;
	t134 = 0.1e1 / t136;
	t133 = t191 * t219;
	t123 = 0.1e1 / t125;
	t121 = (-t184 * t146 - t219 * t221) * t171;
	t119 = t161 * t220 + t222 + (-t145 * t207 - t146 * t161) * t133;
	t118 = -0.2e1 * t191 * t163 * t228 + (-qJD(2) * t183 + 0.2e1 * t163 * t180) * t148;
	t115 = t182 * t198 * t214 + (t182 * t161 * t209 + (-t182 * t202 + ((-qJD(5) * t138 - 0.2e1 * t197) * t170 + (-t126 * t170 + (-qJD(5) * t143 + t127) * t168) * t139) * t171) * t160) * t134;
	t114 = (t119 * t227 - t129 * t161) * t171 * t200 + ((-t129 * t202 + (-t119 * t162 - t117) * t171 * t130) * t161 + (-t129 * t209 - (-t118 * t145 * t169 + t231 * t146 + (-t122 * t220 - t145 * t201 + t146 * t162) * t133) * t196 + (t130 * t202 + t171 * t199) * t119 - ((-t118 + t201) * t146 + ((-t133 * t169 + 0.1e1) * t162 + (t133 - t169) * t122) * t145) * t130 * t212) * t160) * t123;
	t1 = [0, -0.2e1 * t186 * t212 + (-t162 * t183 - t203 * t232) * t148, t118, t118, 0; 0, (t129 * t190 + (-t129 * t213 + (qJD(2) * t121 + t117) * t227) * t123) * t169 + (t130 * t190 * t121 + (-(((t148 * t215 + 0.2e1 * t161 * t228) * t163 + t184 * t122) * t145 + (0.2e1 * t159 * t186 + t213 + (t156 * t158 * t211 + (-t122 * t163 + 0.2e1 * t211) * t161) * t148) * t146) * t196 + (-t130 * t213 + t160 * t199) * t121 + ((t130 * t195 - t129) * t160 - (t159 * t146 * t232 + (t164 * t166 + 0.1e1) * t221) * t148 * t227) * qJD(2)) * t123) * t171, t114, t114, 0; 0, 0.2e1 * (t138 * t181 + t142 * t224) * t226 + (0.2e1 * t142 * t197 - t188 * t138 * t205 + (t162 * t207 + t185) * t225 + (t142 * t126 + t181 * t127 - t185 * t223 - (t188 * t168 + t170 * t215) * t143 * t169) * t139) * t134, t115, t115, t198 + 0.2e1 * (-t126 * t139 * t134 + (-t134 * t229 - t139 * t226) * t143) * t143;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end