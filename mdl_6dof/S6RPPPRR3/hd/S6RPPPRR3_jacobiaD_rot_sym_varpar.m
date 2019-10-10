% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPPPRR3
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
%   Wie in S6RPPPRR3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3,theta4]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:29
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPPPRR3_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR3_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR3_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPPRR3_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR3_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:29:15
	% EndTime: 2019-10-09 23:29:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:29:15
	% EndTime: 2019-10-09 23:29:15
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:29:15
	% EndTime: 2019-10-09 23:29:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:29:15
	% EndTime: 2019-10-09 23:29:15
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (52->9), mult. (138->18), div. (22->4), fcn. (156->4), ass. (0->15)
	t42 = sin(pkin(9));
	t43 = cos(pkin(9));
	t44 = sin(qJ(1));
	t45 = cos(qJ(1));
	t37 = -t44 * t42 - t45 * t43;
	t35 = 0.1e1 / t37 ^ 2;
	t38 = t45 * t42 - t44 * t43;
	t50 = t38 ^ 2 * t35;
	t51 = t35 * t37;
	t32 = t38 * qJD(1);
	t34 = 0.1e1 / t37;
	t49 = t34 * t50;
	t48 = t32 * t51;
	t30 = 0.1e1 + t50;
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0.2e1 * (t34 * t37 + t50) / t30 ^ 2 * (t32 * t49 + t48) + (-0.2e1 * t48 + (t34 - 0.2e1 * t49 - t51) * t32) / t30, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:29:15
	% EndTime: 2019-10-09 23:29:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:29:15
	% EndTime: 2019-10-09 23:29:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:29:15
	% EndTime: 2019-10-09 23:29:16
	% DurationCPUTime: 1.12s
	% Computational Cost: add. (3057->96), mult. (4580->216), div. (480->12), fcn. (5898->11), ass. (0->93)
	t226 = sin(pkin(9));
	t227 = cos(pkin(9));
	t228 = sin(qJ(1));
	t229 = cos(qJ(1));
	t156 = t229 * t226 - t228 * t227;
	t153 = t156 ^ 2;
	t170 = pkin(10) + qJ(5);
	t168 = sin(t170);
	t164 = t168 ^ 2;
	t169 = cos(t170);
	t166 = 0.1e1 / t169 ^ 2;
	t207 = t164 * t166;
	t149 = t153 * t207 + 0.1e1;
	t155 = -t228 * t226 - t229 * t227;
	t151 = t155 * qJD(1);
	t163 = t168 * t164;
	t165 = 0.1e1 / t169;
	t206 = t165 * t168;
	t183 = qJD(5) * (t163 * t165 * t166 + t206);
	t220 = (t156 * t151 * t207 + t153 * t183) / t149 ^ 2;
	t231 = -0.2e1 * t220;
	t192 = 0.1e1 + t207;
	t230 = t156 * t192;
	t209 = t156 * t168;
	t146 = atan2(t209, t169);
	t144 = sin(t146);
	t145 = cos(t146);
	t132 = t144 * t209 + t145 * t169;
	t129 = 0.1e1 / t132;
	t171 = sin(qJ(6));
	t172 = cos(qJ(6));
	t204 = t169 * t172;
	t143 = -t155 * t204 + t156 * t171;
	t137 = 0.1e1 / t143;
	t130 = 0.1e1 / t132 ^ 2;
	t138 = 0.1e1 / t143 ^ 2;
	t154 = t155 ^ 2;
	t210 = t154 * t164;
	t126 = t130 * t210 + 0.1e1;
	t152 = t156 * qJD(1);
	t200 = qJD(5) * t169;
	t147 = 0.1e1 / t149;
	t194 = t156 * t200;
	t202 = qJD(5) * t156;
	t195 = t166 * t202;
	t211 = t151 * t168;
	t121 = ((t194 + t211) * t165 + t164 * t195) * t147;
	t212 = t145 * t168;
	t117 = (t121 * t156 - qJD(5)) * t212 + (t211 + (-t121 + t202) * t169) * t144;
	t224 = t117 * t129 * t130;
	t225 = (-t210 * t224 + (-t152 * t155 * t164 + t154 * t168 * t200) * t130) / t126 ^ 2;
	t128 = t147 * t230;
	t213 = t144 * t169;
	t119 = t156 * t213 - t212 + (t145 * t209 - t213) * t128;
	t223 = t119 * t130;
	t208 = t164 * t165;
	t196 = t156 * t208;
	t189 = t145 * t196;
	t214 = t144 * t168;
	t120 = (t214 + (t189 - t214) * t147) * t155;
	t222 = t120 * t130;
	t201 = qJD(5) * t168;
	t180 = qJD(6) * t156 + t152 * t169 + t155 * t201;
	t199 = qJD(6) * t169;
	t186 = t155 * t199 + t151;
	t123 = t186 * t171 + t180 * t172;
	t221 = t123 * t137 * t138;
	t219 = t130 * t155;
	t122 = t180 * t171 - t186 * t172;
	t205 = t169 * t171;
	t142 = -t155 * t205 - t156 * t172;
	t136 = t142 ^ 2;
	t135 = t136 * t138 + 0.1e1;
	t216 = t138 * t142;
	t218 = 0.1e1 / t135 ^ 2 * (t122 * t216 - t136 * t221);
	t217 = t137 * t171;
	t215 = t142 * t172;
	t203 = t128 - t156;
	t198 = -0.2e1 * t224;
	t197 = 0.2e1 * t218;
	t193 = t128 * t156 - 0.1e1;
	t191 = -0.2e1 * t168 * t225;
	t190 = 0.2e1 * t142 * t221;
	t185 = t156 * t199 + t152;
	t184 = t138 * t215 - t217;
	t182 = t184 * t168;
	t181 = qJD(6) * t155 + t151 * t169 - t156 * t201;
	t141 = t155 * t171 + t156 * t204;
	t140 = -t155 * t172 + t156 * t205;
	t133 = 0.1e1 / t135;
	t124 = 0.1e1 / t126;
	t116 = t230 * t231 + (t192 * t151 + 0.2e1 * t156 * t183) * t147;
	t1 = [t155 * t206 * t231 + (t192 * t155 * qJD(5) - t152 * t206) * t147, 0, 0, 0, t116, 0; t156 * t129 * t191 + (t129 * t194 + (t129 * t151 + (-t117 * t156 - t120 * t152) * t130) * t168) * t124 + (t191 * t222 + (t200 * t222 + (t120 * t198 + ((t144 * t200 + t189 * t231) * t155 + (-t144 * t152 + (t121 * t145 + 0.2e1 * t144 * t220) * t155) * t168 + (((-t121 * t196 - t200) * t155 + t152 * t168) * t144 + (-t152 * t196 + (t163 * t195 + t151 * t208 + (-t121 + 0.2e1 * t202) * t168) * t155) * t145) * t147) * t130) * t168) * t124) * t155, 0, 0, 0, 0.2e1 * (t129 * t169 - t168 * t223) * t155 * t225 + ((t152 * t129 + (qJD(5) * t119 + t117) * t219) * t169 + (-t152 * t223 + (qJD(5) * t129 + t119 * t198 + ((t116 * t156 + t128 * t151) * t212 + (t203 * qJD(5) - t193 * t121) * t214) * t130) * t155 + ((-t116 + t151) * t144 + (t193 * qJD(5) - t203 * t121) * t145) * t169 * t219) * t168) * t124, 0; (-t137 * t140 + t141 * t216) * t197 + (t141 * t190 + t185 * t137 * t172 + t181 * t217 + (t185 * t142 * t171 - t141 * t122 - t140 * t123 - t181 * t215) * t138) * t133, 0, 0, 0, t155 * t182 * t197 + (t152 * t182 + (-t184 * t200 + ((qJD(6) * t137 + t190) * t172 + (-t122 * t172 + (qJD(6) * t142 - t123) * t171) * t138) * t168) * t155) * t133, -0.2e1 * t218 + 0.2e1 * (t122 * t138 * t133 + (-t133 * t221 - t138 * t218) * t142) * t142;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end