% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPPPRR5
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
%   Wie in S6RPPPRR5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta4]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:32
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPPPRR5_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR5_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR5_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPPRR5_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR5_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:32:36
	% EndTime: 2019-10-09 23:32:36
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:32:36
	% EndTime: 2019-10-09 23:32:36
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:32:36
	% EndTime: 2019-10-09 23:32:36
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:32:36
	% EndTime: 2019-10-09 23:32:36
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:32:36
	% EndTime: 2019-10-09 23:32:36
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (46->7), mult. (124->20), div. (18->4), fcn. (140->4), ass. (0->15)
	t36 = cos(pkin(9));
	t37 = cos(qJ(1));
	t40 = sin(pkin(9));
	t44 = sin(qJ(1));
	t33 = t44 * t36 + t37 * t40;
	t29 = 0.1e1 / t33 ^ 2;
	t45 = qJD(1) * t29;
	t32 = t37 * t36 - t44 * t40;
	t31 = t32 ^ 2;
	t26 = t31 * t29 + 0.1e1;
	t41 = t32 / t33 * t45;
	t42 = t33 * t45;
	t43 = (-t31 * t41 - t32 * t42) / t26 ^ 2;
	t24 = 0.1e1 / t26;
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -0.2e1 * t43 + 0.2e1 * (-t24 * t42 + (-t24 * t41 - t29 * t43) * t32) * t32, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:32:36
	% EndTime: 2019-10-09 23:32:36
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:32:36
	% EndTime: 2019-10-09 23:32:37
	% DurationCPUTime: 1.07s
	% Computational Cost: add. (1978->96), mult. (4580->220), div. (480->12), fcn. (5898->11), ass. (0->92)
	t150 = cos(pkin(9));
	t155 = cos(qJ(1));
	t208 = sin(pkin(9));
	t209 = sin(qJ(1));
	t141 = t155 * t150 - t209 * t208;
	t152 = sin(qJ(5));
	t146 = t152 ^ 2;
	t154 = cos(qJ(5));
	t148 = 0.1e1 / t154 ^ 2;
	t188 = t146 * t148;
	t170 = 0.1e1 + t188;
	t210 = t141 * t170;
	t190 = t141 * t152;
	t132 = atan2(t190, -t154);
	t130 = sin(t132);
	t131 = cos(t132);
	t121 = t130 * t190 - t131 * t154;
	t118 = 0.1e1 / t121;
	t142 = t209 * t150 + t155 * t208;
	t151 = sin(qJ(6));
	t153 = cos(qJ(6));
	t185 = t153 * t154;
	t129 = -t141 * t151 + t142 * t185;
	t123 = 0.1e1 / t129;
	t147 = 0.1e1 / t154;
	t119 = 0.1e1 / t121 ^ 2;
	t124 = 0.1e1 / t129 ^ 2;
	t139 = t142 ^ 2;
	t191 = t139 * t146;
	t112 = t119 * t191 + 0.1e1;
	t138 = t141 * qJD(1);
	t181 = qJD(5) * t154;
	t140 = t141 ^ 2;
	t135 = t140 * t188 + 0.1e1;
	t133 = 0.1e1 / t135;
	t172 = t141 * t181;
	t183 = qJD(5) * t141;
	t173 = t148 * t183;
	t137 = t142 * qJD(1);
	t193 = t137 * t152;
	t107 = (-(t172 - t193) * t147 - t146 * t173) * t133;
	t194 = t131 * t152;
	t103 = (t107 * t141 + qJD(5)) * t194 + (-t193 + (t107 + t183) * t154) * t130;
	t206 = t103 * t118 * t119;
	t207 = (-t191 * t206 + (t138 * t142 * t146 + t139 * t152 * t181) * t119) / t112 ^ 2;
	t182 = qJD(5) * t152;
	t163 = -qJD(6) * t141 + t138 * t154 - t142 * t182;
	t180 = qJD(6) * t154;
	t168 = t142 * t180 - t137;
	t108 = t163 * t151 + t168 * t153;
	t186 = t151 * t154;
	t128 = t141 * t153 + t142 * t186;
	t122 = t128 ^ 2;
	t117 = t122 * t124 + 0.1e1;
	t198 = t124 * t128;
	t109 = -t168 * t151 + t163 * t153;
	t202 = t109 * t123 * t124;
	t205 = (t108 * t198 - t122 * t202) / t117 ^ 2;
	t116 = t133 * t210;
	t195 = t130 * t154;
	t105 = t141 * t195 + t194 - (t131 * t190 + t195) * t116;
	t204 = t105 * t119;
	t189 = t146 * t147;
	t174 = t141 * t189;
	t169 = t131 * t174;
	t196 = t130 * t152;
	t106 = (-t196 + (t169 + t196) * t133) * t142;
	t203 = t106 * t119;
	t145 = t152 * t146;
	t187 = t147 * t152;
	t166 = t145 * t147 * t148 + t187;
	t201 = (t166 * t140 * qJD(5) - t141 * t137 * t188) / t135 ^ 2;
	t200 = t119 * t142;
	t199 = t123 * t151;
	t197 = t128 * t153;
	t192 = t138 * t152;
	t184 = -t116 + t141;
	t179 = 0.2e1 * t206;
	t178 = -0.2e1 * t205;
	t177 = -0.2e1 * t201;
	t176 = t152 * t207;
	t175 = t128 * t202;
	t171 = -t116 * t141 + 0.1e1;
	t167 = t141 * t180 - t138;
	t165 = t124 * t197 - t199;
	t164 = qJD(6) * t142 - t137 * t154 - t141 * t182;
	t127 = t141 * t185 + t142 * t151;
	t126 = t141 * t186 - t142 * t153;
	t114 = 0.1e1 / t117;
	t110 = 0.1e1 / t112;
	t102 = 0.2e1 * t201 * t210 + (t170 * t137 - 0.2e1 * t166 * t183) * t133;
	t1 = [t142 * t177 * t187 + (t170 * t142 * qJD(5) + t138 * t187) * t133, 0, 0, 0, t102, 0; -0.2e1 * t141 * t118 * t176 + (t118 * t172 + (-t118 * t137 + (-t103 * t141 - t106 * t138) * t119) * t152) * t110 + (0.2e1 * t176 * t203 + (-t181 * t203 + (t106 * t179 + (-(-t130 * t181 + t169 * t177) * t142 + (t130 * t138 - (-t107 * t131 + t130 * t177) * t142) * t152 + ((-(-t107 * t174 + t181) * t142 - t192) * t130 + (-t138 * t174 + (-t145 * t173 + t137 * t189 - (t107 + 0.2e1 * t183) * t152) * t142) * t131) * t133) * t119) * t152) * t110) * t142, 0, 0, 0, 0.2e1 * (-t118 * t154 + t152 * t204) * t142 * t207 + ((t138 * t118 + (-qJD(5) * t105 - t103) * t200) * t154 + (-t138 * t204 + (-qJD(5) * t118 + t105 * t179 + (-(t102 * t141 + t116 * t137) * t194 - (-t184 * qJD(5) - t171 * t107) * t196) * t119) * t142 - ((t102 - t137) * t130 + (t171 * qJD(5) + t184 * t107) * t131) * t154 * t200) * t152) * t110, 0; 0.2e1 * (-t123 * t126 + t127 * t198) * t205 + (0.2e1 * t127 * t175 + t167 * t123 * t153 + t164 * t199 + (t167 * t128 * t151 - t127 * t108 - t126 * t109 - t164 * t197) * t124) * t114, 0, 0, 0, t165 * t152 * t142 * t178 + (t165 * t192 + (t165 * t181 + ((-qJD(6) * t123 - 0.2e1 * t175) * t153 + (t108 * t153 + (-qJD(6) * t128 + t109) * t151) * t124) * t152) * t142) * t114, t178 + 0.2e1 * (t108 * t124 * t114 + (-t114 * t202 - t124 * t205) * t128) * t128;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end