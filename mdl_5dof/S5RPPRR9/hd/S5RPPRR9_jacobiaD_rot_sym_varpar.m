% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RPPRR9
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
%   Wie in S5RPPRR9_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 16:21
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RPPRR9_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR9_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR9_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPRR9_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR9_jacobiaD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:21:28
	% EndTime: 2019-12-29 16:21:28
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:21:28
	% EndTime: 2019-12-29 16:21:28
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:21:28
	% EndTime: 2019-12-29 16:21:28
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:21:23
	% EndTime: 2019-12-29 16:21:23
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (52->9), mult. (138->18), div. (22->4), fcn. (156->4), ass. (0->15)
	t42 = sin(pkin(8));
	t43 = cos(pkin(8));
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
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0.2e1 * (t34 * t37 + t50) / t30 ^ 2 * (t32 * t49 + t48) + (-0.2e1 * t48 + (t34 - 0.2e1 * t49 - t51) * t32) / t30, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:21:28
	% EndTime: 2019-12-29 16:21:28
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:21:24
	% EndTime: 2019-12-29 16:21:25
	% DurationCPUTime: 1.70s
	% Computational Cost: add. (1789->95), mult. (4580->216), div. (480->12), fcn. (5898->11), ass. (0->92)
	t218 = sin(pkin(8));
	t219 = cos(pkin(8));
	t220 = sin(qJ(1));
	t221 = cos(qJ(1));
	t149 = t221 * t218 - t220 * t219;
	t146 = t149 ^ 2;
	t162 = sin(qJ(4));
	t157 = t162 ^ 2;
	t164 = cos(qJ(4));
	t159 = 0.1e1 / t164 ^ 2;
	t199 = t157 * t159;
	t142 = t146 * t199 + 0.1e1;
	t148 = -t220 * t218 - t221 * t219;
	t144 = t148 * qJD(1);
	t156 = t162 * t157;
	t158 = 0.1e1 / t164;
	t198 = t158 * t162;
	t175 = qJD(4) * (t156 * t158 * t159 + t198);
	t212 = (t149 * t144 * t199 + t146 * t175) / t142 ^ 2;
	t223 = -0.2e1 * t212;
	t184 = 0.1e1 + t199;
	t222 = t149 * t184;
	t201 = t149 * t162;
	t139 = atan2(t201, t164);
	t137 = sin(t139);
	t138 = cos(t139);
	t128 = t137 * t201 + t138 * t164;
	t125 = 0.1e1 / t128;
	t161 = sin(qJ(5));
	t163 = cos(qJ(5));
	t196 = t163 * t164;
	t136 = -t148 * t196 + t149 * t161;
	t130 = 0.1e1 / t136;
	t126 = 0.1e1 / t128 ^ 2;
	t131 = 0.1e1 / t136 ^ 2;
	t147 = t148 ^ 2;
	t202 = t147 * t157;
	t119 = t126 * t202 + 0.1e1;
	t145 = t149 * qJD(1);
	t192 = qJD(4) * t164;
	t140 = 0.1e1 / t142;
	t186 = t149 * t192;
	t194 = qJD(4) * t149;
	t187 = t159 * t194;
	t203 = t144 * t162;
	t114 = ((t186 + t203) * t158 + t157 * t187) * t140;
	t204 = t138 * t162;
	t110 = (t114 * t149 - qJD(4)) * t204 + (t203 + (-t114 + t194) * t164) * t137;
	t216 = t110 * t125 * t126;
	t217 = (-t202 * t216 + (-t145 * t148 * t157 + t147 * t162 * t192) * t126) / t119 ^ 2;
	t123 = t140 * t222;
	t205 = t137 * t164;
	t112 = t149 * t205 - t204 + (t138 * t201 - t205) * t123;
	t215 = t112 * t126;
	t200 = t157 * t158;
	t188 = t149 * t200;
	t181 = t138 * t188;
	t206 = t137 * t162;
	t113 = (t206 + (t181 - t206) * t140) * t148;
	t214 = t113 * t126;
	t193 = qJD(4) * t162;
	t172 = qJD(5) * t149 + t145 * t164 + t148 * t193;
	t191 = qJD(5) * t164;
	t178 = t148 * t191 + t144;
	t116 = t178 * t161 + t172 * t163;
	t213 = t116 * t130 * t131;
	t115 = t172 * t161 - t178 * t163;
	t197 = t161 * t164;
	t135 = -t148 * t197 - t149 * t163;
	t129 = t135 ^ 2;
	t124 = t129 * t131 + 0.1e1;
	t208 = t131 * t135;
	t211 = 0.1e1 / t124 ^ 2 * (t115 * t208 - t129 * t213);
	t210 = t126 * t148;
	t209 = t130 * t161;
	t207 = t135 * t163;
	t195 = t123 - t149;
	t190 = -0.2e1 * t216;
	t189 = 0.2e1 * t211;
	t185 = t123 * t149 - 0.1e1;
	t183 = -0.2e1 * t162 * t217;
	t182 = 0.2e1 * t135 * t213;
	t177 = t149 * t191 + t145;
	t176 = t131 * t207 - t209;
	t174 = t176 * t162;
	t173 = qJD(5) * t148 + t144 * t164 - t149 * t193;
	t134 = t148 * t161 + t149 * t196;
	t133 = -t148 * t163 + t149 * t197;
	t121 = 0.1e1 / t124;
	t117 = 0.1e1 / t119;
	t109 = t222 * t223 + (t184 * t144 + 0.2e1 * t149 * t175) * t140;
	t1 = [t148 * t198 * t223 + (t184 * t148 * qJD(4) - t145 * t198) * t140, 0, 0, t109, 0; t149 * t125 * t183 + (t125 * t186 + (t125 * t144 + (-t110 * t149 - t113 * t145) * t126) * t162) * t117 + (t183 * t214 + (t192 * t214 + (t113 * t190 + ((t137 * t192 + t181 * t223) * t148 + (-t137 * t145 + (t114 * t138 + 0.2e1 * t137 * t212) * t148) * t162 + (((-t114 * t188 - t192) * t148 + t145 * t162) * t137 + (-t145 * t188 + (t156 * t187 + t144 * t200 + (-t114 + 0.2e1 * t194) * t162) * t148) * t138) * t140) * t126) * t162) * t117) * t148, 0, 0, 0.2e1 * (t125 * t164 - t162 * t215) * t148 * t217 + ((t145 * t125 + (qJD(4) * t112 + t110) * t210) * t164 + (-t145 * t215 + (qJD(4) * t125 + t112 * t190 + ((t109 * t149 + t123 * t144) * t204 + (t195 * qJD(4) - t185 * t114) * t206) * t126) * t148 + ((-t109 + t144) * t137 + (t185 * qJD(4) - t195 * t114) * t138) * t164 * t210) * t162) * t117, 0; (-t130 * t133 + t134 * t208) * t189 + (t134 * t182 + t177 * t130 * t163 + t173 * t209 + (t177 * t135 * t161 - t134 * t115 - t133 * t116 - t173 * t207) * t131) * t121, 0, 0, t148 * t174 * t189 + (t145 * t174 + (-t176 * t192 + ((qJD(5) * t130 + t182) * t163 + (-t115 * t163 + (qJD(5) * t135 - t116) * t161) * t131) * t162) * t148) * t121, -0.2e1 * t211 + 0.2e1 * (t115 * t131 * t121 + (-t121 * t213 - t131 * t211) * t135) * t135;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end