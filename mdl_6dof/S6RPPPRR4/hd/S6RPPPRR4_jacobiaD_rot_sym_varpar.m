% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPPPRR4
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
%   Wie in S6RPPPRR4_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:30
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPPPRR4_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR4_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR4_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPPRR4_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR4_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:30:55
	% EndTime: 2019-10-09 23:30:55
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:30:55
	% EndTime: 2019-10-09 23:30:55
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:30:55
	% EndTime: 2019-10-09 23:30:55
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:30:55
	% EndTime: 2019-10-09 23:30:55
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
	% StartTime: 2019-10-09 23:30:55
	% EndTime: 2019-10-09 23:30:55
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:30:55
	% EndTime: 2019-10-09 23:30:55
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:30:56
	% EndTime: 2019-10-09 23:30:57
	% DurationCPUTime: 1.04s
	% Computational Cost: add. (1789->97), mult. (4580->218), div. (480->12), fcn. (5898->11), ass. (0->93)
	t215 = sin(pkin(9));
	t216 = cos(pkin(9));
	t217 = sin(qJ(1));
	t218 = cos(qJ(1));
	t146 = -t217 * t215 - t218 * t216;
	t158 = sin(qJ(5));
	t153 = 0.1e1 / t158 ^ 2;
	t160 = cos(qJ(5));
	t156 = t160 ^ 2;
	t194 = t153 * t156;
	t178 = 0.1e1 + t194;
	t220 = t146 * t178;
	t147 = t218 * t215 - t217 * t216;
	t197 = t146 * t160;
	t137 = atan2(-t197, -t158);
	t135 = sin(t137);
	t136 = cos(t137);
	t126 = -t135 * t197 - t136 * t158;
	t123 = 0.1e1 / t126;
	t157 = sin(qJ(6));
	t159 = cos(qJ(6));
	t192 = t158 * t159;
	t134 = -t146 * t157 + t147 * t192;
	t128 = 0.1e1 / t134;
	t152 = 0.1e1 / t158;
	t124 = 0.1e1 / t126 ^ 2;
	t129 = 0.1e1 / t134 ^ 2;
	t219 = 0.2e1 * t147;
	t144 = t147 ^ 2;
	t198 = t144 * t156;
	t117 = t124 * t198 + 0.1e1;
	t142 = t146 * qJD(1);
	t188 = qJD(5) * t160;
	t145 = t146 ^ 2;
	t140 = t145 * t194 + 0.1e1;
	t138 = 0.1e1 / t140;
	t189 = qJD(5) * t158;
	t180 = t146 * t189;
	t190 = qJD(5) * t146;
	t181 = t153 * t190;
	t143 = t147 * qJD(1);
	t199 = t143 * t160;
	t112 = (-(t180 + t199) * t152 - t156 * t181) * t138;
	t201 = t136 * t160;
	t108 = (-t112 * t146 - qJD(5)) * t201 + (t199 + (t112 + t190) * t158) * t135;
	t213 = t108 * t123 * t124;
	t214 = (-t198 * t213 + (t142 * t147 * t156 - t144 * t158 * t188) * t124) / t117 ^ 2;
	t169 = -qJD(6) * t146 + t142 * t158 + t147 * t188;
	t187 = qJD(6) * t158;
	t172 = t147 * t187 - t143;
	t113 = t169 * t157 + t172 * t159;
	t193 = t157 * t158;
	t133 = t146 * t159 + t147 * t193;
	t127 = t133 ^ 2;
	t122 = t127 * t129 + 0.1e1;
	t205 = t129 * t133;
	t114 = -t172 * t157 + t169 * t159;
	t209 = t114 * t128 * t129;
	t212 = (t113 * t205 - t127 * t209) / t122 ^ 2;
	t121 = t138 * t220;
	t203 = t135 * t158;
	t110 = t146 * t203 - t201 - (-t136 * t197 + t203) * t121;
	t211 = t110 * t124;
	t196 = t152 * t156;
	t182 = t146 * t196;
	t176 = t136 * t182;
	t202 = t135 * t160;
	t111 = (t202 + (t176 - t202) * t138) * t147;
	t210 = t111 * t124;
	t155 = t160 * t156;
	t195 = t152 * t160;
	t171 = t152 * t153 * t155 + t195;
	t208 = (-t171 * t145 * qJD(5) - t146 * t143 * t194) / t140 ^ 2;
	t207 = t124 * t147;
	t206 = t128 * t157;
	t204 = t133 * t159;
	t200 = t142 * t160;
	t191 = -t121 + t146;
	t186 = -0.2e1 * t213;
	t185 = 0.2e1 * t208;
	t184 = t160 * t219;
	t183 = t160 * t214;
	t179 = -t121 * t146 + 0.1e1;
	t177 = 0.2e1 * t133 * t209;
	t173 = t146 * t187 - t142;
	t170 = t129 * t204 - t206;
	t168 = qJD(6) * t147 - t143 * t158 + t146 * t188;
	t132 = t146 * t192 + t147 * t157;
	t131 = t146 * t193 - t147 * t159;
	t119 = 0.1e1 / t122;
	t115 = 0.1e1 / t117;
	t107 = t185 * t220 + (t178 * t143 + 0.2e1 * t171 * t190) * t138;
	t1 = [t152 * t184 * t208 + (t178 * t147 * qJD(5) - t142 * t195) * t138, 0, 0, 0, t107, 0; 0.2e1 * t146 * t123 * t183 + (t123 * t180 + (t123 * t143 + (t108 * t146 + t111 * t142) * t124) * t160) * t115 + (-0.2e1 * t183 * t210 + (-t189 * t210 + (t111 * t186 + ((-t135 * t189 - 0.2e1 * t176 * t208) * t147 + (t135 * t142 + (t112 * t136 + t135 * t185) * t147) * t160 + (((-t112 * t182 + t189) * t147 - t200) * t135 + (t142 * t182 + (-t155 * t181 - t143 * t196 + (-t112 - 0.2e1 * t190) * t160) * t147) * t136) * t138) * t124) * t160) * t115) * t147, 0, 0, 0, (-t123 * t158 - t160 * t211) * t214 * t219 + ((t142 * t123 + (-qJD(5) * t110 - t108) * t207) * t158 + (t142 * t211 + (qJD(5) * t123 + t110 * t186 + ((-t107 * t146 - t121 * t143) * t201 + (t191 * qJD(5) + t179 * t112) * t202) * t124) * t147 + ((t107 - t143) * t135 + (t179 * qJD(5) + t191 * t112) * t136) * t158 * t207) * t160) * t115, 0; 0.2e1 * (-t128 * t131 + t132 * t205) * t212 + (t132 * t177 + t173 * t128 * t159 + t168 * t206 + (t173 * t133 * t157 - t132 * t113 - t131 * t114 - t168 * t204) * t129) * t119, 0, 0, 0, t170 * t184 * t212 + (-t170 * t200 + (t170 * t189 + ((qJD(6) * t128 + t177) * t159 + (-t113 * t159 + (qJD(6) * t133 - t114) * t157) * t129) * t160) * t147) * t119, -0.2e1 * t212 + 0.2e1 * (t113 * t129 * t119 + (-t119 * t209 - t129 * t212) * t133) * t133;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end