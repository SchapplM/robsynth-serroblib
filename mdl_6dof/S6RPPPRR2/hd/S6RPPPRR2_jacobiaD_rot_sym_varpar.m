% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPPPRR2
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
%   Wie in S6RPPPRR2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2,theta4]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:27
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPPPRR2_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR2_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR2_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPPRR2_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR2_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:27:33
	% EndTime: 2019-10-09 23:27:33
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:27:33
	% EndTime: 2019-10-09 23:27:33
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:27:33
	% EndTime: 2019-10-09 23:27:33
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (31->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:27:33
	% EndTime: 2019-10-09 23:27:33
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:27:33
	% EndTime: 2019-10-09 23:27:33
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:27:33
	% EndTime: 2019-10-09 23:27:33
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:27:34
	% EndTime: 2019-10-09 23:27:35
	% DurationCPUTime: 1.02s
	% Computational Cost: add. (3055->92), mult. (2519->202), div. (480->12), fcn. (2968->9), ass. (0->96)
	t142 = pkin(10) + qJ(5);
	t138 = sin(t142);
	t132 = 0.1e1 / t138 ^ 2;
	t140 = cos(t142);
	t136 = t140 ^ 2;
	t190 = t132 * t136;
	t143 = qJ(1) + pkin(9);
	t141 = cos(t143);
	t207 = 0.2e1 * t141;
	t206 = t140 * t190;
	t185 = t141 * t140;
	t126 = atan2(-t185, t138);
	t124 = sin(t126);
	t125 = cos(t126);
	t111 = -t124 * t185 + t125 * t138;
	t108 = 0.1e1 / t111;
	t144 = sin(qJ(6));
	t184 = t141 * t144;
	t139 = sin(t143);
	t145 = cos(qJ(6));
	t186 = t139 * t145;
	t121 = t138 * t186 + t184;
	t117 = 0.1e1 / t121;
	t131 = 0.1e1 / t138;
	t109 = 0.1e1 / t111 ^ 2;
	t118 = 0.1e1 / t121 ^ 2;
	t137 = t141 ^ 2;
	t129 = t137 * t190 + 0.1e1;
	t127 = 0.1e1 / t129;
	t205 = t127 - 0.1e1;
	t134 = t139 ^ 2;
	t189 = t134 * t136;
	t104 = t109 * t189 + 0.1e1;
	t179 = qJD(1) * t141;
	t160 = t136 * t139 * t179;
	t177 = qJD(5) * t140;
	t180 = qJD(1) * t140;
	t170 = t139 * t180;
	t176 = qJD(5) * t141;
	t101 = ((t138 * t176 + t170) * t131 + t176 * t190) * t127;
	t192 = t125 * t140;
	t96 = (-t101 * t141 + qJD(5)) * t192 + (t170 + (-t101 + t176) * t138) * t124;
	t203 = t108 * t109 * t96;
	t204 = 0.1e1 / t104 ^ 2 * (-t189 * t203 + (-t134 * t138 * t177 + t160) * t109);
	t164 = qJD(6) * t138 + qJD(1);
	t154 = t144 * t177 + t164 * t145;
	t163 = qJD(1) * t138 + qJD(6);
	t157 = t141 * t163;
	t105 = t154 * t139 + t144 * t157;
	t183 = t141 * t145;
	t187 = t139 * t144;
	t120 = t138 * t187 - t183;
	t116 = t120 ^ 2;
	t115 = t116 * t118 + 0.1e1;
	t195 = t118 * t120;
	t153 = -t164 * t144 + t145 * t177;
	t106 = t153 * t139 + t145 * t157;
	t200 = t106 * t117 * t118;
	t202 = 0.1e1 / t115 ^ 2 * (t105 * t195 - t116 * t200);
	t201 = t101 * t140;
	t155 = qJD(5) * (-t140 - t206) * t131;
	t199 = (-t132 * t160 + t137 * t155) / t129 ^ 2;
	t198 = t109 * t139;
	t197 = t109 * t140;
	t196 = t117 * t144;
	t194 = t120 * t145;
	t193 = t124 * t138;
	t191 = t131 * t136;
	t188 = t139 * t140;
	t167 = 0.1e1 + t190;
	t112 = t167 * t141 * t127;
	t182 = -t112 + t141;
	t181 = qJD(1) * t139;
	t178 = qJD(5) * t138;
	t175 = -0.2e1 * t203;
	t174 = 0.2e1 * t202;
	t173 = t140 * t204;
	t172 = t140 * t199;
	t171 = t127 * t191;
	t169 = t140 * t179;
	t168 = t112 * t141 - 0.1e1;
	t166 = 0.2e1 * t120 * t200;
	t165 = t199 * t207;
	t162 = t141 * t171;
	t161 = t205 * t140 * t124;
	t159 = t167 * t139;
	t158 = t139 * t163;
	t156 = t118 * t194 - t196;
	t123 = t138 * t183 - t187;
	t122 = t138 * t184 + t186;
	t113 = 0.1e1 / t115;
	t102 = 0.1e1 / t104;
	t100 = (-t125 * t162 - t161) * t139;
	t98 = t141 * t193 + t192 + (-t125 * t185 - t193) * t112;
	t97 = -t167 * t165 + (-qJD(1) * t159 + t155 * t207) * t127;
	t1 = [-0.2e1 * t131 * t139 * t172 + (-qJD(5) * t159 + t131 * t169) * t127, 0, 0, 0, t97, 0; (0.2e1 * t108 * t173 + (t108 * t178 + (qJD(1) * t100 + t96) * t197) * t102) * t141 + (-0.2e1 * t109 * t173 * t100 + (((t101 * t162 + t205 * t178 + 0.2e1 * t172) * t124 + (t165 * t191 + t201 + (-t201 + (0.2e1 * t140 + t206) * t176) * t127) * t125) * t109 * t188 + (-t109 * t178 + t140 * t175) * t100 + (t108 + ((t134 - t137) * t125 * t171 - t141 * t161) * t109) * t180) * t102) * t139, 0, 0, 0, 0.2e1 * (-t108 * t138 - t98 * t197) * t139 * t204 + ((t108 * t179 + (-qJD(5) * t98 - t96) * t198) * t138 + ((qJD(5) * t108 + t98 * t175) * t139 + (t98 * t179 + ((t112 * t181 - t141 * t97) * t125 + (t182 * qJD(5) + t168 * t101) * t124) * t188) * t109 + ((-t97 - t181) * t124 + (t168 * qJD(5) + t182 * t101) * t125) * t138 * t198) * t140) * t102, 0; (-t117 * t122 + t123 * t195) * t174 + (t123 * t166 - t158 * t196 + t154 * t117 * t141 + (-t153 * t120 * t141 - t123 * t105 - t122 * t106 + t158 * t194) * t118) * t113, 0, 0, 0, t156 * t174 * t188 + (-t156 * t169 + (t156 * t178 + ((qJD(6) * t117 + t166) * t145 + (-t105 * t145 + (qJD(6) * t120 - t106) * t144) * t118) * t140) * t139) * t113, -0.2e1 * t202 + 0.2e1 * (t105 * t113 * t118 + (-t113 * t200 - t118 * t202) * t120) * t120;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end