% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RPRPR7
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
%   Wie in S5RPRPR7_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 16:49
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RPRPR7_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR7_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR7_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPR7_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR7_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:49:05
	% EndTime: 2019-12-29 16:49:05
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:49:05
	% EndTime: 2019-12-29 16:49:05
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:49:05
	% EndTime: 2019-12-29 16:49:05
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (31->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:49:05
	% EndTime: 2019-12-29 16:49:05
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:49:10
	% EndTime: 2019-12-29 16:49:10
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:49:06
	% EndTime: 2019-12-29 16:49:07
	% DurationCPUTime: 1.66s
	% Computational Cost: add. (3244->96), mult. (2519->203), div. (480->12), fcn. (2968->9), ass. (0->95)
	t146 = qJ(3) + pkin(9);
	t142 = sin(t146);
	t136 = t142 ^ 2;
	t144 = cos(t146);
	t139 = 0.1e1 / t144 ^ 2;
	t195 = t136 * t139;
	t147 = qJ(1) + pkin(8);
	t143 = sin(t147);
	t213 = 0.2e1 * t143;
	t212 = t142 * t195;
	t145 = cos(t147);
	t149 = cos(qJ(5));
	t187 = t145 * t149;
	t148 = sin(qJ(5));
	t190 = t143 * t148;
	t125 = t144 * t187 + t190;
	t165 = qJD(5) * t144 - qJD(1);
	t184 = qJD(3) * t142;
	t211 = t165 * t148 + t149 * t184;
	t191 = t143 * t142;
	t128 = atan2(-t191, -t144);
	t127 = cos(t128);
	t126 = sin(t128);
	t175 = t126 * t191;
	t112 = -t127 * t144 - t175;
	t109 = 0.1e1 / t112;
	t119 = 0.1e1 / t125;
	t138 = 0.1e1 / t144;
	t110 = 0.1e1 / t112 ^ 2;
	t120 = 0.1e1 / t125 ^ 2;
	t210 = -0.2e1 * t142;
	t137 = t143 ^ 2;
	t131 = t137 * t195 + 0.1e1;
	t129 = 0.1e1 / t131;
	t209 = t129 - 0.1e1;
	t185 = qJD(1) * t145;
	t172 = t142 * t185;
	t182 = qJD(3) * t144;
	t183 = qJD(3) * t143;
	t103 = (-(-t143 * t182 - t172) * t138 + t183 * t195) * t129;
	t197 = t127 * t142;
	t98 = (-t103 * t143 + qJD(3)) * t197 + (-t172 + (t103 - t183) * t144) * t126;
	t208 = t109 * t110 * t98;
	t158 = t144 * t190 + t187;
	t171 = t148 * t184;
	t107 = t158 * qJD(1) - qJD(5) * t125 + t145 * t171;
	t188 = t145 * t148;
	t189 = t143 * t149;
	t124 = t144 * t188 - t189;
	t118 = t124 ^ 2;
	t117 = t118 * t120 + 0.1e1;
	t199 = t120 * t124;
	t164 = -qJD(1) * t144 + qJD(5);
	t160 = t164 * t149;
	t108 = t143 * t160 - t145 * t211;
	t204 = t108 * t119 * t120;
	t207 = (-t107 * t199 - t118 * t204) / t117 ^ 2;
	t206 = t103 * t126;
	t205 = t103 * t142;
	t203 = t110 * t142;
	t202 = t110 * t145;
	t193 = t138 * t142;
	t157 = qJD(3) * (t138 * t212 + t193);
	t162 = t136 * t143 * t185;
	t201 = (t137 * t157 + t139 * t162) / t131 ^ 2;
	t169 = 0.1e1 + t195;
	t114 = t169 * t143 * t129;
	t200 = t114 * t143;
	t198 = t126 * t144;
	t196 = t136 * t138;
	t141 = t145 ^ 2;
	t194 = t136 * t141;
	t192 = t142 * t145;
	t186 = qJD(1) * t143;
	t181 = qJD(3) * t145;
	t106 = t110 * t194 + 0.1e1;
	t180 = 0.2e1 / t106 ^ 2 * (-t194 * t208 + (t141 * t142 * t182 - t162) * t110);
	t179 = 0.2e1 * t208;
	t178 = -0.2e1 * t207;
	t177 = t124 * t204;
	t176 = t110 * t192;
	t174 = t129 * t196;
	t168 = t142 * t180;
	t167 = t201 * t210;
	t166 = t201 * t213;
	t163 = t143 * t174;
	t161 = t169 * t145;
	t159 = -t119 * t148 + t149 * t199;
	t123 = -t144 * t189 + t188;
	t115 = 0.1e1 / t117;
	t104 = 0.1e1 / t106;
	t102 = (t209 * t142 * t126 - t127 * t163) * t145;
	t100 = -t143 * t198 + t197 + (-t127 * t191 + t198) * t114;
	t99 = -t169 * t166 + (qJD(1) * t161 + t157 * t213) * t129;
	t1 = [t138 * t145 * t167 + (qJD(3) * t161 - t186 * t193) * t129, 0, t99, 0, 0; (t109 * t168 + (-t109 * t182 + (qJD(1) * t102 + t98) * t203) * t104) * t143 + (t110 * t168 * t102 + (-((t103 * t163 + t209 * t182 + t167) * t126 + (t166 * t196 - t205 + (t205 + (t210 - t212) * t183) * t129) * t127) * t176 + (-t110 * t182 + t142 * t179) * t102 + (-t109 + ((-t137 + t141) * t127 * t174 + t209 * t175) * t110) * t142 * qJD(1)) * t104) * t145, 0, (t100 * t203 - t109 * t144) * t145 * t180 + ((-t109 * t186 + (-qJD(3) * t100 - t98) * t202) * t144 + (-t109 * t181 - (-t127 * t143 * t99 + t126 * t183 + t200 * t206 - t206 + (-qJD(3) * t126 - t127 * t185) * t114) * t176 + (t110 * t186 + t145 * t179) * t100 - ((t99 - t185) * t126 + ((0.1e1 - t200) * qJD(3) + (t114 - t143) * t103) * t127) * t144 * t202) * t142) * t104, 0, 0; 0.2e1 * (t119 * t158 + t123 * t199) * t207 + (0.2e1 * t123 * t177 + (t123 * t107 + t158 * t108 + (-t143 * t211 - t145 * t160) * t124) * t120 + (t164 * t188 + (-t165 * t149 + t171) * t143) * t119) * t115, 0, t159 * t178 * t192 + (t159 * t144 * t181 + (-t159 * t186 + ((-qJD(5) * t119 - 0.2e1 * t177) * t149 + (-t107 * t149 + (-qJD(5) * t124 + t108) * t148) * t120) * t145) * t142) * t115, 0, t178 + 0.2e1 * (-t107 * t115 * t120 + (-t115 * t204 - t120 * t207) * t124) * t124;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end