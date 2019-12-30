% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RPPRR12
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
%   Wie in S5RPPRR12_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 16:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RPPRR12_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR12_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR12_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPRR12_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR12_jacobiaD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:29:47
	% EndTime: 2019-12-29 16:29:47
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:29:42
	% EndTime: 2019-12-29 16:29:42
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:29:47
	% EndTime: 2019-12-29 16:29:47
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:29:42
	% EndTime: 2019-12-29 16:29:42
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:29:42
	% EndTime: 2019-12-29 16:29:42
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:29:47
	% EndTime: 2019-12-29 16:29:49
	% DurationCPUTime: 1.49s
	% Computational Cost: add. (2081->92), mult. (2519->206), div. (480->12), fcn. (2968->9), ass. (0->95)
	t136 = cos(qJ(1));
	t200 = 0.2e1 * t136;
	t130 = pkin(8) + qJ(4);
	t128 = sin(t130);
	t129 = cos(t130);
	t177 = t136 * t129;
	t118 = atan2(-t177, t128);
	t116 = sin(t118);
	t117 = cos(t118);
	t102 = -t116 * t177 + t117 * t128;
	t99 = 0.1e1 / t102;
	t134 = sin(qJ(1));
	t135 = cos(qJ(5));
	t178 = t134 * t135;
	t133 = sin(qJ(5));
	t180 = t133 * t136;
	t113 = t128 * t178 + t180;
	t109 = 0.1e1 / t113;
	t123 = 0.1e1 / t128;
	t100 = 0.1e1 / t102 ^ 2;
	t110 = 0.1e1 / t113 ^ 2;
	t124 = 0.1e1 / t128 ^ 2;
	t132 = t136 ^ 2;
	t127 = t129 ^ 2;
	t183 = t124 * t127;
	t121 = t132 * t183 + 0.1e1;
	t119 = 0.1e1 / t121;
	t199 = t119 - 0.1e1;
	t131 = t134 ^ 2;
	t174 = qJD(1) * t136;
	t152 = t127 * t134 * t174;
	t172 = qJD(4) * t129;
	t182 = t127 * t131;
	t171 = qJD(4) * t136;
	t162 = t124 * t171;
	t175 = qJD(1) * t134;
	t163 = t129 * t175;
	t93 = ((t128 * t171 + t163) * t123 + t127 * t162) * t119;
	t157 = -t93 + t171;
	t158 = -t136 * t93 + qJD(4);
	t186 = t117 * t129;
	t88 = t158 * t186 + (t157 * t128 + t163) * t116;
	t196 = t99 * t100 * t88;
	t96 = t100 * t182 + 0.1e1;
	t198 = (-t182 * t196 + (-t128 * t131 * t172 + t152) * t100) / t96 ^ 2;
	t94 = 0.1e1 / t96;
	t197 = t100 * t94;
	t176 = t136 * t135;
	t179 = t134 * t133;
	t112 = t128 * t179 - t176;
	t108 = t112 ^ 2;
	t107 = t108 * t110 + 0.1e1;
	t189 = t110 * t112;
	t155 = qJD(1) * t128 + qJD(5);
	t148 = t155 * t136;
	t156 = qJD(5) * t128 + qJD(1);
	t150 = t156 * t133;
	t98 = t135 * t148 + (t135 * t172 - t150) * t134;
	t194 = t109 * t110 * t98;
	t149 = t156 * t135;
	t97 = t134 * t149 + (t134 * t172 + t148) * t133;
	t195 = 0.1e1 / t107 ^ 2 * (-t108 * t194 + t97 * t189);
	t126 = t129 * t127;
	t184 = t123 * t129;
	t146 = qJD(4) * (-t123 * t124 * t126 - t184);
	t191 = (-t124 * t152 + t132 * t146) / t121 ^ 2;
	t190 = t109 * t133;
	t188 = t112 * t135;
	t187 = t116 * t128;
	t185 = t123 * t127;
	t173 = qJD(4) * t128;
	t170 = -0.2e1 * t198;
	t169 = -0.2e1 * t196;
	t168 = 0.2e1 * t195;
	t167 = t99 * t198;
	t166 = t94 * t173;
	t165 = t129 * t191;
	t164 = t119 * t185;
	t161 = 0.1e1 + t183;
	t160 = 0.2e1 * t112 * t194;
	t159 = t191 * t200;
	t154 = t136 * t164;
	t153 = t199 * t129 * t116;
	t151 = t161 * t134;
	t147 = t110 * t188 - t190;
	t145 = t147 * t134;
	t144 = t129 * t171 - t155 * t134;
	t115 = t128 * t176 - t179;
	t114 = t128 * t180 + t178;
	t105 = 0.1e1 / t107;
	t104 = t161 * t136 * t119;
	t92 = (-t117 * t154 - t153) * t134;
	t90 = t136 * t187 + t186 + (-t117 * t177 - t187) * t104;
	t89 = -t161 * t159 + (-qJD(1) * t151 + t146 * t200) * t119;
	t1 = [-0.2e1 * t134 * t123 * t165 + (-qJD(4) * t151 + t174 * t184) * t119, 0, 0, t89, 0; (t99 * t166 + (0.2e1 * t167 + (qJD(1) * t92 + t88) * t197) * t129) * t136 + ((-t92 * t166 + (t92 * t170 + ((t93 * t154 + t199 * t173 + 0.2e1 * t165) * t116 + (t159 * t185 + t93 * t129 + (t126 * t162 + (-t93 + 0.2e1 * t171) * t129) * t119) * t117) * t94 * t134) * t129) * t100 + (t92 * t169 + (t99 + ((t131 - t132) * t117 * t164 - t136 * t153) * t100) * qJD(1)) * t129 * t94) * t134, 0, 0, (t99 * t94 * t174 + (-0.2e1 * t167 + (-qJD(4) * t90 - t88) * t197) * t134) * t128 + (((qJD(4) * t99 + t90 * t169) * t134 + (t90 * t174 + ((t104 * t175 - t136 * t89) * t117 + ((t104 * t136 - 0.1e1) * t93 + (-t104 + t136) * qJD(4)) * t116) * t129 * t134) * t100) * t94 + (t90 * t170 + ((-t89 - t175) * t116 + (t157 * t104 - t158) * t117) * t94 * t128) * t100 * t134) * t129, 0; (-t109 * t114 + t115 * t189) * t168 + (t115 * t160 + t136 * t109 * t149 + t144 * t190 + (t136 * t112 * t150 - t114 * t98 - t115 * t97 - t144 * t188) * t110) * t105, 0, 0, t129 * t145 * t168 + (t145 * t173 + (-t147 * t174 + ((qJD(5) * t109 + t160) * t135 + (-t135 * t97 + (qJD(5) * t112 - t98) * t133) * t110) * t134) * t129) * t105, -0.2e1 * t195 + 0.2e1 * (t105 * t110 * t97 + (-t105 * t194 - t110 * t195) * t112) * t112;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end