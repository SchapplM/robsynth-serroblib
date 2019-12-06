% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PPRRR1
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
%   Wie in S5PPRRR1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5PPRRR1_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR1_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR1_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRRR1_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR1_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:13:11
	% EndTime: 2019-12-05 15:13:11
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:13:11
	% EndTime: 2019-12-05 15:13:11
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:13:11
	% EndTime: 2019-12-05 15:13:11
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:13:11
	% EndTime: 2019-12-05 15:13:11
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:13:11
	% EndTime: 2019-12-05 15:13:11
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:13:11
	% EndTime: 2019-12-05 15:13:11
	% DurationCPUTime: 0.35s
	% Computational Cost: add. (3038->45), mult. (2152->102), div. (440->12), fcn. (2434->9), ass. (0->60)
	t137 = pkin(9) + qJ(3) + qJ(4);
	t135 = sin(t137);
	t130 = t135 ^ 2;
	t136 = cos(t137);
	t166 = t130 / t136 ^ 2;
	t156 = 0.1e1 + t166;
	t140 = qJD(3) + qJD(4);
	t141 = sin(pkin(8));
	t176 = t140 * t141;
	t143 = sin(qJ(5));
	t144 = cos(qJ(5));
	t142 = cos(pkin(8));
	t164 = t136 * t142;
	t122 = t141 * t143 + t144 * t164;
	t138 = t141 ^ 2;
	t126 = t138 * t166 + 0.1e1;
	t124 = 0.1e1 / t126;
	t108 = t156 * t141 * t124;
	t163 = t141 * t135;
	t123 = atan2(-t163, -t136);
	t116 = cos(t123);
	t115 = sin(t123);
	t168 = t115 * t136;
	t154 = t116 * t135 - t141 * t168;
	t159 = t116 * t163;
	t155 = -t159 + t168;
	t97 = t108 * t155 + t154;
	t175 = 0.2e1 * t97;
	t107 = -t115 * t163 - t116 * t136;
	t104 = 0.1e1 / t107;
	t118 = 0.1e1 / t122;
	t105 = 0.1e1 / t107 ^ 2;
	t119 = 0.1e1 / t122 ^ 2;
	t139 = t142 ^ 2;
	t102 = t105 * t130 * t139 + 0.1e1;
	t165 = t136 * t140;
	t170 = t105 * t135;
	t103 = t108 * t140;
	t96 = t103 * t155 + t140 * t154;
	t173 = t104 * t105 * t96;
	t174 = 0.1e1 / t102 ^ 2 * (-t130 * t173 + t165 * t170) * t139;
	t121 = -t141 * t144 + t143 * t164;
	t117 = t121 ^ 2;
	t111 = t117 * t119 + 0.1e1;
	t158 = t135 * t140 * t142;
	t113 = -qJD(5) * t122 + t143 * t158;
	t167 = t119 * t121;
	t161 = qJD(5) * t121;
	t114 = -t144 * t158 - t161;
	t169 = t114 * t118 * t119;
	t172 = 0.1e1 / t111 ^ 2 * (-t113 * t167 - t117 * t169);
	t100 = 0.1e1 / t102;
	t171 = t100 * t105;
	t160 = -0.2e1 * t172;
	t153 = -t118 * t143 + t144 * t167;
	t109 = 0.1e1 / t111;
	t98 = 0.2e1 * (t124 - t156 / t126 ^ 2 * t138) * t156 / t136 * t135 * t176;
	t94 = (t153 * t135 * t160 + (t153 * t165 + ((-qJD(5) * t118 - 0.2e1 * t121 * t169) * t144 + (-t113 * t144 + (t114 - t161) * t143) * t119) * t135) * t109) * t142;
	t93 = ((-0.2e1 * t104 * t174 + (-t140 * t97 - t96) * t171) * t136 + (t105 * t174 * t175 + (t98 * t105 * t159 + t173 * t175 - t140 * t104 - (t176 - t103 + (t103 * t141 - t140) * t108) * t115 * t170) * t100 - (t115 * t98 + ((-t108 * t141 + 0.1e1) * t140 + (t108 - t141) * t103) * t116) * t136 * t171) * t135) * t142;
	t1 = [0, 0, t98, t98, 0; 0, 0, t93, t93, 0; 0, 0, t94, t94, t160 + 0.2e1 * (-t109 * t113 * t119 + (-t109 * t169 - t119 * t172) * t121) * t121;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end