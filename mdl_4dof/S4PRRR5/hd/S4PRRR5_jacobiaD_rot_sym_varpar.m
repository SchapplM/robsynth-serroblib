% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S4PRRR5
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S4PRRR5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% JaD_rot [3x4]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 12:30
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S4PRRR5_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR5_jacobiaD_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR5_jacobiaD_rot_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4PRRR5_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR5_jacobiaD_rot_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:29:51
	% EndTime: 2019-12-29 12:29:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:29:51
	% EndTime: 2019-12-29 12:29:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:29:51
	% EndTime: 2019-12-29 12:29:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:29:56
	% EndTime: 2019-12-29 12:29:56
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:29:56
	% EndTime: 2019-12-29 12:29:57
	% DurationCPUTime: 0.80s
	% Computational Cost: add. (1992->45), mult. (2152->102), div. (440->12), fcn. (2434->9), ass. (0->60)
	t139 = qJ(2) + qJ(3);
	t134 = sin(t139);
	t130 = t134 ^ 2;
	t135 = cos(t139);
	t165 = t130 / t135 ^ 2;
	t155 = 0.1e1 + t165;
	t138 = qJD(2) + qJD(3);
	t140 = sin(pkin(7));
	t175 = t138 * t140;
	t142 = sin(qJ(4));
	t143 = cos(qJ(4));
	t141 = cos(pkin(7));
	t163 = t135 * t141;
	t119 = t140 * t142 + t143 * t163;
	t136 = t140 ^ 2;
	t125 = t136 * t165 + 0.1e1;
	t123 = 0.1e1 / t125;
	t107 = t155 * t140 * t123;
	t162 = t140 * t134;
	t122 = atan2(-t162, -t135);
	t121 = cos(t122);
	t120 = sin(t122);
	t166 = t120 * t135;
	t152 = t121 * t134 - t140 * t166;
	t158 = t121 * t162;
	t153 = -t158 + t166;
	t96 = t153 * t107 + t152;
	t174 = 0.2e1 * t96;
	t106 = -t120 * t162 - t121 * t135;
	t103 = 0.1e1 / t106;
	t115 = 0.1e1 / t119;
	t104 = 0.1e1 / t106 ^ 2;
	t116 = 0.1e1 / t119 ^ 2;
	t137 = t141 ^ 2;
	t101 = t137 * t130 * t104 + 0.1e1;
	t164 = t135 * t138;
	t169 = t104 * t134;
	t102 = t107 * t138;
	t95 = t153 * t102 + t152 * t138;
	t171 = t103 * t104 * t95;
	t173 = 0.1e1 / t101 ^ 2 * (-t130 * t171 + t164 * t169) * t137;
	t99 = 0.1e1 / t101;
	t172 = t104 * t99;
	t118 = -t140 * t143 + t142 * t163;
	t114 = t118 ^ 2;
	t110 = t114 * t116 + 0.1e1;
	t157 = t134 * t138 * t141;
	t111 = -t119 * qJD(4) + t142 * t157;
	t167 = t116 * t118;
	t160 = qJD(4) * t118;
	t112 = -t143 * t157 - t160;
	t168 = t112 * t115 * t116;
	t170 = 0.1e1 / t110 ^ 2 * (-t111 * t167 - t114 * t168);
	t159 = -0.2e1 * t170;
	t154 = -t115 * t142 + t143 * t167;
	t108 = 0.1e1 / t110;
	t97 = 0.2e1 * (t123 - t155 / t125 ^ 2 * t136) * t155 / t135 * t134 * t175;
	t93 = (t154 * t134 * t159 + (t154 * t164 + ((-qJD(4) * t115 - 0.2e1 * t118 * t168) * t143 + (-t111 * t143 + (t112 - t160) * t142) * t116) * t134) * t108) * t141;
	t92 = ((-0.2e1 * t103 * t173 + (-t138 * t96 - t95) * t172) * t135 + (t104 * t173 * t174 + (t97 * t104 * t158 + t171 * t174 - t138 * t103 - (t175 - t102 + (t102 * t140 - t138) * t107) * t120 * t169) * t99 - (t120 * t97 + ((-t107 * t140 + 0.1e1) * t138 + (t107 - t140) * t102) * t121) * t135 * t172) * t134) * t141;
	t1 = [0, t97, t97, 0; 0, t92, t92, 0; 0, t93, t93, t159 + 0.2e1 * (-t108 * t111 * t116 + (-t108 * t168 - t116 * t170) * t118) * t118;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,4);
end