% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PRRRR12
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
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S5PRRRR12_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d2,d3,d4,d5,theta1]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-28 18:09
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5PRRRR12_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR12_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRR12_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR12_jacobia_rot_sym_varpar: pkin has to be [11x1] (double)');
Ja_rot=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-28 18:09:12
	% EndTime: 2024-09-28 18:09:12
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-28 18:09:12
	% EndTime: 2024-09-28 18:09:12
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-28 18:09:12
	% EndTime: 2024-09-28 18:09:12
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (10->0), mult. (27->0), div. (5->0), fcn. (35->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 1, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-28 18:09:12
	% EndTime: 2024-09-28 18:09:12
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (132->0), mult. (84->0), div. (12->0), fcn. (90->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 1, 1, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-28 18:09:12
	% EndTime: 2024-09-28 18:09:12
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (279->0), mult. (126->0), div. (18->0), fcn. (135->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 1, 1, 1, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-28 18:09:12
	% EndTime: 2024-09-28 18:09:12
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (1233->28), mult. (1557->63), div. (96->9), fcn. (2142->13), ass. (0->45)
	t132 = cos(pkin(6));
	t134 = sin(qJ(5));
	t128 = sin(pkin(11));
	t133 = cos(pkin(5));
	t144 = t128 * t133;
	t138 = t134 * t144;
	t131 = cos(pkin(11));
	t135 = cos(qJ(5));
	t141 = t131 * t135;
	t118 = t132 * t138 - t141;
	t137 = t135 * t144;
	t142 = t131 * t134;
	t120 = t132 * t142 + t137;
	t127 = qJ(2) + qJ(3) + qJ(4);
	t125 = sin(t127);
	t126 = cos(t127);
	t129 = sin(pkin(6));
	t130 = sin(pkin(5));
	t143 = t129 * t130;
	t139 = t128 * t143;
	t106 = -t118 * t126 - t120 * t125 + t134 * t139;
	t105 = 0.1e1 / t106 ^ 2;
	t119 = -t132 * t137 - t142;
	t121 = -t132 * t141 + t138;
	t107 = t119 * t126 + t121 * t125 + t135 * t139;
	t146 = t107 ^ 2 * t105;
	t145 = t125 * t131;
	t140 = t125 * t143;
	t136 = t126 * t129 * t133 + t130 * t132;
	t117 = -t126 * t143 + t133 * t132;
	t116 = 0.1e1 / t117 ^ 2;
	t115 = (t126 * t128 + t133 * t145) * t129;
	t114 = -t129 * t128 * t125 + t136 * t131;
	t113 = t136 * t128 + t129 * t145;
	t112 = atan2(t114, t117);
	t110 = cos(t112);
	t109 = sin(t112);
	t104 = 0.1e1 / t106;
	t103 = 0.1e1 / (0.1e1 + t146);
	t102 = t109 * t114 + t110 * t117;
	t101 = 0.1e1 / t102 ^ 2;
	t99 = (-t115 / t117 - t114 * t116 * t140) / (t114 ^ 2 * t116 + 0.1e1);
	t98 = ((t119 * t125 - t121 * t126) * t104 + (t118 * t125 - t120 * t126) * t107 * t105) * t103;
	t97 = ((-t125 * t144 + t126 * t131) * t129 / t102 - ((t114 * t99 + t140) * t110 + (-t117 * t99 - t115) * t109) * t113 * t101) / (t113 ^ 2 * t101 + 0.1e1);
	t1 = [0, t99, t99, t99, 0; 0, t97, t97, t97, 0; 0, t98, t98, t98, (t106 * t104 + t146) * t103;];
	Ja_rot = t1;
end