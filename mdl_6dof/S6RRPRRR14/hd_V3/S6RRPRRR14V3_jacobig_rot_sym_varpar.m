% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRRR14V3
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:13
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRPRRR14V3_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14V3_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR14V3_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S6RRPRRR14V3_jacobig_rot_sym_varpar: pkin has to be [1x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:13:23
	% EndTime: 2019-10-10 11:13:23
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:13:23
	% EndTime: 2019-10-10 11:13:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:13:23
	% EndTime: 2019-10-10 11:13:23
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, sin(qJ(1)), 0, 0, 0, 0; 0, -cos(qJ(1)), 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:13:23
	% EndTime: 2019-10-10 11:13:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, sin(qJ(1)), 0, 0, 0, 0; 0, -cos(qJ(1)), 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:13:23
	% EndTime: 2019-10-10 11:13:23
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (2->2), div. (0->0), fcn. (7->4), ass. (0->4)
	t53 = cos(qJ(1));
	t52 = sin(qJ(1));
	t51 = sin(qJ(2));
	t1 = [0, t52, 0, t53 * t51, 0, 0; 0, -t53, 0, t52 * t51, 0, 0; 1, 0, 0, -cos(qJ(2)), 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:13:23
	% EndTime: 2019-10-10 11:13:23
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (4->4), mult. (9->8), div. (0->0), fcn. (19->6), ass. (0->8)
	t79 = sin(qJ(4));
	t83 = cos(qJ(2));
	t85 = t79 * t83;
	t84 = cos(qJ(1));
	t82 = cos(qJ(4));
	t81 = sin(qJ(1));
	t80 = sin(qJ(2));
	t1 = [0, t81, 0, t84 * t80, -t81 * t82 + t84 * t85, 0; 0, -t84, 0, t81 * t80, t81 * t85 + t84 * t82, 0; 1, 0, 0, -t83, t80 * t79, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:13:23
	% EndTime: 2019-10-10 11:13:23
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (9->9), mult. (24->19), div. (0->0), fcn. (42->8), ass. (0->14)
	t120 = sin(qJ(2));
	t121 = sin(qJ(1));
	t130 = t121 * t120;
	t124 = cos(qJ(2));
	t129 = t121 * t124;
	t119 = sin(qJ(4));
	t125 = cos(qJ(1));
	t128 = t125 * t119;
	t127 = t125 * t120;
	t123 = cos(qJ(4));
	t126 = t125 * t123;
	t122 = cos(qJ(5));
	t118 = sin(qJ(5));
	t1 = [0, t121, 0, t127, -t121 * t123 + t124 * t128, (t121 * t119 + t124 * t126) * t118 - t122 * t127; 0, -t125, 0, t130, t119 * t129 + t126, (t123 * t129 - t128) * t118 - t122 * t130; 1, 0, 0, -t124, t120 * t119, t120 * t123 * t118 + t124 * t122;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end