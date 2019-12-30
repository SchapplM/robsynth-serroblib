% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S5RRRRR10
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% Jg_rot [3x5]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 21:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S5RRRRR10_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR10_jacobig_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR10_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR10_jacobig_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 21:01:26
	% EndTime: 2019-12-29 21:01:26
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 21:01:26
	% EndTime: 2019-12-29 21:01:26
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 21:01:26
	% EndTime: 2019-12-29 21:01:26
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t43 = sin(pkin(5));
	t1 = [0, sin(qJ(1)) * t43, 0, 0, 0; 0, -cos(qJ(1)) * t43, 0, 0, 0; 1, cos(pkin(5)), 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 21:01:26
	% EndTime: 2019-12-29 21:01:26
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (4->4), mult. (9->8), div. (0->0), fcn. (17->6), ass. (0->8)
	t70 = cos(pkin(5));
	t73 = cos(qJ(2));
	t75 = t70 * t73;
	t74 = cos(qJ(1));
	t72 = sin(qJ(1));
	t71 = sin(qJ(2));
	t69 = sin(pkin(5));
	t1 = [0, t72 * t69, t74 * t71 + t72 * t75, 0, 0; 0, -t74 * t69, t72 * t71 - t74 * t75, 0, 0; 1, t70, -t69 * t73, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 21:01:26
	% EndTime: 2019-12-29 21:01:26
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (7->5), mult. (16->8), div. (0->0), fcn. (29->6), ass. (0->11)
	t86 = sin(pkin(5));
	t90 = cos(qJ(2));
	t93 = t86 * t90;
	t87 = cos(pkin(5));
	t92 = t87 * t90;
	t91 = cos(qJ(1));
	t89 = sin(qJ(1));
	t88 = sin(qJ(2));
	t85 = t91 * t88 + t89 * t92;
	t84 = t89 * t88 - t91 * t92;
	t1 = [0, t89 * t86, t85, t85, 0; 0, -t91 * t86, t84, t84, 0; 1, t87, -t93, -t93, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 21:01:27
	% EndTime: 2019-12-29 21:01:27
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (18->11), mult. (31->18), div. (0->0), fcn. (52->8), ass. (0->19)
	t130 = sin(pkin(5));
	t134 = cos(qJ(2));
	t142 = t130 * t134;
	t133 = sin(qJ(1));
	t141 = t133 * t130;
	t132 = sin(qJ(2));
	t140 = t133 * t132;
	t139 = t133 * t134;
	t135 = cos(qJ(1));
	t138 = t135 * t130;
	t137 = t135 * t132;
	t136 = t135 * t134;
	t131 = cos(pkin(5));
	t129 = qJ(3) + qJ(4);
	t128 = cos(t129);
	t127 = sin(t129);
	t126 = t131 * t139 + t137;
	t125 = -t131 * t136 + t140;
	t1 = [0, t141, t126, t126, (-t131 * t140 + t136) * t127 - t128 * t141; 0, -t138, t125, t125, (t131 * t137 + t139) * t127 + t128 * t138; 1, t131, -t142, -t142, t130 * t132 * t127 - t131 * t128;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,5);
end