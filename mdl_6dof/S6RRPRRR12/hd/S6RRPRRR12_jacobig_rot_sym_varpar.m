% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRRR12
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:07
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRPRRR12_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR12_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR12_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR12_jacobig_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:07:01
	% EndTime: 2019-10-10 11:07:01
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:07:01
	% EndTime: 2019-10-10 11:07:01
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:07:02
	% EndTime: 2019-10-10 11:07:02
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t43 = sin(pkin(6));
	t1 = [0, sin(qJ(1)) * t43, 0, 0, 0, 0; 0, -cos(qJ(1)) * t43, 0, 0, 0, 0; 1, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:07:02
	% EndTime: 2019-10-10 11:07:02
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t55 = sin(pkin(6));
	t1 = [0, sin(qJ(1)) * t55, 0, 0, 0, 0; 0, -cos(qJ(1)) * t55, 0, 0, 0, 0; 1, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:07:02
	% EndTime: 2019-10-10 11:07:02
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (3->3), mult. (9->8), div. (0->0), fcn. (17->6), ass. (0->8)
	t71 = cos(pkin(6));
	t72 = sin(qJ(2));
	t76 = t71 * t72;
	t75 = cos(qJ(1));
	t74 = cos(qJ(2));
	t73 = sin(qJ(1));
	t70 = sin(pkin(6));
	t1 = [0, t73 * t70, 0, -t73 * t76 + t75 * t74, 0, 0; 0, -t75 * t70, 0, t73 * t74 + t75 * t76, 0, 0; 1, t71, 0, t70 * t72, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:07:02
	% EndTime: 2019-10-10 11:07:02
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (5->3), mult. (16->8), div. (0->0), fcn. (29->6), ass. (0->11)
	t91 = cos(pkin(6));
	t92 = sin(qJ(2));
	t96 = t91 * t92;
	t95 = cos(qJ(1));
	t94 = cos(qJ(2));
	t93 = sin(qJ(1));
	t90 = sin(pkin(6));
	t89 = t90 * t92;
	t88 = -t93 * t96 + t95 * t94;
	t87 = t93 * t94 + t95 * t96;
	t1 = [0, t93 * t90, 0, t88, t88, 0; 0, -t95 * t90, 0, t87, t87, 0; 1, t91, 0, t89, t89, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:07:02
	% EndTime: 2019-10-10 11:07:02
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (16->9), mult. (31->18), div. (0->0), fcn. (52->8), ass. (0->19)
	t129 = sin(pkin(6));
	t132 = sin(qJ(1));
	t140 = t132 * t129;
	t131 = sin(qJ(2));
	t139 = t132 * t131;
	t133 = cos(qJ(2));
	t138 = t132 * t133;
	t134 = cos(qJ(1));
	t137 = t134 * t129;
	t136 = t134 * t131;
	t135 = t134 * t133;
	t130 = cos(pkin(6));
	t128 = qJ(4) + qJ(5);
	t127 = cos(t128);
	t126 = sin(t128);
	t125 = t129 * t131;
	t124 = -t130 * t139 + t135;
	t123 = t130 * t136 + t138;
	t1 = [0, t140, 0, t124, t124, t126 * t140 - (t130 * t138 + t136) * t127; 0, -t137, 0, t123, t123, -t126 * t137 - (-t130 * t135 + t139) * t127; 1, t130, 0, t125, t125, t129 * t133 * t127 + t130 * t126;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end