% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRRRPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:46
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRRPR10_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR10_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPR10_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR10_jacobig_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:46:05
	% EndTime: 2019-10-10 12:46:05
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:46:05
	% EndTime: 2019-10-10 12:46:05
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:46:05
	% EndTime: 2019-10-10 12:46:05
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t43 = sin(pkin(6));
	t1 = [0, sin(qJ(1)) * t43, 0, 0, 0, 0; 0, -cos(qJ(1)) * t43, 0, 0, 0, 0; 1, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:46:05
	% EndTime: 2019-10-10 12:46:05
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (4->4), mult. (9->8), div. (0->0), fcn. (17->6), ass. (0->8)
	t70 = cos(pkin(6));
	t73 = cos(qJ(2));
	t75 = t70 * t73;
	t74 = cos(qJ(1));
	t72 = sin(qJ(1));
	t71 = sin(qJ(2));
	t69 = sin(pkin(6));
	t1 = [0, t72 * t69, t74 * t71 + t72 * t75, 0, 0, 0; 0, -t74 * t69, t72 * t71 - t74 * t75, 0, 0, 0; 1, t70, -t69 * t73, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:46:05
	% EndTime: 2019-10-10 12:46:05
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (7->5), mult. (16->8), div. (0->0), fcn. (29->6), ass. (0->11)
	t86 = sin(pkin(6));
	t90 = cos(qJ(2));
	t93 = t86 * t90;
	t87 = cos(pkin(6));
	t92 = t87 * t90;
	t91 = cos(qJ(1));
	t89 = sin(qJ(1));
	t88 = sin(qJ(2));
	t85 = t91 * t88 + t89 * t92;
	t84 = t89 * t88 - t91 * t92;
	t1 = [0, t89 * t86, t85, t85, 0, 0; 0, -t91 * t86, t84, t84, 0, 0; 1, t87, -t93, -t93, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:46:05
	% EndTime: 2019-10-10 12:46:05
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (7->5), mult. (16->8), div. (0->0), fcn. (29->6), ass. (0->11)
	t109 = sin(pkin(6));
	t113 = cos(qJ(2));
	t116 = t109 * t113;
	t110 = cos(pkin(6));
	t115 = t110 * t113;
	t114 = cos(qJ(1));
	t112 = sin(qJ(1));
	t111 = sin(qJ(2));
	t108 = t114 * t111 + t112 * t115;
	t107 = t112 * t111 - t114 * t115;
	t1 = [0, t112 * t109, t108, t108, 0, 0; 0, -t114 * t109, t107, t107, 0, 0; 1, t110, -t116, -t116, 0, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:46:06
	% EndTime: 2019-10-10 12:46:06
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (18->11), mult. (31->18), div. (0->0), fcn. (52->8), ass. (0->19)
	t129 = sin(pkin(6));
	t133 = cos(qJ(2));
	t141 = t129 * t133;
	t132 = sin(qJ(1));
	t140 = t132 * t129;
	t131 = sin(qJ(2));
	t139 = t132 * t131;
	t138 = t132 * t133;
	t134 = cos(qJ(1));
	t137 = t134 * t129;
	t136 = t134 * t131;
	t135 = t134 * t133;
	t130 = cos(pkin(6));
	t128 = qJ(3) + qJ(4);
	t127 = cos(t128);
	t126 = sin(t128);
	t125 = t130 * t138 + t136;
	t124 = -t130 * t135 + t139;
	t1 = [0, t140, t125, t125, 0, (-t130 * t139 + t135) * t127 + t126 * t140; 0, -t137, t124, t124, 0, (t130 * t136 + t138) * t127 - t126 * t137; 1, t130, -t141, -t141, 0, t129 * t131 * t127 + t130 * t126;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end