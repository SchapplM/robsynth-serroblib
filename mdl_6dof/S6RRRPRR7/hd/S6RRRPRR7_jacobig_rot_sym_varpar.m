% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRRPRR7
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:04
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRPRR7_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR7_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR7_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR7_jacobig_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:04:12
	% EndTime: 2019-10-10 12:04:12
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:04:12
	% EndTime: 2019-10-10 12:04:12
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:04:12
	% EndTime: 2019-10-10 12:04:12
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t43 = sin(pkin(6));
	t1 = [0, sin(qJ(1)) * t43, 0, 0, 0, 0; 0, -cos(qJ(1)) * t43, 0, 0, 0, 0; 1, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:04:12
	% EndTime: 2019-10-10 12:04:12
	% DurationCPUTime: 0.03s
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
	% StartTime: 2019-10-10 12:04:12
	% EndTime: 2019-10-10 12:04:12
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (4->4), mult. (9->8), div. (0->0), fcn. (17->6), ass. (0->8)
	t75 = cos(pkin(6));
	t78 = cos(qJ(2));
	t80 = t75 * t78;
	t79 = cos(qJ(1));
	t77 = sin(qJ(1));
	t76 = sin(qJ(2));
	t74 = sin(pkin(6));
	t1 = [0, t77 * t74, t79 * t76 + t77 * t80, 0, 0, 0; 0, -t79 * t74, t77 * t76 - t79 * t80, 0, 0, 0; 1, t75, -t74 * t78, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:04:12
	% EndTime: 2019-10-10 12:04:12
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (7->5), mult. (16->8), div. (0->0), fcn. (29->6), ass. (0->11)
	t89 = sin(pkin(6));
	t93 = cos(qJ(2));
	t96 = t89 * t93;
	t90 = cos(pkin(6));
	t95 = t90 * t93;
	t94 = cos(qJ(1));
	t92 = sin(qJ(1));
	t91 = sin(qJ(2));
	t88 = t94 * t91 + t92 * t95;
	t87 = t92 * t91 - t94 * t95;
	t1 = [0, t92 * t89, t88, 0, t88, 0; 0, -t94 * t89, t87, 0, t87, 0; 1, t90, -t96, 0, -t96, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:04:12
	% EndTime: 2019-10-10 12:04:12
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (24->11), mult. (31->18), div. (0->0), fcn. (52->8), ass. (0->19)
	t133 = sin(pkin(6));
	t137 = cos(qJ(2));
	t145 = t133 * t137;
	t136 = sin(qJ(1));
	t144 = t136 * t133;
	t135 = sin(qJ(2));
	t143 = t136 * t135;
	t142 = t136 * t137;
	t138 = cos(qJ(1));
	t141 = t138 * t133;
	t140 = t138 * t135;
	t139 = t138 * t137;
	t134 = cos(pkin(6));
	t132 = qJ(3) + pkin(12) + qJ(5);
	t131 = cos(t132);
	t130 = sin(t132);
	t129 = t134 * t142 + t140;
	t128 = -t134 * t139 + t143;
	t1 = [0, t144, t129, 0, t129, (-t134 * t143 + t139) * t130 - t131 * t144; 0, -t141, t128, 0, t128, (t134 * t140 + t142) * t130 + t131 * t141; 1, t134, -t145, 0, -t145, t133 * t135 * t130 - t134 * t131;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end