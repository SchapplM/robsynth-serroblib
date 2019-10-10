% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRRRPR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:42
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRRPR7_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR7_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPR7_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR7_jacobig_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:42:14
	% EndTime: 2019-10-10 12:42:14
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:42:14
	% EndTime: 2019-10-10 12:42:14
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:42:14
	% EndTime: 2019-10-10 12:42:14
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t43 = sin(pkin(6));
	t1 = [0, sin(qJ(1)) * t43, 0, 0, 0, 0; 0, -cos(qJ(1)) * t43, 0, 0, 0, 0; 1, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:42:14
	% EndTime: 2019-10-10 12:42:14
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
	% StartTime: 2019-10-10 12:42:14
	% EndTime: 2019-10-10 12:42:14
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
	% StartTime: 2019-10-10 12:42:14
	% EndTime: 2019-10-10 12:42:14
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (7->5), mult. (16->8), div. (0->0), fcn. (29->6), ass. (0->11)
	t91 = sin(pkin(6));
	t95 = cos(qJ(2));
	t98 = t91 * t95;
	t92 = cos(pkin(6));
	t97 = t92 * t95;
	t96 = cos(qJ(1));
	t94 = sin(qJ(1));
	t93 = sin(qJ(2));
	t90 = t96 * t93 + t94 * t97;
	t89 = t94 * t93 - t96 * t97;
	t1 = [0, t94 * t91, t90, t90, 0, 0; 0, -t96 * t91, t89, t89, 0, 0; 1, t92, -t98, -t98, 0, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:42:15
	% EndTime: 2019-10-10 12:42:15
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (24->11), mult. (31->18), div. (0->0), fcn. (52->8), ass. (0->19)
	t135 = sin(pkin(6));
	t139 = cos(qJ(2));
	t147 = t135 * t139;
	t138 = sin(qJ(1));
	t146 = t138 * t135;
	t137 = sin(qJ(2));
	t145 = t138 * t137;
	t144 = t138 * t139;
	t140 = cos(qJ(1));
	t143 = t140 * t135;
	t142 = t140 * t137;
	t141 = t140 * t139;
	t136 = cos(pkin(6));
	t134 = qJ(3) + qJ(4) + pkin(12);
	t133 = cos(t134);
	t132 = sin(t134);
	t131 = t136 * t144 + t142;
	t130 = -t136 * t141 + t145;
	t1 = [0, t146, t131, t131, 0, (-t136 * t145 + t141) * t132 - t133 * t146; 0, -t143, t130, t130, 0, (t136 * t142 + t144) * t132 + t133 * t143; 1, t136, -t147, -t147, 0, t135 * t137 * t132 - t136 * t133;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end