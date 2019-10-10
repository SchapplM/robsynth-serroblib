% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:03
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRPRRR9_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR9_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR9_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR9_jacobig_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:03:13
	% EndTime: 2019-10-10 11:03:13
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:03:13
	% EndTime: 2019-10-10 11:03:13
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:03:13
	% EndTime: 2019-10-10 11:03:13
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t43 = sin(pkin(6));
	t1 = [0, sin(qJ(1)) * t43, 0, 0, 0, 0; 0, -cos(qJ(1)) * t43, 0, 0, 0, 0; 1, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:03:13
	% EndTime: 2019-10-10 11:03:13
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t60 = sin(pkin(6));
	t1 = [0, sin(qJ(1)) * t60, 0, 0, 0, 0; 0, -cos(qJ(1)) * t60, 0, 0, 0, 0; 1, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:03:13
	% EndTime: 2019-10-10 11:03:13
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (4->4), mult. (9->8), div. (0->0), fcn. (17->6), ass. (0->8)
	t74 = cos(pkin(6));
	t77 = cos(qJ(2));
	t79 = t74 * t77;
	t78 = cos(qJ(1));
	t76 = sin(qJ(1));
	t75 = sin(qJ(2));
	t73 = sin(pkin(6));
	t1 = [0, t76 * t73, 0, t78 * t75 + t76 * t79, 0, 0; 0, -t78 * t73, 0, t76 * t75 - t78 * t79, 0, 0; 1, t74, 0, -t73 * t77, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:03:13
	% EndTime: 2019-10-10 11:03:13
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (7->5), mult. (16->8), div. (0->0), fcn. (29->6), ass. (0->11)
	t88 = sin(pkin(6));
	t92 = cos(qJ(2));
	t95 = t88 * t92;
	t89 = cos(pkin(6));
	t94 = t89 * t92;
	t93 = cos(qJ(1));
	t91 = sin(qJ(1));
	t90 = sin(qJ(2));
	t87 = t93 * t90 + t91 * t94;
	t86 = t91 * t90 - t93 * t94;
	t1 = [0, t91 * t88, 0, t87, t87, 0; 0, -t93 * t88, 0, t86, t86, 0; 1, t89, 0, -t95, -t95, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:03:13
	% EndTime: 2019-10-10 11:03:13
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (24->11), mult. (31->18), div. (0->0), fcn. (52->8), ass. (0->19)
	t132 = sin(pkin(6));
	t136 = cos(qJ(2));
	t144 = t132 * t136;
	t135 = sin(qJ(1));
	t143 = t135 * t132;
	t134 = sin(qJ(2));
	t142 = t135 * t134;
	t141 = t135 * t136;
	t137 = cos(qJ(1));
	t140 = t137 * t132;
	t139 = t137 * t134;
	t138 = t137 * t136;
	t133 = cos(pkin(6));
	t131 = pkin(12) + qJ(4) + qJ(5);
	t130 = cos(t131);
	t129 = sin(t131);
	t128 = t133 * t141 + t139;
	t127 = -t133 * t138 + t142;
	t1 = [0, t143, 0, t128, t128, (-t133 * t142 + t138) * t129 - t130 * t143; 0, -t140, 0, t127, t127, (t133 * t139 + t141) * t129 + t130 * t140; 1, t133, 0, -t144, -t144, t132 * t134 * t129 - t133 * t130;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end