% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRRP6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:37
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRPRRP6_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP6_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRP6_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP6_jacobig_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:37:18
	% EndTime: 2019-10-10 10:37:18
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:37:18
	% EndTime: 2019-10-10 10:37:18
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:37:19
	% EndTime: 2019-10-10 10:37:19
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t43 = sin(pkin(6));
	t1 = [0, sin(qJ(1)) * t43, 0, 0, 0, 0; 0, -cos(qJ(1)) * t43, 0, 0, 0, 0; 1, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:37:19
	% EndTime: 2019-10-10 10:37:19
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t59 = sin(pkin(6));
	t1 = [0, sin(qJ(1)) * t59, 0, 0, 0, 0; 0, -cos(qJ(1)) * t59, 0, 0, 0, 0; 1, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:37:19
	% EndTime: 2019-10-10 10:37:19
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (8->5), mult. (22->12), div. (0->0), fcn. (35->8), ass. (0->12)
	t91 = sin(pkin(11));
	t93 = cos(pkin(11));
	t95 = sin(qJ(2));
	t97 = cos(qJ(2));
	t99 = t91 * t95 - t93 * t97;
	t98 = cos(qJ(1));
	t96 = sin(qJ(1));
	t94 = cos(pkin(6));
	t92 = sin(pkin(6));
	t90 = -t97 * t91 - t95 * t93;
	t89 = t99 * t94;
	t1 = [0, t96 * t92, 0, -t96 * t89 - t98 * t90, 0, 0; 0, -t98 * t92, 0, t98 * t89 - t96 * t90, 0, 0; 1, t94, 0, t99 * t92, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:37:19
	% EndTime: 2019-10-10 10:37:19
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (18->10), mult. (50->24), div. (0->0), fcn. (76->10), ass. (0->17)
	t131 = sin(pkin(6));
	t136 = sin(qJ(1));
	t143 = t136 * t131;
	t139 = cos(qJ(1));
	t142 = t139 * t131;
	t130 = sin(pkin(11));
	t132 = cos(pkin(11));
	t135 = sin(qJ(2));
	t138 = cos(qJ(2));
	t141 = t138 * t130 + t135 * t132;
	t140 = t135 * t130 - t138 * t132;
	t137 = cos(qJ(4));
	t134 = sin(qJ(4));
	t133 = cos(pkin(6));
	t127 = t141 * t133;
	t126 = t140 * t133;
	t1 = [0, t143, 0, -t136 * t126 + t139 * t141, (-t136 * t127 - t139 * t140) * t134 - t137 * t143, 0; 0, -t142, 0, t139 * t126 + t136 * t141, (t139 * t127 - t136 * t140) * t134 + t137 * t142, 0; 1, t133, 0, t140 * t131, t141 * t134 * t131 - t133 * t137, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:37:19
	% EndTime: 2019-10-10 10:37:19
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (18->10), mult. (50->24), div. (0->0), fcn. (76->10), ass. (0->17)
	t160 = sin(pkin(6));
	t165 = sin(qJ(1));
	t172 = t165 * t160;
	t168 = cos(qJ(1));
	t171 = t168 * t160;
	t159 = sin(pkin(11));
	t161 = cos(pkin(11));
	t164 = sin(qJ(2));
	t167 = cos(qJ(2));
	t170 = t167 * t159 + t164 * t161;
	t169 = t164 * t159 - t167 * t161;
	t166 = cos(qJ(4));
	t163 = sin(qJ(4));
	t162 = cos(pkin(6));
	t156 = t170 * t162;
	t155 = t169 * t162;
	t1 = [0, t172, 0, -t165 * t155 + t168 * t170, (-t165 * t156 - t168 * t169) * t163 - t166 * t172, 0; 0, -t171, 0, t168 * t155 + t165 * t170, (t168 * t156 - t165 * t169) * t163 + t166 * t171, 0; 1, t162, 0, t169 * t160, t170 * t163 * t160 - t162 * t166, 0;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end