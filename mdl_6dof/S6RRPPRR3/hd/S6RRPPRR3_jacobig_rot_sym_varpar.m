% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3,theta4]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:39
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRPPRR3_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR3_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_jacobig_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:39:18
	% EndTime: 2019-10-10 09:39:18
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:39:18
	% EndTime: 2019-10-10 09:39:18
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:39:18
	% EndTime: 2019-10-10 09:39:18
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t43 = sin(pkin(6));
	t1 = [0, sin(qJ(1)) * t43, 0, 0, 0, 0; 0, -cos(qJ(1)) * t43, 0, 0, 0, 0; 1, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:39:19
	% EndTime: 2019-10-10 09:39:19
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t59 = sin(pkin(6));
	t1 = [0, sin(qJ(1)) * t59, 0, 0, 0, 0; 0, -cos(qJ(1)) * t59, 0, 0, 0, 0; 1, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:39:19
	% EndTime: 2019-10-10 09:39:19
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t87 = sin(pkin(6));
	t1 = [0, sin(qJ(1)) * t87, 0, 0, 0, 0; 0, -cos(qJ(1)) * t87, 0, 0, 0, 0; 1, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:39:19
	% EndTime: 2019-10-10 09:39:19
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (8->6), mult. (22->12), div. (0->0), fcn. (35->8), ass. (0->12)
	t101 = cos(pkin(11));
	t103 = sin(qJ(2));
	t105 = cos(qJ(2));
	t99 = sin(pkin(11));
	t107 = t101 * t105 - t103 * t99;
	t106 = cos(qJ(1));
	t104 = sin(qJ(1));
	t102 = cos(pkin(6));
	t100 = sin(pkin(6));
	t98 = -t103 * t101 - t105 * t99;
	t97 = t107 * t102;
	t1 = [0, t104 * t100, 0, 0, t104 * t97 - t106 * t98, 0; 0, -t106 * t100, 0, 0, -t104 * t98 - t106 * t97, 0; 1, t102, 0, 0, -t107 * t100, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:39:19
	% EndTime: 2019-10-10 09:39:19
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (24->11), mult. (50->24), div. (0->0), fcn. (76->10), ass. (0->18)
	t139 = sin(pkin(6));
	t143 = sin(qJ(1));
	t149 = t143 * t139;
	t145 = cos(qJ(1));
	t148 = t145 * t139;
	t138 = sin(pkin(11));
	t140 = cos(pkin(11));
	t142 = sin(qJ(2));
	t144 = cos(qJ(2));
	t147 = t144 * t138 + t142 * t140;
	t146 = t142 * t138 - t144 * t140;
	t141 = cos(pkin(6));
	t137 = pkin(12) + qJ(5);
	t136 = cos(t137);
	t135 = sin(t137);
	t132 = t147 * t141;
	t131 = t146 * t141;
	t1 = [0, t149, 0, 0, -t143 * t131 + t145 * t147, (-t143 * t132 - t145 * t146) * t135 - t136 * t149; 0, -t148, 0, 0, t145 * t131 + t143 * t147, (t145 * t132 - t143 * t146) * t135 + t136 * t148; 1, t141, 0, 0, t146 * t139, t147 * t135 * t139 - t141 * t136;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end