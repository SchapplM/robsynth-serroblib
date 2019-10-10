% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRRR10
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
% Datum: 2019-10-10 11:05
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRPRRR10_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR10_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR10_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR10_jacobig_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:05:08
	% EndTime: 2019-10-10 11:05:08
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:05:08
	% EndTime: 2019-10-10 11:05:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:05:08
	% EndTime: 2019-10-10 11:05:08
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t43 = sin(pkin(6));
	t1 = [0, sin(qJ(1)) * t43, 0, 0, 0, 0; 0, -cos(qJ(1)) * t43, 0, 0, 0, 0; 1, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:05:08
	% EndTime: 2019-10-10 11:05:08
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t60 = sin(pkin(6));
	t1 = [0, sin(qJ(1)) * t60, 0, 0, 0, 0; 0, -cos(qJ(1)) * t60, 0, 0, 0, 0; 1, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:05:08
	% EndTime: 2019-10-10 11:05:08
	% DurationCPUTime: 0.04s
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
	% StartTime: 2019-10-10 11:05:08
	% EndTime: 2019-10-10 11:05:08
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (15->10), mult. (24->18), div. (0->0), fcn. (40->8), ass. (0->16)
	t111 = sin(pkin(6));
	t114 = sin(qJ(1));
	t122 = t114 * t111;
	t113 = sin(qJ(2));
	t121 = t114 * t113;
	t115 = cos(qJ(2));
	t120 = t114 * t115;
	t116 = cos(qJ(1));
	t119 = t116 * t111;
	t118 = t116 * t113;
	t117 = t116 * t115;
	t112 = cos(pkin(6));
	t110 = pkin(12) + qJ(4);
	t109 = cos(t110);
	t108 = sin(t110);
	t1 = [0, t122, 0, t112 * t120 + t118, (-t112 * t121 + t117) * t108 - t109 * t122, 0; 0, -t119, 0, -t112 * t117 + t121, (t112 * t118 + t120) * t108 + t109 * t119, 0; 1, t112, 0, -t111 * t115, t111 * t113 * t108 - t112 * t109, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:05:08
	% EndTime: 2019-10-10 11:05:08
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (26->10), mult. (39->18), div. (0->0), fcn. (63->8), ass. (0->19)
	t127 = sin(pkin(6));
	t130 = sin(qJ(1));
	t138 = t130 * t127;
	t129 = sin(qJ(2));
	t137 = t130 * t129;
	t131 = cos(qJ(2));
	t136 = t130 * t131;
	t132 = cos(qJ(1));
	t135 = t132 * t127;
	t134 = t132 * t129;
	t133 = t132 * t131;
	t128 = cos(pkin(6));
	t126 = pkin(12) + qJ(4);
	t125 = cos(t126);
	t124 = sin(t126);
	t123 = t127 * t129 * t124 - t128 * t125;
	t122 = (-t128 * t137 + t133) * t124 - t125 * t138;
	t121 = (t128 * t134 + t136) * t124 + t125 * t135;
	t1 = [0, t138, 0, t128 * t136 + t134, t122, t122; 0, -t135, 0, -t128 * t133 + t137, t121, t121; 1, t128, 0, -t127 * t131, t123, t123;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end