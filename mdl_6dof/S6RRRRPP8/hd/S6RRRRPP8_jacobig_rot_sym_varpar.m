% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRRRPP8
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:31
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRRPP8_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP8_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPP8_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP8_jacobig_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:31:15
	% EndTime: 2019-10-10 12:31:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:31:15
	% EndTime: 2019-10-10 12:31:16
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:31:16
	% EndTime: 2019-10-10 12:31:16
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t43 = sin(pkin(6));
	t1 = [0, sin(qJ(1)) * t43, 0, 0, 0, 0; 0, -cos(qJ(1)) * t43, 0, 0, 0, 0; 1, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:31:16
	% EndTime: 2019-10-10 12:31:16
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
	% StartTime: 2019-10-10 12:31:16
	% EndTime: 2019-10-10 12:31:16
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->9), mult. (24->18), div. (0->0), fcn. (40->8), ass. (0->15)
	t101 = sin(pkin(6));
	t105 = sin(qJ(1));
	t114 = t105 * t101;
	t104 = sin(qJ(2));
	t113 = t105 * t104;
	t107 = cos(qJ(2));
	t112 = t105 * t107;
	t108 = cos(qJ(1));
	t111 = t108 * t101;
	t110 = t108 * t104;
	t109 = t108 * t107;
	t106 = cos(qJ(3));
	t103 = sin(qJ(3));
	t102 = cos(pkin(6));
	t1 = [0, t114, t102 * t112 + t110, (-t102 * t113 + t109) * t103 - t106 * t114, 0, 0; 0, -t111, -t102 * t109 + t113, (t102 * t110 + t112) * t103 + t106 * t111, 0, 0; 1, t102, -t101 * t107, t101 * t104 * t103 - t102 * t106, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:31:16
	% EndTime: 2019-10-10 12:31:16
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (9->9), mult. (24->18), div. (0->0), fcn. (40->8), ass. (0->15)
	t134 = sin(pkin(6));
	t138 = sin(qJ(1));
	t147 = t138 * t134;
	t137 = sin(qJ(2));
	t146 = t138 * t137;
	t140 = cos(qJ(2));
	t145 = t138 * t140;
	t141 = cos(qJ(1));
	t144 = t141 * t134;
	t143 = t141 * t137;
	t142 = t141 * t140;
	t139 = cos(qJ(3));
	t136 = sin(qJ(3));
	t135 = cos(pkin(6));
	t1 = [0, t147, t135 * t145 + t143, (-t135 * t146 + t142) * t136 - t139 * t147, 0, 0; 0, -t144, -t135 * t142 + t146, (t135 * t143 + t145) * t136 + t139 * t144, 0, 0; 1, t135, -t134 * t140, t134 * t136 * t137 - t135 * t139, 0, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:31:16
	% EndTime: 2019-10-10 12:31:16
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->9), mult. (24->18), div. (0->0), fcn. (40->8), ass. (0->15)
	t120 = sin(pkin(6));
	t124 = sin(qJ(1));
	t133 = t124 * t120;
	t123 = sin(qJ(2));
	t132 = t124 * t123;
	t126 = cos(qJ(2));
	t131 = t124 * t126;
	t127 = cos(qJ(1));
	t130 = t127 * t120;
	t129 = t127 * t123;
	t128 = t127 * t126;
	t125 = cos(qJ(3));
	t122 = sin(qJ(3));
	t121 = cos(pkin(6));
	t1 = [0, t133, t121 * t131 + t129, (-t121 * t132 + t128) * t122 - t125 * t133, 0, 0; 0, -t130, -t121 * t128 + t132, (t121 * t129 + t131) * t122 + t125 * t130, 0, 0; 1, t121, -t120 * t126, t120 * t123 * t122 - t121 * t125, 0, 0;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end