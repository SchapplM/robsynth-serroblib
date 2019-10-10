% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRRPPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4,theta5]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:24
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRPPR5_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR5_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPPR5_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPPR5_jacobig_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:24:03
	% EndTime: 2019-10-10 11:24:03
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:24:03
	% EndTime: 2019-10-10 11:24:03
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:24:03
	% EndTime: 2019-10-10 11:24:03
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t43 = sin(pkin(6));
	t1 = [0, sin(qJ(1)) * t43, 0, 0, 0, 0; 0, -cos(qJ(1)) * t43, 0, 0, 0, 0; 1, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:24:03
	% EndTime: 2019-10-10 11:24:03
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
	% StartTime: 2019-10-10 11:24:03
	% EndTime: 2019-10-10 11:24:03
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
	% StartTime: 2019-10-10 11:24:03
	% EndTime: 2019-10-10 11:24:03
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (4->4), mult. (9->8), div. (0->0), fcn. (17->6), ass. (0->8)
	t103 = cos(pkin(6));
	t106 = cos(qJ(2));
	t108 = t103 * t106;
	t107 = cos(qJ(1));
	t105 = sin(qJ(1));
	t104 = sin(qJ(2));
	t102 = sin(pkin(6));
	t1 = [0, t105 * t102, t104 * t107 + t105 * t108, 0, 0, 0; 0, -t107 * t102, t104 * t105 - t107 * t108, 0, 0, 0; 1, t103, -t102 * t106, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:24:03
	% EndTime: 2019-10-10 11:24:03
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (15->10), mult. (24->18), div. (0->0), fcn. (40->8), ass. (0->16)
	t113 = sin(pkin(6));
	t116 = sin(qJ(1));
	t124 = t116 * t113;
	t115 = sin(qJ(2));
	t123 = t116 * t115;
	t117 = cos(qJ(2));
	t122 = t116 * t117;
	t118 = cos(qJ(1));
	t121 = t118 * t113;
	t120 = t118 * t115;
	t119 = t118 * t117;
	t114 = cos(pkin(6));
	t112 = qJ(3) + pkin(11);
	t111 = cos(t112);
	t110 = sin(t112);
	t1 = [0, t124, t114 * t122 + t120, 0, 0, (-t114 * t123 + t119) * t110 - t111 * t124; 0, -t121, -t114 * t119 + t123, 0, 0, (t114 * t120 + t122) * t110 + t111 * t121; 1, t114, -t113 * t117, 0, 0, t113 * t115 * t110 - t114 * t111;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end