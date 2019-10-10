% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRRPRP6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:44
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRPRP6_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP6_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRP6_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP6_jacobig_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:43:59
	% EndTime: 2019-10-10 11:43:59
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:43:59
	% EndTime: 2019-10-10 11:43:59
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:43:59
	% EndTime: 2019-10-10 11:43:59
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t43 = sin(pkin(6));
	t1 = [0, sin(qJ(1)) * t43, 0, 0, 0, 0; 0, -cos(qJ(1)) * t43, 0, 0, 0, 0; 1, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:43:59
	% EndTime: 2019-10-10 11:43:59
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
	% StartTime: 2019-10-10 11:43:59
	% EndTime: 2019-10-10 11:43:59
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
	% StartTime: 2019-10-10 11:44:00
	% EndTime: 2019-10-10 11:44:00
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (15->10), mult. (24->18), div. (0->0), fcn. (40->8), ass. (0->16)
	t112 = sin(pkin(6));
	t115 = sin(qJ(1));
	t123 = t115 * t112;
	t114 = sin(qJ(2));
	t122 = t115 * t114;
	t116 = cos(qJ(2));
	t121 = t115 * t116;
	t117 = cos(qJ(1));
	t120 = t117 * t112;
	t119 = t117 * t114;
	t118 = t117 * t116;
	t113 = cos(pkin(6));
	t111 = qJ(3) + pkin(11);
	t110 = cos(t111);
	t109 = sin(t111);
	t1 = [0, t123, t113 * t121 + t119, 0, (-t113 * t122 + t118) * t109 - t110 * t123, 0; 0, -t120, -t113 * t118 + t122, 0, (t113 * t119 + t121) * t109 + t110 * t120, 0; 1, t113, -t112 * t116, 0, t112 * t114 * t109 - t113 * t110, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:44:00
	% EndTime: 2019-10-10 11:44:00
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (15->10), mult. (24->18), div. (0->0), fcn. (40->8), ass. (0->16)
	t114 = sin(pkin(6));
	t117 = sin(qJ(1));
	t125 = t117 * t114;
	t116 = sin(qJ(2));
	t124 = t117 * t116;
	t118 = cos(qJ(2));
	t123 = t117 * t118;
	t119 = cos(qJ(1));
	t122 = t119 * t114;
	t121 = t119 * t116;
	t120 = t119 * t118;
	t115 = cos(pkin(6));
	t113 = qJ(3) + pkin(11);
	t112 = cos(t113);
	t111 = sin(t113);
	t1 = [0, t125, t115 * t123 + t121, 0, (-t115 * t124 + t120) * t111 - t112 * t125, 0; 0, -t122, -t115 * t120 + t124, 0, (t115 * t121 + t123) * t111 + t112 * t122, 0; 1, t115, -t114 * t118, 0, t114 * t116 * t111 - t115 * t112, 0;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end