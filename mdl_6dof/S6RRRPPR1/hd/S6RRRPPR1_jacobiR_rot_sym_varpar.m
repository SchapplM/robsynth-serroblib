% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4,theta5]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:17
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPPR1_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR1_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPPR1_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR1_jacobiR_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:17:04
	% EndTime: 2019-10-10 11:17:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:17:04
	% EndTime: 2019-10-10 11:17:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8, 0, 0, 0, 0, 0; t9, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0, 0; -t8, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:17:04
	% EndTime: 2019-10-10 11:17:04
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (7->7), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->9)
	t10 = sin(qJ(1));
	t9 = sin(qJ(2));
	t16 = t10 * t9;
	t12 = cos(qJ(1));
	t15 = t12 * t9;
	t11 = cos(qJ(2));
	t14 = t10 * t11;
	t13 = t12 * t11;
	t1 = [-t14, -t15, 0, 0, 0, 0; t13, -t16, 0, 0, 0, 0; 0, t11, 0, 0, 0, 0; t16, -t13, 0, 0, 0, 0; -t15, -t14, 0, 0, 0, 0; 0, -t9, 0, 0, 0, 0; t12, 0, 0, 0, 0, 0; t10, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:17:04
	% EndTime: 2019-10-10 11:17:04
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (28->13), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->10)
	t22 = qJ(2) + qJ(3);
	t20 = sin(t22);
	t23 = sin(qJ(1));
	t28 = t23 * t20;
	t21 = cos(t22);
	t27 = t23 * t21;
	t24 = cos(qJ(1));
	t26 = t24 * t20;
	t25 = t24 * t21;
	t1 = [-t27, -t26, -t26, 0, 0, 0; t25, -t28, -t28, 0, 0, 0; 0, t21, t21, 0, 0, 0; t28, -t25, -t25, 0, 0, 0; -t26, -t27, -t27, 0, 0, 0; 0, -t20, -t20, 0, 0, 0; t24, 0, 0, 0, 0, 0; t23, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:17:04
	% EndTime: 2019-10-10 11:17:04
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (44->13), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->10)
	t27 = qJ(2) + qJ(3) + pkin(10);
	t25 = sin(t27);
	t28 = sin(qJ(1));
	t33 = t28 * t25;
	t26 = cos(t27);
	t32 = t28 * t26;
	t29 = cos(qJ(1));
	t31 = t29 * t25;
	t30 = t29 * t26;
	t1 = [-t32, -t31, -t31, 0, 0, 0; t30, -t33, -t33, 0, 0, 0; 0, t26, t26, 0, 0, 0; t33, -t30, -t30, 0, 0, 0; -t31, -t32, -t32, 0, 0, 0; 0, -t25, -t25, 0, 0, 0; t29, 0, 0, 0, 0, 0; t28, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:17:04
	% EndTime: 2019-10-10 11:17:04
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (59->12), mult. (38->18), div. (0->0), fcn. (66->6), ass. (0->20)
	t92 = qJ(2) + qJ(3) + pkin(10);
	t91 = cos(t92);
	t93 = sin(pkin(11));
	t103 = t91 * t93;
	t95 = sin(qJ(1));
	t102 = t95 * t93;
	t94 = cos(pkin(11));
	t101 = t95 * t94;
	t96 = cos(qJ(1));
	t100 = t96 * t93;
	t99 = t96 * t94;
	t90 = sin(t92);
	t98 = t90 * t101;
	t97 = t90 * t99;
	t89 = t96 * t91;
	t88 = t95 * t91;
	t87 = t91 * t94;
	t86 = t90 * t100;
	t85 = t90 * t102;
	t1 = [-t91 * t101 + t100, -t97, -t97, 0, 0, 0; t91 * t99 + t102, -t98, -t98, 0, 0, 0; 0, t87, t87, 0, 0, 0; t91 * t102 + t99, t86, t86, 0, 0, 0; -t91 * t100 + t101, t85, t85, 0, 0, 0; 0, -t103, -t103, 0, 0, 0; -t95 * t90, t89, t89, 0, 0, 0; t96 * t90, t88, t88, 0, 0, 0; 0, t90, t90, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:17:04
	% EndTime: 2019-10-10 11:17:04
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (107->17), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->25)
	t108 = qJ(2) + qJ(3) + pkin(10);
	t105 = cos(t108);
	t109 = pkin(11) + qJ(6);
	t106 = sin(t109);
	t118 = t105 * t106;
	t110 = sin(qJ(1));
	t117 = t110 * t106;
	t107 = cos(t109);
	t116 = t110 * t107;
	t111 = cos(qJ(1));
	t115 = t111 * t106;
	t114 = t111 * t107;
	t104 = sin(t108);
	t113 = t104 * t116;
	t112 = t104 * t114;
	t103 = t111 * t105;
	t102 = t110 * t105;
	t101 = t105 * t107;
	t100 = t104 * t115;
	t99 = t104 * t117;
	t98 = t105 * t114 + t117;
	t97 = -t105 * t115 + t116;
	t96 = -t105 * t116 + t115;
	t95 = t105 * t117 + t114;
	t1 = [t96, -t112, -t112, 0, 0, t97; t98, -t113, -t113, 0, 0, -t95; 0, t101, t101, 0, 0, -t104 * t106; t95, t100, t100, 0, 0, -t98; t97, t99, t99, 0, 0, t96; 0, -t118, -t118, 0, 0, -t104 * t107; -t110 * t104, t103, t103, 0, 0, 0; t111 * t104, t102, t102, 0, 0, 0; 0, t104, t104, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end