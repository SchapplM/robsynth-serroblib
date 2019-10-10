% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRPRR7
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:56
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRPRR7_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR7_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR7_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR7_jacobiR_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:56:39
	% EndTime: 2019-10-10 00:56:39
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:56:39
	% EndTime: 2019-10-10 00:56:39
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
	% StartTime: 2019-10-10 00:56:39
	% EndTime: 2019-10-10 00:56:39
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t6 = cos(qJ(1));
	t5 = sin(qJ(1));
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t5, 0, 0, 0, 0, 0; -t6, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t6, 0, 0, 0, 0, 0; t5, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:56:39
	% EndTime: 2019-10-10 00:56:39
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (5->5), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->9)
	t11 = sin(qJ(3));
	t12 = sin(qJ(1));
	t16 = t12 * t11;
	t13 = cos(qJ(3));
	t14 = cos(qJ(1));
	t15 = t14 * t13;
	t10 = t14 * t11;
	t9 = t12 * t13;
	t1 = [t10, 0, t9, 0, 0, 0; t16, 0, -t15, 0, 0, 0; 0, 0, -t11, 0, 0, 0; t15, 0, -t16, 0, 0, 0; t9, 0, t10, 0, 0, 0; 0, 0, -t13, 0, 0, 0; -t12, 0, 0, 0, 0, 0; t14, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:56:39
	% EndTime: 2019-10-10 00:56:39
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (15->6), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->10)
	t17 = qJ(3) + pkin(10);
	t15 = sin(t17);
	t18 = sin(qJ(1));
	t21 = t18 * t15;
	t16 = cos(t17);
	t19 = cos(qJ(1));
	t20 = t19 * t16;
	t14 = t19 * t15;
	t13 = t18 * t16;
	t1 = [t14, 0, t13, 0, 0, 0; t21, 0, -t20, 0, 0, 0; 0, 0, -t15, 0, 0, 0; t20, 0, -t21, 0, 0, 0; t13, 0, t14, 0, 0, 0; 0, 0, -t16, 0, 0, 0; -t18, 0, 0, 0, 0, 0; t19, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:56:39
	% EndTime: 2019-10-10 00:56:39
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (41->10), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->10)
	t26 = qJ(3) + pkin(10) + qJ(5);
	t24 = sin(t26);
	t27 = sin(qJ(1));
	t30 = t27 * t24;
	t25 = cos(t26);
	t28 = cos(qJ(1));
	t29 = t28 * t25;
	t23 = t28 * t24;
	t22 = t27 * t25;
	t1 = [t23, 0, t22, 0, t22, 0; t30, 0, -t29, 0, -t29, 0; 0, 0, -t24, 0, -t24, 0; t29, 0, -t30, 0, -t30, 0; t22, 0, t23, 0, t23, 0; 0, 0, -t25, 0, -t25, 0; -t27, 0, 0, 0, 0, 0; t28, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:56:39
	% EndTime: 2019-10-10 00:56:39
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (80->19), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->24)
	t95 = qJ(3) + pkin(10) + qJ(5);
	t93 = sin(t95);
	t98 = cos(qJ(6));
	t107 = t93 * t98;
	t96 = sin(qJ(6));
	t97 = sin(qJ(1));
	t106 = t97 * t96;
	t105 = t97 * t98;
	t99 = cos(qJ(1));
	t104 = t99 * t93;
	t103 = t99 * t96;
	t102 = t99 * t98;
	t94 = cos(t95);
	t101 = t94 * t106;
	t100 = t94 * t102;
	t92 = t97 * t93;
	t91 = t93 * t96;
	t90 = t94 * t103;
	t89 = t94 * t105;
	t88 = t93 * t102 - t106;
	t87 = t93 * t103 + t105;
	t86 = t93 * t105 + t103;
	t85 = -t93 * t106 + t102;
	t1 = [t88, 0, t89, 0, t89, t85; t86, 0, -t100, 0, -t100, t87; 0, 0, -t107, 0, -t107, -t94 * t96; -t87, 0, -t101, 0, -t101, -t86; t85, 0, t90, 0, t90, t88; 0, 0, t91, 0, t91, -t94 * t98; -t99 * t94, 0, t92, 0, t92, 0; -t97 * t94, 0, -t104, 0, -t104, 0; 0, 0, t94, 0, t94, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end