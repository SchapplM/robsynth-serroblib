% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPPRRR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:11
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPPRRR7_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR7_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRR7_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR7_jacobiR_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:11:46
	% EndTime: 2019-10-10 00:11:46
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:11:46
	% EndTime: 2019-10-10 00:11:46
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
	% StartTime: 2019-10-10 00:11:46
	% EndTime: 2019-10-10 00:11:46
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
	% StartTime: 2019-10-10 00:11:46
	% EndTime: 2019-10-10 00:11:46
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (1->1), mult. (4->4), div. (0->0), fcn. (10->4), ass. (0->5)
	t10 = cos(qJ(1));
	t9 = sin(qJ(1));
	t8 = cos(pkin(10));
	t7 = sin(pkin(10));
	t1 = [t10 * t7, 0, 0, 0, 0, 0; t9 * t7, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t10 * t8, 0, 0, 0, 0, 0; t9 * t8, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0, 0; t10, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:11:46
	% EndTime: 2019-10-10 00:11:46
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (15->6), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->10)
	t16 = pkin(10) + qJ(4);
	t14 = sin(t16);
	t17 = sin(qJ(1));
	t20 = t17 * t14;
	t15 = cos(t16);
	t18 = cos(qJ(1));
	t19 = t18 * t15;
	t13 = t18 * t14;
	t12 = t17 * t15;
	t1 = [t13, 0, 0, t12, 0, 0; t20, 0, 0, -t19, 0, 0; 0, 0, 0, -t14, 0, 0; t19, 0, 0, -t20, 0, 0; t12, 0, 0, t13, 0, 0; 0, 0, 0, -t15, 0, 0; -t17, 0, 0, 0, 0, 0; t18, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:11:46
	% EndTime: 2019-10-10 00:11:46
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (41->10), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->10)
	t27 = pkin(10) + qJ(4) + qJ(5);
	t25 = sin(t27);
	t28 = sin(qJ(1));
	t31 = t28 * t25;
	t26 = cos(t27);
	t29 = cos(qJ(1));
	t30 = t29 * t26;
	t24 = t29 * t25;
	t23 = t28 * t26;
	t1 = [t24, 0, 0, t23, t23, 0; t31, 0, 0, -t30, -t30, 0; 0, 0, 0, -t25, -t25, 0; t30, 0, 0, -t31, -t31, 0; t23, 0, 0, t24, t24, 0; 0, 0, 0, -t26, -t26, 0; -t28, 0, 0, 0, 0, 0; t29, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:11:46
	% EndTime: 2019-10-10 00:11:46
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (80->19), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->24)
	t96 = pkin(10) + qJ(4) + qJ(5);
	t94 = sin(t96);
	t99 = cos(qJ(6));
	t108 = t94 * t99;
	t97 = sin(qJ(6));
	t98 = sin(qJ(1));
	t107 = t98 * t97;
	t106 = t98 * t99;
	t100 = cos(qJ(1));
	t105 = t100 * t94;
	t104 = t100 * t97;
	t103 = t100 * t99;
	t95 = cos(t96);
	t102 = t95 * t107;
	t101 = t95 * t103;
	t93 = t98 * t94;
	t92 = t94 * t97;
	t91 = t95 * t104;
	t90 = t95 * t106;
	t89 = t94 * t103 - t107;
	t88 = t94 * t104 + t106;
	t87 = t94 * t106 + t104;
	t86 = -t94 * t107 + t103;
	t1 = [t89, 0, 0, t90, t90, t86; t87, 0, 0, -t101, -t101, t88; 0, 0, 0, -t108, -t108, -t95 * t97; -t88, 0, 0, -t102, -t102, -t87; t86, 0, 0, t91, t91, t89; 0, 0, 0, t92, t92, -t95 * t99; -t100 * t95, 0, 0, t93, t93, 0; -t98 * t95, 0, 0, -t105, -t105, 0; 0, 0, 0, t95, t95, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end