% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:02
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RRRRR2_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR2_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR2_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5RRRRR2_jacobiR_rot_sym_varpar: pkin has to be [2x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:02:47
	% EndTime: 2019-10-09 21:02:47
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:02:47
	% EndTime: 2019-10-09 21:02:47
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8, 0, 0, 0, 0; t9, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0; -t8, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:02:47
	% EndTime: 2019-10-09 21:02:47
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (14->7), mult. (0->0), div. (0->0), fcn. (8->2), ass. (0->4)
	t16 = qJ(1) + qJ(2);
	t15 = cos(t16);
	t14 = sin(t16);
	t1 = [-t14, -t14, 0, 0, 0; t15, t15, 0, 0, 0; 0, 0, 0, 0, 0; -t15, -t15, 0, 0, 0; -t14, -t14, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:02:47
	% EndTime: 2019-10-09 21:02:47
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (25->10), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->10)
	t21 = qJ(1) + qJ(2);
	t19 = sin(t21);
	t23 = cos(qJ(3));
	t25 = t19 * t23;
	t20 = cos(t21);
	t22 = sin(qJ(3));
	t24 = t20 * t22;
	t18 = t20 * t23;
	t17 = t19 * t22;
	t1 = [-t25, -t25, -t24, 0, 0; t18, t18, -t17, 0, 0; 0, 0, t23, 0, 0; t17, t17, -t18, 0, 0; -t24, -t24, -t25, 0, 0; 0, 0, -t22, 0, 0; t20, t20, 0, 0, 0; t19, t19, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:02:47
	% EndTime: 2019-10-09 21:02:47
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (54->16), mult. (16->4), div. (0->0), fcn. (40->4), ass. (0->11)
	t34 = qJ(1) + qJ(2);
	t30 = sin(t34);
	t33 = qJ(3) + qJ(4);
	t31 = cos(t33);
	t36 = t30 * t31;
	t29 = sin(t33);
	t32 = cos(t34);
	t35 = t32 * t29;
	t28 = t32 * t31;
	t27 = t30 * t29;
	t1 = [-t36, -t36, -t35, -t35, 0; t28, t28, -t27, -t27, 0; 0, 0, t31, t31, 0; t27, t27, -t28, -t28, 0; -t35, -t35, -t36, -t36, 0; 0, 0, -t29, -t29, 0; t32, t32, 0, 0, 0; t30, t30, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:02:47
	% EndTime: 2019-10-09 21:02:47
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (98->18), mult. (66->20), div. (0->0), fcn. (114->6), ass. (0->27)
	t104 = qJ(3) + qJ(4);
	t100 = sin(t104);
	t105 = qJ(1) + qJ(2);
	t101 = sin(t105);
	t115 = t101 * t100;
	t106 = sin(qJ(5));
	t114 = t101 * t106;
	t107 = cos(qJ(5));
	t113 = t101 * t107;
	t102 = cos(t104);
	t112 = t102 * t106;
	t99 = t102 * t107;
	t103 = cos(t105);
	t111 = t103 * t106;
	t110 = t103 * t107;
	t109 = t100 * t113;
	t108 = t100 * t110;
	t98 = t103 * t102;
	t97 = t103 * t100;
	t96 = t101 * t102;
	t95 = t100 * t111;
	t94 = t100 * t114;
	t93 = t102 * t110 + t114;
	t92 = -t102 * t111 + t113;
	t91 = -t101 * t99 + t111;
	t90 = t101 * t112 + t110;
	t1 = [t91, t91, -t108, -t108, t92; t93, t93, -t109, -t109, -t90; 0, 0, t99, t99, -t100 * t106; t90, t90, t95, t95, -t93; t92, t92, t94, t94, t91; 0, 0, -t112, -t112, -t100 * t107; -t115, -t115, t98, t98, 0; t97, t97, t96, t96, 0; 0, 0, t100, t100, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end