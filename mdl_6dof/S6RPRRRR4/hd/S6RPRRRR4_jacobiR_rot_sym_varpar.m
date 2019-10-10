% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRRRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:04
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRRR4_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR4_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRR4_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR4_jacobiR_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:04:21
	% EndTime: 2019-10-10 09:04:21
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:04:21
	% EndTime: 2019-10-10 09:04:21
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
	% StartTime: 2019-10-10 09:04:21
	% EndTime: 2019-10-10 09:04:21
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (10->4), ass. (0->5)
	t10 = cos(qJ(1));
	t9 = sin(qJ(1));
	t8 = cos(pkin(11));
	t7 = sin(pkin(11));
	t1 = [-t9 * t8, 0, 0, 0, 0, 0; t10 * t8, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t9 * t7, 0, 0, 0, 0, 0; -t10 * t7, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t10, 0, 0, 0, 0, 0; t9, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:04:21
	% EndTime: 2019-10-10 09:04:21
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (17->8), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->10)
	t14 = pkin(11) + qJ(3);
	t12 = sin(t14);
	t15 = sin(qJ(1));
	t20 = t15 * t12;
	t13 = cos(t14);
	t19 = t15 * t13;
	t16 = cos(qJ(1));
	t18 = t16 * t12;
	t17 = t16 * t13;
	t1 = [-t19, 0, -t18, 0, 0, 0; t17, 0, -t20, 0, 0, 0; 0, 0, t13, 0, 0, 0; t20, 0, -t17, 0, 0, 0; -t18, 0, -t19, 0, 0, 0; 0, 0, -t12, 0, 0, 0; t16, 0, 0, 0, 0, 0; t15, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:04:21
	% EndTime: 2019-10-10 09:04:21
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (44->13), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->10)
	t23 = pkin(11) + qJ(3) + qJ(4);
	t21 = sin(t23);
	t24 = sin(qJ(1));
	t29 = t24 * t21;
	t22 = cos(t23);
	t28 = t24 * t22;
	t25 = cos(qJ(1));
	t27 = t25 * t21;
	t26 = t25 * t22;
	t1 = [-t28, 0, -t27, -t27, 0, 0; t26, 0, -t29, -t29, 0, 0; 0, 0, t22, t22, 0, 0; t29, 0, -t26, -t26, 0, 0; -t27, 0, -t28, -t28, 0, 0; 0, 0, -t21, -t21, 0, 0; t25, 0, 0, 0, 0, 0; t24, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:04:21
	% EndTime: 2019-10-10 09:04:21
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (83->18), mult. (16->4), div. (0->0), fcn. (40->4), ass. (0->10)
	t29 = pkin(11) + qJ(3) + qJ(4) + qJ(5);
	t27 = sin(t29);
	t30 = sin(qJ(1));
	t35 = t30 * t27;
	t28 = cos(t29);
	t34 = t30 * t28;
	t31 = cos(qJ(1));
	t33 = t31 * t27;
	t32 = t31 * t28;
	t1 = [-t34, 0, -t33, -t33, -t33, 0; t32, 0, -t35, -t35, -t35, 0; 0, 0, t28, t28, t28, 0; t35, 0, -t32, -t32, -t32, 0; -t33, 0, -t34, -t34, -t34, 0; 0, 0, -t27, -t27, -t27, 0; t31, 0, 0, 0, 0, 0; t30, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:04:22
	% EndTime: 2019-10-10 09:04:22
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (137->19), mult. (64->20), div. (0->0), fcn. (111->6), ass. (0->24)
	t107 = pkin(11) + qJ(3) + qJ(4) + qJ(5);
	t106 = cos(t107);
	t108 = sin(qJ(6));
	t118 = t106 * t108;
	t109 = sin(qJ(1));
	t117 = t109 * t108;
	t110 = cos(qJ(6));
	t116 = t109 * t110;
	t111 = cos(qJ(1));
	t115 = t111 * t108;
	t114 = t111 * t110;
	t105 = sin(t107);
	t113 = t105 * t116;
	t112 = t105 * t114;
	t104 = t111 * t106;
	t103 = t106 * t110;
	t102 = t109 * t106;
	t101 = t105 * t115;
	t100 = t105 * t117;
	t99 = t106 * t114 + t117;
	t98 = -t106 * t115 + t116;
	t97 = -t106 * t116 + t115;
	t96 = t106 * t117 + t114;
	t1 = [t97, 0, -t112, -t112, -t112, t98; t99, 0, -t113, -t113, -t113, -t96; 0, 0, t103, t103, t103, -t105 * t108; t96, 0, t101, t101, t101, -t99; t98, 0, t100, t100, t100, t97; 0, 0, -t118, -t118, -t118, -t105 * t110; -t109 * t105, 0, t104, t104, t104, 0; t111 * t105, 0, t102, t102, t102, 0; 0, 0, t105, t105, t105, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end