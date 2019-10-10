% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:46
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRPRR1_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR1_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR1_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR1_jacobiR_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:46:14
	% EndTime: 2019-10-10 00:46:14
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:46:14
	% EndTime: 2019-10-10 00:46:14
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
	% StartTime: 2019-10-10 00:46:14
	% EndTime: 2019-10-10 00:46:14
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (7->4), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->4)
	t12 = qJ(1) + pkin(10);
	t11 = cos(t12);
	t10 = sin(t12);
	t1 = [-t10, 0, 0, 0, 0, 0; t11, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t11, 0, 0, 0, 0, 0; -t10, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:46:14
	% EndTime: 2019-10-10 00:46:14
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (17->8), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->10)
	t13 = qJ(1) + pkin(10);
	t11 = sin(t13);
	t14 = sin(qJ(3));
	t19 = t11 * t14;
	t15 = cos(qJ(3));
	t18 = t11 * t15;
	t12 = cos(t13);
	t17 = t12 * t14;
	t16 = t12 * t15;
	t1 = [-t18, 0, -t17, 0, 0, 0; t16, 0, -t19, 0, 0, 0; 0, 0, t15, 0, 0, 0; t19, 0, -t16, 0, 0, 0; -t17, 0, -t18, 0, 0, 0; 0, 0, -t14, 0, 0, 0; t12, 0, 0, 0, 0, 0; t11, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:46:14
	% EndTime: 2019-10-10 00:46:14
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (27->9), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->11)
	t21 = qJ(3) + pkin(11);
	t17 = sin(t21);
	t22 = qJ(1) + pkin(10);
	t18 = sin(t22);
	t26 = t18 * t17;
	t19 = cos(t21);
	t25 = t18 * t19;
	t20 = cos(t22);
	t24 = t20 * t17;
	t23 = t20 * t19;
	t1 = [-t25, 0, -t24, 0, 0, 0; t23, 0, -t26, 0, 0, 0; 0, 0, t19, 0, 0, 0; t26, 0, -t23, 0, 0, 0; -t24, 0, -t25, 0, 0, 0; 0, 0, -t17, 0, 0, 0; t20, 0, 0, 0, 0, 0; t18, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:46:14
	% EndTime: 2019-10-10 00:46:14
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (58->14), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->11)
	t29 = qJ(3) + pkin(11) + qJ(5);
	t25 = sin(t29);
	t30 = qJ(1) + pkin(10);
	t27 = sin(t30);
	t34 = t27 * t25;
	t26 = cos(t29);
	t33 = t27 * t26;
	t28 = cos(t30);
	t32 = t28 * t25;
	t31 = t28 * t26;
	t1 = [-t33, 0, -t32, 0, -t32, 0; t31, 0, -t34, 0, -t34, 0; 0, 0, t26, 0, t26, 0; t34, 0, -t31, 0, -t31, 0; -t32, 0, -t33, 0, -t33, 0; 0, 0, -t25, 0, -t25, 0; t28, 0, 0, 0, 0, 0; t27, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:46:14
	% EndTime: 2019-10-10 00:46:14
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (107->17), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->25)
	t108 = qJ(3) + pkin(11) + qJ(5);
	t105 = cos(t108);
	t110 = sin(qJ(6));
	t118 = t105 * t110;
	t109 = qJ(1) + pkin(10);
	t106 = sin(t109);
	t117 = t106 * t110;
	t111 = cos(qJ(6));
	t116 = t106 * t111;
	t107 = cos(t109);
	t115 = t107 * t110;
	t114 = t107 * t111;
	t104 = sin(t108);
	t113 = t104 * t116;
	t112 = t104 * t114;
	t103 = t105 * t111;
	t102 = t107 * t105;
	t101 = t106 * t105;
	t100 = t104 * t115;
	t99 = t104 * t117;
	t98 = t105 * t114 + t117;
	t97 = -t105 * t115 + t116;
	t96 = -t105 * t116 + t115;
	t95 = t105 * t117 + t114;
	t1 = [t96, 0, -t112, 0, -t112, t97; t98, 0, -t113, 0, -t113, -t95; 0, 0, t103, 0, t103, -t104 * t110; t95, 0, t100, 0, t100, -t98; t97, 0, t99, 0, t99, t96; 0, 0, -t118, 0, -t118, -t104 * t111; -t106 * t104, 0, t102, 0, t102, 0; t107 * t104, 0, t101, 0, t101, 0; 0, 0, t104, 0, t104, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end