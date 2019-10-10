% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRRRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:50
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRRP4_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP4_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRP4_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP4_jacobiR_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:50:39
	% EndTime: 2019-10-10 01:50:39
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:50:39
	% EndTime: 2019-10-10 01:50:39
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
	% StartTime: 2019-10-10 01:50:39
	% EndTime: 2019-10-10 01:50:39
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (10->4), ass. (0->5)
	t10 = cos(qJ(1));
	t9 = sin(qJ(1));
	t8 = cos(pkin(10));
	t7 = sin(pkin(10));
	t1 = [-t9 * t8, 0, 0, 0, 0, 0; t10 * t8, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t9 * t7, 0, 0, 0, 0, 0; -t10 * t7, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t10, 0, 0, 0, 0, 0; t9, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:50:39
	% EndTime: 2019-10-10 01:50:39
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (17->8), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->10)
	t14 = pkin(10) + qJ(3);
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
	% StartTime: 2019-10-10 01:50:39
	% EndTime: 2019-10-10 01:50:39
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (44->13), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->10)
	t23 = pkin(10) + qJ(3) + qJ(4);
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
	% StartTime: 2019-10-10 01:50:39
	% EndTime: 2019-10-10 01:50:39
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (77->16), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->24)
	t101 = pkin(10) + qJ(3) + qJ(4);
	t100 = cos(t101);
	t102 = sin(qJ(5));
	t112 = t100 * t102;
	t103 = sin(qJ(1));
	t111 = t103 * t102;
	t104 = cos(qJ(5));
	t110 = t103 * t104;
	t105 = cos(qJ(1));
	t109 = t105 * t102;
	t108 = t105 * t104;
	t99 = sin(t101);
	t107 = t99 * t110;
	t106 = t99 * t108;
	t98 = t105 * t100;
	t97 = t100 * t104;
	t96 = t103 * t100;
	t95 = t99 * t109;
	t94 = t99 * t111;
	t93 = t100 * t108 + t111;
	t92 = -t100 * t109 + t110;
	t91 = -t100 * t110 + t109;
	t90 = t100 * t111 + t108;
	t1 = [t91, 0, -t106, -t106, t92, 0; t93, 0, -t107, -t107, -t90, 0; 0, 0, t97, t97, -t99 * t102, 0; t90, 0, t95, t95, -t93, 0; t92, 0, t94, t94, t91, 0; 0, 0, -t112, -t112, -t99 * t104, 0; -t103 * t99, 0, t98, t98, 0, 0; t105 * t99, 0, t96, t96, 0, 0; 0, 0, t99, t99, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:50:39
	% EndTime: 2019-10-10 01:50:39
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (77->16), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->24)
	t100 = sin(qJ(5));
	t99 = pkin(10) + qJ(3) + qJ(4);
	t98 = cos(t99);
	t110 = t98 * t100;
	t101 = sin(qJ(1));
	t109 = t101 * t100;
	t102 = cos(qJ(5));
	t108 = t101 * t102;
	t103 = cos(qJ(1));
	t107 = t103 * t100;
	t106 = t103 * t102;
	t97 = sin(t99);
	t105 = t97 * t108;
	t104 = t97 * t106;
	t96 = t103 * t98;
	t95 = t98 * t102;
	t94 = t101 * t98;
	t93 = t97 * t107;
	t92 = t97 * t109;
	t91 = t98 * t106 + t109;
	t90 = -t98 * t107 + t108;
	t89 = -t98 * t108 + t107;
	t88 = t98 * t109 + t106;
	t1 = [t89, 0, -t104, -t104, t90, 0; t91, 0, -t105, -t105, -t88, 0; 0, 0, t95, t95, -t97 * t100, 0; t88, 0, t93, t93, -t91, 0; t90, 0, t92, t92, t89, 0; 0, 0, -t110, -t110, -t97 * t102, 0; -t101 * t97, 0, t96, t96, 0, 0; t103 * t97, 0, t94, t94, 0, 0; 0, 0, t97, t97, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end