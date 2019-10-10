% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRRPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:33
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRPR7_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR7_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR7_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR7_jacobiR_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:33:56
	% EndTime: 2019-10-10 01:33:56
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:33:56
	% EndTime: 2019-10-10 01:33:56
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
	% StartTime: 2019-10-10 01:33:56
	% EndTime: 2019-10-10 01:33:56
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
	% StartTime: 2019-10-10 01:33:56
	% EndTime: 2019-10-10 01:33:56
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
	% StartTime: 2019-10-10 01:33:56
	% EndTime: 2019-10-10 01:33:56
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (25->10), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->10)
	t24 = qJ(3) + qJ(4);
	t22 = sin(t24);
	t25 = sin(qJ(1));
	t28 = t25 * t22;
	t23 = cos(t24);
	t26 = cos(qJ(1));
	t27 = t26 * t23;
	t21 = t26 * t22;
	t20 = t25 * t23;
	t1 = [t21, 0, t20, t20, 0, 0; t28, 0, -t27, -t27, 0, 0; 0, 0, -t22, -t22, 0, 0; t27, 0, -t28, -t28, 0, 0; t20, 0, t21, t21, 0, 0; 0, 0, -t23, -t23, 0, 0; -t25, 0, 0, 0, 0, 0; t26, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:33:56
	% EndTime: 2019-10-10 01:33:56
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (41->10), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->10)
	t31 = qJ(3) + qJ(4) + pkin(10);
	t29 = sin(t31);
	t32 = sin(qJ(1));
	t35 = t32 * t29;
	t30 = cos(t31);
	t33 = cos(qJ(1));
	t34 = t33 * t30;
	t28 = t33 * t29;
	t27 = t32 * t30;
	t1 = [t28, 0, t27, t27, 0, 0; t35, 0, -t34, -t34, 0, 0; 0, 0, -t29, -t29, 0, 0; t34, 0, -t35, -t35, 0, 0; t27, 0, t28, t28, 0, 0; 0, 0, -t30, -t30, 0, 0; -t32, 0, 0, 0, 0, 0; t33, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:33:56
	% EndTime: 2019-10-10 01:33:56
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (80->19), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->24)
	t103 = cos(qJ(1));
	t99 = qJ(3) + qJ(4) + pkin(10);
	t97 = sin(t99);
	t111 = t103 * t97;
	t102 = cos(qJ(6));
	t110 = t97 * t102;
	t100 = sin(qJ(6));
	t101 = sin(qJ(1));
	t109 = t101 * t100;
	t108 = t101 * t102;
	t107 = t103 * t100;
	t106 = t103 * t102;
	t98 = cos(t99);
	t105 = t98 * t109;
	t104 = t98 * t106;
	t96 = t101 * t97;
	t95 = t97 * t100;
	t94 = t98 * t107;
	t93 = t98 * t108;
	t92 = t97 * t106 - t109;
	t91 = t97 * t107 + t108;
	t90 = t97 * t108 + t107;
	t89 = -t97 * t109 + t106;
	t1 = [t92, 0, t93, t93, 0, t89; t90, 0, -t104, -t104, 0, t91; 0, 0, -t110, -t110, 0, -t98 * t100; -t91, 0, -t105, -t105, 0, -t90; t89, 0, t94, t94, 0, t92; 0, 0, t95, t95, 0, -t98 * t102; -t103 * t98, 0, t96, t96, 0, 0; -t101 * t98, 0, -t111, -t111, 0, 0; 0, 0, t98, t98, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end