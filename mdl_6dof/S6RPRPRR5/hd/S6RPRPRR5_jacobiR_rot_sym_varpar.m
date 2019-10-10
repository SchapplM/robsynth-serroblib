% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRPRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:53
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRPRR5_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR5_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR5_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR5_jacobiR_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:53:10
	% EndTime: 2019-10-10 00:53:10
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:53:10
	% EndTime: 2019-10-10 00:53:10
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
	% StartTime: 2019-10-10 00:53:10
	% EndTime: 2019-10-10 00:53:10
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
	% StartTime: 2019-10-10 00:53:10
	% EndTime: 2019-10-10 00:53:10
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
	% StartTime: 2019-10-10 00:53:10
	% EndTime: 2019-10-10 00:53:10
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (14->5), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->10)
	t44 = pkin(10) + qJ(3);
	t42 = sin(t44);
	t45 = sin(qJ(1));
	t49 = t45 * t42;
	t43 = cos(t44);
	t48 = t45 * t43;
	t46 = cos(qJ(1));
	t47 = t46 * t42;
	t41 = t46 * t43;
	t1 = [-t48, 0, -t47, 0, 0, 0; t41, 0, -t49, 0, 0, 0; 0, 0, t43, 0, 0, 0; t46, 0, 0, 0, 0, 0; t45, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t49, 0, t41, 0, 0, 0; t47, 0, t48, 0, 0, 0; 0, 0, t42, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:53:10
	% EndTime: 2019-10-10 00:53:10
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (50->13), mult. (56->8), div. (0->0), fcn. (90->6), ass. (0->14)
	t41 = pkin(10) + qJ(3);
	t39 = sin(t41);
	t40 = cos(t41);
	t42 = sin(qJ(5));
	t44 = cos(qJ(5));
	t49 = -t39 * t44 + t40 * t42;
	t46 = t39 * t42 + t40 * t44;
	t45 = cos(qJ(1));
	t43 = sin(qJ(1));
	t35 = t46 * t45;
	t34 = t49 * t45;
	t33 = t46 * t43;
	t32 = t49 * t43;
	t1 = [-t33, 0, t34, 0, -t34, 0; t35, 0, t32, 0, -t32, 0; 0, 0, t46, 0, -t46, 0; t32, 0, t35, 0, -t35, 0; -t34, 0, t33, 0, -t33, 0; 0, 0, -t49, 0, t49, 0; -t45, 0, 0, 0, 0, 0; -t43, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:53:10
	% EndTime: 2019-10-10 00:53:10
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (106->21), mult. (134->24), div. (0->0), fcn. (202->8), ass. (0->26)
	t120 = pkin(10) + qJ(3);
	t119 = cos(t120);
	t122 = sin(qJ(5));
	t138 = sin(t120);
	t139 = cos(qJ(5));
	t114 = -t119 * t122 + t138 * t139;
	t113 = t119 * t139 + t138 * t122;
	t123 = sin(qJ(1));
	t109 = t114 * t123;
	t121 = sin(qJ(6));
	t137 = t109 * t121;
	t124 = cos(qJ(6));
	t136 = t109 * t124;
	t125 = cos(qJ(1));
	t112 = t114 * t125;
	t135 = t112 * t121;
	t134 = t112 * t124;
	t133 = t113 * t121;
	t132 = t113 * t124;
	t110 = t113 * t123;
	t127 = -t110 * t124 - t125 * t121;
	t126 = t110 * t121 - t125 * t124;
	t111 = t113 * t125;
	t108 = t111 * t124 - t123 * t121;
	t107 = -t111 * t121 - t123 * t124;
	t1 = [t127, 0, -t134, 0, t134, t107; t108, 0, -t136, 0, t136, -t126; 0, 0, t132, 0, -t132, -t114 * t121; t126, 0, t135, 0, -t135, -t108; t107, 0, t137, 0, -t137, t127; 0, 0, -t133, 0, t133, -t114 * t124; t109, 0, -t111, 0, t111, 0; -t112, 0, -t110, 0, t110, 0; 0, 0, -t114, 0, t114, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end