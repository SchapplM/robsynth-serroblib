% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PRPPR5
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta4]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5PRPPR5_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR5_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPPR5_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR5_jacobiR_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:38:50
	% EndTime: 2019-12-31 17:38:50
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:38:50
	% EndTime: 2019-12-31 17:38:50
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:38:50
	% EndTime: 2019-12-31 17:38:50
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (5->5), mult. (4->4), div. (0->0), fcn. (10->4), ass. (0->5)
	t9 = cos(qJ(2));
	t8 = sin(qJ(2));
	t7 = cos(pkin(7));
	t6 = sin(pkin(7));
	t1 = [0, -t7 * t8, 0, 0, 0; 0, -t6 * t8, 0, 0, 0; 0, t9, 0, 0, 0; 0, -t7 * t9, 0, 0, 0; 0, -t6 * t9, 0, 0, 0; 0, -t8, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:38:50
	% EndTime: 2019-12-31 17:38:50
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (10->4), ass. (0->5)
	t30 = cos(qJ(2));
	t29 = sin(qJ(2));
	t28 = cos(pkin(7));
	t27 = sin(pkin(7));
	t1 = [0, -t28 * t29, 0, 0, 0; 0, -t27 * t29, 0, 0, 0; 0, t30, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, t28 * t30, 0, 0, 0; 0, t27 * t30, 0, 0, 0; 0, t29, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:38:50
	% EndTime: 2019-12-31 17:38:50
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (6->3), mult. (20->8), div. (0->0), fcn. (32->6), ass. (0->9)
	t10 = sin(pkin(8));
	t12 = cos(pkin(8));
	t14 = sin(qJ(2));
	t15 = cos(qJ(2));
	t17 = t15 * t10 - t14 * t12;
	t16 = t14 * t10 + t15 * t12;
	t13 = cos(pkin(7));
	t11 = sin(pkin(7));
	t1 = [0, t17 * t13, 0, 0, 0; 0, t17 * t11, 0, 0, 0; 0, t16, 0, 0, 0; 0, t16 * t13, 0, 0, 0; 0, t16 * t11, 0, 0, 0; 0, -t17, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:38:50
	% EndTime: 2019-12-31 17:38:50
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (24->11), mult. (66->24), div. (0->0), fcn. (100->8), ass. (0->15)
	t82 = sin(pkin(8));
	t84 = cos(pkin(8));
	t87 = sin(qJ(2));
	t89 = cos(qJ(2));
	t76 = t87 * t82 + t89 * t84;
	t90 = t89 * t82 - t87 * t84;
	t88 = cos(qJ(5));
	t86 = sin(qJ(5));
	t85 = cos(pkin(7));
	t83 = sin(pkin(7));
	t75 = t90 * t85;
	t74 = t76 * t85;
	t73 = t90 * t83;
	t72 = t76 * t83;
	t1 = [0, t75 * t88, 0, 0, -t74 * t86 - t83 * t88; 0, t73 * t88, 0, 0, -t72 * t86 + t85 * t88; 0, t76 * t88, 0, 0, t90 * t86; 0, -t75 * t86, 0, 0, -t74 * t88 + t83 * t86; 0, -t73 * t86, 0, 0, -t72 * t88 - t85 * t86; 0, -t76 * t86, 0, 0, t90 * t88; 0, -t74, 0, 0, 0; 0, -t72, 0, 0, 0; 0, t90, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end