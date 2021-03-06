% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PRPRR7
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5PRPRR7_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR7_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPRR7_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR7_jacobiR_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:01:10
	% EndTime: 2019-12-05 16:01:10
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:01:10
	% EndTime: 2019-12-05 16:01:10
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:01:10
	% EndTime: 2019-12-05 16:01:10
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (5->5), mult. (4->4), div. (0->0), fcn. (10->4), ass. (0->5)
	t9 = cos(qJ(2));
	t8 = sin(qJ(2));
	t7 = cos(pkin(8));
	t6 = sin(pkin(8));
	t1 = [0, -t7 * t8, 0, 0, 0; 0, -t6 * t8, 0, 0, 0; 0, t9, 0, 0, 0; 0, -t7 * t9, 0, 0, 0; 0, -t6 * t9, 0, 0, 0; 0, -t8, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:01:11
	% EndTime: 2019-12-05 16:01:11
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (4->4), div. (0->0), fcn. (10->4), ass. (0->5)
	t31 = cos(qJ(2));
	t30 = sin(qJ(2));
	t29 = cos(pkin(8));
	t28 = sin(pkin(8));
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, t29 * t30, 0, 0, 0; 0, t28 * t30, 0, 0, 0; 0, -t31, 0, 0, 0; 0, t29 * t31, 0, 0, 0; 0, t28 * t31, 0, 0, 0; 0, t30, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:01:11
	% EndTime: 2019-12-05 16:01:11
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (7->7), mult. (26->18), div. (0->0), fcn. (45->6), ass. (0->11)
	t42 = sin(qJ(4));
	t43 = sin(qJ(2));
	t49 = t43 * t42;
	t44 = cos(qJ(4));
	t48 = t43 * t44;
	t45 = cos(qJ(2));
	t47 = t45 * t42;
	t46 = t45 * t44;
	t41 = cos(pkin(8));
	t40 = sin(pkin(8));
	t1 = [0, t41 * t47, 0, -t40 * t42 + t41 * t48, 0; 0, t40 * t47, 0, t40 * t48 + t41 * t42, 0; 0, t49, 0, -t46, 0; 0, t41 * t46, 0, -t40 * t44 - t41 * t49, 0; 0, t40 * t46, 0, -t40 * t49 + t41 * t44, 0; 0, t48, 0, t47, 0; 0, -t41 * t43, 0, 0, 0; 0, -t40 * t43, 0, 0, 0; 0, t45, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:01:11
	% EndTime: 2019-12-05 16:01:11
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (38->9), mult. (40->18), div. (0->0), fcn. (69->6), ass. (0->16)
	t65 = qJ(4) + qJ(5);
	t63 = sin(t65);
	t68 = sin(qJ(2));
	t72 = t68 * t63;
	t64 = cos(t65);
	t71 = t68 * t64;
	t69 = cos(qJ(2));
	t62 = t69 * t63;
	t70 = t69 * t64;
	t67 = cos(pkin(8));
	t66 = sin(pkin(8));
	t61 = t67 * t64 - t66 * t72;
	t60 = t67 * t63 + t66 * t71;
	t59 = -t66 * t64 - t67 * t72;
	t58 = -t66 * t63 + t67 * t71;
	t1 = [0, t67 * t62, 0, t58, t58; 0, t66 * t62, 0, t60, t60; 0, t72, 0, -t70, -t70; 0, t67 * t70, 0, t59, t59; 0, t66 * t70, 0, t61, t61; 0, t71, 0, t62, t62; 0, -t67 * t68, 0, 0, 0; 0, -t66 * t68, 0, 0, 0; 0, t69, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end