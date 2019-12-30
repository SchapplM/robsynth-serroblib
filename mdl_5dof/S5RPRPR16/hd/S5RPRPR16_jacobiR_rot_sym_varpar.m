% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPRPR16
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 17:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RPRPR16_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR16_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPR16_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPR16_jacobiR_rot_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:13:28
	% EndTime: 2019-12-29 17:13:28
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:13:28
	% EndTime: 2019-12-29 17:13:28
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8, 0, 0, 0, 0; t9, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0; -t8, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:13:33
	% EndTime: 2019-12-29 17:13:33
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t6 = cos(qJ(1));
	t5 = sin(qJ(1));
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; t5, 0, 0, 0, 0; -t6, 0, 0, 0, 0; 0, 0, 0, 0, 0; t6, 0, 0, 0, 0; t5, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:13:28
	% EndTime: 2019-12-29 17:13:28
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (5->5), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->9)
	t11 = sin(qJ(3));
	t12 = sin(qJ(1));
	t16 = t12 * t11;
	t13 = cos(qJ(3));
	t14 = cos(qJ(1));
	t15 = t14 * t13;
	t10 = t14 * t11;
	t9 = t12 * t13;
	t1 = [t10, 0, t9, 0, 0; t16, 0, -t15, 0, 0; 0, 0, -t11, 0, 0; t15, 0, -t16, 0, 0; t9, 0, t10, 0, 0; 0, 0, -t13, 0, 0; -t12, 0, 0, 0, 0; t14, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:13:28
	% EndTime: 2019-12-29 17:13:28
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (7->7), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->9)
	t39 = sin(qJ(3));
	t40 = sin(qJ(1));
	t46 = t40 * t39;
	t41 = cos(qJ(3));
	t45 = t40 * t41;
	t42 = cos(qJ(1));
	t44 = t42 * t39;
	t43 = t42 * t41;
	t1 = [-t40, 0, 0, 0, 0; t42, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t44, 0, -t45, 0, 0; -t46, 0, t43, 0, 0; 0, 0, t39, 0, 0; -t43, 0, t46, 0, 0; -t45, 0, -t44, 0, 0; 0, 0, t41, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:13:28
	% EndTime: 2019-12-29 17:13:28
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (13->11), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->18)
	t61 = sin(qJ(3));
	t62 = sin(qJ(1));
	t72 = t62 * t61;
	t63 = cos(qJ(5));
	t71 = t62 * t63;
	t60 = sin(qJ(5));
	t64 = cos(qJ(3));
	t70 = t64 * t60;
	t69 = t64 * t63;
	t65 = cos(qJ(1));
	t68 = t65 * t61;
	t67 = t65 * t63;
	t66 = t65 * t64;
	t59 = -t62 * t70 + t67;
	t58 = -t65 * t60 - t62 * t69;
	t57 = -t60 * t66 - t71;
	t56 = t62 * t60 - t63 * t66;
	t1 = [t57, 0, t60 * t72, 0, t58; t59, 0, -t60 * t68, 0, -t56; 0, 0, t70, 0, t61 * t63; t56, 0, t61 * t71, 0, -t59; t58, 0, -t61 * t67, 0, t57; 0, 0, t69, 0, -t61 * t60; t68, 0, t62 * t64, 0, 0; t72, 0, -t66, 0, 0; 0, 0, -t61, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end