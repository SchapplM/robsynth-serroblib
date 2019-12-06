% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRPRR6
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RRPRR6_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRR6_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_jacobiR_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:36:54
	% EndTime: 2019-12-05 18:36:54
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:36:54
	% EndTime: 2019-12-05 18:36:55
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t4 = cos(qJ(1));
	t3 = sin(qJ(1));
	t1 = [0, 0, 0, 0, 0; -t4, 0, 0, 0, 0; -t3, 0, 0, 0, 0; 0, 0, 0, 0, 0; t3, 0, 0, 0, 0; -t4, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:36:54
	% EndTime: 2019-12-05 18:36:54
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (14->7), mult. (0->0), div. (0->0), fcn. (8->2), ass. (0->4)
	t9 = qJ(1) + qJ(2);
	t8 = cos(t9);
	t7 = sin(t9);
	t1 = [0, 0, 0, 0, 0; -t8, -t8, 0, 0, 0; -t7, -t7, 0, 0, 0; 0, 0, 0, 0, 0; t7, t7, 0, 0, 0; -t8, -t8, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:36:55
	% EndTime: 2019-12-05 18:36:55
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (18->7), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->10)
	t27 = qJ(1) + qJ(2);
	t25 = sin(t27);
	t29 = cos(pkin(9));
	t31 = t25 * t29;
	t26 = cos(t27);
	t30 = t26 * t29;
	t28 = sin(pkin(9));
	t24 = t26 * t28;
	t23 = t25 * t28;
	t1 = [0, 0, 0, 0, 0; -t30, -t30, 0, 0, 0; -t31, -t31, 0, 0, 0; 0, 0, 0, 0, 0; t24, t24, 0, 0, 0; t23, t23, 0, 0, 0; 0, 0, 0, 0, 0; -t25, -t25, 0, 0, 0; t26, t26, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:36:55
	% EndTime: 2019-12-05 18:36:55
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (46->14), mult. (42->14), div. (0->0), fcn. (72->6), ass. (0->16)
	t51 = qJ(1) + qJ(2);
	t49 = sin(t51);
	t52 = sin(pkin(9));
	t59 = t49 * t52;
	t50 = cos(t51);
	t58 = t50 * t52;
	t53 = cos(pkin(9));
	t54 = sin(qJ(4));
	t57 = t53 * t54;
	t55 = cos(qJ(4));
	t56 = t53 * t55;
	t47 = -t49 * t54 - t50 * t56;
	t46 = -t49 * t55 + t50 * t57;
	t45 = t49 * t56 - t50 * t54;
	t44 = t49 * t57 + t50 * t55;
	t1 = [0, 0, 0, -t52 * t54, 0; t47, t47, 0, t44, 0; -t45, -t45, 0, -t46, 0; 0, 0, 0, -t52 * t55, 0; t46, t46, 0, t45, 0; t44, t44, 0, t47, 0; 0, 0, 0, 0, 0; -t58, -t58, 0, 0, 0; -t59, -t59, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:36:55
	% EndTime: 2019-12-05 18:36:55
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (96->18), mult. (56->14), div. (0->0), fcn. (96->6), ass. (0->19)
	t81 = qJ(1) + qJ(2);
	t77 = sin(t81);
	t82 = sin(pkin(9));
	t89 = t77 * t82;
	t83 = cos(pkin(9));
	t88 = t77 * t83;
	t79 = cos(t81);
	t87 = t79 * t82;
	t86 = t79 * t83;
	t80 = qJ(4) + qJ(5);
	t76 = sin(t80);
	t85 = t82 * t76;
	t78 = cos(t80);
	t84 = t82 * t78;
	t73 = -t77 * t76 - t78 * t86;
	t72 = t76 * t86 - t77 * t78;
	t71 = -t79 * t76 + t78 * t88;
	t70 = t76 * t88 + t79 * t78;
	t1 = [0, 0, 0, -t85, -t85; t73, t73, 0, t70, t70; -t71, -t71, 0, -t72, -t72; 0, 0, 0, -t84, -t84; t72, t72, 0, t71, t71; t70, t70, 0, t73, t73; 0, 0, 0, 0, 0; -t87, -t87, 0, 0, 0; -t89, -t89, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end