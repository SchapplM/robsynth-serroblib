% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPPPR2
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RPPPR2_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPPR2_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_jacobiR_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:23:58
	% EndTime: 2020-01-03 11:23:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:23:57
	% EndTime: 2020-01-03 11:23:57
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t4 = cos(qJ(1));
	t3 = sin(qJ(1));
	t1 = [0, 0, 0, 0, 0; t4, 0, 0, 0, 0; t3, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t3, 0, 0, 0, 0; t4, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:23:57
	% EndTime: 2020-01-03 11:23:57
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->3), mult. (4->4), div. (0->0), fcn. (10->4), ass. (0->5)
	t15 = cos(qJ(1));
	t14 = sin(qJ(1));
	t13 = cos(pkin(7));
	t12 = sin(pkin(7));
	t1 = [0, 0, 0, 0, 0; t15 * t13, 0, 0, 0, 0; t14 * t13, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t15 * t12, 0, 0, 0, 0; -t14 * t12, 0, 0, 0, 0; 0, 0, 0, 0, 0; t14, 0, 0, 0, 0; -t15, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:23:58
	% EndTime: 2020-01-03 11:23:58
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (4->4), mult. (14->10), div. (0->0), fcn. (24->6), ass. (0->11)
	t22 = sin(pkin(8));
	t26 = sin(qJ(1));
	t31 = t26 * t22;
	t24 = cos(pkin(8));
	t30 = t26 * t24;
	t27 = cos(qJ(1));
	t29 = t27 * t22;
	t28 = t27 * t24;
	t25 = cos(pkin(7));
	t23 = sin(pkin(7));
	t1 = [0, 0, 0, 0, 0; t25 * t28 + t31, 0, 0, 0, 0; t25 * t30 - t29, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t25 * t29 + t30, 0, 0, 0, 0; -t25 * t31 - t28, 0, 0, 0, 0; 0, 0, 0, 0, 0; t27 * t23, 0, 0, 0, 0; t26 * t23, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:23:58
	% EndTime: 2020-01-03 11:23:58
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (10->8), mult. (30->18), div. (0->0), fcn. (46->8), ass. (0->17)
	t44 = sin(pkin(7));
	t48 = sin(qJ(1));
	t55 = t44 * t48;
	t49 = cos(qJ(1));
	t54 = t44 * t49;
	t43 = sin(pkin(8));
	t53 = t48 * t43;
	t46 = cos(pkin(8));
	t52 = t48 * t46;
	t51 = t49 * t43;
	t50 = t49 * t46;
	t47 = cos(pkin(7));
	t45 = cos(pkin(9));
	t42 = sin(pkin(9));
	t41 = t47 * t50 + t53;
	t40 = t47 * t52 - t51;
	t1 = [0, 0, 0, 0, 0; t41 * t45 + t42 * t54, 0, 0, 0, 0; t40 * t45 + t42 * t55, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t41 * t42 + t45 * t54, 0, 0, 0, 0; -t40 * t42 + t45 * t55, 0, 0, 0, 0; 0, 0, 0, 0, 0; t47 * t51 - t52, 0, 0, 0, 0; t47 * t53 + t50, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:23:58
	% EndTime: 2020-01-03 11:23:58
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (40->17), mult. (112->36), div. (0->0), fcn. (164->10), ass. (0->27)
	t101 = cos(qJ(5));
	t102 = cos(qJ(1));
	t95 = sin(pkin(7));
	t105 = t102 * t95;
	t100 = sin(qJ(1));
	t98 = cos(pkin(7));
	t104 = t102 * t98;
	t94 = sin(pkin(8));
	t97 = cos(pkin(8));
	t91 = t100 * t94 + t97 * t104;
	t93 = sin(pkin(9));
	t96 = cos(pkin(9));
	t85 = t93 * t105 + t91 * t96;
	t90 = -t100 * t97 + t94 * t104;
	t99 = sin(qJ(5));
	t111 = -t90 * t101 + t85 * t99;
	t110 = t85 * t101 + t90 * t99;
	t108 = t94 * t95;
	t107 = t100 * t95;
	t106 = t100 * t98;
	t89 = -t102 * t94 + t97 * t106;
	t88 = t102 * t97 + t94 * t106;
	t87 = t95 * t97 * t96 - t98 * t93;
	t84 = t93 * t107 + t89 * t96;
	t83 = t84 * t101 + t88 * t99;
	t82 = t88 * t101 - t84 * t99;
	t1 = [0, 0, 0, 0, t101 * t108 - t87 * t99; t110, 0, 0, 0, t82; t83, 0, 0, 0, t111; 0, 0, 0, 0, -t87 * t101 - t99 * t108; -t111, 0, 0, 0, -t83; t82, 0, 0, 0, t110; 0, 0, 0, 0, 0; -t96 * t105 + t91 * t93, 0, 0, 0, 0; -t96 * t107 + t89 * t93, 0, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end