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
% Datum: 2019-10-24 10:39
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
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
	% StartTime: 2019-10-24 10:39:19
	% EndTime: 2019-10-24 10:39:20
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:39:19
	% EndTime: 2019-10-24 10:39:20
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
	% StartTime: 2019-10-24 10:39:19
	% EndTime: 2019-10-24 10:39:19
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->3), mult. (4->4), div. (0->0), fcn. (10->4), ass. (0->5)
	t15 = cos(qJ(1));
	t14 = sin(qJ(1));
	t13 = cos(pkin(7));
	t12 = sin(pkin(7));
	t1 = [0, 0, 0, 0, 0; -t15 * t13, 0, 0, 0, 0; -t14 * t13, 0, 0, 0, 0; 0, 0, 0, 0, 0; t15 * t12, 0, 0, 0, 0; t14 * t12, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t14, 0, 0, 0, 0; t15, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:39:20
	% EndTime: 2019-10-24 10:39:20
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (6->6), mult. (14->10), div. (0->0), fcn. (24->6), ass. (0->11)
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
	t1 = [0, 0, 0, 0, 0; -t25 * t28 - t31, 0, 0, 0, 0; -t25 * t30 + t29, 0, 0, 0, 0; 0, 0, 0, 0, 0; t25 * t29 - t30, 0, 0, 0, 0; t25 * t31 + t28, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t27 * t23, 0, 0, 0, 0; -t26 * t23, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:39:20
	% EndTime: 2019-10-24 10:39:20
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (10->8), mult. (30->18), div. (0->0), fcn. (46->8), ass. (0->17)
	t39 = sin(pkin(7));
	t43 = sin(qJ(1));
	t50 = t39 * t43;
	t44 = cos(qJ(1));
	t49 = t39 * t44;
	t38 = sin(pkin(8));
	t48 = t43 * t38;
	t41 = cos(pkin(8));
	t47 = t43 * t41;
	t46 = t44 * t38;
	t45 = t44 * t41;
	t42 = cos(pkin(7));
	t40 = cos(pkin(9));
	t37 = sin(pkin(9));
	t36 = -t42 * t45 - t48;
	t35 = -t42 * t47 + t46;
	t1 = [0, 0, 0, 0, 0; t36 * t40 - t37 * t49, 0, 0, 0, 0; t35 * t40 - t37 * t50, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t36 * t37 - t40 * t49, 0, 0, 0, 0; -t35 * t37 - t40 * t50, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t42 * t46 + t47, 0, 0, 0, 0; -t42 * t48 - t45, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:39:20
	% EndTime: 2019-10-24 10:39:20
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (40->17), mult. (112->35), div. (0->0), fcn. (164->10), ass. (0->28)
	t100 = cos(qJ(1));
	t93 = sin(pkin(7));
	t102 = t100 * t93;
	t96 = cos(pkin(7));
	t101 = t100 * t96;
	t92 = sin(pkin(8));
	t98 = sin(qJ(1));
	t104 = t98 * t92;
	t95 = cos(pkin(8));
	t89 = t95 * t101 + t104;
	t91 = sin(pkin(9));
	t94 = cos(pkin(9));
	t83 = t91 * t102 + t89 * t94;
	t103 = t98 * t95;
	t88 = t92 * t101 - t103;
	t97 = sin(qJ(5));
	t99 = cos(qJ(5));
	t110 = t83 * t97 - t88 * t99;
	t109 = -t83 * t99 - t88 * t97;
	t106 = t92 * t93;
	t105 = t93 * t98;
	t87 = t100 * t92 - t96 * t103;
	t86 = -t100 * t95 - t96 * t104;
	t85 = t93 * t95 * t94 - t96 * t91;
	t82 = -t91 * t105 + t87 * t94;
	t81 = t82 * t99 + t86 * t97;
	t80 = -t82 * t97 + t86 * t99;
	t1 = [0, 0, 0, 0, t99 * t106 - t85 * t97; t109, 0, 0, 0, t80; t81, 0, 0, 0, -t110; 0, 0, 0, 0, -t97 * t106 - t85 * t99; t110, 0, 0, 0, -t81; t80, 0, 0, 0, t109; 0, 0, 0, 0, 0; t94 * t102 - t89 * t91, 0, 0, 0, 0; t94 * t105 + t87 * t91, 0, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end