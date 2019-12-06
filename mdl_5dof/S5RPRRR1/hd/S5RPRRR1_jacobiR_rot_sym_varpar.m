% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPRRR1
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
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RPRRR1_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRR1_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_jacobiR_rot_sym_varpar: pkin has to be [1x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:10:18
	% EndTime: 2019-12-05 18:10:18
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:10:18
	% EndTime: 2019-12-05 18:10:18
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8, 0, 0, 0, 0; t9, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0; -t8, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:10:18
	% EndTime: 2019-12-05 18:10:18
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t5 = cos(qJ(1));
	t4 = sin(qJ(1));
	t1 = [-t4, 0, 0, 0, 0; t5, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; t5, 0, 0, 0, 0; t4, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:10:18
	% EndTime: 2019-12-05 18:10:18
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (7->7), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->9)
	t8 = sin(qJ(3));
	t9 = sin(qJ(1));
	t15 = t9 * t8;
	t11 = cos(qJ(1));
	t14 = t11 * t8;
	t10 = cos(qJ(3));
	t13 = t9 * t10;
	t12 = t11 * t10;
	t1 = [-t13, 0, -t14, 0, 0; t12, 0, -t15, 0, 0; 0, 0, t10, 0, 0; t15, 0, -t12, 0, 0; -t14, 0, -t13, 0, 0; 0, 0, -t8, 0, 0; t11, 0, 0, 0, 0; t9, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:10:18
	% EndTime: 2019-12-05 18:10:18
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (14->12), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->18)
	t57 = sin(qJ(3));
	t58 = sin(qJ(1));
	t68 = t58 * t57;
	t59 = cos(qJ(4));
	t67 = t58 * t59;
	t56 = sin(qJ(4));
	t60 = cos(qJ(3));
	t66 = t60 * t56;
	t65 = t60 * t59;
	t61 = cos(qJ(1));
	t64 = t61 * t57;
	t63 = t61 * t59;
	t62 = t61 * t60;
	t55 = t58 * t56 + t59 * t62;
	t54 = -t56 * t62 + t67;
	t53 = t61 * t56 - t58 * t65;
	t52 = t58 * t66 + t63;
	t1 = [t53, 0, -t57 * t63, t54, 0; t55, 0, -t57 * t67, -t52, 0; 0, 0, t65, -t57 * t56, 0; t52, 0, t56 * t64, -t55, 0; t54, 0, t56 * t68, t53, 0; 0, 0, -t66, -t57 * t59, 0; -t68, 0, t62, 0, 0; t64, 0, t58 * t60, 0, 0; 0, 0, t57, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:10:18
	% EndTime: 2019-12-05 18:10:18
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (37->20), mult. (118->41), div. (0->0), fcn. (180->8), ass. (0->27)
	t94 = sin(qJ(5));
	t96 = sin(qJ(3));
	t113 = t96 * t94;
	t98 = cos(qJ(5));
	t112 = t96 * t98;
	t99 = cos(qJ(4));
	t111 = t96 * t99;
	t95 = sin(qJ(4));
	t97 = sin(qJ(1));
	t110 = t97 * t95;
	t100 = cos(qJ(3));
	t109 = t100 * t95;
	t108 = t100 * t99;
	t101 = cos(qJ(1));
	t107 = t101 * t96;
	t106 = t100 * t101;
	t89 = -t101 * t95 + t97 * t108;
	t105 = t97 * t112 - t89 * t94;
	t104 = -t97 * t113 - t89 * t98;
	t103 = t100 * t94 - t98 * t111;
	t102 = t100 * t98 + t94 * t111;
	t91 = t99 * t106 + t110;
	t90 = t95 * t106 - t97 * t99;
	t88 = -t101 * t99 - t97 * t109;
	t87 = t94 * t107 + t91 * t98;
	t86 = t98 * t107 - t91 * t94;
	t1 = [t104, 0, t103 * t101, -t90 * t98, t86; t87, 0, t103 * t97, t88 * t98, t105; 0, 0, t98 * t108 + t113, -t95 * t112, -t102; -t105, 0, t102 * t101, t90 * t94, -t87; t86, 0, t102 * t97, -t88 * t94, t104; 0, 0, -t94 * t108 + t112, t95 * t113, t103; t88, 0, -t95 * t107, t91, 0; t90, 0, -t96 * t110, t89, 0; 0, 0, t109, t111, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end