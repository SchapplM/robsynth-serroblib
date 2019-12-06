% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PPRRP3
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
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:11
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5PPRRP3_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP3_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRRP3_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP3_jacobiR_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:11:38
	% EndTime: 2019-12-05 15:11:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:11:38
	% EndTime: 2019-12-05 15:11:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:11:38
	% EndTime: 2019-12-05 15:11:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:11:38
	% EndTime: 2019-12-05 15:11:38
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (6->6), mult. (14->10), div. (0->0), fcn. (24->6), ass. (0->11)
	t22 = sin(pkin(7));
	t25 = sin(qJ(3));
	t30 = t22 * t25;
	t26 = cos(qJ(3));
	t29 = t22 * t26;
	t24 = cos(pkin(7));
	t28 = t24 * t25;
	t27 = t24 * t26;
	t23 = cos(pkin(8));
	t21 = sin(pkin(8));
	t1 = [0, 0, -t23 * t28 + t29, 0, 0; 0, 0, -t23 * t30 - t27, 0, 0; 0, 0, -t21 * t25, 0, 0; 0, 0, -t23 * t27 - t30, 0, 0; 0, 0, -t23 * t29 + t28, 0, 0; 0, 0, -t21 * t26, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:11:38
	% EndTime: 2019-12-05 15:11:38
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (19->13), mult. (57->29), div. (0->0), fcn. (88->8), ass. (0->20)
	t68 = sin(pkin(8));
	t72 = sin(qJ(4));
	t82 = t68 * t72;
	t74 = cos(qJ(4));
	t81 = t68 * t74;
	t75 = cos(qJ(3));
	t80 = t68 * t75;
	t69 = sin(pkin(7));
	t73 = sin(qJ(3));
	t79 = t69 * t73;
	t78 = t69 * t75;
	t71 = cos(pkin(7));
	t77 = t71 * t73;
	t76 = t71 * t75;
	t70 = cos(pkin(8));
	t67 = t70 * t76 + t79;
	t66 = -t70 * t77 + t78;
	t65 = t70 * t78 - t77;
	t64 = -t70 * t79 - t76;
	t1 = [0, 0, t66 * t74, -t67 * t72 + t71 * t81, 0; 0, 0, t64 * t74, -t65 * t72 + t69 * t81, 0; 0, 0, -t73 * t81, -t70 * t74 - t72 * t80, 0; 0, 0, -t66 * t72, -t67 * t74 - t71 * t82, 0; 0, 0, -t64 * t72, -t65 * t74 - t69 * t82, 0; 0, 0, t73 * t82, t70 * t72 - t74 * t80, 0; 0, 0, t67, 0, 0; 0, 0, t65, 0, 0; 0, 0, t80, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:11:38
	% EndTime: 2019-12-05 15:11:38
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (18->12), mult. (57->29), div. (0->0), fcn. (88->8), ass. (0->20)
	t90 = sin(pkin(8));
	t94 = sin(qJ(4));
	t104 = t90 * t94;
	t96 = cos(qJ(4));
	t103 = t90 * t96;
	t97 = cos(qJ(3));
	t102 = t90 * t97;
	t91 = sin(pkin(7));
	t95 = sin(qJ(3));
	t101 = t91 * t95;
	t100 = t91 * t97;
	t93 = cos(pkin(7));
	t99 = t93 * t95;
	t98 = t93 * t97;
	t92 = cos(pkin(8));
	t89 = t92 * t98 + t101;
	t88 = -t92 * t99 + t100;
	t87 = t92 * t100 - t99;
	t86 = -t92 * t101 - t98;
	t1 = [0, 0, t88 * t96, t93 * t103 - t89 * t94, 0; 0, 0, t86 * t96, t91 * t103 - t87 * t94, 0; 0, 0, -t95 * t103, -t94 * t102 - t92 * t96, 0; 0, 0, t89, 0, 0; 0, 0, t87, 0, 0; 0, 0, t102, 0, 0; 0, 0, t88 * t94, t93 * t104 + t89 * t96, 0; 0, 0, t86 * t94, t91 * t104 + t87 * t96, 0; 0, 0, -t95 * t104, t96 * t102 - t92 * t94, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end