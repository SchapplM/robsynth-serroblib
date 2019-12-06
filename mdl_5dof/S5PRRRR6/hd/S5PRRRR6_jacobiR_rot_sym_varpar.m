% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PRRRR6
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5PRRRR6_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR6_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRR6_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR6_jacobiR_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:10:46
	% EndTime: 2019-12-05 17:10:46
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:10:46
	% EndTime: 2019-12-05 17:10:46
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:10:47
	% EndTime: 2019-12-05 17:10:47
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (5->5), mult. (4->4), div. (0->0), fcn. (10->4), ass. (0->5)
	t9 = cos(qJ(2));
	t8 = sin(qJ(2));
	t7 = cos(pkin(9));
	t6 = sin(pkin(9));
	t1 = [0, -t7 * t8, 0, 0, 0; 0, -t6 * t8, 0, 0, 0; 0, t9, 0, 0, 0; 0, -t7 * t9, 0, 0, 0; 0, -t6 * t9, 0, 0, 0; 0, -t8, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:10:47
	% EndTime: 2019-12-05 17:10:47
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (22->11), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->10)
	t15 = qJ(2) + qJ(3);
	t13 = sin(t15);
	t16 = sin(pkin(9));
	t21 = t16 * t13;
	t14 = cos(t15);
	t20 = t16 * t14;
	t17 = cos(pkin(9));
	t19 = t17 * t13;
	t18 = t17 * t14;
	t1 = [0, -t19, -t19, 0, 0; 0, -t21, -t21, 0, 0; 0, t14, t14, 0, 0; 0, -t18, -t18, 0, 0; 0, -t20, -t20, 0, 0; 0, -t13, -t13, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:10:47
	% EndTime: 2019-12-05 17:10:47
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (36->13), mult. (38->18), div. (0->0), fcn. (66->6), ass. (0->20)
	t76 = qJ(2) + qJ(3);
	t75 = cos(t76);
	t79 = sin(qJ(4));
	t87 = t75 * t79;
	t77 = sin(pkin(9));
	t86 = t77 * t79;
	t80 = cos(qJ(4));
	t85 = t77 * t80;
	t78 = cos(pkin(9));
	t84 = t78 * t79;
	t83 = t78 * t80;
	t74 = sin(t76);
	t82 = t74 * t85;
	t81 = t74 * t83;
	t73 = t75 * t80;
	t72 = t78 * t75;
	t71 = t77 * t75;
	t70 = t74 * t84;
	t69 = t74 * t86;
	t1 = [0, -t81, -t81, -t75 * t84 + t85, 0; 0, -t82, -t82, -t75 * t86 - t83, 0; 0, t73, t73, -t74 * t79, 0; 0, t70, t70, -t75 * t83 - t86, 0; 0, t69, t69, -t75 * t85 + t84, 0; 0, -t87, -t87, -t74 * t80, 0; 0, t72, t72, 0, 0; 0, t71, t71, 0, 0; 0, t74, t74, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:10:47
	% EndTime: 2019-12-05 17:10:47
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (80->16), mult. (52->18), div. (0->0), fcn. (90->6), ass. (0->27)
	t90 = qJ(4) + qJ(5);
	t86 = sin(t90);
	t91 = qJ(2) + qJ(3);
	t87 = sin(t91);
	t102 = t87 * t86;
	t88 = cos(t90);
	t101 = t87 * t88;
	t89 = cos(t91);
	t100 = t89 * t86;
	t92 = sin(pkin(9));
	t99 = t92 * t86;
	t98 = t92 * t88;
	t84 = t92 * t89;
	t93 = cos(pkin(9));
	t97 = t93 * t86;
	t96 = t93 * t88;
	t85 = t93 * t89;
	t95 = t87 * t98;
	t94 = t87 * t96;
	t83 = t89 * t88;
	t82 = t87 * t97;
	t81 = t87 * t99;
	t80 = -t88 * t85 - t99;
	t79 = -t86 * t85 + t98;
	t78 = -t88 * t84 + t97;
	t77 = -t86 * t84 - t96;
	t1 = [0, -t94, -t94, t79, t79; 0, -t95, -t95, t77, t77; 0, t83, t83, -t102, -t102; 0, t82, t82, t80, t80; 0, t81, t81, t78, t78; 0, -t100, -t100, -t101, -t101; 0, t85, t85, 0, 0; 0, t84, t84, 0, 0; 0, t87, t87, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end