% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S4PRRR7
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% JR_rot [9x4]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 12:35
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S4PRRR7_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR7_jacobiR_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4PRRR7_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4PRRR7_jacobiR_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:34:58
	% EndTime: 2019-12-29 12:34:58
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:34:58
	% EndTime: 2019-12-29 12:34:58
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:34:58
	% EndTime: 2019-12-29 12:34:58
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (5->5), mult. (14->12), div. (0->0), fcn. (24->6), ass. (0->9)
	t22 = cos(pkin(4));
	t23 = sin(qJ(2));
	t26 = t22 * t23;
	t24 = cos(qJ(2));
	t25 = t22 * t24;
	t21 = cos(pkin(8));
	t20 = sin(pkin(4));
	t19 = sin(pkin(8));
	t1 = [0, -t19 * t25 - t21 * t23, 0, 0; 0, -t19 * t23 + t21 * t25, 0, 0; 0, t20 * t24, 0, 0; 0, t19 * t26 - t21 * t24, 0, 0; 0, -t19 * t24 - t21 * t26, 0, 0; 0, -t20 * t23, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:34:58
	% EndTime: 2019-12-29 12:34:59
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (19->13), mult. (57->32), div. (0->0), fcn. (88->8), ass. (0->19)
	t65 = sin(pkin(4));
	t68 = sin(qJ(3));
	t77 = t65 * t68;
	t69 = sin(qJ(2));
	t76 = t65 * t69;
	t70 = cos(qJ(3));
	t75 = t65 * t70;
	t71 = cos(qJ(2));
	t74 = t65 * t71;
	t67 = cos(pkin(4));
	t73 = t67 * t69;
	t72 = t67 * t71;
	t66 = cos(pkin(8));
	t64 = sin(pkin(8));
	t63 = -t64 * t73 + t66 * t71;
	t62 = -t64 * t72 - t66 * t69;
	t61 = t64 * t71 + t66 * t73;
	t60 = -t64 * t69 + t66 * t72;
	t1 = [0, t62 * t70, -t63 * t68 + t64 * t75, 0; 0, t60 * t70, -t61 * t68 - t66 * t75, 0; 0, t70 * t74, t67 * t70 - t68 * t76, 0; 0, -t62 * t68, -t63 * t70 - t64 * t77, 0; 0, -t60 * t68, -t61 * t70 + t66 * t77, 0; 0, -t68 * t74, -t67 * t68 - t69 * t75, 0; 0, t63, 0, 0; 0, t61, 0, 0; 0, t76, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:34:59
	% EndTime: 2019-12-29 12:34:59
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (54->27), mult. (163->63), div. (0->0), fcn. (238->10), ass. (0->29)
	t107 = sin(pkin(4));
	t111 = sin(qJ(3));
	t123 = t107 * t111;
	t114 = cos(qJ(3));
	t122 = t107 * t114;
	t115 = cos(qJ(2));
	t121 = t107 * t115;
	t109 = cos(pkin(4));
	t112 = sin(qJ(2));
	t120 = t109 * t112;
	t119 = t109 * t115;
	t110 = sin(qJ(4));
	t118 = t110 * t114;
	t113 = cos(qJ(4));
	t117 = t113 * t114;
	t116 = t114 * t115;
	t108 = cos(pkin(8));
	t106 = sin(pkin(8));
	t104 = t109 * t111 + t112 * t122;
	t103 = t109 * t114 - t112 * t123;
	t102 = -t106 * t120 + t108 * t115;
	t101 = t106 * t119 + t108 * t112;
	t100 = t106 * t115 + t108 * t120;
	t99 = t106 * t112 - t108 * t119;
	t98 = t102 * t114 + t106 * t123;
	t97 = -t102 * t111 + t106 * t122;
	t96 = t100 * t114 - t108 * t123;
	t95 = -t100 * t111 - t108 * t122;
	t1 = [0, -t101 * t117 + t102 * t110, t97 * t113, t101 * t113 - t98 * t110; 0, t100 * t110 - t99 * t117, t95 * t113, -t96 * t110 + t99 * t113; 0, (t110 * t112 + t113 * t116) * t107, t103 * t113, -t104 * t110 - t113 * t121; 0, t101 * t118 + t102 * t113, -t97 * t110, -t101 * t110 - t98 * t113; 0, t100 * t113 + t99 * t118, -t95 * t110, -t99 * t110 - t96 * t113; 0, (-t110 * t116 + t112 * t113) * t107, -t103 * t110, -t104 * t113 + t110 * t121; 0, -t101 * t111, t98, 0; 0, -t99 * t111, t96, 0; 0, t111 * t121, t104, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,4);
end