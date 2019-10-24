% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PRPRR8
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:28
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5PRPRR8_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR8_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPRR8_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR8_jacobiR_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:28:21
	% EndTime: 2019-10-24 10:28:22
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:28:21
	% EndTime: 2019-10-24 10:28:22
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:28:21
	% EndTime: 2019-10-24 10:28:22
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (5->5), mult. (14->12), div. (0->0), fcn. (24->6), ass. (0->9)
	t22 = cos(pkin(5));
	t23 = sin(qJ(2));
	t26 = t22 * t23;
	t24 = cos(qJ(2));
	t25 = t22 * t24;
	t21 = cos(pkin(9));
	t20 = sin(pkin(5));
	t19 = sin(pkin(9));
	t1 = [0, -t19 * t25 - t21 * t23, 0, 0, 0; 0, -t19 * t23 + t21 * t25, 0, 0, 0; 0, t20 * t24, 0, 0, 0; 0, t19 * t26 - t21 * t24, 0, 0, 0; 0, -t19 * t24 - t21 * t26, 0, 0, 0; 0, -t20 * t23, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:28:22
	% EndTime: 2019-10-24 10:28:22
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (5->5), mult. (14->12), div. (0->0), fcn. (24->6), ass. (0->9)
	t46 = cos(pkin(5));
	t47 = sin(qJ(2));
	t50 = t46 * t47;
	t48 = cos(qJ(2));
	t49 = t46 * t48;
	t45 = cos(pkin(9));
	t44 = sin(pkin(5));
	t43 = sin(pkin(9));
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, t43 * t49 + t45 * t47, 0, 0, 0; 0, t43 * t47 - t45 * t49, 0, 0, 0; 0, -t44 * t48, 0, 0, 0; 0, -t43 * t50 + t45 * t48, 0, 0, 0; 0, t43 * t48 + t45 * t50, 0, 0, 0; 0, t44 * t47, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:28:22
	% EndTime: 2019-10-24 10:28:22
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (16->12), mult. (57->31), div. (0->0), fcn. (88->8), ass. (0->18)
	t66 = sin(pkin(5));
	t69 = sin(qJ(4));
	t77 = t66 * t69;
	t71 = cos(qJ(4));
	t76 = t66 * t71;
	t72 = cos(qJ(2));
	t75 = t66 * t72;
	t68 = cos(pkin(5));
	t70 = sin(qJ(2));
	t74 = t68 * t70;
	t73 = t68 * t72;
	t67 = cos(pkin(9));
	t65 = sin(pkin(9));
	t64 = -t65 * t74 + t67 * t72;
	t63 = t65 * t73 + t67 * t70;
	t62 = t65 * t72 + t67 * t74;
	t61 = t65 * t70 - t67 * t73;
	t1 = [0, t64 * t69, 0, t63 * t71 - t65 * t77, 0; 0, t62 * t69, 0, t61 * t71 + t67 * t77, 0; 0, t70 * t77, 0, -t68 * t69 - t71 * t75, 0; 0, t64 * t71, 0, -t63 * t69 - t65 * t76, 0; 0, t62 * t71, 0, -t61 * t69 + t67 * t76, 0; 0, t70 * t76, 0, -t68 * t71 + t69 * t75, 0; 0, -t63, 0, 0, 0; 0, -t61, 0, 0, 0; 0, t75, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:28:22
	% EndTime: 2019-10-24 10:28:22
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (57->28), mult. (163->64), div. (0->0), fcn. (238->10), ass. (0->30)
	t108 = sin(pkin(5));
	t112 = sin(qJ(4));
	t125 = t108 * t112;
	t115 = cos(qJ(4));
	t124 = t108 * t115;
	t116 = cos(qJ(2));
	t123 = t108 * t116;
	t110 = cos(pkin(5));
	t113 = sin(qJ(2));
	t122 = t110 * t113;
	t121 = t110 * t116;
	t111 = sin(qJ(5));
	t120 = t111 * t112;
	t119 = t111 * t113;
	t114 = cos(qJ(5));
	t118 = t112 * t114;
	t117 = t113 * t114;
	t109 = cos(pkin(9));
	t107 = sin(pkin(9));
	t105 = t110 * t115 - t112 * t123;
	t104 = -t110 * t112 - t115 * t123;
	t103 = -t107 * t122 + t109 * t116;
	t102 = t107 * t121 + t109 * t113;
	t101 = t107 * t116 + t109 * t122;
	t100 = t107 * t113 - t109 * t121;
	t99 = t100 * t112 - t109 * t124;
	t98 = t100 * t115 + t109 * t125;
	t97 = t102 * t112 + t107 * t124;
	t96 = t102 * t115 - t107 * t125;
	t1 = [0, -t102 * t111 + t103 * t118, 0, t96 * t114, t103 * t114 - t97 * t111; 0, -t100 * t111 + t101 * t118, 0, t98 * t114, t101 * t114 - t99 * t111; 0, (t111 * t116 + t112 * t117) * t108, 0, t104 * t114, -t105 * t111 + t108 * t117; 0, -t102 * t114 - t103 * t120, 0, -t96 * t111, -t103 * t111 - t97 * t114; 0, -t100 * t114 - t101 * t120, 0, -t98 * t111, -t101 * t111 - t99 * t114; 0, (-t112 * t119 + t114 * t116) * t108, 0, -t104 * t111, -t105 * t114 - t108 * t119; 0, -t103 * t115, 0, t97, 0; 0, -t101 * t115, 0, t99, 0; 0, -t113 * t124, 0, t105, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end