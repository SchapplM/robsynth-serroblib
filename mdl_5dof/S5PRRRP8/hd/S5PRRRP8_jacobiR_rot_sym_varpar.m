% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PRRRP8
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:35
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5PRRRP8_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP8_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRP8_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP8_jacobiR_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:35:19
	% EndTime: 2019-10-24 10:35:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:35:19
	% EndTime: 2019-10-24 10:35:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:35:19
	% EndTime: 2019-10-24 10:35:19
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
	% StartTime: 2019-10-24 10:35:19
	% EndTime: 2019-10-24 10:35:19
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (19->13), mult. (57->32), div. (0->0), fcn. (88->8), ass. (0->19)
	t66 = sin(pkin(5));
	t69 = sin(qJ(3));
	t78 = t66 * t69;
	t70 = sin(qJ(2));
	t77 = t66 * t70;
	t71 = cos(qJ(3));
	t76 = t66 * t71;
	t72 = cos(qJ(2));
	t75 = t66 * t72;
	t68 = cos(pkin(5));
	t74 = t68 * t70;
	t73 = t68 * t72;
	t67 = cos(pkin(9));
	t65 = sin(pkin(9));
	t64 = -t65 * t74 + t67 * t72;
	t63 = -t65 * t73 - t67 * t70;
	t62 = t65 * t72 + t67 * t74;
	t61 = -t65 * t70 + t67 * t73;
	t1 = [0, t63 * t71, -t64 * t69 + t65 * t76, 0, 0; 0, t61 * t71, -t62 * t69 - t67 * t76, 0, 0; 0, t71 * t75, t68 * t71 - t69 * t77, 0, 0; 0, -t63 * t69, -t64 * t71 - t65 * t78, 0, 0; 0, -t61 * t69, -t62 * t71 + t67 * t78, 0, 0; 0, -t69 * t75, -t68 * t69 - t70 * t76, 0, 0; 0, t64, 0, 0, 0; 0, t62, 0, 0, 0; 0, t77, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:35:19
	% EndTime: 2019-10-24 10:35:19
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (54->27), mult. (163->63), div. (0->0), fcn. (238->10), ass. (0->29)
	t107 = sin(pkin(5));
	t111 = sin(qJ(3));
	t123 = t107 * t111;
	t114 = cos(qJ(3));
	t122 = t107 * t114;
	t115 = cos(qJ(2));
	t121 = t107 * t115;
	t109 = cos(pkin(5));
	t112 = sin(qJ(2));
	t120 = t109 * t112;
	t119 = t109 * t115;
	t110 = sin(qJ(4));
	t118 = t110 * t114;
	t113 = cos(qJ(4));
	t117 = t113 * t114;
	t116 = t114 * t115;
	t108 = cos(pkin(9));
	t106 = sin(pkin(9));
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
	t1 = [0, -t101 * t117 + t102 * t110, t97 * t113, t101 * t113 - t98 * t110, 0; 0, t100 * t110 - t99 * t117, t95 * t113, -t96 * t110 + t99 * t113, 0; 0, (t110 * t112 + t113 * t116) * t107, t103 * t113, -t104 * t110 - t113 * t121, 0; 0, t101 * t118 + t102 * t113, -t97 * t110, -t101 * t110 - t98 * t113, 0; 0, t100 * t113 + t99 * t118, -t95 * t110, -t99 * t110 - t96 * t113, 0; 0, (-t110 * t116 + t112 * t113) * t107, -t103 * t110, -t104 * t113 + t110 * t121, 0; 0, -t101 * t111, t98, 0, 0; 0, -t99 * t111, t96, 0, 0; 0, t111 * t121, t104, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:35:20
	% EndTime: 2019-10-24 10:35:20
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (51->24), mult. (163->63), div. (0->0), fcn. (238->10), ass. (0->29)
	t135 = sin(pkin(5));
	t139 = sin(qJ(3));
	t151 = t135 * t139;
	t142 = cos(qJ(3));
	t150 = t135 * t142;
	t143 = cos(qJ(2));
	t149 = t135 * t143;
	t137 = cos(pkin(5));
	t140 = sin(qJ(2));
	t148 = t137 * t140;
	t147 = t137 * t143;
	t138 = sin(qJ(4));
	t146 = t138 * t142;
	t141 = cos(qJ(4));
	t145 = t141 * t142;
	t144 = t142 * t143;
	t136 = cos(pkin(9));
	t134 = sin(pkin(9));
	t132 = t137 * t139 + t140 * t150;
	t131 = t137 * t142 - t140 * t151;
	t130 = -t134 * t148 + t136 * t143;
	t129 = t134 * t147 + t136 * t140;
	t128 = t134 * t143 + t136 * t148;
	t127 = t134 * t140 - t136 * t147;
	t126 = t130 * t142 + t134 * t151;
	t125 = -t130 * t139 + t134 * t150;
	t124 = t128 * t142 - t136 * t151;
	t123 = -t128 * t139 - t136 * t150;
	t1 = [0, -t129 * t145 + t130 * t138, t125 * t141, -t126 * t138 + t129 * t141, 0; 0, -t127 * t145 + t128 * t138, t123 * t141, -t124 * t138 + t127 * t141, 0; 0, (t138 * t140 + t141 * t144) * t135, t131 * t141, -t132 * t138 - t141 * t149, 0; 0, -t129 * t139, t126, 0, 0; 0, -t127 * t139, t124, 0, 0; 0, t139 * t149, t132, 0, 0; 0, -t129 * t146 - t130 * t141, t125 * t138, t126 * t141 + t129 * t138, 0; 0, -t127 * t146 - t128 * t141, t123 * t138, t124 * t141 + t127 * t138, 0; 0, (t138 * t144 - t140 * t141) * t135, t131 * t138, t132 * t141 - t138 * t149, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end