% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:12
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRPPR4_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR4_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPPR4_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR4_jacobiR_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:12:44
	% EndTime: 2019-10-09 22:12:44
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:12:44
	% EndTime: 2019-10-09 22:12:44
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:12:44
	% EndTime: 2019-10-09 22:12:44
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (5->5), mult. (14->12), div. (0->0), fcn. (24->6), ass. (0->9)
	t22 = cos(pkin(6));
	t23 = sin(qJ(2));
	t26 = t22 * t23;
	t24 = cos(qJ(2));
	t25 = t22 * t24;
	t21 = cos(pkin(10));
	t20 = sin(pkin(6));
	t19 = sin(pkin(10));
	t1 = [0, -t19 * t25 - t21 * t23, 0, 0, 0, 0; 0, -t19 * t23 + t21 * t25, 0, 0, 0, 0; 0, t20 * t24, 0, 0, 0, 0; 0, t19 * t26 - t21 * t24, 0, 0, 0, 0; 0, -t19 * t24 - t21 * t26, 0, 0, 0, 0; 0, -t20 * t23, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:12:44
	% EndTime: 2019-10-09 22:12:44
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (19->13), mult. (57->32), div. (0->0), fcn. (88->8), ass. (0->19)
	t66 = sin(pkin(6));
	t69 = sin(qJ(3));
	t78 = t66 * t69;
	t70 = sin(qJ(2));
	t77 = t66 * t70;
	t71 = cos(qJ(3));
	t76 = t66 * t71;
	t72 = cos(qJ(2));
	t75 = t66 * t72;
	t68 = cos(pkin(6));
	t74 = t68 * t70;
	t73 = t68 * t72;
	t67 = cos(pkin(10));
	t65 = sin(pkin(10));
	t64 = -t65 * t74 + t67 * t72;
	t63 = -t65 * t73 - t67 * t70;
	t62 = t65 * t72 + t67 * t74;
	t61 = -t65 * t70 + t67 * t73;
	t1 = [0, t63 * t71, -t64 * t69 + t65 * t76, 0, 0, 0; 0, t61 * t71, -t62 * t69 - t67 * t76, 0, 0, 0; 0, t71 * t75, t68 * t71 - t69 * t77, 0, 0, 0; 0, -t63 * t69, -t64 * t71 - t65 * t78, 0, 0, 0; 0, -t61 * t69, -t62 * t71 + t67 * t78, 0, 0, 0; 0, -t69 * t75, -t68 * t69 - t70 * t76, 0, 0, 0; 0, t64, 0, 0, 0, 0; 0, t62, 0, 0, 0, 0; 0, t77, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:12:44
	% EndTime: 2019-10-09 22:12:44
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (34->19), mult. (107->48), div. (0->0), fcn. (156->10), ass. (0->27)
	t100 = sin(qJ(3));
	t96 = sin(pkin(6));
	t112 = t100 * t96;
	t102 = cos(qJ(3));
	t94 = sin(pkin(11));
	t111 = t102 * t94;
	t110 = t102 * t96;
	t97 = cos(pkin(11));
	t109 = t102 * t97;
	t101 = sin(qJ(2));
	t95 = sin(pkin(10));
	t108 = t95 * t101;
	t103 = cos(qJ(2));
	t107 = t95 * t103;
	t98 = cos(pkin(10));
	t106 = t98 * t101;
	t105 = t98 * t103;
	t104 = t102 * t103;
	t99 = cos(pkin(6));
	t93 = -t101 * t112 + t99 * t102;
	t92 = -t99 * t108 + t105;
	t91 = -t99 * t107 - t106;
	t90 = t99 * t106 + t107;
	t89 = t99 * t105 - t108;
	t88 = -t92 * t100 + t95 * t110;
	t87 = -t90 * t100 - t98 * t110;
	t1 = [0, t91 * t109 + t92 * t94, t88 * t97, 0, 0, 0; 0, t89 * t109 + t90 * t94, t87 * t97, 0, 0, 0; 0, (t101 * t94 + t97 * t104) * t96, t93 * t97, 0, 0, 0; 0, -t91 * t111 + t92 * t97, -t88 * t94, 0, 0, 0; 0, -t89 * t111 + t90 * t97, -t87 * t94, 0, 0, 0; 0, (t101 * t97 - t94 * t104) * t96, -t93 * t94, 0, 0, 0; 0, t91 * t100, t92 * t102 + t95 * t112, 0, 0, 0; 0, t89 * t100, t90 * t102 - t98 * t112, 0, 0, 0; 0, t103 * t112, t99 * t100 + t101 * t110, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:12:44
	% EndTime: 2019-10-09 22:12:44
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (31->16), mult. (107->50), div. (0->0), fcn. (156->10), ass. (0->25)
	t111 = sin(pkin(11));
	t119 = cos(qJ(3));
	t127 = t111 * t119;
	t113 = sin(pkin(6));
	t117 = sin(qJ(3));
	t126 = t113 * t117;
	t125 = t113 * t119;
	t114 = cos(pkin(11));
	t124 = t114 * t119;
	t116 = cos(pkin(6));
	t118 = sin(qJ(2));
	t123 = t116 * t118;
	t120 = cos(qJ(2));
	t122 = t116 * t120;
	t121 = t119 * t120;
	t115 = cos(pkin(10));
	t112 = sin(pkin(10));
	t110 = t116 * t119 - t118 * t126;
	t109 = -t112 * t123 + t115 * t120;
	t108 = -t112 * t122 - t115 * t118;
	t107 = t112 * t120 + t115 * t123;
	t106 = -t112 * t118 + t115 * t122;
	t105 = -t109 * t117 + t112 * t125;
	t104 = -t107 * t117 - t115 * t125;
	t1 = [0, t108 * t124 + t109 * t111, t105 * t114, 0, 0, 0; 0, t106 * t124 + t107 * t111, t104 * t114, 0, 0, 0; 0, (t111 * t118 + t114 * t121) * t113, t110 * t114, 0, 0, 0; 0, t108 * t117, t109 * t119 + t112 * t126, 0, 0, 0; 0, t106 * t117, t107 * t119 - t115 * t126, 0, 0, 0; 0, t120 * t126, t116 * t117 + t118 * t125, 0, 0, 0; 0, t108 * t127 - t109 * t114, t105 * t111, 0, 0, 0; 0, t106 * t127 - t107 * t114, t104 * t111, 0, 0, 0; 0, (t111 * t121 - t114 * t118) * t113, t110 * t111, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:12:44
	% EndTime: 2019-10-09 22:12:44
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (116->40), mult. (337->91), div. (0->0), fcn. (474->12), ass. (0->45)
	t144 = sin(pkin(11));
	t154 = cos(qJ(3));
	t165 = t144 * t154;
	t146 = sin(pkin(6));
	t151 = sin(qJ(3));
	t164 = t146 * t151;
	t163 = t146 * t154;
	t155 = cos(qJ(2));
	t162 = t146 * t155;
	t147 = cos(pkin(11));
	t161 = t147 * t154;
	t149 = cos(pkin(6));
	t152 = sin(qJ(2));
	t160 = t149 * t152;
	t159 = t149 * t155;
	t158 = t154 * t155;
	t150 = sin(qJ(6));
	t153 = cos(qJ(6));
	t157 = t144 * t153 - t147 * t150;
	t156 = t144 * t150 + t147 * t153;
	t148 = cos(pkin(10));
	t145 = sin(pkin(10));
	t142 = t149 * t151 + t152 * t163;
	t141 = t149 * t154 - t152 * t164;
	t140 = -t145 * t160 + t148 * t155;
	t139 = t145 * t159 + t148 * t152;
	t138 = t145 * t155 + t148 * t160;
	t137 = t145 * t152 - t148 * t159;
	t136 = (t144 * t152 + t147 * t158) * t146;
	t135 = (t144 * t158 - t147 * t152) * t146;
	t134 = t140 * t154 + t145 * t164;
	t133 = -t140 * t151 + t145 * t163;
	t132 = t138 * t154 - t148 * t164;
	t131 = -t138 * t151 - t148 * t163;
	t130 = t142 * t147 - t144 * t162;
	t129 = t142 * t144 + t147 * t162;
	t128 = -t139 * t161 + t140 * t144;
	t127 = -t139 * t165 - t140 * t147;
	t126 = -t137 * t161 + t138 * t144;
	t125 = -t137 * t165 - t138 * t147;
	t124 = t134 * t147 + t139 * t144;
	t123 = t134 * t144 - t139 * t147;
	t122 = t132 * t147 + t137 * t144;
	t121 = t132 * t144 - t137 * t147;
	t1 = [0, t127 * t150 + t128 * t153, t156 * t133, 0, 0, t123 * t153 - t124 * t150; 0, t125 * t150 + t126 * t153, t156 * t131, 0, 0, t121 * t153 - t122 * t150; 0, t135 * t150 + t136 * t153, t156 * t141, 0, 0, t129 * t153 - t130 * t150; 0, t127 * t153 - t128 * t150, t157 * t133, 0, 0, -t123 * t150 - t124 * t153; 0, t125 * t153 - t126 * t150, t157 * t131, 0, 0, -t121 * t150 - t122 * t153; 0, t135 * t153 - t136 * t150, t157 * t141, 0, 0, -t129 * t150 - t130 * t153; 0, t139 * t151, -t134, 0, 0, 0; 0, t137 * t151, -t132, 0, 0, 0; 0, -t151 * t162, -t142, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end