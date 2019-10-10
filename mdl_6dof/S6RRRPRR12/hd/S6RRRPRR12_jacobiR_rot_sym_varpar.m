% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRPRR12
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:12
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPRR12_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR12_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR12_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR12_jacobiR_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:12:18
	% EndTime: 2019-10-10 12:12:18
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:12:18
	% EndTime: 2019-10-10 12:12:18
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8, 0, 0, 0, 0, 0; t9, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0, 0; -t8, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:12:18
	% EndTime: 2019-10-10 12:12:18
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->7), mult. (28->12), div. (0->0), fcn. (48->6), ass. (0->15)
	t50 = sin(qJ(2));
	t51 = sin(qJ(1));
	t57 = t51 * t50;
	t52 = cos(qJ(2));
	t56 = t51 * t52;
	t53 = cos(qJ(1));
	t55 = t53 * t50;
	t54 = t53 * t52;
	t49 = cos(pkin(6));
	t48 = sin(pkin(6));
	t47 = -t49 * t57 + t54;
	t46 = -t49 * t56 - t55;
	t45 = -t49 * t55 - t56;
	t44 = -t49 * t54 + t57;
	t1 = [t45, t46, 0, 0, 0, 0; t47, -t44, 0, 0, 0, 0; 0, t48 * t52, 0, 0, 0, 0; t44, -t47, 0, 0, 0, 0; t46, t45, 0, 0, 0, 0; 0, -t48 * t50, 0, 0, 0, 0; t53 * t48, 0, 0, 0, 0, 0; t51 * t48, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:12:18
	% EndTime: 2019-10-10 12:12:19
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (29->15), mult. (87->31), div. (0->0), fcn. (134->8), ass. (0->25)
	t83 = sin(pkin(6));
	t86 = sin(qJ(2));
	t100 = t83 * t86;
	t88 = cos(qJ(3));
	t99 = t83 * t88;
	t89 = cos(qJ(2));
	t98 = t83 * t89;
	t90 = cos(qJ(1));
	t97 = t83 * t90;
	t87 = sin(qJ(1));
	t96 = t87 * t86;
	t95 = t87 * t89;
	t94 = t90 * t86;
	t93 = t90 * t89;
	t84 = cos(pkin(6));
	t79 = t84 * t94 + t95;
	t85 = sin(qJ(3));
	t92 = -t79 * t88 + t85 * t97;
	t91 = t79 * t85 + t88 * t97;
	t81 = -t84 * t96 + t93;
	t80 = t84 * t95 + t94;
	t78 = t84 * t93 - t96;
	t77 = t87 * t83 * t85 + t81 * t88;
	t76 = -t81 * t85 + t87 * t99;
	t1 = [t92, -t80 * t88, t76, 0, 0, 0; t77, t78 * t88, -t91, 0, 0, 0; 0, t88 * t98, -t85 * t100 + t84 * t88, 0, 0, 0; t91, t80 * t85, -t77, 0, 0, 0; t76, -t78 * t85, t92, 0, 0, 0; 0, -t85 * t98, -t84 * t85 - t86 * t99, 0, 0, 0; t78, t81, 0, 0, 0, 0; t80, t79, 0, 0, 0, 0; 0, t100, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:12:19
	% EndTime: 2019-10-10 12:12:19
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (54->25), mult. (163->57), div. (0->0), fcn. (238->10), ass. (0->30)
	t111 = sin(pkin(12));
	t118 = cos(qJ(3));
	t130 = t111 * t118;
	t112 = sin(pkin(6));
	t115 = sin(qJ(3));
	t129 = t112 * t115;
	t128 = t112 * t118;
	t120 = cos(qJ(1));
	t127 = t112 * t120;
	t113 = cos(pkin(12));
	t126 = t113 * t118;
	t116 = sin(qJ(2));
	t117 = sin(qJ(1));
	t125 = t117 * t116;
	t119 = cos(qJ(2));
	t124 = t117 * t119;
	t123 = t118 * t119;
	t122 = t120 * t116;
	t121 = t120 * t119;
	t114 = cos(pkin(6));
	t107 = t114 * t122 + t124;
	t101 = -t107 * t115 - t118 * t127;
	t102 = -t107 * t118 + t115 * t127;
	t109 = -t114 * t125 + t121;
	t108 = t114 * t124 + t122;
	t106 = t114 * t121 - t125;
	t105 = t114 * t118 - t116 * t129;
	t104 = t109 * t118 + t117 * t129;
	t103 = t109 * t115 - t117 * t128;
	t1 = [t102 * t113 + t106 * t111, -t108 * t126 + t109 * t111, -t103 * t113, 0, 0, 0; t104 * t113 + t108 * t111, t106 * t126 + t107 * t111, t101 * t113, 0, 0, 0; 0, (t111 * t116 + t113 * t123) * t112, t105 * t113, 0, 0, 0; -t102 * t111 + t106 * t113, t108 * t130 + t109 * t113, t103 * t111, 0, 0, 0; -t104 * t111 + t108 * t113, -t106 * t130 + t107 * t113, -t101 * t111, 0, 0, 0; 0, (-t111 * t123 + t113 * t116) * t112, -t105 * t111, 0, 0, 0; t101, -t108 * t115, t104, 0, 0, 0; t103, t106 * t115, -t102, 0, 0, 0; 0, t119 * t129, t114 * t115 + t116 * t128, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:12:19
	% EndTime: 2019-10-10 12:12:19
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (112->31), mult. (219->62), div. (0->0), fcn. (320->10), ass. (0->37)
	t136 = cos(pkin(6));
	t138 = sin(qJ(2));
	t142 = cos(qJ(1));
	t144 = t142 * t138;
	t139 = sin(qJ(1));
	t141 = cos(qJ(2));
	t146 = t139 * t141;
	t127 = t136 * t144 + t146;
	t137 = sin(qJ(3));
	t140 = cos(qJ(3));
	t135 = sin(pkin(6));
	t148 = t135 * t142;
	t121 = -t127 * t140 + t137 * t148;
	t143 = t142 * t141;
	t147 = t139 * t138;
	t126 = -t136 * t143 + t147;
	t134 = pkin(12) + qJ(5);
	t132 = sin(t134);
	t133 = cos(t134);
	t157 = t121 * t132 + t126 * t133;
	t156 = t121 * t133 - t126 * t132;
	t153 = t132 * t140;
	t152 = t133 * t140;
	t151 = t135 * t137;
	t150 = t135 * t140;
	t149 = t135 * t141;
	t145 = t140 * t141;
	t119 = -t127 * t137 - t140 * t148;
	t129 = -t136 * t147 + t143;
	t128 = t136 * t146 + t144;
	t125 = t136 * t137 + t138 * t150;
	t124 = t136 * t140 - t138 * t151;
	t123 = t129 * t140 + t139 * t151;
	t122 = t129 * t137 - t139 * t150;
	t118 = t123 * t133 + t128 * t132;
	t117 = -t123 * t132 + t128 * t133;
	t1 = [t156, -t128 * t152 + t129 * t132, -t122 * t133, 0, t117, 0; t118, -t126 * t152 + t127 * t132, t119 * t133, 0, t157, 0; 0, (t132 * t138 + t133 * t145) * t135, t124 * t133, 0, -t125 * t132 - t133 * t149, 0; -t157, t128 * t153 + t129 * t133, t122 * t132, 0, -t118, 0; t117, t126 * t153 + t127 * t133, -t119 * t132, 0, t156, 0; 0, (-t132 * t145 + t133 * t138) * t135, -t124 * t132, 0, -t125 * t133 + t132 * t149, 0; t119, -t128 * t137, t123, 0, 0, 0; t122, -t126 * t137, -t121, 0, 0, 0; 0, t137 * t149, t125, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:12:19
	% EndTime: 2019-10-10 12:12:19
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (194->32), mult. (275->62), div. (0->0), fcn. (402->10), ass. (0->39)
	t156 = cos(pkin(6));
	t158 = sin(qJ(2));
	t162 = cos(qJ(1));
	t164 = t162 * t158;
	t159 = sin(qJ(1));
	t161 = cos(qJ(2));
	t166 = t159 * t161;
	t147 = t156 * t164 + t166;
	t157 = sin(qJ(3));
	t160 = cos(qJ(3));
	t155 = sin(pkin(6));
	t168 = t155 * t162;
	t141 = -t147 * t160 + t157 * t168;
	t163 = t162 * t161;
	t167 = t159 * t158;
	t146 = -t156 * t163 + t167;
	t154 = pkin(12) + qJ(5) + qJ(6);
	t152 = sin(t154);
	t153 = cos(t154);
	t133 = t141 * t152 + t146 * t153;
	t134 = t141 * t153 - t146 * t152;
	t173 = t152 * t160;
	t172 = t153 * t160;
	t171 = t155 * t157;
	t170 = t155 * t160;
	t169 = t155 * t161;
	t165 = t160 * t161;
	t139 = -t147 * t157 - t160 * t168;
	t149 = -t156 * t167 + t163;
	t148 = t156 * t166 + t164;
	t145 = t156 * t157 + t158 * t170;
	t144 = t156 * t160 - t158 * t171;
	t143 = t149 * t160 + t159 * t171;
	t142 = t149 * t157 - t159 * t170;
	t138 = -t145 * t153 + t152 * t169;
	t137 = -t145 * t152 - t153 * t169;
	t136 = t143 * t153 + t148 * t152;
	t135 = -t143 * t152 + t148 * t153;
	t1 = [t134, -t148 * t172 + t149 * t152, -t142 * t153, 0, t135, t135; t136, -t146 * t172 + t147 * t152, t139 * t153, 0, t133, t133; 0, (t152 * t158 + t153 * t165) * t155, t144 * t153, 0, t137, t137; -t133, t148 * t173 + t149 * t153, t142 * t152, 0, -t136, -t136; t135, t146 * t173 + t147 * t153, -t139 * t152, 0, t134, t134; 0, (-t152 * t165 + t153 * t158) * t155, -t144 * t152, 0, t138, t138; t139, -t148 * t157, t143, 0, 0, 0; t142, -t146 * t157, -t141, 0, 0, 0; 0, t157 * t169, t145, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end