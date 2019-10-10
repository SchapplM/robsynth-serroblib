% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPRRR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:07
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRRR12_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR12_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR12_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR12_jacobiR_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:07:01
	% EndTime: 2019-10-10 11:07:01
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:07:01
	% EndTime: 2019-10-10 11:07:01
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
	% StartTime: 2019-10-10 11:07:02
	% EndTime: 2019-10-10 11:07:02
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
	% StartTime: 2019-10-10 11:07:02
	% EndTime: 2019-10-10 11:07:02
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->7), mult. (28->12), div. (0->0), fcn. (48->6), ass. (0->15)
	t62 = sin(qJ(2));
	t63 = sin(qJ(1));
	t69 = t63 * t62;
	t64 = cos(qJ(2));
	t68 = t63 * t64;
	t65 = cos(qJ(1));
	t67 = t65 * t62;
	t66 = t65 * t64;
	t61 = cos(pkin(6));
	t60 = sin(pkin(6));
	t59 = -t61 * t69 + t66;
	t58 = t61 * t68 + t67;
	t57 = t61 * t67 + t68;
	t56 = -t61 * t66 + t69;
	t1 = [t65 * t60, 0, 0, 0, 0, 0; t63 * t60, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t57, t58, 0, 0, 0, 0; -t59, t56, 0, 0, 0, 0; 0, -t60 * t64, 0, 0, 0, 0; -t56, t59, 0, 0, 0, 0; t58, t57, 0, 0, 0, 0; 0, t60 * t62, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:07:02
	% EndTime: 2019-10-10 11:07:02
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (26->15), mult. (87->30), div. (0->0), fcn. (134->8), ass. (0->25)
	t85 = sin(pkin(6));
	t87 = sin(qJ(4));
	t102 = t85 * t87;
	t90 = cos(qJ(4));
	t101 = t85 * t90;
	t91 = cos(qJ(2));
	t100 = t85 * t91;
	t92 = cos(qJ(1));
	t99 = t85 * t92;
	t88 = sin(qJ(2));
	t89 = sin(qJ(1));
	t98 = t89 * t88;
	t97 = t89 * t91;
	t96 = t92 * t88;
	t95 = t92 * t91;
	t86 = cos(pkin(6));
	t79 = -t86 * t95 + t98;
	t94 = -t79 * t87 + t90 * t99;
	t93 = t79 * t90 + t87 * t99;
	t82 = -t86 * t98 + t95;
	t81 = t86 * t97 + t96;
	t80 = t86 * t96 + t97;
	t78 = t101 * t89 + t81 * t87;
	t77 = -t102 * t89 + t81 * t90;
	t1 = [t94, t82 * t87, 0, t77, 0, 0; t78, t80 * t87, 0, t93, 0, 0; 0, t88 * t102, 0, -t100 * t90 - t86 * t87, 0, 0; -t93, t82 * t90, 0, -t78, 0, 0; t77, t80 * t90, 0, t94, 0, 0; 0, t88 * t101, 0, t100 * t87 - t86 * t90, 0, 0; -t80, -t81, 0, 0, 0, 0; t82, -t79, 0, 0, 0, 0; 0, t100, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:07:02
	% EndTime: 2019-10-10 11:07:02
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (74->17), mult. (117->30), div. (0->0), fcn. (180->8), ass. (0->28)
	t112 = sin(pkin(6));
	t114 = sin(qJ(2));
	t125 = t112 * t114;
	t115 = sin(qJ(1));
	t124 = t112 * t115;
	t116 = cos(qJ(2));
	t123 = t112 * t116;
	t117 = cos(qJ(1));
	t122 = t112 * t117;
	t121 = t115 * t114;
	t120 = t115 * t116;
	t119 = t117 * t114;
	t118 = t117 * t116;
	t113 = cos(pkin(6));
	t103 = -t113 * t118 + t121;
	t111 = qJ(4) + qJ(5);
	t109 = sin(t111);
	t110 = cos(t111);
	t100 = -t103 * t109 + t110 * t122;
	t99 = t103 * t110 + t109 * t122;
	t106 = -t113 * t121 + t118;
	t105 = t113 * t120 + t119;
	t104 = t113 * t119 + t120;
	t102 = t109 * t123 - t113 * t110;
	t101 = -t113 * t109 - t110 * t123;
	t98 = t105 * t109 + t110 * t124;
	t97 = t105 * t110 - t109 * t124;
	t1 = [t100, t106 * t109, 0, t97, t97, 0; t98, t104 * t109, 0, t99, t99, 0; 0, t109 * t125, 0, t101, t101, 0; -t99, t106 * t110, 0, -t98, -t98, 0; t97, t104 * t110, 0, t100, t100, 0; 0, t110 * t125, 0, t102, t102, 0; -t104, -t105, 0, 0, 0, 0; t106, -t103, 0, 0, 0, 0; 0, t123, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:07:02
	% EndTime: 2019-10-10 11:07:02
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (164->36), mult. (270->63), div. (0->0), fcn. (395->10), ass. (0->43)
	t164 = cos(pkin(6));
	t169 = cos(qJ(2));
	t170 = cos(qJ(1));
	t171 = t170 * t169;
	t166 = sin(qJ(2));
	t167 = sin(qJ(1));
	t174 = t167 * t166;
	t154 = -t164 * t171 + t174;
	t162 = qJ(4) + qJ(5);
	t160 = sin(t162);
	t161 = cos(t162);
	t163 = sin(pkin(6));
	t177 = t163 * t170;
	t148 = -t154 * t160 + t161 * t177;
	t172 = t170 * t166;
	t173 = t167 * t169;
	t155 = t164 * t172 + t173;
	t165 = sin(qJ(6));
	t168 = cos(qJ(6));
	t188 = t148 * t165 + t155 * t168;
	t187 = t148 * t168 - t155 * t165;
	t156 = t164 * t173 + t172;
	t179 = t163 * t167;
	t145 = -t156 * t161 + t160 * t179;
	t186 = t145 * t165;
	t147 = t154 * t161 + t160 * t177;
	t185 = t147 * t165;
	t178 = t163 * t169;
	t152 = -t164 * t160 - t161 * t178;
	t184 = t152 * t165;
	t181 = t160 * t165;
	t180 = t160 * t168;
	t176 = t165 * t166;
	t175 = t166 * t168;
	t157 = -t164 * t174 + t171;
	t153 = -t160 * t178 + t164 * t161;
	t150 = t152 * t168;
	t146 = t156 * t160 + t161 * t179;
	t144 = t147 * t168;
	t143 = t145 * t168;
	t142 = t146 * t168 + t157 * t165;
	t141 = -t146 * t165 + t157 * t168;
	t1 = [t187, -t156 * t165 + t157 * t180, 0, -t143, -t143, t141; t142, -t154 * t165 + t155 * t180, 0, t144, t144, t188; 0, (t160 * t175 + t165 * t169) * t163, 0, t150, t150, -t153 * t165 + t163 * t175; -t188, -t156 * t168 - t157 * t181, 0, t186, t186, -t142; t141, -t154 * t168 - t155 * t181, 0, -t185, -t185, t187; 0, (-t160 * t176 + t168 * t169) * t163, 0, -t184, -t184, -t153 * t168 - t163 * t176; t147, -t157 * t161, 0, t146, t146, 0; t145, -t155 * t161, 0, -t148, -t148, 0; 0, -t163 * t166 * t161, 0, t153, t153, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end