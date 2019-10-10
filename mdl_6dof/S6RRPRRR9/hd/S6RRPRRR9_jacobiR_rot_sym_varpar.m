% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPRRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:03
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRRR9_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR9_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR9_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR9_jacobiR_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:03:13
	% EndTime: 2019-10-10 11:03:13
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:03:13
	% EndTime: 2019-10-10 11:03:13
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
	% StartTime: 2019-10-10 11:03:13
	% EndTime: 2019-10-10 11:03:13
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
	% StartTime: 2019-10-10 11:03:13
	% EndTime: 2019-10-10 11:03:13
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (19->11), mult. (57->26), div. (0->0), fcn. (88->8), ass. (0->20)
	t66 = sin(pkin(6));
	t70 = sin(qJ(1));
	t79 = t66 * t70;
	t71 = cos(qJ(2));
	t78 = t66 * t71;
	t72 = cos(qJ(1));
	t77 = t66 * t72;
	t69 = sin(qJ(2));
	t76 = t70 * t69;
	t75 = t70 * t71;
	t74 = t72 * t69;
	t73 = t72 * t71;
	t68 = cos(pkin(6));
	t67 = cos(pkin(12));
	t65 = sin(pkin(12));
	t64 = -t68 * t76 + t73;
	t63 = t68 * t75 + t74;
	t62 = t68 * t74 + t75;
	t61 = t68 * t73 - t76;
	t1 = [-t62 * t67 + t65 * t77, -t63 * t67, 0, 0, 0, 0; t64 * t67 + t65 * t79, t61 * t67, 0, 0, 0, 0; 0, t67 * t78, 0, 0, 0, 0; t62 * t65 + t67 * t77, t63 * t65, 0, 0, 0, 0; -t64 * t65 + t67 * t79, -t61 * t65, 0, 0, 0, 0; 0, -t65 * t78, 0, 0, 0, 0; t61, t64, 0, 0, 0, 0; t63, t62, 0, 0, 0, 0; 0, t66 * t69, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:03:13
	% EndTime: 2019-10-10 11:03:13
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (55->16), mult. (87->30), div. (0->0), fcn. (134->8), ass. (0->26)
	t90 = sin(pkin(6));
	t92 = sin(qJ(2));
	t105 = t90 * t92;
	t93 = sin(qJ(1));
	t104 = t90 * t93;
	t94 = cos(qJ(2));
	t103 = t90 * t94;
	t95 = cos(qJ(1));
	t102 = t90 * t95;
	t101 = t93 * t92;
	t100 = t93 * t94;
	t99 = t95 * t92;
	t98 = t95 * t94;
	t91 = cos(pkin(6));
	t83 = t91 * t99 + t100;
	t89 = pkin(12) + qJ(4);
	t87 = sin(t89);
	t88 = cos(t89);
	t97 = t87 * t102 - t83 * t88;
	t96 = t88 * t102 + t83 * t87;
	t85 = -t91 * t101 + t98;
	t84 = t91 * t100 + t99;
	t82 = t91 * t98 - t101;
	t81 = t87 * t104 + t85 * t88;
	t80 = t88 * t104 - t85 * t87;
	t1 = [t97, -t84 * t88, 0, t80, 0, 0; t81, t82 * t88, 0, -t96, 0, 0; 0, t88 * t103, 0, -t87 * t105 + t91 * t88, 0, 0; t96, t84 * t87, 0, -t81, 0, 0; t80, -t82 * t87, 0, t97, 0, 0; 0, -t87 * t103, 0, -t88 * t105 - t91 * t87, 0, 0; t82, t85, 0, 0, 0, 0; t84, t83, 0, 0, 0, 0; 0, t105, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:03:13
	% EndTime: 2019-10-10 11:03:13
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (115->18), mult. (117->30), div. (0->0), fcn. (180->8), ass. (0->28)
	t110 = sin(pkin(6));
	t112 = sin(qJ(2));
	t124 = t110 * t112;
	t113 = sin(qJ(1));
	t123 = t110 * t113;
	t114 = cos(qJ(2));
	t122 = t110 * t114;
	t115 = cos(qJ(1));
	t121 = t110 * t115;
	t120 = t113 * t112;
	t119 = t113 * t114;
	t118 = t115 * t112;
	t117 = t115 * t114;
	t111 = cos(pkin(6));
	t103 = t111 * t118 + t119;
	t109 = pkin(12) + qJ(4) + qJ(5);
	t107 = sin(t109);
	t108 = cos(t109);
	t97 = -t103 * t108 + t107 * t121;
	t116 = t103 * t107 + t108 * t121;
	t105 = -t111 * t120 + t117;
	t104 = t111 * t119 + t118;
	t102 = t111 * t117 - t120;
	t101 = -t111 * t107 - t108 * t124;
	t100 = -t107 * t124 + t111 * t108;
	t99 = t105 * t108 + t107 * t123;
	t98 = -t105 * t107 + t108 * t123;
	t1 = [t97, -t104 * t108, 0, t98, t98, 0; t99, t102 * t108, 0, -t116, -t116, 0; 0, t108 * t122, 0, t100, t100, 0; t116, t104 * t107, 0, -t99, -t99, 0; t98, -t102 * t107, 0, t97, t97, 0; 0, -t107 * t122, 0, t101, t101, 0; t102, t105, 0, 0, 0, 0; t104, t103, 0, 0, 0, 0; 0, t124, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:03:13
	% EndTime: 2019-10-10 11:03:14
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (230->35), mult. (270->63), div. (0->0), fcn. (395->10), ass. (0->43)
	t167 = cos(pkin(6));
	t169 = sin(qJ(2));
	t173 = cos(qJ(1));
	t175 = t173 * t169;
	t170 = sin(qJ(1));
	t172 = cos(qJ(2));
	t177 = t170 * t172;
	t158 = t167 * t175 + t177;
	t165 = pkin(12) + qJ(4) + qJ(5);
	t163 = sin(t165);
	t164 = cos(t165);
	t166 = sin(pkin(6));
	t180 = t166 * t173;
	t151 = -t158 * t164 + t163 * t180;
	t174 = t173 * t172;
	t178 = t170 * t169;
	t157 = -t167 * t174 + t178;
	t168 = sin(qJ(6));
	t171 = cos(qJ(6));
	t191 = t151 * t168 + t157 * t171;
	t190 = t151 * t171 - t157 * t168;
	t149 = -t158 * t163 - t164 * t180;
	t189 = t149 * t168;
	t160 = -t167 * t178 + t174;
	t181 = t166 * t170;
	t152 = t160 * t163 - t164 * t181;
	t188 = t152 * t168;
	t182 = t166 * t169;
	t155 = -t163 * t182 + t167 * t164;
	t187 = t155 * t168;
	t184 = t164 * t168;
	t183 = t164 * t171;
	t179 = t168 * t172;
	t176 = t171 * t172;
	t159 = t167 * t177 + t175;
	t156 = t167 * t163 + t164 * t182;
	t154 = t155 * t171;
	t153 = t160 * t164 + t163 * t181;
	t148 = t152 * t171;
	t147 = t149 * t171;
	t146 = t153 * t171 + t159 * t168;
	t145 = -t153 * t168 + t159 * t171;
	t1 = [t190, -t159 * t183 + t160 * t168, 0, -t148, -t148, t145; t146, -t157 * t183 + t158 * t168, 0, t147, t147, t191; 0, (t164 * t176 + t168 * t169) * t166, 0, t154, t154, -t156 * t168 - t166 * t176; -t191, t159 * t184 + t160 * t171, 0, t188, t188, -t146; t145, t157 * t184 + t158 * t171, 0, -t189, -t189, t190; 0, (-t164 * t179 + t169 * t171) * t166, 0, -t187, -t187, -t156 * t171 + t166 * t179; t149, -t159 * t163, 0, t153, t153, 0; t152, -t157 * t163, 0, -t151, -t151, 0; 0, t166 * t172 * t163, 0, t156, t156, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end