% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3,theta4]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:39
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPPRR3_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR3_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_jacobiR_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:39:18
	% EndTime: 2019-10-10 09:39:18
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:39:18
	% EndTime: 2019-10-10 09:39:18
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
	% StartTime: 2019-10-10 09:39:18
	% EndTime: 2019-10-10 09:39:18
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
	% StartTime: 2019-10-10 09:39:19
	% EndTime: 2019-10-10 09:39:19
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (26->10), mult. (74->18), div. (0->0), fcn. (112->8), ass. (0->17)
	t70 = cos(pkin(6));
	t67 = sin(pkin(11));
	t69 = cos(pkin(11));
	t71 = sin(qJ(2));
	t73 = cos(qJ(2));
	t75 = t73 * t67 + t71 * t69;
	t63 = t75 * t70;
	t64 = t71 * t67 - t73 * t69;
	t72 = sin(qJ(1));
	t74 = cos(qJ(1));
	t77 = -t74 * t63 + t72 * t64;
	t76 = t72 * t63 + t74 * t64;
	t68 = sin(pkin(6));
	t62 = t64 * t70;
	t61 = t72 * t62 - t74 * t75;
	t60 = -t74 * t62 - t72 * t75;
	t1 = [t77, t61, 0, 0, 0, 0; -t76, t60, 0, 0, 0, 0; 0, -t64 * t68, 0, 0, 0, 0; -t60, t76, 0, 0, 0, 0; t61, t77, 0, 0, 0, 0; 0, -t75 * t68, 0, 0, 0, 0; t74 * t68, 0, 0, 0, 0, 0; t72 * t68, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:39:19
	% EndTime: 2019-10-10 09:39:19
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (46->15), mult. (126->32), div. (0->0), fcn. (184->10), ass. (0->22)
	t100 = sin(pkin(6));
	t105 = sin(qJ(1));
	t111 = t100 * t105;
	t107 = cos(qJ(1));
	t110 = t100 * t107;
	t102 = cos(pkin(11));
	t104 = sin(qJ(2));
	t106 = cos(qJ(2));
	t99 = sin(pkin(11));
	t109 = t106 * t102 - t104 * t99;
	t103 = cos(pkin(6));
	t108 = t104 * t102 + t106 * t99;
	t94 = t108 * t103;
	t89 = -t105 * t109 - t107 * t94;
	t91 = -t105 * t94 + t107 * t109;
	t101 = cos(pkin(12));
	t98 = sin(pkin(12));
	t93 = t109 * t103;
	t92 = t109 * t100;
	t90 = -t105 * t93 - t107 * t108;
	t88 = -t105 * t108 + t107 * t93;
	t1 = [t89 * t101 + t98 * t110, t90 * t101, 0, 0, 0, 0; t91 * t101 + t98 * t111, t88 * t101, 0, 0, 0, 0; 0, t92 * t101, 0, 0, 0, 0; t101 * t110 - t89 * t98, -t90 * t98, 0, 0, 0, 0; t101 * t111 - t91 * t98, -t88 * t98, 0, 0, 0, 0; 0, -t92 * t98, 0, 0, 0, 0; t88, t91, 0, 0, 0, 0; -t90, -t89, 0, 0, 0, 0; 0, t108 * t100, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:39:19
	% EndTime: 2019-10-10 09:39:19
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (92->19), mult. (182->36), div. (0->0), fcn. (266->10), ass. (0->28)
	t126 = sin(pkin(6));
	t130 = sin(qJ(1));
	t137 = t126 * t130;
	t132 = cos(qJ(1));
	t136 = t126 * t132;
	t128 = cos(pkin(6));
	t125 = sin(pkin(11));
	t127 = cos(pkin(11));
	t129 = sin(qJ(2));
	t131 = cos(qJ(2));
	t134 = t131 * t125 + t129 * t127;
	t117 = t134 * t128;
	t118 = t129 * t125 - t131 * t127;
	t111 = t132 * t117 - t130 * t118;
	t124 = pkin(12) + qJ(5);
	t122 = sin(t124);
	t123 = cos(t124);
	t135 = -t111 * t123 + t122 * t136;
	t113 = -t130 * t117 - t132 * t118;
	t133 = t111 * t122 + t123 * t136;
	t116 = t118 * t128;
	t115 = t134 * t126;
	t114 = t118 * t126;
	t112 = t130 * t116 - t132 * t134;
	t110 = -t132 * t116 - t130 * t134;
	t109 = t113 * t123 + t122 * t137;
	t108 = -t113 * t122 + t123 * t137;
	t1 = [t135, t112 * t123, 0, 0, t108, 0; t109, t110 * t123, 0, 0, -t133, 0; 0, -t114 * t123, 0, 0, -t115 * t122 + t128 * t123, 0; t133, -t112 * t122, 0, 0, -t109, 0; t108, -t110 * t122, 0, 0, t135, 0; 0, t114 * t122, 0, 0, -t115 * t123 - t128 * t122, 0; t110, t113, 0, 0, 0, 0; -t112, t111, 0, 0, 0, 0; 0, t115, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:39:19
	% EndTime: 2019-10-10 09:39:19
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (205->32), mult. (425->65), div. (0->0), fcn. (606->12), ass. (0->38)
	t176 = pkin(12) + qJ(5);
	t174 = sin(t176);
	t175 = cos(t176);
	t180 = cos(pkin(6));
	t177 = sin(pkin(11));
	t179 = cos(pkin(11));
	t182 = sin(qJ(2));
	t185 = cos(qJ(2));
	t188 = t185 * t177 + t182 * t179;
	t168 = t188 * t180;
	t169 = t182 * t177 - t185 * t179;
	t183 = sin(qJ(1));
	t186 = cos(qJ(1));
	t190 = t186 * t168 - t183 * t169;
	t178 = sin(pkin(6));
	t193 = t178 * t186;
	t154 = t174 * t193 - t175 * t190;
	t187 = t169 * t180;
	t158 = -t183 * t188 - t186 * t187;
	t181 = sin(qJ(6));
	t184 = cos(qJ(6));
	t200 = t154 * t181 - t158 * t184;
	t199 = t154 * t184 + t158 * t181;
	t196 = t175 * t181;
	t195 = t175 * t184;
	t194 = t178 * t183;
	t189 = -t183 * t168 - t186 * t169;
	t152 = -t174 * t190 - t175 * t193;
	t167 = t188 * t178;
	t166 = t169 * t178;
	t164 = t167 * t175 + t180 * t174;
	t163 = -t167 * t174 + t180 * t175;
	t161 = t183 * t187 - t186 * t188;
	t156 = t174 * t194 + t175 * t189;
	t155 = t174 * t189 - t175 * t194;
	t151 = t156 * t184 - t161 * t181;
	t150 = -t156 * t181 - t161 * t184;
	t1 = [t199, t161 * t195 + t181 * t189, 0, 0, -t155 * t184, t150; t151, t158 * t195 + t181 * t190, 0, 0, t152 * t184, t200; 0, -t166 * t195 + t167 * t181, 0, 0, t163 * t184, -t164 * t181 + t166 * t184; -t200, -t161 * t196 + t184 * t189, 0, 0, t155 * t181, -t151; t150, -t158 * t196 + t184 * t190, 0, 0, -t152 * t181, t199; 0, t166 * t196 + t167 * t184, 0, 0, -t163 * t181, -t164 * t184 - t166 * t181; t152, t161 * t174, 0, 0, t156, 0; t155, t158 * t174, 0, 0, -t154, 0; 0, -t166 * t174, 0, 0, t164, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end