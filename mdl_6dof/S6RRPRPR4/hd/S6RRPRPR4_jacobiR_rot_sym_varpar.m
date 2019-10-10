% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPRPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:09
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRPR4_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR4_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR4_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR4_jacobiR_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:09:41
	% EndTime: 2019-10-10 10:09:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:09:41
	% EndTime: 2019-10-10 10:09:41
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
	% StartTime: 2019-10-10 10:09:41
	% EndTime: 2019-10-10 10:09:41
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
	% StartTime: 2019-10-10 10:09:41
	% EndTime: 2019-10-10 10:09:42
	% DurationCPUTime: 0.05s
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
	% StartTime: 2019-10-10 10:09:42
	% EndTime: 2019-10-10 10:09:42
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (66->18), mult. (182->36), div. (0->0), fcn. (266->10), ass. (0->27)
	t115 = sin(pkin(6));
	t120 = sin(qJ(1));
	t128 = t115 * t120;
	t123 = cos(qJ(1));
	t127 = t115 * t123;
	t117 = cos(pkin(6));
	t114 = sin(pkin(11));
	t116 = cos(pkin(11));
	t119 = sin(qJ(2));
	t122 = cos(qJ(2));
	t125 = t122 * t114 + t119 * t116;
	t109 = t125 * t117;
	t110 = t119 * t114 - t122 * t116;
	t103 = t123 * t109 - t120 * t110;
	t118 = sin(qJ(4));
	t121 = cos(qJ(4));
	t126 = -t103 * t121 + t118 * t127;
	t105 = -t120 * t109 - t123 * t110;
	t124 = t103 * t118 + t121 * t127;
	t108 = t110 * t117;
	t107 = t125 * t115;
	t106 = t110 * t115;
	t104 = t120 * t108 - t123 * t125;
	t102 = -t123 * t108 - t120 * t125;
	t101 = t105 * t121 + t118 * t128;
	t100 = -t105 * t118 + t121 * t128;
	t1 = [t126, t104 * t121, 0, t100, 0, 0; t101, t102 * t121, 0, -t124, 0, 0; 0, -t106 * t121, 0, -t107 * t118 + t117 * t121, 0, 0; t124, -t104 * t118, 0, -t101, 0, 0; t100, -t102 * t118, 0, t126, 0, 0; 0, t106 * t118, 0, -t107 * t121 - t117 * t118, 0, 0; t102, t105, 0, 0, 0, 0; -t104, t103, 0, 0, 0, 0; 0, t107, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:09:42
	% EndTime: 2019-10-10 10:09:42
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (92->19), mult. (182->36), div. (0->0), fcn. (266->10), ass. (0->28)
	t127 = sin(pkin(6));
	t131 = sin(qJ(1));
	t138 = t127 * t131;
	t133 = cos(qJ(1));
	t137 = t127 * t133;
	t129 = cos(pkin(6));
	t126 = sin(pkin(11));
	t128 = cos(pkin(11));
	t130 = sin(qJ(2));
	t132 = cos(qJ(2));
	t135 = t132 * t126 + t130 * t128;
	t118 = t135 * t129;
	t119 = t130 * t126 - t132 * t128;
	t112 = t133 * t118 - t131 * t119;
	t125 = qJ(4) + pkin(12);
	t123 = sin(t125);
	t124 = cos(t125);
	t136 = -t112 * t124 + t123 * t137;
	t114 = -t131 * t118 - t133 * t119;
	t134 = t112 * t123 + t124 * t137;
	t117 = t119 * t129;
	t116 = t135 * t127;
	t115 = t119 * t127;
	t113 = t131 * t117 - t133 * t135;
	t111 = -t133 * t117 - t131 * t135;
	t110 = t114 * t124 + t123 * t138;
	t109 = -t114 * t123 + t124 * t138;
	t1 = [t136, t113 * t124, 0, t109, 0, 0; t110, t111 * t124, 0, -t134, 0, 0; 0, -t115 * t124, 0, -t116 * t123 + t129 * t124, 0, 0; t134, -t113 * t123, 0, -t110, 0, 0; t109, -t111 * t123, 0, t136, 0, 0; 0, t115 * t123, 0, -t116 * t124 - t129 * t123, 0, 0; t111, t114, 0, 0, 0, 0; -t113, t112, 0, 0, 0, 0; 0, t116, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:09:42
	% EndTime: 2019-10-10 10:09:42
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (205->32), mult. (425->65), div. (0->0), fcn. (606->12), ass. (0->38)
	t177 = qJ(4) + pkin(12);
	t175 = sin(t177);
	t176 = cos(t177);
	t181 = cos(pkin(6));
	t178 = sin(pkin(11));
	t180 = cos(pkin(11));
	t183 = sin(qJ(2));
	t186 = cos(qJ(2));
	t189 = t186 * t178 + t183 * t180;
	t169 = t189 * t181;
	t170 = t183 * t178 - t186 * t180;
	t184 = sin(qJ(1));
	t187 = cos(qJ(1));
	t191 = t187 * t169 - t184 * t170;
	t179 = sin(pkin(6));
	t194 = t179 * t187;
	t155 = t175 * t194 - t176 * t191;
	t188 = t170 * t181;
	t159 = -t184 * t189 - t187 * t188;
	t182 = sin(qJ(6));
	t185 = cos(qJ(6));
	t201 = t155 * t182 - t159 * t185;
	t200 = t155 * t185 + t159 * t182;
	t197 = t176 * t182;
	t196 = t176 * t185;
	t195 = t179 * t184;
	t190 = -t184 * t169 - t187 * t170;
	t153 = -t175 * t191 - t176 * t194;
	t168 = t189 * t179;
	t167 = t170 * t179;
	t165 = t168 * t176 + t181 * t175;
	t164 = -t168 * t175 + t181 * t176;
	t162 = t184 * t188 - t187 * t189;
	t157 = t175 * t195 + t176 * t190;
	t156 = t175 * t190 - t176 * t195;
	t152 = t157 * t185 - t162 * t182;
	t151 = -t157 * t182 - t162 * t185;
	t1 = [t200, t162 * t196 + t182 * t190, 0, -t156 * t185, 0, t151; t152, t159 * t196 + t182 * t191, 0, t153 * t185, 0, t201; 0, -t167 * t196 + t168 * t182, 0, t164 * t185, 0, -t165 * t182 + t167 * t185; -t201, -t162 * t197 + t185 * t190, 0, t156 * t182, 0, -t152; t151, -t159 * t197 + t185 * t191, 0, -t153 * t182, 0, t200; 0, t167 * t197 + t168 * t185, 0, -t164 * t182, 0, -t165 * t185 - t167 * t182; t153, t162 * t175, 0, t157, 0, 0; t156, t159 * t175, 0, -t155, 0, 0; 0, -t167 * t175, 0, t165, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end