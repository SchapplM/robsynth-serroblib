% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRPRR7
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
% Datum: 2019-10-10 12:04
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPRR7_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR7_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR7_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR7_jacobiR_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:04:12
	% EndTime: 2019-10-10 12:04:12
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:04:12
	% EndTime: 2019-10-10 12:04:12
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
	% StartTime: 2019-10-10 12:04:12
	% EndTime: 2019-10-10 12:04:12
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
	% StartTime: 2019-10-10 12:04:12
	% EndTime: 2019-10-10 12:04:12
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
	% StartTime: 2019-10-10 12:04:12
	% EndTime: 2019-10-10 12:04:12
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (55->16), mult. (87->30), div. (0->0), fcn. (134->8), ass. (0->26)
	t91 = sin(pkin(6));
	t93 = sin(qJ(2));
	t106 = t91 * t93;
	t94 = sin(qJ(1));
	t105 = t91 * t94;
	t95 = cos(qJ(2));
	t104 = t91 * t95;
	t96 = cos(qJ(1));
	t103 = t91 * t96;
	t102 = t94 * t93;
	t101 = t94 * t95;
	t100 = t96 * t93;
	t99 = t96 * t95;
	t92 = cos(pkin(6));
	t84 = t92 * t100 + t101;
	t90 = qJ(3) + pkin(12);
	t88 = sin(t90);
	t89 = cos(t90);
	t98 = t88 * t103 - t84 * t89;
	t97 = t89 * t103 + t84 * t88;
	t86 = -t92 * t102 + t99;
	t85 = t92 * t101 + t100;
	t83 = t92 * t99 - t102;
	t82 = t88 * t105 + t86 * t89;
	t81 = t89 * t105 - t86 * t88;
	t1 = [t98, -t85 * t89, t81, 0, 0, 0; t82, t83 * t89, -t97, 0, 0, 0; 0, t89 * t104, -t88 * t106 + t92 * t89, 0, 0, 0; t97, t85 * t88, -t82, 0, 0, 0; t81, -t83 * t88, t98, 0, 0, 0; 0, -t88 * t104, -t89 * t106 - t92 * t88, 0, 0, 0; t83, t86, 0, 0, 0, 0; t85, t84, 0, 0, 0, 0; 0, t106, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:04:12
	% EndTime: 2019-10-10 12:04:12
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (115->18), mult. (117->30), div. (0->0), fcn. (180->8), ass. (0->28)
	t111 = sin(pkin(6));
	t113 = sin(qJ(2));
	t125 = t111 * t113;
	t114 = sin(qJ(1));
	t124 = t111 * t114;
	t115 = cos(qJ(2));
	t123 = t111 * t115;
	t116 = cos(qJ(1));
	t122 = t111 * t116;
	t121 = t114 * t113;
	t120 = t114 * t115;
	t119 = t116 * t113;
	t118 = t116 * t115;
	t112 = cos(pkin(6));
	t104 = t112 * t119 + t120;
	t110 = qJ(3) + pkin(12) + qJ(5);
	t108 = sin(t110);
	t109 = cos(t110);
	t98 = -t104 * t109 + t108 * t122;
	t117 = t104 * t108 + t109 * t122;
	t106 = -t112 * t121 + t118;
	t105 = t112 * t120 + t119;
	t103 = t112 * t118 - t121;
	t102 = -t112 * t108 - t109 * t125;
	t101 = -t108 * t125 + t112 * t109;
	t100 = t106 * t109 + t108 * t124;
	t99 = -t106 * t108 + t109 * t124;
	t1 = [t98, -t105 * t109, t99, 0, t99, 0; t100, t103 * t109, -t117, 0, -t117, 0; 0, t109 * t123, t101, 0, t101, 0; t117, t105 * t108, -t100, 0, -t100, 0; t99, -t103 * t108, t98, 0, t98, 0; 0, -t108 * t123, t102, 0, t102, 0; t103, t106, 0, 0, 0, 0; t105, t104, 0, 0, 0, 0; 0, t125, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:04:12
	% EndTime: 2019-10-10 12:04:13
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (230->35), mult. (270->63), div. (0->0), fcn. (395->10), ass. (0->43)
	t168 = cos(pkin(6));
	t170 = sin(qJ(2));
	t174 = cos(qJ(1));
	t176 = t174 * t170;
	t171 = sin(qJ(1));
	t173 = cos(qJ(2));
	t178 = t171 * t173;
	t159 = t168 * t176 + t178;
	t166 = qJ(3) + pkin(12) + qJ(5);
	t164 = sin(t166);
	t165 = cos(t166);
	t167 = sin(pkin(6));
	t181 = t167 * t174;
	t152 = -t159 * t165 + t164 * t181;
	t175 = t174 * t173;
	t179 = t171 * t170;
	t158 = -t168 * t175 + t179;
	t169 = sin(qJ(6));
	t172 = cos(qJ(6));
	t192 = t152 * t169 + t158 * t172;
	t191 = t152 * t172 - t158 * t169;
	t150 = -t159 * t164 - t165 * t181;
	t190 = t150 * t169;
	t161 = -t168 * t179 + t175;
	t182 = t167 * t171;
	t153 = t161 * t164 - t165 * t182;
	t189 = t153 * t169;
	t183 = t167 * t170;
	t156 = -t164 * t183 + t168 * t165;
	t188 = t156 * t169;
	t185 = t165 * t169;
	t184 = t165 * t172;
	t180 = t169 * t173;
	t177 = t172 * t173;
	t160 = t168 * t178 + t176;
	t157 = t168 * t164 + t165 * t183;
	t155 = t156 * t172;
	t154 = t161 * t165 + t164 * t182;
	t149 = t153 * t172;
	t148 = t150 * t172;
	t147 = t154 * t172 + t160 * t169;
	t146 = -t154 * t169 + t160 * t172;
	t1 = [t191, -t160 * t184 + t161 * t169, -t149, 0, -t149, t146; t147, -t158 * t184 + t159 * t169, t148, 0, t148, t192; 0, (t165 * t177 + t169 * t170) * t167, t155, 0, t155, -t157 * t169 - t167 * t177; -t192, t160 * t185 + t161 * t172, t189, 0, t189, -t147; t146, t158 * t185 + t159 * t172, -t190, 0, -t190, t191; 0, (-t165 * t180 + t170 * t172) * t167, -t188, 0, -t188, -t157 * t172 + t167 * t180; t150, -t160 * t164, t154, 0, t154, 0; t153, -t158 * t164, -t152, 0, -t152, 0; 0, t167 * t173 * t164, t157, 0, t157, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end