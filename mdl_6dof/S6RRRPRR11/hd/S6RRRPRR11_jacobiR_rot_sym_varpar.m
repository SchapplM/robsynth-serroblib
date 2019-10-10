% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:10
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPRR11_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR11_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR11_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR11_jacobiR_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:10:18
	% EndTime: 2019-10-10 12:10:18
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:10:18
	% EndTime: 2019-10-10 12:10:18
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
	% StartTime: 2019-10-10 12:10:18
	% EndTime: 2019-10-10 12:10:18
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
	% StartTime: 2019-10-10 12:10:18
	% EndTime: 2019-10-10 12:10:19
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
	% StartTime: 2019-10-10 12:10:19
	% EndTime: 2019-10-10 12:10:19
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (26->14), mult. (87->31), div. (0->0), fcn. (134->8), ass. (0->25)
	t102 = sin(pkin(6));
	t105 = sin(qJ(2));
	t119 = t102 * t105;
	t107 = cos(qJ(3));
	t118 = t102 * t107;
	t108 = cos(qJ(2));
	t117 = t102 * t108;
	t109 = cos(qJ(1));
	t116 = t102 * t109;
	t106 = sin(qJ(1));
	t115 = t106 * t105;
	t114 = t106 * t108;
	t113 = t109 * t105;
	t112 = t109 * t108;
	t104 = sin(qJ(3));
	t103 = cos(pkin(6));
	t99 = t103 * t113 + t114;
	t111 = -t99 * t104 - t107 * t116;
	t110 = t104 * t116 - t99 * t107;
	t101 = -t103 * t115 + t112;
	t100 = t103 * t114 + t113;
	t98 = t103 * t112 - t115;
	t97 = t106 * t102 * t104 + t101 * t107;
	t96 = t101 * t104 - t106 * t118;
	t1 = [t110, -t100 * t107, -t96, 0, 0, 0; t97, t98 * t107, t111, 0, 0, 0; 0, t107 * t117, t103 * t107 - t104 * t119, 0, 0, 0; t98, t101, 0, 0, 0, 0; t100, t99, 0, 0, 0, 0; 0, t119, 0, 0, 0, 0; t111, -t100 * t104, t97, 0, 0, 0; t96, t98 * t104, -t110, 0, 0, 0; 0, t104 * t117, t103 * t104 + t105 * t118, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:10:19
	% EndTime: 2019-10-10 12:10:19
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (91->32), mult. (265->47), div. (0->0), fcn. (382->10), ass. (0->37)
	t120 = cos(pkin(6));
	t123 = sin(qJ(2));
	t128 = cos(qJ(1));
	t136 = t128 * t123;
	t124 = sin(qJ(1));
	t127 = cos(qJ(2));
	t137 = t124 * t127;
	t113 = t120 * t136 + t137;
	t122 = sin(qJ(3));
	t126 = cos(qJ(3));
	t119 = sin(pkin(6));
	t139 = t119 * t128;
	t104 = t113 * t122 + t126 * t139;
	t107 = -t113 * t126 + t122 * t139;
	t121 = sin(qJ(5));
	t125 = cos(qJ(5));
	t134 = t104 * t125 + t107 * t121;
	t133 = t104 * t121 - t107 * t125;
	t142 = t119 * t123;
	t141 = t119 * t126;
	t140 = t119 * t127;
	t138 = t124 * t123;
	t135 = t128 * t127;
	t115 = -t120 * t138 + t135;
	t108 = t115 * t122 - t124 * t141;
	t109 = t124 * t119 * t122 + t115 * t126;
	t102 = t108 * t125 - t109 * t121;
	t103 = t108 * t121 + t109 * t125;
	t110 = -t120 * t126 + t122 * t142;
	t111 = t120 * t122 + t123 * t141;
	t132 = t110 * t125 - t111 * t121;
	t131 = t110 * t121 + t111 * t125;
	t130 = t121 * t126 - t122 * t125;
	t129 = t121 * t122 + t125 * t126;
	t114 = -t120 * t137 - t136;
	t112 = -t120 * t135 + t138;
	t1 = [-t133, t129 * t114, -t102, 0, t102, 0; t103, -t129 * t112, -t134, 0, t134, 0; 0, t129 * t140, -t132, 0, t132, 0; -t134, -t130 * t114, t103, 0, -t103, 0; t102, t130 * t112, t133, 0, -t133, 0; 0, -t130 * t140, t131, 0, -t131, 0; t112, -t115, 0, 0, 0, 0; t114, -t113, 0, 0, 0, 0; 0, -t142, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:10:19
	% EndTime: 2019-10-10 12:10:20
	% DurationCPUTime: 0.28s
	% Computational Cost: add. (189->42), mult. (534->79), div. (0->0), fcn. (756->12), ass. (0->50)
	t231 = sin(qJ(2));
	t232 = sin(qJ(1));
	t236 = cos(qJ(2));
	t237 = cos(qJ(1));
	t249 = cos(pkin(6));
	t241 = t237 * t249;
	t220 = t231 * t241 + t232 * t236;
	t230 = sin(qJ(3));
	t235 = cos(qJ(3));
	t227 = sin(pkin(6));
	t243 = t227 * t237;
	t210 = t220 * t230 + t235 * t243;
	t213 = -t220 * t235 + t230 * t243;
	t229 = sin(qJ(5));
	t234 = cos(qJ(5));
	t199 = t210 * t229 - t213 * t234;
	t219 = t232 * t231 - t236 * t241;
	t228 = sin(qJ(6));
	t233 = cos(qJ(6));
	t257 = t199 * t228 + t219 * t233;
	t256 = -t199 * t233 + t219 * t228;
	t198 = t210 * t234 + t213 * t229;
	t255 = t198 * t228;
	t254 = t198 * t233;
	t246 = t227 * t231;
	t217 = t230 * t246 - t249 * t235;
	t245 = t227 * t235;
	t218 = t249 * t230 + t231 * t245;
	t207 = t217 * t234 - t218 * t229;
	t253 = t207 * t228;
	t252 = t207 * t233;
	t242 = t232 * t249;
	t222 = -t231 * t242 + t237 * t236;
	t214 = t222 * t230 - t232 * t245;
	t215 = t232 * t227 * t230 + t222 * t235;
	t240 = -t214 * t234 + t215 * t229;
	t251 = t240 * t228;
	t250 = t240 * t233;
	t244 = t227 * t236;
	t203 = t214 * t229 + t215 * t234;
	t208 = t217 * t229 + t218 * t234;
	t239 = t229 * t235 - t230 * t234;
	t238 = t229 * t230 + t234 * t235;
	t221 = -t237 * t231 - t236 * t242;
	t216 = t238 * t244;
	t205 = t238 * t221;
	t204 = t238 * t219;
	t196 = t203 * t233 + t221 * t228;
	t195 = -t203 * t228 + t221 * t233;
	t1 = [t256, t205 * t233 - t222 * t228, t250, 0, -t250, t195; t196, -t204 * t233 - t220 * t228, -t254, 0, t254, -t257; 0, t216 * t233 - t228 * t246, -t252, 0, t252, -t208 * t228 + t233 * t244; t257, -t205 * t228 - t222 * t233, -t251, 0, t251, -t196; t195, t204 * t228 - t220 * t233, t255, 0, -t255, t256; 0, -t216 * t228 - t233 * t246, t253, 0, -t253, -t208 * t233 - t228 * t244; t198, t239 * t221, -t203, 0, t203, 0; t240, -t239 * t219, -t199, 0, t199, 0; 0, t239 * t244, -t208, 0, t208, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end