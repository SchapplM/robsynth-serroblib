% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRRPRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:31
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRPRR4_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRR4_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_jacobiR_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:31:29
	% EndTime: 2019-10-09 22:31:29
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:31:29
	% EndTime: 2019-10-09 22:31:29
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:31:29
	% EndTime: 2019-10-09 22:31:29
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (5->5), mult. (14->12), div. (0->0), fcn. (24->6), ass. (0->9)
	t22 = cos(pkin(6));
	t23 = sin(qJ(2));
	t26 = t22 * t23;
	t24 = cos(qJ(2));
	t25 = t22 * t24;
	t21 = cos(pkin(11));
	t20 = sin(pkin(6));
	t19 = sin(pkin(11));
	t1 = [0, -t19 * t25 - t21 * t23, 0, 0, 0, 0; 0, -t19 * t23 + t21 * t25, 0, 0, 0, 0; 0, t20 * t24, 0, 0, 0, 0; 0, t19 * t26 - t21 * t24, 0, 0, 0, 0; 0, -t19 * t24 - t21 * t26, 0, 0, 0, 0; 0, -t20 * t23, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:31:29
	% EndTime: 2019-10-09 22:31:29
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
	t67 = cos(pkin(11));
	t65 = sin(pkin(11));
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
	% StartTime: 2019-10-09 22:31:29
	% EndTime: 2019-10-09 22:31:29
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (16->10), mult. (57->32), div. (0->0), fcn. (88->8), ass. (0->19)
	t84 = sin(pkin(6));
	t87 = sin(qJ(3));
	t96 = t84 * t87;
	t88 = sin(qJ(2));
	t95 = t84 * t88;
	t89 = cos(qJ(3));
	t94 = t84 * t89;
	t90 = cos(qJ(2));
	t93 = t84 * t90;
	t86 = cos(pkin(6));
	t92 = t86 * t88;
	t91 = t86 * t90;
	t85 = cos(pkin(11));
	t83 = sin(pkin(11));
	t82 = -t83 * t92 + t85 * t90;
	t81 = -t83 * t91 - t85 * t88;
	t80 = t83 * t90 + t85 * t92;
	t79 = -t83 * t88 + t85 * t91;
	t1 = [0, t81 * t89, -t82 * t87 + t83 * t94, 0, 0, 0; 0, t79 * t89, -t80 * t87 - t85 * t94, 0, 0, 0; 0, t89 * t93, t86 * t89 - t87 * t95, 0, 0, 0; 0, t82, 0, 0, 0, 0; 0, t80, 0, 0, 0, 0; 0, t95, 0, 0, 0, 0; 0, t81 * t87, t82 * t89 + t83 * t96, 0, 0, 0; 0, t79 * t87, t80 * t89 - t85 * t96, 0, 0, 0; 0, t87 * t93, t86 * t87 + t88 * t94, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:31:29
	% EndTime: 2019-10-09 22:31:29
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (69->27), mult. (203->47), div. (0->0), fcn. (292->10), ass. (0->36)
	t94 = sin(pkin(6));
	t98 = sin(qJ(3));
	t117 = t94 * t98;
	t99 = sin(qJ(2));
	t116 = t94 * t99;
	t96 = cos(pkin(6));
	t115 = t96 * t99;
	t101 = cos(qJ(3));
	t114 = t101 * t94;
	t102 = cos(qJ(2));
	t113 = t102 * t94;
	t93 = sin(pkin(11));
	t112 = t93 * t102;
	t95 = cos(pkin(11));
	t111 = t95 * t102;
	t100 = cos(qJ(5));
	t86 = t95 * t115 + t112;
	t81 = t95 * t114 + t86 * t98;
	t82 = t86 * t101 - t95 * t117;
	t97 = sin(qJ(5));
	t110 = t82 * t100 + t81 * t97;
	t109 = t81 * t100 - t82 * t97;
	t88 = -t93 * t115 + t111;
	t83 = -t93 * t114 + t88 * t98;
	t84 = t88 * t101 + t93 * t117;
	t108 = t84 * t100 + t83 * t97;
	t107 = t83 * t100 - t84 * t97;
	t89 = -t96 * t101 + t98 * t116;
	t90 = t99 * t114 + t96 * t98;
	t106 = t90 * t100 + t89 * t97;
	t105 = t89 * t100 - t90 * t97;
	t104 = t100 * t98 - t101 * t97;
	t103 = t100 * t101 + t97 * t98;
	t87 = -t96 * t112 - t95 * t99;
	t85 = t96 * t111 - t93 * t99;
	t1 = [0, t103 * t87, -t107, 0, t107, 0; 0, t103 * t85, -t109, 0, t109, 0; 0, t103 * t113, -t105, 0, t105, 0; 0, t104 * t87, t108, 0, -t108, 0; 0, t104 * t85, t110, 0, -t110, 0; 0, t104 * t113, t106, 0, -t106, 0; 0, -t88, 0, 0, 0, 0; 0, -t86, 0, 0, 0, 0; 0, -t116, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:31:30
	% EndTime: 2019-10-09 22:31:30
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (151->39), mult. (430->78), div. (0->0), fcn. (608->12), ass. (0->46)
	t192 = sin(pkin(11));
	t194 = cos(pkin(11));
	t203 = cos(qJ(2));
	t195 = cos(pkin(6));
	t199 = sin(qJ(2));
	t207 = t195 * t199;
	t184 = t192 * t203 + t194 * t207;
	t198 = sin(qJ(3));
	t193 = sin(pkin(6));
	t202 = cos(qJ(3));
	t209 = t193 * t202;
	t178 = t184 * t198 + t194 * t209;
	t211 = t193 * t198;
	t179 = t184 * t202 - t194 * t211;
	t197 = sin(qJ(5));
	t201 = cos(qJ(5));
	t168 = t178 * t201 - t179 * t197;
	t196 = sin(qJ(6));
	t217 = t168 * t196;
	t200 = cos(qJ(6));
	t216 = t168 * t200;
	t186 = -t192 * t207 + t194 * t203;
	t180 = t186 * t198 - t192 * t209;
	t181 = t186 * t202 + t192 * t211;
	t171 = t180 * t201 - t181 * t197;
	t215 = t171 * t196;
	t214 = t171 * t200;
	t210 = t193 * t199;
	t187 = -t195 * t202 + t198 * t210;
	t188 = t195 * t198 + t199 * t209;
	t176 = t187 * t201 - t188 * t197;
	t213 = t176 * t196;
	t212 = t176 * t200;
	t208 = t193 * t203;
	t206 = t195 * t203;
	t169 = t178 * t197 + t179 * t201;
	t172 = t180 * t197 + t181 * t201;
	t177 = t187 * t197 + t188 * t201;
	t205 = t197 * t202 - t198 * t201;
	t204 = t197 * t198 + t201 * t202;
	t185 = -t192 * t206 - t194 * t199;
	t183 = -t192 * t199 + t194 * t206;
	t182 = t204 * t208;
	t174 = t204 * t185;
	t173 = t204 * t183;
	t1 = [0, t174 * t200 - t186 * t196, -t214, 0, t214, -t172 * t196 + t185 * t200; 0, t173 * t200 - t184 * t196, -t216, 0, t216, -t169 * t196 + t183 * t200; 0, t182 * t200 - t196 * t210, -t212, 0, t212, -t177 * t196 + t200 * t208; 0, -t174 * t196 - t186 * t200, t215, 0, -t215, -t172 * t200 - t185 * t196; 0, -t173 * t196 - t184 * t200, t217, 0, -t217, -t169 * t200 - t183 * t196; 0, -t182 * t196 - t200 * t210, t213, 0, -t213, -t177 * t200 - t196 * t208; 0, t205 * t185, -t172, 0, t172, 0; 0, t205 * t183, -t169, 0, t169, 0; 0, t205 * t208, -t177, 0, t177, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end