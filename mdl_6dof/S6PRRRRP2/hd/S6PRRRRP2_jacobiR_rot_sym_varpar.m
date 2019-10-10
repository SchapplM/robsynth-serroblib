% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRRRRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:04
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRRRP2_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP2_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRRP2_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP2_jacobiR_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:03:56
	% EndTime: 2019-10-09 23:03:56
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:03:56
	% EndTime: 2019-10-09 23:03:56
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:03:56
	% EndTime: 2019-10-09 23:03:56
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
	% StartTime: 2019-10-09 23:03:56
	% EndTime: 2019-10-09 23:03:56
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
	% StartTime: 2019-10-09 23:03:56
	% EndTime: 2019-10-09 23:03:56
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (59->14), mult. (87->32), div. (0->0), fcn. (134->8), ass. (0->26)
	t93 = sin(pkin(11));
	t94 = sin(pkin(6));
	t104 = t93 * t94;
	t95 = cos(pkin(11));
	t103 = t94 * t95;
	t97 = sin(qJ(2));
	t102 = t94 * t97;
	t98 = cos(qJ(2));
	t101 = t94 * t98;
	t96 = cos(pkin(6));
	t100 = t96 * t97;
	t99 = t96 * t98;
	t92 = qJ(3) + qJ(4);
	t91 = cos(t92);
	t90 = sin(t92);
	t89 = -t93 * t100 + t95 * t98;
	t88 = -t93 * t99 - t95 * t97;
	t87 = t95 * t100 + t93 * t98;
	t86 = -t93 * t97 + t95 * t99;
	t85 = -t91 * t102 - t96 * t90;
	t84 = -t90 * t102 + t96 * t91;
	t83 = -t90 * t104 - t89 * t91;
	t82 = t91 * t104 - t89 * t90;
	t81 = t90 * t103 - t87 * t91;
	t80 = -t91 * t103 - t87 * t90;
	t1 = [0, t88 * t91, t82, t82, 0, 0; 0, t86 * t91, t80, t80, 0, 0; 0, t91 * t101, t84, t84, 0, 0; 0, -t88 * t90, t83, t83, 0, 0; 0, -t86 * t90, t81, t81, 0, 0; 0, -t90 * t101, t85, t85, 0, 0; 0, t89, 0, 0, 0, 0; 0, t87, 0, 0, 0, 0; 0, t102, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:03:56
	% EndTime: 2019-10-09 23:03:57
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (129->31), mult. (214->65), div. (0->0), fcn. (313->10), ass. (0->37)
	t145 = sin(pkin(11));
	t147 = cos(pkin(11));
	t152 = cos(qJ(2));
	t148 = cos(pkin(6));
	t150 = sin(qJ(2));
	t156 = t148 * t150;
	t138 = t145 * t152 + t147 * t156;
	t144 = qJ(3) + qJ(4);
	t142 = sin(t144);
	t143 = cos(t144);
	t146 = sin(pkin(6));
	t158 = t146 * t147;
	t130 = -t138 * t142 - t143 * t158;
	t149 = sin(qJ(5));
	t164 = t130 * t149;
	t140 = -t145 * t156 + t147 * t152;
	t159 = t145 * t146;
	t132 = -t140 * t142 + t143 * t159;
	t163 = t132 * t149;
	t157 = t146 * t150;
	t135 = -t142 * t157 + t148 * t143;
	t162 = t135 * t149;
	t161 = t143 * t149;
	t151 = cos(qJ(5));
	t160 = t143 * t151;
	t155 = t148 * t152;
	t154 = t149 * t152;
	t153 = t151 * t152;
	t139 = t145 * t155 + t147 * t150;
	t137 = t145 * t150 - t147 * t155;
	t136 = t148 * t142 + t143 * t157;
	t134 = t135 * t151;
	t133 = t140 * t143 + t142 * t159;
	t131 = t138 * t143 - t142 * t158;
	t129 = t132 * t151;
	t128 = t130 * t151;
	t1 = [0, -t139 * t160 + t140 * t149, t129, t129, -t133 * t149 + t139 * t151, 0; 0, -t137 * t160 + t138 * t149, t128, t128, -t131 * t149 + t137 * t151, 0; 0, (t143 * t153 + t149 * t150) * t146, t134, t134, -t136 * t149 - t146 * t153, 0; 0, t139 * t161 + t140 * t151, -t163, -t163, -t133 * t151 - t139 * t149, 0; 0, t137 * t161 + t138 * t151, -t164, -t164, -t131 * t151 - t137 * t149, 0; 0, (-t143 * t154 + t150 * t151) * t146, -t162, -t162, -t136 * t151 + t146 * t154, 0; 0, -t139 * t142, t133, t133, 0, 0; 0, -t137 * t142, t131, t131, 0, 0; 0, t146 * t152 * t142, t136, t136, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:03:57
	% EndTime: 2019-10-09 23:03:57
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (123->25), mult. (214->65), div. (0->0), fcn. (313->10), ass. (0->37)
	t186 = qJ(3) + qJ(4);
	t185 = cos(t186);
	t191 = sin(qJ(5));
	t203 = t185 * t191;
	t193 = cos(qJ(5));
	t202 = t185 * t193;
	t187 = sin(pkin(11));
	t188 = sin(pkin(6));
	t201 = t187 * t188;
	t189 = cos(pkin(11));
	t200 = t188 * t189;
	t192 = sin(qJ(2));
	t199 = t188 * t192;
	t190 = cos(pkin(6));
	t198 = t190 * t192;
	t194 = cos(qJ(2));
	t197 = t190 * t194;
	t196 = t191 * t194;
	t195 = t193 * t194;
	t184 = sin(t186);
	t182 = -t187 * t198 + t189 * t194;
	t181 = t187 * t197 + t189 * t192;
	t180 = t187 * t194 + t189 * t198;
	t179 = t187 * t192 - t189 * t197;
	t178 = t190 * t184 + t185 * t199;
	t177 = -t184 * t199 + t190 * t185;
	t176 = t177 * t193;
	t175 = t177 * t191;
	t174 = t182 * t185 + t184 * t201;
	t173 = -t182 * t184 + t185 * t201;
	t172 = t180 * t185 - t184 * t200;
	t171 = -t180 * t184 - t185 * t200;
	t170 = t173 * t193;
	t169 = t173 * t191;
	t168 = t171 * t193;
	t167 = t171 * t191;
	t1 = [0, -t181 * t202 + t182 * t191, t170, t170, -t174 * t191 + t181 * t193, 0; 0, -t179 * t202 + t180 * t191, t168, t168, -t172 * t191 + t179 * t193, 0; 0, (t185 * t195 + t191 * t192) * t188, t176, t176, -t178 * t191 - t188 * t195, 0; 0, -t181 * t184, t174, t174, 0, 0; 0, -t179 * t184, t172, t172, 0, 0; 0, t188 * t194 * t184, t178, t178, 0, 0; 0, -t181 * t203 - t182 * t193, t169, t169, t174 * t193 + t181 * t191, 0; 0, -t179 * t203 - t180 * t193, t167, t167, t172 * t193 + t179 * t191, 0; 0, (t185 * t196 - t192 * t193) * t188, t175, t175, t178 * t193 - t188 * t196, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end