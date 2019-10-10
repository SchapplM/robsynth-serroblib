% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRRRPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:48
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRRPR2_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR2_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPR2_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR2_jacobiR_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:48:33
	% EndTime: 2019-10-09 22:48:33
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:48:33
	% EndTime: 2019-10-09 22:48:33
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:48:33
	% EndTime: 2019-10-09 22:48:34
	% DurationCPUTime: 0.04s
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
	% StartTime: 2019-10-09 22:48:34
	% EndTime: 2019-10-09 22:48:34
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
	% StartTime: 2019-10-09 22:48:34
	% EndTime: 2019-10-09 22:48:34
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
	% StartTime: 2019-10-09 22:48:34
	% EndTime: 2019-10-09 22:48:34
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (97->23), mult. (158->52), div. (0->0), fcn. (231->10), ass. (0->36)
	t134 = sin(pkin(11));
	t137 = cos(pkin(11));
	t140 = cos(qJ(2));
	t138 = cos(pkin(6));
	t139 = sin(qJ(2));
	t142 = t138 * t139;
	t127 = t134 * t140 + t137 * t142;
	t132 = qJ(3) + qJ(4);
	t130 = sin(t132);
	t131 = cos(t132);
	t135 = sin(pkin(6));
	t144 = t135 * t137;
	t119 = -t127 * t130 - t131 * t144;
	t133 = sin(pkin(12));
	t151 = t119 * t133;
	t129 = -t134 * t142 + t137 * t140;
	t145 = t134 * t135;
	t121 = -t129 * t130 + t131 * t145;
	t150 = t121 * t133;
	t143 = t135 * t139;
	t124 = -t130 * t143 + t138 * t131;
	t149 = t124 * t133;
	t148 = t131 * t133;
	t136 = cos(pkin(12));
	t147 = t131 * t136;
	t146 = t131 * t140;
	t141 = t138 * t140;
	t128 = -t134 * t141 - t137 * t139;
	t126 = -t134 * t139 + t137 * t141;
	t125 = t138 * t130 + t131 * t143;
	t123 = t124 * t136;
	t122 = t129 * t131 + t130 * t145;
	t120 = t127 * t131 - t130 * t144;
	t118 = t121 * t136;
	t117 = t119 * t136;
	t1 = [0, t128 * t147 + t129 * t133, t118, t118, 0, 0; 0, t126 * t147 + t127 * t133, t117, t117, 0, 0; 0, (t133 * t139 + t136 * t146) * t135, t123, t123, 0, 0; 0, -t128 * t148 + t129 * t136, -t150, -t150, 0, 0; 0, -t126 * t148 + t127 * t136, -t151, -t151, 0, 0; 0, (-t133 * t146 + t136 * t139) * t135, -t149, -t149, 0, 0; 0, t128 * t130, t122, t122, 0, 0; 0, t126 * t130, t120, t120, 0, 0; 0, t135 * t140 * t130, t125, t125, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:48:34
	% EndTime: 2019-10-09 22:48:34
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (165->32), mult. (214->64), div. (0->0), fcn. (313->10), ass. (0->38)
	t155 = sin(pkin(11));
	t157 = cos(pkin(11));
	t160 = cos(qJ(2));
	t158 = cos(pkin(6));
	t159 = sin(qJ(2));
	t162 = t158 * t159;
	t145 = t155 * t160 + t157 * t162;
	t154 = qJ(3) + qJ(4);
	t151 = sin(t154);
	t152 = cos(t154);
	t156 = sin(pkin(6));
	t165 = t156 * t157;
	t137 = -t145 * t151 - t152 * t165;
	t153 = pkin(12) + qJ(6);
	t149 = sin(t153);
	t172 = t137 * t149;
	t147 = -t155 * t162 + t157 * t160;
	t166 = t155 * t156;
	t139 = -t147 * t151 + t152 * t166;
	t171 = t139 * t149;
	t164 = t156 * t159;
	t142 = -t151 * t164 + t152 * t158;
	t170 = t142 * t149;
	t169 = t149 * t152;
	t150 = cos(t153);
	t168 = t150 * t152;
	t167 = t152 * t160;
	t163 = t156 * t160;
	t161 = t158 * t160;
	t146 = t155 * t161 + t157 * t159;
	t144 = t155 * t159 - t157 * t161;
	t143 = t151 * t158 + t152 * t164;
	t141 = t142 * t150;
	t140 = t147 * t152 + t151 * t166;
	t138 = t145 * t152 - t151 * t165;
	t136 = t139 * t150;
	t135 = t137 * t150;
	t1 = [0, -t146 * t168 + t147 * t149, t136, t136, 0, -t140 * t149 + t146 * t150; 0, -t144 * t168 + t145 * t149, t135, t135, 0, -t138 * t149 + t144 * t150; 0, (t149 * t159 + t150 * t167) * t156, t141, t141, 0, -t143 * t149 - t150 * t163; 0, t146 * t169 + t147 * t150, -t171, -t171, 0, -t140 * t150 - t146 * t149; 0, t144 * t169 + t145 * t150, -t172, -t172, 0, -t138 * t150 - t144 * t149; 0, (-t149 * t167 + t150 * t159) * t156, -t170, -t170, 0, -t143 * t150 + t149 * t163; 0, -t146 * t151, t140, t140, 0, 0; 0, -t144 * t151, t138, t138, 0, 0; 0, t151 * t163, t143, t143, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end