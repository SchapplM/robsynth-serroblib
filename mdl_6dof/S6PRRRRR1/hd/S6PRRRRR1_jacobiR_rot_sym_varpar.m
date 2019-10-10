% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRRRRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:13
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRRRR1_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR1_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRRR1_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR1_jacobiR_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:13:44
	% EndTime: 2019-10-09 23:13:44
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:13:44
	% EndTime: 2019-10-09 23:13:44
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:13:44
	% EndTime: 2019-10-09 23:13:44
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (5->5), mult. (14->12), div. (0->0), fcn. (24->6), ass. (0->9)
	t22 = cos(pkin(6));
	t23 = sin(qJ(2));
	t26 = t22 * t23;
	t24 = cos(qJ(2));
	t25 = t22 * t24;
	t21 = cos(pkin(12));
	t20 = sin(pkin(6));
	t19 = sin(pkin(12));
	t1 = [0, -t19 * t25 - t21 * t23, 0, 0, 0, 0; 0, -t19 * t23 + t21 * t25, 0, 0, 0, 0; 0, t20 * t24, 0, 0, 0, 0; 0, t19 * t26 - t21 * t24, 0, 0, 0, 0; 0, -t19 * t24 - t21 * t26, 0, 0, 0, 0; 0, -t20 * t23, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:13:45
	% EndTime: 2019-10-09 23:13:45
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
	t67 = cos(pkin(12));
	t65 = sin(pkin(12));
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
	% StartTime: 2019-10-09 23:13:45
	% EndTime: 2019-10-09 23:13:45
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (59->14), mult. (87->32), div. (0->0), fcn. (134->8), ass. (0->26)
	t93 = sin(pkin(12));
	t94 = sin(pkin(6));
	t104 = t93 * t94;
	t95 = cos(pkin(12));
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
	% StartTime: 2019-10-09 23:13:45
	% EndTime: 2019-10-09 23:13:45
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (123->14), mult. (117->32), div. (0->0), fcn. (180->8), ass. (0->26)
	t97 = sin(pkin(12));
	t98 = sin(pkin(6));
	t108 = t97 * t98;
	t99 = cos(pkin(12));
	t107 = t98 * t99;
	t102 = cos(qJ(2));
	t106 = t102 * t98;
	t101 = sin(qJ(2));
	t105 = t98 * t101;
	t100 = cos(pkin(6));
	t104 = t100 * t101;
	t103 = t100 * t102;
	t96 = qJ(3) + qJ(4) + qJ(5);
	t95 = cos(t96);
	t94 = sin(t96);
	t93 = t99 * t102 - t97 * t104;
	t92 = -t99 * t101 - t97 * t103;
	t91 = t97 * t102 + t99 * t104;
	t90 = -t97 * t101 + t99 * t103;
	t89 = -t100 * t94 - t95 * t105;
	t88 = t100 * t95 - t94 * t105;
	t87 = -t94 * t108 - t93 * t95;
	t86 = t95 * t108 - t93 * t94;
	t85 = t94 * t107 - t91 * t95;
	t84 = -t95 * t107 - t91 * t94;
	t1 = [0, t92 * t95, t86, t86, t86, 0; 0, t90 * t95, t84, t84, t84, 0; 0, t95 * t106, t88, t88, t88, 0; 0, -t92 * t94, t87, t87, t87, 0; 0, -t90 * t94, t85, t85, t85, 0; 0, -t94 * t106, t89, t89, t89, 0; 0, t93, 0, 0, 0, 0; 0, t91, 0, 0, 0, 0; 0, t105, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:13:45
	% EndTime: 2019-10-09 23:13:45
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (240->34), mult. (265->65), div. (0->0), fcn. (388->10), ass. (0->37)
	t152 = sin(pkin(12));
	t154 = cos(pkin(12));
	t159 = cos(qJ(2));
	t155 = cos(pkin(6));
	t157 = sin(qJ(2));
	t163 = t155 * t157;
	t145 = t152 * t159 + t154 * t163;
	t151 = qJ(3) + qJ(4) + qJ(5);
	t149 = sin(t151);
	t150 = cos(t151);
	t153 = sin(pkin(6));
	t165 = t153 * t154;
	t137 = -t145 * t149 - t150 * t165;
	t156 = sin(qJ(6));
	t171 = t137 * t156;
	t147 = -t152 * t163 + t154 * t159;
	t166 = t152 * t153;
	t139 = -t147 * t149 + t150 * t166;
	t170 = t139 * t156;
	t164 = t153 * t157;
	t142 = -t149 * t164 + t155 * t150;
	t169 = t142 * t156;
	t168 = t150 * t156;
	t158 = cos(qJ(6));
	t167 = t150 * t158;
	t162 = t155 * t159;
	t161 = t156 * t159;
	t160 = t158 * t159;
	t146 = t152 * t162 + t154 * t157;
	t144 = t152 * t157 - t154 * t162;
	t143 = t155 * t149 + t150 * t164;
	t141 = t142 * t158;
	t140 = t147 * t150 + t149 * t166;
	t138 = t145 * t150 - t149 * t165;
	t136 = t139 * t158;
	t135 = t137 * t158;
	t1 = [0, -t146 * t167 + t147 * t156, t136, t136, t136, -t140 * t156 + t146 * t158; 0, -t144 * t167 + t145 * t156, t135, t135, t135, -t138 * t156 + t144 * t158; 0, (t150 * t160 + t156 * t157) * t153, t141, t141, t141, -t143 * t156 - t153 * t160; 0, t146 * t168 + t147 * t158, -t170, -t170, -t170, -t140 * t158 - t146 * t156; 0, t144 * t168 + t145 * t158, -t171, -t171, -t171, -t138 * t158 - t144 * t156; 0, (-t150 * t161 + t157 * t158) * t153, -t169, -t169, -t169, -t143 * t158 + t153 * t161; 0, -t146 * t149, t140, t140, t140, 0; 0, -t144 * t149, t138, t138, t138, 0; 0, t153 * t159 * t149, t143, t143, t143, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end