% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPRRR14V3
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
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:13
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRRR14V3_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14V3_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR14V3_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S6RRPRRR14V3_jacobiR_rot_sym_varpar: pkin has to be [1x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:13:23
	% EndTime: 2019-10-10 11:13:23
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:13:23
	% EndTime: 2019-10-10 11:13:23
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
	% StartTime: 2019-10-10 11:13:23
	% EndTime: 2019-10-10 11:13:23
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (7->7), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->9)
	t7 = sin(qJ(2));
	t8 = sin(qJ(1));
	t14 = t8 * t7;
	t9 = cos(qJ(2));
	t13 = t8 * t9;
	t10 = cos(qJ(1));
	t12 = t10 * t7;
	t11 = t10 * t9;
	t1 = [-t13, -t12, 0, 0, 0, 0; t11, -t14, 0, 0, 0, 0; 0, t9, 0, 0, 0, 0; t14, -t11, 0, 0, 0, 0; -t12, -t13, 0, 0, 0, 0; 0, -t7, 0, 0, 0, 0; t10, 0, 0, 0, 0, 0; t8, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:13:23
	% EndTime: 2019-10-10 11:13:23
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (4->4), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->9)
	t36 = sin(qJ(2));
	t37 = sin(qJ(1));
	t42 = t37 * t36;
	t38 = cos(qJ(2));
	t41 = t37 * t38;
	t39 = cos(qJ(1));
	t40 = t39 * t36;
	t35 = t39 * t38;
	t1 = [-t41, -t40, 0, 0, 0, 0; t35, -t42, 0, 0, 0, 0; 0, t38, 0, 0, 0, 0; t39, 0, 0, 0, 0, 0; t37, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t42, t35, 0, 0, 0, 0; t40, t41, 0, 0, 0, 0; 0, t36, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:13:23
	% EndTime: 2019-10-10 11:13:23
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (14->12), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->18)
	t59 = sin(qJ(2));
	t60 = sin(qJ(1));
	t70 = t60 * t59;
	t61 = cos(qJ(4));
	t69 = t60 * t61;
	t58 = sin(qJ(4));
	t62 = cos(qJ(2));
	t68 = t62 * t58;
	t67 = t62 * t61;
	t63 = cos(qJ(1));
	t66 = t63 * t59;
	t65 = t63 * t61;
	t64 = t63 * t62;
	t57 = t60 * t58 + t61 * t64;
	t56 = -t58 * t64 + t69;
	t55 = t63 * t58 - t60 * t67;
	t54 = t60 * t68 + t65;
	t1 = [t55, -t59 * t65, 0, t56, 0, 0; t57, -t59 * t69, 0, -t54, 0, 0; 0, t67, 0, -t59 * t58, 0, 0; t54, t58 * t66, 0, -t57, 0, 0; t56, t58 * t70, 0, t55, 0, 0; 0, -t68, 0, -t59 * t61, 0, 0; -t70, t64, 0, 0, 0, 0; t66, t60 * t62, 0, 0, 0, 0; 0, t59, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:13:23
	% EndTime: 2019-10-10 11:13:23
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (37->20), mult. (118->41), div. (0->0), fcn. (180->8), ass. (0->27)
	t94 = sin(qJ(5));
	t96 = sin(qJ(2));
	t113 = t96 * t94;
	t98 = cos(qJ(5));
	t112 = t96 * t98;
	t99 = cos(qJ(4));
	t111 = t96 * t99;
	t95 = sin(qJ(4));
	t97 = sin(qJ(1));
	t110 = t97 * t95;
	t100 = cos(qJ(2));
	t109 = t100 * t95;
	t108 = t100 * t99;
	t101 = cos(qJ(1));
	t107 = t101 * t96;
	t106 = t100 * t101;
	t89 = -t101 * t95 + t97 * t108;
	t105 = t97 * t112 - t89 * t94;
	t104 = -t97 * t113 - t89 * t98;
	t103 = t100 * t94 - t98 * t111;
	t102 = t100 * t98 + t94 * t111;
	t91 = t99 * t106 + t110;
	t90 = t95 * t106 - t97 * t99;
	t88 = -t101 * t99 - t97 * t109;
	t87 = t94 * t107 + t91 * t98;
	t86 = t98 * t107 - t91 * t94;
	t1 = [t104, t103 * t101, 0, -t90 * t98, t86, 0; t87, t103 * t97, 0, t88 * t98, t105, 0; 0, t98 * t108 + t113, 0, -t95 * t112, -t102, 0; -t105, t102 * t101, 0, t90 * t94, -t87, 0; t86, t102 * t97, 0, -t88 * t94, t104, 0; 0, -t94 * t108 + t112, 0, t95 * t113, t103, 0; t88, -t95 * t107, 0, t91, 0, 0; t90, -t96 * t110, 0, t89, 0, 0; 0, t109, 0, t111, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:13:23
	% EndTime: 2019-10-10 11:13:23
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (90->38), mult. (280->83), div. (0->0), fcn. (405->10), ass. (0->43)
	t155 = sin(qJ(1));
	t153 = sin(qJ(4));
	t160 = cos(qJ(1));
	t164 = t160 * t153;
	t158 = cos(qJ(4));
	t159 = cos(qJ(2));
	t166 = t158 * t159;
	t143 = t155 * t166 - t164;
	t157 = cos(qJ(5));
	t152 = sin(qJ(5));
	t154 = sin(qJ(2));
	t174 = t152 * t154;
	t134 = t143 * t157 + t155 * t174;
	t165 = t158 * t160;
	t168 = t155 * t153;
	t142 = t159 * t168 + t165;
	t151 = sin(qJ(6));
	t156 = cos(qJ(6));
	t179 = t134 * t151 - t142 * t156;
	t178 = -t134 * t156 - t142 * t151;
	t175 = t151 * t157;
	t173 = t153 * t154;
	t172 = t153 * t159;
	t171 = t154 * t157;
	t170 = t154 * t158;
	t169 = t154 * t160;
	t167 = t156 * t157;
	t163 = t151 * t173;
	t162 = t156 * t173;
	t161 = t154 * t164;
	t133 = -t143 * t152 + t155 * t171;
	t141 = -t152 * t159 + t157 * t170;
	t140 = -t152 * t170 - t157 * t159;
	t146 = t159 * t165 + t168;
	t145 = -t155 * t158 + t159 * t164;
	t144 = t157 * t166 + t174;
	t139 = t141 * t160;
	t138 = t141 * t155;
	t137 = t146 * t157 + t152 * t169;
	t136 = t146 * t152 - t157 * t169;
	t132 = t137 * t156 + t145 * t151;
	t131 = -t137 * t151 + t145 * t156;
	t1 = [t178, -t139 * t156 - t151 * t161, 0, -t145 * t167 + t146 * t151, -t136 * t156, t131; t132, -t138 * t156 - t155 * t163, 0, -t142 * t167 + t143 * t151, t133 * t156, -t179; 0, t144 * t156 + t151 * t172, 0, (t151 * t158 - t153 * t167) * t154, t140 * t156, -t141 * t151 + t162; t179, t139 * t151 - t156 * t161, 0, t145 * t175 + t146 * t156, t136 * t151, -t132; t131, t138 * t151 - t155 * t162, 0, t142 * t175 + t143 * t156, -t133 * t151, t178; 0, -t144 * t151 + t156 * t172, 0, (t153 * t175 + t156 * t158) * t154, -t140 * t151, -t141 * t156 - t163; t133, t140 * t160, 0, -t145 * t152, t137, 0; t136, t140 * t155, 0, -t142 * t152, t134, 0; 0, t152 * t166 - t171, 0, -t152 * t173, t141, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end