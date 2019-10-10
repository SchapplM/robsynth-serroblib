% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PPRRPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:12
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PPRRPR2_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPRRPR2_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_jacobiR_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:12:33
	% EndTime: 2019-10-09 21:12:33
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:12:33
	% EndTime: 2019-10-09 21:12:33
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:12:33
	% EndTime: 2019-10-09 21:12:33
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:12:33
	% EndTime: 2019-10-09 21:12:34
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (20->14), mult. (62->33), div. (0->0), fcn. (88->10), ass. (0->20)
	t53 = sin(pkin(11));
	t59 = cos(pkin(6));
	t68 = t53 * t59;
	t54 = sin(pkin(7));
	t55 = sin(pkin(6));
	t67 = t54 * t55;
	t66 = t54 * t59;
	t56 = cos(pkin(12));
	t58 = cos(pkin(7));
	t65 = t56 * t58;
	t57 = cos(pkin(11));
	t64 = t57 * t59;
	t52 = sin(pkin(12));
	t63 = -(-t53 * t52 + t56 * t64) * t58 + t57 * t67;
	t62 = (-t57 * t52 - t56 * t68) * t58 + t53 * t67;
	t61 = cos(qJ(3));
	t60 = sin(qJ(3));
	t51 = -t52 * t68 + t57 * t56;
	t49 = t52 * t64 + t53 * t56;
	t1 = [0, 0, -t51 * t60 + t62 * t61, 0, 0, 0; 0, 0, -t49 * t60 - t63 * t61, 0, 0, 0; 0, 0, t61 * t66 + (-t52 * t60 + t61 * t65) * t55, 0, 0, 0; 0, 0, -t51 * t61 - t62 * t60, 0, 0, 0; 0, 0, -t49 * t61 + t63 * t60, 0, 0, 0; 0, 0, -t60 * t66 + (-t52 * t61 - t60 * t65) * t55, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:12:34
	% EndTime: 2019-10-09 21:12:34
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (69->26), mult. (203->58), div. (0->0), fcn. (284->12), ass. (0->34)
	t116 = sin(pkin(11));
	t122 = cos(pkin(6));
	t134 = t116 * t122;
	t117 = sin(pkin(7));
	t118 = sin(pkin(6));
	t133 = t117 * t118;
	t132 = t117 * t122;
	t121 = cos(pkin(7));
	t131 = t118 * t121;
	t119 = cos(pkin(12));
	t130 = t119 * t121;
	t120 = cos(pkin(11));
	t129 = t120 * t122;
	t115 = sin(pkin(12));
	t111 = -t116 * t115 + t119 * t129;
	t128 = t111 * t121 - t120 * t133;
	t113 = -t120 * t115 - t119 * t134;
	t127 = t113 * t121 + t116 * t133;
	t126 = cos(qJ(3));
	t125 = cos(qJ(4));
	t124 = sin(qJ(3));
	t123 = sin(qJ(4));
	t114 = -t115 * t134 + t120 * t119;
	t112 = t115 * t129 + t116 * t119;
	t110 = -t119 * t133 + t122 * t121;
	t109 = -t113 * t117 + t116 * t131;
	t108 = -t111 * t117 - t120 * t131;
	t107 = t124 * t132 + (t115 * t126 + t124 * t130) * t118;
	t106 = t126 * t132 + (-t115 * t124 + t126 * t130) * t118;
	t105 = t114 * t126 + t127 * t124;
	t104 = -t114 * t124 + t127 * t126;
	t103 = t112 * t126 + t128 * t124;
	t102 = -t112 * t124 + t128 * t126;
	t1 = [0, 0, t104 * t125, -t105 * t123 + t109 * t125, 0, 0; 0, 0, t102 * t125, -t103 * t123 + t108 * t125, 0, 0; 0, 0, t106 * t125, -t107 * t123 + t110 * t125, 0, 0; 0, 0, -t104 * t123, -t105 * t125 - t109 * t123, 0, 0; 0, 0, -t102 * t123, -t103 * t125 - t108 * t123, 0, 0; 0, 0, -t106 * t123, -t107 * t125 - t110 * t123, 0, 0; 0, 0, t105, 0, 0, 0; 0, 0, t103, 0, 0, 0; 0, 0, t107, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:12:34
	% EndTime: 2019-10-09 21:12:34
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (69->26), mult. (203->58), div. (0->0), fcn. (284->12), ass. (0->34)
	t136 = sin(pkin(11));
	t142 = cos(pkin(6));
	t154 = t136 * t142;
	t137 = sin(pkin(7));
	t138 = sin(pkin(6));
	t153 = t137 * t138;
	t152 = t137 * t142;
	t141 = cos(pkin(7));
	t151 = t138 * t141;
	t139 = cos(pkin(12));
	t150 = t139 * t141;
	t140 = cos(pkin(11));
	t149 = t140 * t142;
	t135 = sin(pkin(12));
	t131 = -t136 * t135 + t139 * t149;
	t148 = t131 * t141 - t140 * t153;
	t133 = -t140 * t135 - t139 * t154;
	t147 = t133 * t141 + t136 * t153;
	t146 = cos(qJ(3));
	t145 = cos(qJ(4));
	t144 = sin(qJ(3));
	t143 = sin(qJ(4));
	t134 = -t135 * t154 + t140 * t139;
	t132 = t135 * t149 + t136 * t139;
	t130 = -t139 * t153 + t142 * t141;
	t129 = -t133 * t137 + t136 * t151;
	t128 = -t131 * t137 - t140 * t151;
	t127 = t144 * t152 + (t135 * t146 + t144 * t150) * t138;
	t126 = t146 * t152 + (-t135 * t144 + t146 * t150) * t138;
	t125 = t134 * t146 + t147 * t144;
	t124 = -t134 * t144 + t147 * t146;
	t123 = t132 * t146 + t148 * t144;
	t122 = -t132 * t144 + t148 * t146;
	t1 = [0, 0, t125, 0, 0, 0; 0, 0, t123, 0, 0, 0; 0, 0, t127, 0, 0, 0; 0, 0, -t124 * t145, t125 * t143 - t129 * t145, 0, 0; 0, 0, -t122 * t145, t123 * t143 - t128 * t145, 0, 0; 0, 0, -t126 * t145, t127 * t143 - t130 * t145, 0, 0; 0, 0, t124 * t143, t125 * t145 + t129 * t143, 0, 0; 0, 0, t122 * t143, t123 * t145 + t128 * t143, 0, 0; 0, 0, t126 * t143, t127 * t145 + t130 * t143, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:12:34
	% EndTime: 2019-10-09 21:12:34
	% DurationCPUTime: 0.25s
	% Computational Cost: add. (172->43), mult. (516->91), div. (0->0), fcn. (712->14), ass. (0->46)
	t209 = cos(qJ(3));
	t183 = sin(pkin(11));
	t189 = cos(pkin(6));
	t208 = t183 * t189;
	t184 = sin(pkin(7));
	t207 = t184 * t189;
	t185 = sin(pkin(6));
	t206 = t185 * t184;
	t188 = cos(pkin(7));
	t205 = t185 * t188;
	t186 = cos(pkin(12));
	t204 = t186 * t188;
	t187 = cos(pkin(11));
	t203 = t187 * t189;
	t190 = sin(qJ(6));
	t191 = sin(qJ(4));
	t202 = t190 * t191;
	t193 = cos(qJ(6));
	t201 = t191 * t193;
	t200 = t185 * t209;
	t199 = t184 * t200;
	t182 = sin(pkin(12));
	t198 = -t183 * t182 + t186 * t203;
	t197 = t187 * t182 + t186 * t208;
	t196 = t198 * t188;
	t195 = t197 * t188;
	t194 = cos(qJ(4));
	t192 = sin(qJ(3));
	t178 = -t182 * t208 + t187 * t186;
	t177 = t182 * t203 + t183 * t186;
	t176 = -t186 * t206 + t189 * t188;
	t173 = t183 * t205 + t197 * t184;
	t172 = -t198 * t184 - t187 * t205;
	t171 = t192 * t207 + (t209 * t182 + t192 * t204) * t185;
	t170 = t185 * t182 * t192 - t200 * t204 - t209 * t207;
	t169 = t171 * t194 + t176 * t191;
	t168 = t171 * t191 - t176 * t194;
	t167 = t178 * t209 + (t183 * t206 - t195) * t192;
	t166 = t178 * t192 - t183 * t199 + t209 * t195;
	t165 = t177 * t209 + (-t187 * t206 + t196) * t192;
	t164 = t177 * t192 + t187 * t199 - t209 * t196;
	t163 = t167 * t194 + t173 * t191;
	t162 = t167 * t191 - t173 * t194;
	t161 = t165 * t194 + t172 * t191;
	t160 = t165 * t191 - t172 * t194;
	t1 = [0, 0, -t166 * t202 + t167 * t193, t163 * t190, 0, t162 * t193 - t166 * t190; 0, 0, -t164 * t202 + t165 * t193, t161 * t190, 0, t160 * t193 - t164 * t190; 0, 0, -t170 * t202 + t171 * t193, t169 * t190, 0, t168 * t193 - t170 * t190; 0, 0, -t166 * t201 - t167 * t190, t163 * t193, 0, -t162 * t190 - t166 * t193; 0, 0, -t164 * t201 - t165 * t190, t161 * t193, 0, -t160 * t190 - t164 * t193; 0, 0, -t170 * t201 - t171 * t190, t169 * t193, 0, -t168 * t190 - t170 * t193; 0, 0, -t166 * t194, -t162, 0, 0; 0, 0, -t164 * t194, -t160, 0, 0; 0, 0, -t170 * t194, -t168, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end