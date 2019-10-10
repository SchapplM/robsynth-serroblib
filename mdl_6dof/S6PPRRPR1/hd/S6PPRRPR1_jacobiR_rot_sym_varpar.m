% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PPRRPR1
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:10
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PPRRPR1_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPRRPR1_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_jacobiR_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:10:40
	% EndTime: 2019-10-09 21:10:40
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:10:40
	% EndTime: 2019-10-09 21:10:40
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:10:40
	% EndTime: 2019-10-09 21:10:40
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:10:40
	% EndTime: 2019-10-09 21:10:40
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
	% StartTime: 2019-10-09 21:10:40
	% EndTime: 2019-10-09 21:10:40
	% DurationCPUTime: 0.10s
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
	% StartTime: 2019-10-09 21:10:40
	% EndTime: 2019-10-09 21:10:40
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (113->32), mult. (338->75), div. (0->0), fcn. (466->14), ass. (0->41)
	t152 = sin(pkin(13));
	t164 = cos(qJ(4));
	t175 = t152 * t164;
	t154 = sin(pkin(11));
	t161 = cos(pkin(6));
	t174 = t154 * t161;
	t155 = sin(pkin(7));
	t156 = sin(pkin(6));
	t173 = t155 * t156;
	t172 = t155 * t161;
	t160 = cos(pkin(7));
	t171 = t156 * t160;
	t157 = cos(pkin(13));
	t170 = t157 * t164;
	t158 = cos(pkin(12));
	t169 = t158 * t160;
	t159 = cos(pkin(11));
	t168 = t159 * t161;
	t153 = sin(pkin(12));
	t148 = -t153 * t154 + t158 * t168;
	t167 = t148 * t160 - t159 * t173;
	t150 = -t153 * t159 - t158 * t174;
	t166 = t150 * t160 + t154 * t173;
	t165 = cos(qJ(3));
	t163 = sin(qJ(3));
	t162 = sin(qJ(4));
	t151 = -t153 * t174 + t158 * t159;
	t149 = t153 * t168 + t154 * t158;
	t147 = -t158 * t173 + t160 * t161;
	t146 = -t150 * t155 + t154 * t171;
	t145 = -t148 * t155 - t159 * t171;
	t144 = t163 * t172 + (t153 * t165 + t163 * t169) * t156;
	t143 = t165 * t172 + (-t153 * t163 + t165 * t169) * t156;
	t142 = -t144 * t162 + t147 * t164;
	t141 = t151 * t165 + t163 * t166;
	t140 = -t151 * t163 + t165 * t166;
	t139 = t149 * t165 + t163 * t167;
	t138 = -t149 * t163 + t165 * t167;
	t137 = -t141 * t162 + t146 * t164;
	t136 = -t139 * t162 + t145 * t164;
	t1 = [0, 0, t140 * t170 + t141 * t152, t137 * t157, 0, 0; 0, 0, t138 * t170 + t139 * t152, t136 * t157, 0, 0; 0, 0, t143 * t170 + t144 * t152, t142 * t157, 0, 0; 0, 0, -t140 * t175 + t141 * t157, -t137 * t152, 0, 0; 0, 0, -t138 * t175 + t139 * t157, -t136 * t152, 0, 0; 0, 0, -t143 * t175 + t144 * t157, -t142 * t152, 0, 0; 0, 0, t140 * t162, t141 * t164 + t146 * t162, 0, 0; 0, 0, t138 * t162, t139 * t164 + t145 * t162, 0, 0; 0, 0, t143 * t162, t144 * t164 + t147 * t162, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:10:40
	% EndTime: 2019-10-09 21:10:41
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (205->44), mult. (516->91), div. (0->0), fcn. (712->14), ass. (0->47)
	t217 = cos(qJ(3));
	t191 = pkin(13) + qJ(6);
	t189 = sin(t191);
	t202 = cos(qJ(4));
	t216 = t189 * t202;
	t190 = cos(t191);
	t215 = t190 * t202;
	t193 = sin(pkin(11));
	t199 = cos(pkin(6));
	t214 = t193 * t199;
	t194 = sin(pkin(7));
	t213 = t194 * t199;
	t195 = sin(pkin(6));
	t212 = t195 * t194;
	t198 = cos(pkin(7));
	t211 = t195 * t198;
	t196 = cos(pkin(12));
	t210 = t196 * t198;
	t197 = cos(pkin(11));
	t209 = t197 * t199;
	t208 = t195 * t217;
	t207 = t194 * t208;
	t192 = sin(pkin(12));
	t206 = -t192 * t193 + t196 * t209;
	t205 = t192 * t197 + t196 * t214;
	t204 = t206 * t198;
	t203 = t205 * t198;
	t201 = sin(qJ(3));
	t200 = sin(qJ(4));
	t185 = -t192 * t214 + t196 * t197;
	t184 = t192 * t209 + t193 * t196;
	t183 = -t196 * t212 + t198 * t199;
	t180 = t193 * t211 + t194 * t205;
	t179 = -t194 * t206 - t197 * t211;
	t178 = t201 * t213 + (t192 * t217 + t201 * t210) * t195;
	t177 = t195 * t192 * t201 - t208 * t210 - t213 * t217;
	t176 = t178 * t202 + t183 * t200;
	t175 = -t178 * t200 + t183 * t202;
	t174 = t185 * t217 + (t193 * t212 - t203) * t201;
	t173 = t185 * t201 - t193 * t207 + t203 * t217;
	t172 = t184 * t217 + (-t197 * t212 + t204) * t201;
	t171 = t184 * t201 + t197 * t207 - t204 * t217;
	t170 = t174 * t202 + t180 * t200;
	t169 = -t174 * t200 + t180 * t202;
	t168 = t172 * t202 + t179 * t200;
	t167 = -t172 * t200 + t179 * t202;
	t1 = [0, 0, -t173 * t215 + t174 * t189, t169 * t190, 0, -t170 * t189 + t173 * t190; 0, 0, -t171 * t215 + t172 * t189, t167 * t190, 0, -t168 * t189 + t171 * t190; 0, 0, -t177 * t215 + t178 * t189, t175 * t190, 0, -t176 * t189 + t177 * t190; 0, 0, t173 * t216 + t174 * t190, -t169 * t189, 0, -t170 * t190 - t173 * t189; 0, 0, t171 * t216 + t172 * t190, -t167 * t189, 0, -t168 * t190 - t171 * t189; 0, 0, t177 * t216 + t178 * t190, -t175 * t189, 0, -t176 * t190 - t177 * t189; 0, 0, -t173 * t200, t170, 0, 0; 0, 0, -t171 * t200, t168, 0, 0; 0, 0, -t177 * t200, t176, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end