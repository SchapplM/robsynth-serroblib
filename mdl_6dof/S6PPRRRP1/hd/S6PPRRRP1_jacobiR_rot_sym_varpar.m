% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PPRRRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:14
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PPRRRP1_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPRRRP1_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_jacobiR_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:14:28
	% EndTime: 2019-10-09 21:14:28
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:14:28
	% EndTime: 2019-10-09 21:14:28
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:14:28
	% EndTime: 2019-10-09 21:14:28
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:14:28
	% EndTime: 2019-10-09 21:14:28
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
	% StartTime: 2019-10-09 21:14:28
	% EndTime: 2019-10-09 21:14:28
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
	% StartTime: 2019-10-09 21:14:29
	% EndTime: 2019-10-09 21:14:29
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (175->43), mult. (516->91), div. (0->0), fcn. (712->14), ass. (0->46)
	t206 = cos(qJ(3));
	t180 = sin(pkin(11));
	t186 = cos(pkin(6));
	t205 = t180 * t186;
	t181 = sin(pkin(7));
	t204 = t181 * t186;
	t182 = sin(pkin(6));
	t203 = t182 * t181;
	t185 = cos(pkin(7));
	t202 = t182 * t185;
	t183 = cos(pkin(12));
	t201 = t183 * t185;
	t184 = cos(pkin(11));
	t200 = t184 * t186;
	t187 = sin(qJ(5));
	t191 = cos(qJ(4));
	t199 = t187 * t191;
	t190 = cos(qJ(5));
	t198 = t190 * t191;
	t197 = t182 * t206;
	t196 = t181 * t197;
	t179 = sin(pkin(12));
	t195 = -t180 * t179 + t183 * t200;
	t194 = t184 * t179 + t183 * t205;
	t193 = t195 * t185;
	t192 = t194 * t185;
	t189 = sin(qJ(3));
	t188 = sin(qJ(4));
	t175 = -t179 * t205 + t184 * t183;
	t174 = t179 * t200 + t180 * t183;
	t173 = -t183 * t203 + t186 * t185;
	t170 = t180 * t202 + t194 * t181;
	t169 = -t195 * t181 - t184 * t202;
	t168 = t189 * t204 + (t206 * t179 + t189 * t201) * t182;
	t167 = t182 * t179 * t189 - t197 * t201 - t206 * t204;
	t166 = t168 * t191 + t173 * t188;
	t165 = -t168 * t188 + t173 * t191;
	t164 = t175 * t206 + (t180 * t203 - t192) * t189;
	t163 = t175 * t189 - t180 * t196 + t206 * t192;
	t162 = t174 * t206 + (-t184 * t203 + t193) * t189;
	t161 = t174 * t189 + t184 * t196 - t206 * t193;
	t160 = t164 * t191 + t170 * t188;
	t159 = -t164 * t188 + t170 * t191;
	t158 = t162 * t191 + t169 * t188;
	t157 = -t162 * t188 + t169 * t191;
	t1 = [0, 0, -t163 * t198 + t164 * t187, t159 * t190, -t160 * t187 + t163 * t190, 0; 0, 0, -t161 * t198 + t162 * t187, t157 * t190, -t158 * t187 + t161 * t190, 0; 0, 0, -t167 * t198 + t168 * t187, t165 * t190, -t166 * t187 + t167 * t190, 0; 0, 0, t163 * t199 + t164 * t190, -t159 * t187, -t160 * t190 - t163 * t187, 0; 0, 0, t161 * t199 + t162 * t190, -t157 * t187, -t158 * t190 - t161 * t187, 0; 0, 0, t167 * t199 + t168 * t190, -t165 * t187, -t166 * t190 - t167 * t187, 0; 0, 0, -t163 * t188, t160, 0, 0; 0, 0, -t161 * t188, t158, 0, 0; 0, 0, -t167 * t188, t166, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:14:29
	% EndTime: 2019-10-09 21:14:29
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (175->43), mult. (516->91), div. (0->0), fcn. (712->14), ass. (0->46)
	t217 = cos(qJ(3));
	t191 = sin(pkin(11));
	t197 = cos(pkin(6));
	t216 = t191 * t197;
	t192 = sin(pkin(7));
	t215 = t192 * t197;
	t193 = sin(pkin(6));
	t214 = t193 * t192;
	t196 = cos(pkin(7));
	t213 = t193 * t196;
	t194 = cos(pkin(12));
	t212 = t194 * t196;
	t195 = cos(pkin(11));
	t211 = t195 * t197;
	t198 = sin(qJ(5));
	t202 = cos(qJ(4));
	t210 = t198 * t202;
	t201 = cos(qJ(5));
	t209 = t201 * t202;
	t208 = t193 * t217;
	t207 = t192 * t208;
	t190 = sin(pkin(12));
	t206 = -t191 * t190 + t194 * t211;
	t205 = t195 * t190 + t194 * t216;
	t204 = t206 * t196;
	t203 = t205 * t196;
	t200 = sin(qJ(3));
	t199 = sin(qJ(4));
	t186 = -t190 * t216 + t195 * t194;
	t185 = t190 * t211 + t191 * t194;
	t184 = -t194 * t214 + t197 * t196;
	t181 = t191 * t213 + t205 * t192;
	t180 = -t206 * t192 - t195 * t213;
	t179 = t200 * t215 + (t217 * t190 + t200 * t212) * t193;
	t178 = t193 * t190 * t200 - t208 * t212 - t217 * t215;
	t177 = t179 * t202 + t184 * t199;
	t176 = -t179 * t199 + t184 * t202;
	t175 = t186 * t217 + (t191 * t214 - t203) * t200;
	t174 = t186 * t200 - t191 * t207 + t217 * t203;
	t173 = t185 * t217 + (-t195 * t214 + t204) * t200;
	t172 = t185 * t200 + t195 * t207 - t217 * t204;
	t171 = t175 * t202 + t181 * t199;
	t170 = -t175 * t199 + t181 * t202;
	t169 = t173 * t202 + t180 * t199;
	t168 = -t173 * t199 + t180 * t202;
	t1 = [0, 0, -t174 * t209 + t175 * t198, t170 * t201, -t171 * t198 + t174 * t201, 0; 0, 0, -t172 * t209 + t173 * t198, t168 * t201, -t169 * t198 + t172 * t201, 0; 0, 0, -t178 * t209 + t179 * t198, t176 * t201, -t177 * t198 + t178 * t201, 0; 0, 0, t174 * t210 + t175 * t201, -t170 * t198, -t171 * t201 - t174 * t198, 0; 0, 0, t172 * t210 + t173 * t201, -t168 * t198, -t169 * t201 - t172 * t198, 0; 0, 0, t178 * t210 + t179 * t201, -t176 * t198, -t177 * t201 - t178 * t198, 0; 0, 0, -t174 * t199, t171, 0, 0; 0, 0, -t172 * t199, t169, 0, 0; 0, 0, -t178 * t199, t177, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end