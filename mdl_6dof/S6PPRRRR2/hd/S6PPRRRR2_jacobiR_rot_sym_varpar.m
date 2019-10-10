% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PPRRRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:20
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PPRRRR2_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR2_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPRRRR2_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR2_jacobiR_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:20:19
	% EndTime: 2019-10-09 21:20:19
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:20:19
	% EndTime: 2019-10-09 21:20:19
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:20:19
	% EndTime: 2019-10-09 21:20:19
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:20:19
	% EndTime: 2019-10-09 21:20:19
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (20->14), mult. (62->33), div. (0->0), fcn. (88->10), ass. (0->20)
	t53 = sin(pkin(12));
	t59 = cos(pkin(6));
	t68 = t53 * t59;
	t54 = sin(pkin(7));
	t55 = sin(pkin(6));
	t67 = t54 * t55;
	t66 = t54 * t59;
	t56 = cos(pkin(13));
	t58 = cos(pkin(7));
	t65 = t56 * t58;
	t57 = cos(pkin(12));
	t64 = t57 * t59;
	t52 = sin(pkin(13));
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
	% StartTime: 2019-10-09 21:20:19
	% EndTime: 2019-10-09 21:20:19
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (69->26), mult. (203->58), div. (0->0), fcn. (284->12), ass. (0->34)
	t116 = sin(pkin(12));
	t122 = cos(pkin(6));
	t134 = t116 * t122;
	t117 = sin(pkin(7));
	t118 = sin(pkin(6));
	t133 = t117 * t118;
	t132 = t117 * t122;
	t121 = cos(pkin(7));
	t131 = t118 * t121;
	t119 = cos(pkin(13));
	t130 = t119 * t121;
	t120 = cos(pkin(12));
	t129 = t120 * t122;
	t115 = sin(pkin(13));
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
	% StartTime: 2019-10-09 21:20:19
	% EndTime: 2019-10-09 21:20:20
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (175->43), mult. (516->91), div. (0->0), fcn. (712->14), ass. (0->46)
	t206 = cos(qJ(3));
	t180 = sin(pkin(12));
	t186 = cos(pkin(6));
	t205 = t180 * t186;
	t181 = sin(pkin(7));
	t204 = t181 * t186;
	t182 = sin(pkin(6));
	t203 = t182 * t181;
	t185 = cos(pkin(7));
	t202 = t182 * t185;
	t183 = cos(pkin(13));
	t201 = t183 * t185;
	t184 = cos(pkin(12));
	t200 = t184 * t186;
	t187 = sin(qJ(5));
	t191 = cos(qJ(4));
	t199 = t187 * t191;
	t190 = cos(qJ(5));
	t198 = t190 * t191;
	t197 = t182 * t206;
	t196 = t181 * t197;
	t179 = sin(pkin(13));
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
	% StartTime: 2019-10-09 21:20:19
	% EndTime: 2019-10-09 21:20:20
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (279->44), mult. (694->91), div. (0->0), fcn. (958->14), ass. (0->53)
	t232 = cos(qJ(3));
	t206 = qJ(5) + qJ(6);
	t204 = sin(t206);
	t217 = cos(qJ(4));
	t231 = t204 * t217;
	t205 = cos(t206);
	t230 = t205 * t217;
	t208 = sin(pkin(12));
	t214 = cos(pkin(6));
	t229 = t208 * t214;
	t209 = sin(pkin(7));
	t228 = t209 * t214;
	t210 = sin(pkin(6));
	t227 = t210 * t209;
	t213 = cos(pkin(7));
	t226 = t210 * t213;
	t211 = cos(pkin(13));
	t225 = t211 * t213;
	t212 = cos(pkin(12));
	t224 = t212 * t214;
	t223 = t210 * t232;
	t222 = t209 * t223;
	t207 = sin(pkin(13));
	t221 = -t208 * t207 + t211 * t224;
	t220 = t212 * t207 + t211 * t229;
	t219 = t221 * t213;
	t218 = t220 * t213;
	t216 = sin(qJ(3));
	t215 = sin(qJ(4));
	t200 = -t207 * t229 + t212 * t211;
	t199 = t207 * t224 + t208 * t211;
	t198 = -t211 * t227 + t214 * t213;
	t195 = t208 * t226 + t220 * t209;
	t194 = -t221 * t209 - t212 * t226;
	t193 = t216 * t228 + (t232 * t207 + t216 * t225) * t210;
	t192 = t210 * t207 * t216 - t223 * t225 - t232 * t228;
	t191 = t193 * t217 + t198 * t215;
	t190 = -t193 * t215 + t198 * t217;
	t189 = t200 * t232 + (t208 * t227 - t218) * t216;
	t188 = t200 * t216 - t208 * t222 + t232 * t218;
	t187 = t199 * t232 + (-t212 * t227 + t219) * t216;
	t186 = t199 * t216 + t212 * t222 - t232 * t219;
	t185 = t189 * t217 + t195 * t215;
	t184 = -t189 * t215 + t195 * t217;
	t183 = t187 * t217 + t194 * t215;
	t182 = -t187 * t215 + t194 * t217;
	t181 = -t191 * t205 - t192 * t204;
	t180 = -t191 * t204 + t192 * t205;
	t179 = -t185 * t205 - t188 * t204;
	t178 = -t185 * t204 + t188 * t205;
	t177 = -t183 * t205 - t186 * t204;
	t176 = -t183 * t204 + t186 * t205;
	t1 = [0, 0, -t188 * t230 + t189 * t204, t184 * t205, t178, t178; 0, 0, -t186 * t230 + t187 * t204, t182 * t205, t176, t176; 0, 0, -t192 * t230 + t193 * t204, t190 * t205, t180, t180; 0, 0, t188 * t231 + t189 * t205, -t184 * t204, t179, t179; 0, 0, t186 * t231 + t187 * t205, -t182 * t204, t177, t177; 0, 0, t192 * t231 + t193 * t205, -t190 * t204, t181, t181; 0, 0, -t188 * t215, t185, 0, 0; 0, 0, -t186 * t215, t183, 0, 0; 0, 0, -t192 * t215, t191, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end