% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRPRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:00
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRPRR9_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR9_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_jacobiR_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:00:09
	% EndTime: 2019-10-10 01:00:09
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:00:09
	% EndTime: 2019-10-10 01:00:09
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
	% StartTime: 2019-10-10 01:00:10
	% EndTime: 2019-10-10 01:00:10
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (4->4), mult. (14->10), div. (0->0), fcn. (24->6), ass. (0->11)
	t38 = sin(pkin(12));
	t42 = sin(qJ(1));
	t47 = t42 * t38;
	t40 = cos(pkin(12));
	t46 = t42 * t40;
	t43 = cos(qJ(1));
	t45 = t43 * t38;
	t44 = t43 * t40;
	t41 = cos(pkin(6));
	t39 = sin(pkin(6));
	t1 = [-t41 * t45 - t46, 0, 0, 0, 0, 0; -t41 * t47 + t44, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t41 * t44 + t47, 0, 0, 0, 0, 0; -t41 * t46 - t45, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t43 * t39, 0, 0, 0, 0, 0; t42 * t39, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:00:10
	% EndTime: 2019-10-10 01:00:10
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (40->17), mult. (122->36), div. (0->0), fcn. (174->10), ass. (0->27)
	t103 = cos(qJ(1));
	t96 = sin(pkin(6));
	t111 = t103 * t96;
	t101 = sin(qJ(1));
	t99 = cos(pkin(6));
	t110 = t103 * t99;
	t94 = sin(pkin(12));
	t97 = cos(pkin(12));
	t88 = t101 * t94 - t97 * t110;
	t95 = sin(pkin(7));
	t98 = cos(pkin(7));
	t117 = t95 * t111 + t88 * t98;
	t100 = sin(qJ(3));
	t102 = cos(qJ(3));
	t89 = t101 * t97 + t94 * t110;
	t116 = t89 * t100 + t117 * t102;
	t114 = t94 * t96;
	t113 = t101 * t96;
	t112 = t101 * t99;
	t107 = t96 * t97 * t98 + t95 * t99;
	t90 = -t103 * t94 - t97 * t112;
	t106 = t95 * t113 + t90 * t98;
	t104 = t117 * t100 - t89 * t102;
	t91 = t103 * t97 - t94 * t112;
	t87 = t106 * t100 + t91 * t102;
	t86 = -t91 * t100 + t106 * t102;
	t1 = [t104, 0, t86, 0, 0, 0; t87, 0, -t116, 0, 0, 0; 0, 0, -t100 * t114 + t107 * t102, 0, 0, 0; t116, 0, -t87, 0, 0, 0; t86, 0, t104, 0, 0, 0; 0, 0, -t107 * t100 - t102 * t114, 0, 0, 0; t98 * t111 - t88 * t95, 0, 0, 0, 0, 0; t98 * t113 - t90 * t95, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:00:10
	% EndTime: 2019-10-10 01:00:10
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (70->22), mult. (202->42), div. (0->0), fcn. (284->12), ass. (0->33)
	t119 = sin(pkin(6));
	t125 = sin(qJ(1));
	t137 = t119 * t125;
	t127 = cos(qJ(1));
	t136 = t119 * t127;
	t117 = sin(pkin(12));
	t135 = t125 * t117;
	t121 = cos(pkin(12));
	t134 = t125 * t121;
	t133 = t127 * t117;
	t132 = t127 * t121;
	t116 = sin(pkin(13));
	t120 = cos(pkin(13));
	t124 = sin(qJ(3));
	t126 = cos(qJ(3));
	t131 = t126 * t116 + t124 * t120;
	t112 = t124 * t116 - t126 * t120;
	t118 = sin(pkin(7));
	t104 = t112 * t118;
	t122 = cos(pkin(7));
	t106 = t112 * t122;
	t123 = cos(pkin(6));
	t108 = -t123 * t132 + t135;
	t109 = t123 * t133 + t134;
	t130 = -t104 * t136 - t108 * t106 + t109 * t131;
	t105 = t131 * t118;
	t107 = t131 * t122;
	t129 = t105 * t136 + t108 * t107 + t109 * t112;
	t110 = -t123 * t134 - t133;
	t111 = -t123 * t135 + t132;
	t128 = t105 * t137 + t110 * t107 - t111 * t112;
	t103 = -t104 * t137 - t110 * t106 - t111 * t131;
	t1 = [t129, 0, t103, 0, 0, 0; t128, 0, -t130, 0, 0, 0; 0, 0, -t123 * t104 + (-t106 * t121 - t117 * t131) * t119, 0, 0, 0; t130, 0, -t128, 0, 0, 0; t103, 0, t129, 0, 0, 0; 0, 0, -t123 * t105 + (-t107 * t121 + t112 * t117) * t119, 0, 0, 0; -t108 * t118 + t122 * t136, 0, 0, 0, 0, 0; -t110 * t118 + t122 * t137, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:00:10
	% EndTime: 2019-10-10 01:00:10
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (168->34), mult. (475->63), div. (0->0), fcn. (661->14), ass. (0->44)
	t180 = sin(pkin(7));
	t178 = sin(pkin(13));
	t182 = cos(pkin(13));
	t187 = sin(qJ(3));
	t190 = cos(qJ(3));
	t192 = t190 * t178 + t187 * t182;
	t166 = t192 * t180;
	t184 = cos(pkin(7));
	t168 = t192 * t184;
	t185 = cos(pkin(6));
	t183 = cos(pkin(12));
	t191 = cos(qJ(1));
	t193 = t191 * t183;
	t179 = sin(pkin(12));
	t188 = sin(qJ(1));
	t196 = t188 * t179;
	t170 = -t185 * t193 + t196;
	t194 = t191 * t179;
	t195 = t188 * t183;
	t171 = t185 * t194 + t195;
	t174 = t187 * t178 - t190 * t182;
	t181 = sin(pkin(6));
	t197 = t181 * t191;
	t157 = t166 * t197 + t170 * t168 + t171 * t174;
	t162 = -t170 * t180 + t184 * t197;
	t186 = sin(qJ(5));
	t189 = cos(qJ(5));
	t200 = t157 * t189 + t162 * t186;
	t199 = t157 * t186 - t162 * t189;
	t198 = t181 * t188;
	t165 = t174 * t180;
	t167 = t174 * t184;
	t155 = t165 * t197 + t170 * t167 - t171 * t192;
	t172 = -t185 * t195 - t194;
	t173 = -t185 * t196 + t193;
	t159 = t166 * t198 + t172 * t168 - t173 * t174;
	t161 = t185 * t166 + (t168 * t183 - t174 * t179) * t181;
	t169 = -t181 * t183 * t180 + t185 * t184;
	t164 = -t172 * t180 + t184 * t198;
	t160 = -t185 * t165 + (-t167 * t183 - t179 * t192) * t181;
	t158 = -t165 * t198 - t172 * t167 - t173 * t192;
	t154 = t159 * t189 + t164 * t186;
	t153 = -t159 * t186 + t164 * t189;
	t1 = [t200, 0, t158 * t189, 0, t153, 0; t154, 0, t155 * t189, 0, t199, 0; 0, 0, t160 * t189, 0, -t161 * t186 + t169 * t189, 0; -t199, 0, -t158 * t186, 0, -t154, 0; t153, 0, -t155 * t186, 0, t200, 0; 0, 0, -t160 * t186, 0, -t161 * t189 - t169 * t186, 0; t155, 0, t159, 0, 0, 0; -t158, 0, -t157, 0, 0, 0; 0, 0, t161, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:00:11
	% EndTime: 2019-10-10 01:00:11
	% DurationCPUTime: 0.41s
	% Computational Cost: add. (375->45), mult. (1060->92), div. (0->0), fcn. (1462->16), ass. (0->54)
	t239 = sin(pkin(7));
	t237 = sin(pkin(13));
	t241 = cos(pkin(13));
	t247 = sin(qJ(3));
	t251 = cos(qJ(3));
	t257 = t251 * t237 + t247 * t241;
	t225 = t257 * t239;
	t243 = cos(pkin(7));
	t227 = t257 * t243;
	t244 = cos(pkin(6));
	t242 = cos(pkin(12));
	t252 = cos(qJ(1));
	t258 = t252 * t242;
	t238 = sin(pkin(12));
	t248 = sin(qJ(1));
	t262 = t248 * t238;
	t229 = -t244 * t258 + t262;
	t259 = t252 * t238;
	t261 = t248 * t242;
	t230 = t244 * t259 + t261;
	t233 = t247 * t237 - t251 * t241;
	t240 = sin(pkin(6));
	t264 = t240 * t252;
	t213 = t225 * t264 + t229 * t227 + t230 * t233;
	t221 = -t229 * t239 + t243 * t264;
	t246 = sin(qJ(5));
	t250 = cos(qJ(5));
	t204 = t213 * t250 + t221 * t246;
	t245 = sin(qJ(6));
	t249 = cos(qJ(6));
	t224 = t233 * t239;
	t226 = t233 * t243;
	t255 = t224 * t264 + t229 * t226 - t230 * t257;
	t269 = t204 * t245 - t249 * t255;
	t268 = t204 * t249 + t245 * t255;
	t202 = t213 * t246 - t221 * t250;
	t265 = t240 * t248;
	t263 = t245 * t250;
	t260 = t249 * t250;
	t231 = -t244 * t261 - t259;
	t256 = -t231 * t239 + t243 * t265;
	t232 = -t244 * t262 + t258;
	t254 = t225 * t265 + t231 * t227 - t232 * t233;
	t253 = t244 * t225 + (t227 * t242 - t233 * t238) * t240;
	t228 = -t240 * t242 * t239 + t244 * t243;
	t218 = -t244 * t224 + (-t226 * t242 - t238 * t257) * t240;
	t215 = -t224 * t265 - t231 * t226 - t232 * t257;
	t208 = t228 * t246 + t250 * t253;
	t207 = t228 * t250 - t246 * t253;
	t206 = t246 * t256 + t250 * t254;
	t205 = t246 * t254 - t250 * t256;
	t201 = t206 * t249 - t215 * t245;
	t200 = -t206 * t245 - t215 * t249;
	t1 = [t268, 0, t215 * t260 + t245 * t254, 0, -t205 * t249, t200; t201, 0, -t213 * t245 + t255 * t260, 0, t202 * t249, t269; 0, 0, t218 * t260 + t245 * t253, 0, t207 * t249, -t208 * t245 - t218 * t249; -t269, 0, -t215 * t263 + t249 * t254, 0, t205 * t245, -t201; t200, 0, -t213 * t249 - t255 * t263, 0, -t202 * t245, t268; 0, 0, -t218 * t263 + t249 * t253, 0, -t207 * t245, -t208 * t249 + t218 * t245; t202, 0, t215 * t246, 0, t206, 0; t205, 0, t255 * t246, 0, -t204, 0; 0, 0, t218 * t246, 0, t208, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end