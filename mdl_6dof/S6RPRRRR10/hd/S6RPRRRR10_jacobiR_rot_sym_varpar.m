% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:11
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRRR10_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR10_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRR10_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR10_jacobiR_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:11:29
	% EndTime: 2019-10-10 09:11:29
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:11:29
	% EndTime: 2019-10-10 09:11:29
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
	% StartTime: 2019-10-10 09:11:29
	% EndTime: 2019-10-10 09:11:29
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (4->4), mult. (14->10), div. (0->0), fcn. (24->6), ass. (0->11)
	t38 = sin(pkin(13));
	t42 = sin(qJ(1));
	t47 = t42 * t38;
	t40 = cos(pkin(13));
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
	% StartTime: 2019-10-10 09:11:30
	% EndTime: 2019-10-10 09:11:30
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (40->17), mult. (122->36), div. (0->0), fcn. (174->10), ass. (0->27)
	t103 = cos(qJ(1));
	t96 = sin(pkin(6));
	t111 = t103 * t96;
	t101 = sin(qJ(1));
	t99 = cos(pkin(6));
	t110 = t103 * t99;
	t94 = sin(pkin(13));
	t97 = cos(pkin(13));
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
	% StartTime: 2019-10-10 09:11:30
	% EndTime: 2019-10-10 09:11:30
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (105->30), mult. (307->59), div. (0->0), fcn. (430->12), ass. (0->40)
	t160 = cos(pkin(6));
	t155 = sin(pkin(13));
	t166 = cos(qJ(1));
	t170 = t166 * t155;
	t158 = cos(pkin(13));
	t163 = sin(qJ(1));
	t171 = t163 * t158;
	t150 = t160 * t170 + t171;
	t162 = sin(qJ(3));
	t165 = cos(qJ(3));
	t169 = t166 * t158;
	t172 = t163 * t155;
	t149 = -t160 * t169 + t172;
	t156 = sin(pkin(7));
	t159 = cos(pkin(7));
	t157 = sin(pkin(6));
	t174 = t157 * t166;
	t167 = t149 * t159 + t156 * t174;
	t139 = -t150 * t165 + t167 * t162;
	t144 = -t149 * t156 + t159 * t174;
	t161 = sin(qJ(4));
	t164 = cos(qJ(4));
	t181 = t139 * t161 - t144 * t164;
	t180 = t139 * t164 + t144 * t161;
	t176 = t156 * t160;
	t175 = t157 * t163;
	t173 = t159 * t165;
	t168 = t156 * t175;
	t137 = -t150 * t162 - t167 * t165;
	t152 = -t160 * t172 + t169;
	t151 = -t160 * t171 - t170;
	t148 = -t157 * t158 * t156 + t160 * t159;
	t146 = -t151 * t156 + t159 * t175;
	t143 = t162 * t176 + (t158 * t159 * t162 + t155 * t165) * t157;
	t142 = t165 * t176 + (-t155 * t162 + t158 * t173) * t157;
	t141 = t152 * t165 + (t151 * t159 + t168) * t162;
	t140 = -t151 * t173 + t152 * t162 - t165 * t168;
	t136 = t141 * t164 + t146 * t161;
	t135 = -t141 * t161 + t146 * t164;
	t1 = [t180, 0, -t140 * t164, t135, 0, 0; t136, 0, t137 * t164, t181, 0, 0; 0, 0, t142 * t164, -t143 * t161 + t148 * t164, 0, 0; -t181, 0, t140 * t161, -t136, 0, 0; t135, 0, -t137 * t161, t180, 0, 0; 0, 0, -t142 * t161, -t143 * t164 - t148 * t161, 0, 0; t137, 0, t141, 0, 0, 0; t140, 0, -t139, 0, 0, 0; 0, 0, t143, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:11:30
	% EndTime: 2019-10-10 09:11:30
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (179->32), mult. (411->59), div. (0->0), fcn. (576->12), ass. (0->43)
	t186 = cos(pkin(6));
	t181 = sin(pkin(13));
	t190 = cos(qJ(1));
	t194 = t190 * t181;
	t184 = cos(pkin(13));
	t188 = sin(qJ(1));
	t195 = t188 * t184;
	t173 = t186 * t194 + t195;
	t187 = sin(qJ(3));
	t189 = cos(qJ(3));
	t193 = t190 * t184;
	t196 = t188 * t181;
	t172 = -t186 * t193 + t196;
	t182 = sin(pkin(7));
	t185 = cos(pkin(7));
	t183 = sin(pkin(6));
	t198 = t183 * t190;
	t191 = t172 * t185 + t182 * t198;
	t162 = -t173 * t189 + t187 * t191;
	t167 = -t172 * t182 + t185 * t198;
	t180 = qJ(4) + qJ(5);
	t178 = sin(t180);
	t179 = cos(t180);
	t154 = t162 * t178 - t167 * t179;
	t155 = t162 * t179 + t167 * t178;
	t200 = t182 * t186;
	t199 = t183 * t188;
	t197 = t185 * t189;
	t192 = t182 * t199;
	t160 = -t173 * t187 - t189 * t191;
	t175 = -t186 * t196 + t193;
	t174 = -t186 * t195 - t194;
	t171 = -t183 * t184 * t182 + t186 * t185;
	t169 = -t174 * t182 + t185 * t199;
	t166 = t187 * t200 + (t184 * t185 * t187 + t181 * t189) * t183;
	t165 = t189 * t200 + (-t181 * t187 + t184 * t197) * t183;
	t164 = t175 * t189 + (t174 * t185 + t192) * t187;
	t163 = -t174 * t197 + t175 * t187 - t189 * t192;
	t159 = -t166 * t179 - t171 * t178;
	t158 = -t166 * t178 + t171 * t179;
	t157 = t164 * t179 + t169 * t178;
	t156 = -t164 * t178 + t169 * t179;
	t1 = [t155, 0, -t163 * t179, t156, t156, 0; t157, 0, t160 * t179, t154, t154, 0; 0, 0, t165 * t179, t158, t158, 0; -t154, 0, t163 * t178, -t157, -t157, 0; t156, 0, -t160 * t178, t155, t155, 0; 0, 0, -t165 * t178, t159, t159, 0; t160, 0, t164, 0, 0, 0; t163, 0, -t162, 0, 0, 0; 0, 0, t166, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:11:31
	% EndTime: 2019-10-10 09:11:31
	% DurationCPUTime: 0.38s
	% Computational Cost: add. (363->52), mult. (854->91), div. (0->0), fcn. (1181->14), ass. (0->61)
	t251 = cos(pkin(6));
	t246 = sin(pkin(13));
	t257 = cos(qJ(1));
	t265 = t257 * t246;
	t249 = cos(pkin(13));
	t254 = sin(qJ(1));
	t266 = t254 * t249;
	t238 = t251 * t265 + t266;
	t253 = sin(qJ(3));
	t256 = cos(qJ(3));
	t264 = t257 * t249;
	t267 = t254 * t246;
	t237 = -t251 * t264 + t267;
	t250 = cos(pkin(7));
	t247 = sin(pkin(7));
	t248 = sin(pkin(6));
	t269 = t248 * t257;
	t263 = t247 * t269;
	t261 = t237 * t250 + t263;
	t226 = -t238 * t256 + t261 * t253;
	t232 = -t237 * t247 + t250 * t269;
	t245 = qJ(4) + qJ(5);
	t243 = sin(t245);
	t244 = cos(t245);
	t217 = t226 * t244 + t232 * t243;
	t252 = sin(qJ(6));
	t284 = t217 * t252;
	t255 = cos(qJ(6));
	t283 = t217 * t255;
	t215 = t226 * t243 - t232 * t244;
	t260 = t251 * t266 + t265;
	t270 = t248 * t254;
	t280 = -t247 * t270 + t260 * t250;
	t279 = t215 * t252;
	t239 = -t251 * t267 + t264;
	t228 = t239 * t256 - t280 * t253;
	t258 = t260 * t247 + t250 * t270;
	t218 = t228 * t243 - t258 * t244;
	t278 = t218 * t252;
	t272 = t247 * t251;
	t231 = t253 * t272 + (t249 * t250 * t253 + t246 * t256) * t248;
	t271 = t248 * t249;
	t236 = -t247 * t271 + t251 * t250;
	t221 = -t231 * t243 + t236 * t244;
	t277 = t221 * t252;
	t276 = t238 * t253;
	t274 = t244 * t252;
	t273 = t244 * t255;
	t268 = t250 * t256;
	t230 = t248 * t246 * t253 - t256 * t272 - t268 * t271;
	t227 = t239 * t253 + t280 * t256;
	t225 = -t261 * t256 - t276;
	t223 = t237 * t268 + t256 * t263 + t276;
	t222 = t231 * t244 + t236 * t243;
	t220 = t221 * t255;
	t219 = t228 * t244 + t258 * t243;
	t214 = t218 * t255;
	t213 = t215 * t255;
	t212 = t219 * t255 + t227 * t252;
	t211 = -t219 * t252 + t227 * t255;
	t1 = [t225 * t252 + t283, 0, -t227 * t273 + t228 * t252, -t214, -t214, t211; t212, 0, -t223 * t273 - t226 * t252, t213, t213, t223 * t255 + t284; 0, 0, -t230 * t273 + t231 * t252, t220, t220, -t222 * t252 + t230 * t255; t225 * t255 - t284, 0, t227 * t274 + t228 * t255, t278, t278, -t212; t211, 0, t223 * t274 - t226 * t255, -t279, -t279, -t223 * t252 + t283; 0, 0, t230 * t274 + t231 * t255, -t277, -t277, -t222 * t255 - t230 * t252; t215, 0, -t227 * t243, t219, t219, 0; t218, 0, -t223 * t243, -t217, -t217, 0; 0, 0, -t230 * t243, t222, t222, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end