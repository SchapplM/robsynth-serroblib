% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PPPRRR1
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 08:49
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PPPRRR1_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPPRRR1_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_jacobiR_rot_sym_varpar: pkin has to be [14x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:49:21
	% EndTime: 2019-10-10 08:49:21
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:49:21
	% EndTime: 2019-10-10 08:49:21
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:49:21
	% EndTime: 2019-10-10 08:49:21
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:49:21
	% EndTime: 2019-10-10 08:49:21
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:49:21
	% EndTime: 2019-10-10 08:49:21
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (62->25), mult. (184->54), div. (0->0), fcn. (252->14), ass. (0->34)
	t95 = sin(pkin(13));
	t99 = sin(pkin(6));
	t119 = t95 * t99;
	t96 = sin(pkin(12));
	t118 = t96 * t95;
	t98 = sin(pkin(7));
	t117 = t98 * t99;
	t104 = cos(pkin(7));
	t116 = t104 * t99;
	t101 = cos(pkin(13));
	t115 = t96 * t101;
	t102 = cos(pkin(12));
	t105 = cos(pkin(6));
	t114 = t102 * t105;
	t100 = cos(pkin(14));
	t103 = cos(pkin(8));
	t90 = t101 * t114 - t118;
	t109 = -t102 * t117 + t104 * t90;
	t91 = t95 * t114 + t115;
	t94 = sin(pkin(14));
	t97 = sin(pkin(8));
	t113 = t103 * (t109 * t100 - t91 * t94) + (-t102 * t116 - t90 * t98) * t97;
	t92 = -t102 * t95 - t105 * t115;
	t110 = t104 * t92 + t96 * t117;
	t93 = t102 * t101 - t105 * t118;
	t112 = t103 * (t110 * t100 - t93 * t94) + (t96 * t116 - t92 * t98) * t97;
	t108 = t101 * t116 + t105 * t98;
	t111 = t103 * (t108 * t100 - t94 * t119) + (-t101 * t117 + t105 * t104) * t97;
	t107 = cos(qJ(4));
	t106 = sin(qJ(4));
	t86 = t100 * t119 + t108 * t94;
	t84 = t93 * t100 + t110 * t94;
	t82 = t91 * t100 + t109 * t94;
	t1 = [0, 0, 0, -t84 * t106 + t112 * t107, 0, 0; 0, 0, 0, -t82 * t106 + t113 * t107, 0, 0; 0, 0, 0, -t86 * t106 + t111 * t107, 0, 0; 0, 0, 0, -t112 * t106 - t84 * t107, 0, 0; 0, 0, 0, -t113 * t106 - t82 * t107, 0, 0; 0, 0, 0, -t111 * t106 - t86 * t107, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:49:22
	% EndTime: 2019-10-10 08:49:22
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (200->38), mult. (582->82), div. (0->0), fcn. (794->16), ass. (0->50)
	t177 = sin(pkin(12));
	t186 = cos(pkin(6));
	t201 = t177 * t186;
	t179 = sin(pkin(7));
	t180 = sin(pkin(6));
	t200 = t179 * t180;
	t199 = t179 * t186;
	t185 = cos(pkin(7));
	t198 = t180 * t185;
	t182 = cos(pkin(13));
	t197 = t182 * t185;
	t183 = cos(pkin(12));
	t196 = t183 * t186;
	t176 = sin(pkin(13));
	t172 = t176 * t196 + t177 * t182;
	t175 = sin(pkin(14));
	t181 = cos(pkin(14));
	t171 = -t177 * t176 + t182 * t196;
	t192 = t171 * t185 - t183 * t200;
	t161 = -t172 * t175 + t192 * t181;
	t168 = -t171 * t179 - t183 * t198;
	t178 = sin(pkin(8));
	t184 = cos(pkin(8));
	t195 = t161 * t184 + t168 * t178;
	t174 = -t176 * t201 + t183 * t182;
	t173 = -t183 * t176 - t182 * t201;
	t191 = t173 * t185 + t177 * t200;
	t163 = -t174 * t175 + t191 * t181;
	t169 = -t173 * t179 + t177 * t198;
	t194 = t163 * t184 + t169 * t178;
	t166 = t181 * t199 + (-t175 * t176 + t181 * t197) * t180;
	t170 = -t182 * t200 + t186 * t185;
	t193 = t166 * t184 + t170 * t178;
	t190 = cos(qJ(4));
	t189 = cos(qJ(5));
	t188 = sin(qJ(4));
	t187 = sin(qJ(5));
	t167 = t180 * t176 * t181 + (t180 * t197 + t199) * t175;
	t165 = -t166 * t178 + t170 * t184;
	t164 = t174 * t181 + t191 * t175;
	t162 = t172 * t181 + t192 * t175;
	t160 = -t163 * t178 + t169 * t184;
	t159 = -t161 * t178 + t168 * t184;
	t158 = t167 * t190 + t193 * t188;
	t157 = -t167 * t188 + t193 * t190;
	t156 = t164 * t190 + t194 * t188;
	t155 = -t164 * t188 + t194 * t190;
	t154 = t162 * t190 + t195 * t188;
	t153 = -t162 * t188 + t195 * t190;
	t1 = [0, 0, 0, t155 * t189, -t156 * t187 + t160 * t189, 0; 0, 0, 0, t153 * t189, -t154 * t187 + t159 * t189, 0; 0, 0, 0, t157 * t189, -t158 * t187 + t165 * t189, 0; 0, 0, 0, -t155 * t187, -t156 * t189 - t160 * t187, 0; 0, 0, 0, -t153 * t187, -t154 * t189 - t159 * t187, 0; 0, 0, 0, -t157 * t187, -t158 * t189 - t165 * t187, 0; 0, 0, 0, t156, 0, 0; 0, 0, 0, t154, 0, 0; 0, 0, 0, t158, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:49:22
	% EndTime: 2019-10-10 08:49:23
	% DurationCPUTime: 0.33s
	% Computational Cost: add. (492->53), mult. (1433->111), div. (0->0), fcn. (1946->18), ass. (0->60)
	t242 = sin(pkin(8));
	t248 = cos(pkin(8));
	t239 = sin(pkin(14));
	t240 = sin(pkin(13));
	t244 = sin(pkin(6));
	t245 = cos(pkin(14));
	t246 = cos(pkin(13));
	t249 = cos(pkin(7));
	t273 = t246 * t249;
	t243 = sin(pkin(7));
	t250 = cos(pkin(6));
	t275 = t243 * t250;
	t259 = t245 * t275 + (-t239 * t240 + t245 * t273) * t244;
	t276 = t243 * t244;
	t265 = -t246 * t276 + t250 * t249;
	t281 = t265 * t242 + t259 * t248;
	t247 = cos(pkin(12));
	t241 = sin(pkin(12));
	t277 = t241 * t250;
	t238 = -t240 * t277 + t247 * t246;
	t237 = -t247 * t240 - t246 * t277;
	t266 = t237 * t249 + t241 * t276;
	t260 = t238 * t239 - t266 * t245;
	t274 = t244 * t249;
	t267 = -t237 * t243 + t241 * t274;
	t280 = -t267 * t242 + t260 * t248;
	t272 = t247 * t250;
	t236 = t240 * t272 + t241 * t246;
	t235 = -t241 * t240 + t246 * t272;
	t268 = t235 * t249 - t247 * t276;
	t261 = t236 * t239 - t268 * t245;
	t269 = -t235 * t243 - t247 * t274;
	t279 = -t269 * t242 + t261 * t248;
	t278 = cos(qJ(4));
	t251 = sin(qJ(6));
	t255 = cos(qJ(5));
	t271 = t251 * t255;
	t254 = cos(qJ(6));
	t270 = t254 * t255;
	t253 = sin(qJ(4));
	t252 = sin(qJ(5));
	t233 = t244 * t240 * t245 + (t244 * t273 + t275) * t239;
	t229 = -t259 * t242 + t265 * t248;
	t228 = t238 * t245 + t266 * t239;
	t227 = t236 * t245 + t268 * t239;
	t224 = t260 * t242 + t267 * t248;
	t223 = t261 * t242 + t269 * t248;
	t222 = t233 * t278 + t281 * t253;
	t221 = t233 * t253 - t281 * t278;
	t220 = t228 * t278 - t280 * t253;
	t219 = t228 * t253 + t280 * t278;
	t218 = t227 * t278 - t279 * t253;
	t217 = t227 * t253 + t279 * t278;
	t216 = t222 * t255 + t229 * t252;
	t215 = -t222 * t252 + t229 * t255;
	t214 = t220 * t255 + t224 * t252;
	t213 = -t220 * t252 + t224 * t255;
	t212 = t218 * t255 + t223 * t252;
	t211 = -t218 * t252 + t223 * t255;
	t1 = [0, 0, 0, -t219 * t270 + t220 * t251, t213 * t254, -t214 * t251 + t219 * t254; 0, 0, 0, -t217 * t270 + t218 * t251, t211 * t254, -t212 * t251 + t217 * t254; 0, 0, 0, -t221 * t270 + t222 * t251, t215 * t254, -t216 * t251 + t221 * t254; 0, 0, 0, t219 * t271 + t220 * t254, -t213 * t251, -t214 * t254 - t219 * t251; 0, 0, 0, t217 * t271 + t218 * t254, -t211 * t251, -t212 * t254 - t217 * t251; 0, 0, 0, t221 * t271 + t222 * t254, -t215 * t251, -t216 * t254 - t221 * t251; 0, 0, 0, -t219 * t252, t214, 0; 0, 0, 0, -t217 * t252, t212, 0; 0, 0, 0, -t221 * t252, t216, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end