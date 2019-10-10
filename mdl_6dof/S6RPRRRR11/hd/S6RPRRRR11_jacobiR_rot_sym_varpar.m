% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRRRR11
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
% Datum: 2019-10-10 09:13
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRRR11_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR11_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRR11_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR11_jacobiR_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:13:37
	% EndTime: 2019-10-10 09:13:37
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:13:37
	% EndTime: 2019-10-10 09:13:37
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
	% StartTime: 2019-10-10 09:13:37
	% EndTime: 2019-10-10 09:13:37
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
	% StartTime: 2019-10-10 09:13:37
	% EndTime: 2019-10-10 09:13:37
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
	% StartTime: 2019-10-10 09:13:37
	% EndTime: 2019-10-10 09:13:37
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
	% StartTime: 2019-10-10 09:13:38
	% EndTime: 2019-10-10 09:13:38
	% DurationCPUTime: 0.39s
	% Computational Cost: add. (237->47), mult. (692->91), div. (0->0), fcn. (956->14), ass. (0->54)
	t218 = cos(pkin(6));
	t213 = sin(pkin(13));
	t226 = cos(qJ(1));
	t234 = t226 * t213;
	t216 = cos(pkin(13));
	t222 = sin(qJ(1));
	t236 = t222 * t216;
	t208 = t218 * t234 + t236;
	t221 = sin(qJ(3));
	t225 = cos(qJ(3));
	t233 = t226 * t216;
	t237 = t222 * t213;
	t207 = -t218 * t233 + t237;
	t217 = cos(pkin(7));
	t214 = sin(pkin(7));
	t215 = sin(pkin(6));
	t240 = t215 * t226;
	t232 = t214 * t240;
	t230 = t207 * t217 + t232;
	t196 = -t208 * t225 + t230 * t221;
	t202 = -t207 * t214 + t217 * t240;
	t220 = sin(qJ(4));
	t224 = cos(qJ(4));
	t188 = t196 * t224 + t202 * t220;
	t219 = sin(qJ(5));
	t250 = t188 * t219;
	t223 = cos(qJ(5));
	t249 = t188 * t223;
	t186 = t196 * t220 - t202 * t224;
	t229 = t218 * t236 + t234;
	t241 = t215 * t222;
	t246 = -t214 * t241 + t229 * t217;
	t245 = t208 * t221;
	t243 = t214 * t218;
	t242 = t215 * t216;
	t239 = t217 * t225;
	t238 = t219 * t224;
	t235 = t223 * t224;
	t227 = t229 * t214 + t217 * t241;
	t209 = -t218 * t237 + t233;
	t206 = -t214 * t242 + t218 * t217;
	t201 = t221 * t243 + (t216 * t217 * t221 + t213 * t225) * t215;
	t200 = t215 * t213 * t221 - t225 * t243 - t239 * t242;
	t198 = t209 * t225 - t246 * t221;
	t197 = t209 * t221 + t246 * t225;
	t195 = -t230 * t225 - t245;
	t193 = t207 * t239 + t225 * t232 + t245;
	t192 = t201 * t224 + t206 * t220;
	t191 = -t201 * t220 + t206 * t224;
	t190 = t198 * t224 + t227 * t220;
	t189 = t198 * t220 - t227 * t224;
	t185 = t190 * t223 + t197 * t219;
	t184 = -t190 * t219 + t197 * t223;
	t1 = [t195 * t219 + t249, 0, -t197 * t235 + t198 * t219, -t189 * t223, t184, 0; t185, 0, -t193 * t235 - t196 * t219, t186 * t223, t193 * t223 + t250, 0; 0, 0, -t200 * t235 + t201 * t219, t191 * t223, -t192 * t219 + t200 * t223, 0; t195 * t223 - t250, 0, t197 * t238 + t198 * t223, t189 * t219, -t185, 0; t184, 0, t193 * t238 - t196 * t223, -t186 * t219, -t193 * t219 + t249, 0; 0, 0, t200 * t238 + t201 * t223, -t191 * t219, -t192 * t223 - t200 * t219, 0; t186, 0, -t197 * t220, t190, 0, 0; t189, 0, -t193 * t220, -t188, 0, 0; 0, 0, -t200 * t220, t192, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:13:38
	% EndTime: 2019-10-10 09:13:38
	% DurationCPUTime: 0.34s
	% Computational Cost: add. (349->49), mult. (870->91), div. (0->0), fcn. (1202->14), ass. (0->59)
	t241 = cos(pkin(6));
	t236 = sin(pkin(13));
	t247 = cos(qJ(1));
	t255 = t247 * t236;
	t239 = cos(pkin(13));
	t244 = sin(qJ(1));
	t256 = t244 * t239;
	t228 = t241 * t255 + t256;
	t243 = sin(qJ(3));
	t246 = cos(qJ(3));
	t254 = t247 * t239;
	t257 = t244 * t236;
	t227 = -t241 * t254 + t257;
	t240 = cos(pkin(7));
	t237 = sin(pkin(7));
	t238 = sin(pkin(6));
	t259 = t238 * t247;
	t253 = t237 * t259;
	t251 = t227 * t240 + t253;
	t216 = -t228 * t246 + t251 * t243;
	t222 = -t227 * t237 + t240 * t259;
	t242 = sin(qJ(4));
	t245 = cos(qJ(4));
	t208 = t216 * t245 + t222 * t242;
	t235 = qJ(5) + qJ(6);
	t233 = sin(t235);
	t271 = t208 * t233;
	t234 = cos(t235);
	t270 = t208 * t234;
	t206 = t216 * t242 - t222 * t245;
	t250 = t241 * t256 + t255;
	t260 = t238 * t244;
	t267 = -t237 * t260 + t250 * t240;
	t266 = t228 * t243;
	t264 = t233 * t245;
	t263 = t234 * t245;
	t262 = t237 * t241;
	t261 = t238 * t239;
	t258 = t240 * t246;
	t248 = t250 * t237 + t240 * t260;
	t229 = -t241 * t257 + t254;
	t226 = -t237 * t261 + t241 * t240;
	t221 = t243 * t262 + (t239 * t240 * t243 + t236 * t246) * t238;
	t220 = t238 * t236 * t243 - t246 * t262 - t258 * t261;
	t218 = t229 * t246 - t267 * t243;
	t217 = t229 * t243 + t267 * t246;
	t215 = -t251 * t246 - t266;
	t213 = t227 * t258 + t246 * t253 + t266;
	t212 = t221 * t245 + t226 * t242;
	t211 = -t221 * t242 + t226 * t245;
	t210 = t218 * t245 + t248 * t242;
	t209 = t218 * t242 - t248 * t245;
	t205 = -t212 * t234 - t220 * t233;
	t204 = -t212 * t233 + t220 * t234;
	t203 = t210 * t234 + t217 * t233;
	t202 = -t210 * t233 + t217 * t234;
	t201 = -t213 * t233 + t270;
	t200 = t213 * t234 + t271;
	t1 = [t215 * t233 + t270, 0, -t217 * t263 + t218 * t233, -t209 * t234, t202, t202; t203, 0, -t213 * t263 - t216 * t233, t206 * t234, t200, t200; 0, 0, -t220 * t263 + t221 * t233, t211 * t234, t204, t204; t215 * t234 - t271, 0, t217 * t264 + t218 * t234, t209 * t233, -t203, -t203; t202, 0, t213 * t264 - t216 * t234, -t206 * t233, t201, t201; 0, 0, t220 * t264 + t221 * t234, -t211 * t233, t205, t205; t206, 0, -t217 * t242, t210, 0, 0; t209, 0, -t213 * t242, -t208, 0, 0; 0, 0, -t220 * t242, t212, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end