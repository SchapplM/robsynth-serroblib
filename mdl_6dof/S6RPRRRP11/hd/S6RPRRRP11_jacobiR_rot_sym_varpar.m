% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRRRP11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 08:56
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRRP11_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP11_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRP11_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP11_jacobiR_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:56:10
	% EndTime: 2019-10-10 08:56:10
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:56:10
	% EndTime: 2019-10-10 08:56:10
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
	% StartTime: 2019-10-10 08:56:10
	% EndTime: 2019-10-10 08:56:10
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
	% StartTime: 2019-10-10 08:56:10
	% EndTime: 2019-10-10 08:56:10
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
	% StartTime: 2019-10-10 08:56:10
	% EndTime: 2019-10-10 08:56:10
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (105->30), mult. (307->59), div. (0->0), fcn. (430->12), ass. (0->40)
	t160 = cos(pkin(6));
	t155 = sin(pkin(12));
	t166 = cos(qJ(1));
	t170 = t166 * t155;
	t158 = cos(pkin(12));
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
	% StartTime: 2019-10-10 08:56:11
	% EndTime: 2019-10-10 08:56:11
	% DurationCPUTime: 0.40s
	% Computational Cost: add. (237->47), mult. (692->91), div. (0->0), fcn. (956->14), ass. (0->54)
	t218 = cos(pkin(6));
	t213 = sin(pkin(12));
	t226 = cos(qJ(1));
	t234 = t226 * t213;
	t216 = cos(pkin(12));
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
	% StartTime: 2019-10-10 08:56:11
	% EndTime: 2019-10-10 08:56:11
	% DurationCPUTime: 0.36s
	% Computational Cost: add. (237->47), mult. (692->91), div. (0->0), fcn. (956->14), ass. (0->54)
	t229 = cos(pkin(6));
	t224 = sin(pkin(12));
	t237 = cos(qJ(1));
	t245 = t237 * t224;
	t227 = cos(pkin(12));
	t233 = sin(qJ(1));
	t247 = t233 * t227;
	t219 = t229 * t245 + t247;
	t232 = sin(qJ(3));
	t236 = cos(qJ(3));
	t244 = t237 * t227;
	t248 = t233 * t224;
	t218 = -t229 * t244 + t248;
	t228 = cos(pkin(7));
	t225 = sin(pkin(7));
	t226 = sin(pkin(6));
	t251 = t226 * t237;
	t243 = t225 * t251;
	t241 = t218 * t228 + t243;
	t207 = -t219 * t236 + t241 * t232;
	t213 = -t218 * t225 + t228 * t251;
	t231 = sin(qJ(4));
	t235 = cos(qJ(4));
	t199 = t207 * t235 + t213 * t231;
	t230 = sin(qJ(5));
	t261 = t199 * t230;
	t234 = cos(qJ(5));
	t260 = t199 * t234;
	t197 = t207 * t231 - t213 * t235;
	t240 = t229 * t247 + t245;
	t252 = t226 * t233;
	t257 = -t225 * t252 + t240 * t228;
	t256 = t219 * t232;
	t254 = t225 * t229;
	t253 = t226 * t227;
	t250 = t228 * t236;
	t249 = t230 * t235;
	t246 = t234 * t235;
	t238 = t240 * t225 + t228 * t252;
	t220 = -t229 * t248 + t244;
	t217 = -t225 * t253 + t229 * t228;
	t212 = t232 * t254 + (t227 * t228 * t232 + t224 * t236) * t226;
	t211 = t226 * t224 * t232 - t236 * t254 - t250 * t253;
	t209 = t220 * t236 - t257 * t232;
	t208 = t220 * t232 + t257 * t236;
	t206 = -t241 * t236 - t256;
	t204 = t218 * t250 + t236 * t243 + t256;
	t203 = t212 * t235 + t217 * t231;
	t202 = -t212 * t231 + t217 * t235;
	t201 = t209 * t235 + t238 * t231;
	t200 = t209 * t231 - t238 * t235;
	t196 = t201 * t234 + t208 * t230;
	t195 = -t201 * t230 + t208 * t234;
	t1 = [t206 * t230 + t260, 0, -t208 * t246 + t209 * t230, -t200 * t234, t195, 0; t196, 0, -t204 * t246 - t207 * t230, t197 * t234, t204 * t234 + t261, 0; 0, 0, -t211 * t246 + t212 * t230, t202 * t234, -t203 * t230 + t211 * t234, 0; t206 * t234 - t261, 0, t208 * t249 + t209 * t234, t200 * t230, -t196, 0; t195, 0, t204 * t249 - t207 * t234, -t197 * t230, -t204 * t230 + t260, 0; 0, 0, t211 * t249 + t212 * t234, -t202 * t230, -t203 * t234 - t211 * t230, 0; t197, 0, -t208 * t231, t201, 0, 0; t200, 0, -t204 * t231, -t199, 0, 0; 0, 0, -t211 * t231, t203, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end