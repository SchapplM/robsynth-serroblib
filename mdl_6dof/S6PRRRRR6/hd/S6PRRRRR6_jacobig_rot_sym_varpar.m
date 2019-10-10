% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PRRRRR6
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:23
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRRRR6_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR6_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRRR6_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRRRRR6_jacobig_rot_sym_varpar: pkin has to be [14x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:23:35
	% EndTime: 2019-10-09 23:23:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:23:35
	% EndTime: 2019-10-09 23:23:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:23:35
	% EndTime: 2019-10-09 23:23:35
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t18 = sin(pkin(6));
	t1 = [0, sin(pkin(14)) * t18, 0, 0, 0, 0; 0, -cos(pkin(14)) * t18, 0, 0, 0, 0; 0, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:23:36
	% EndTime: 2019-10-09 23:23:36
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (6->6), mult. (17->14), div. (0->0), fcn. (28->8), ass. (0->12)
	t64 = sin(pkin(14));
	t66 = sin(pkin(6));
	t74 = t64 * t66;
	t67 = cos(pkin(14));
	t73 = t67 * t66;
	t69 = cos(pkin(6));
	t71 = cos(qJ(2));
	t72 = t69 * t71;
	t70 = sin(qJ(2));
	t68 = cos(pkin(7));
	t65 = sin(pkin(7));
	t1 = [0, t74, -(-t64 * t72 - t67 * t70) * t65 + t68 * t74, 0, 0, 0; 0, -t73, -(-t64 * t70 + t67 * t72) * t65 - t68 * t73, 0, 0, 0; 0, t69, -t66 * t71 * t65 + t69 * t68, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:23:36
	% EndTime: 2019-10-09 23:23:36
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (24->17), mult. (69->39), div. (0->0), fcn. (101->12), ass. (0->23)
	t142 = sin(pkin(14));
	t145 = sin(pkin(6));
	t158 = t142 * t145;
	t144 = sin(pkin(7));
	t157 = t144 * t145;
	t146 = cos(pkin(14));
	t156 = t146 * t145;
	t149 = cos(pkin(6));
	t151 = sin(qJ(2));
	t155 = t149 * t151;
	t153 = cos(qJ(2));
	t154 = t149 * t153;
	t152 = cos(qJ(3));
	t150 = sin(qJ(3));
	t148 = cos(pkin(7));
	t147 = cos(pkin(8));
	t143 = sin(pkin(8));
	t141 = -t142 * t154 - t146 * t151;
	t140 = -t142 * t151 + t146 * t154;
	t139 = t149 * t148 - t153 * t157;
	t138 = -t141 * t144 + t148 * t158;
	t137 = -t140 * t144 - t148 * t156;
	t1 = [0, t158, t138, -(-(-t142 * t155 + t146 * t153) * t150 + (t141 * t148 + t142 * t157) * t152) * t143 + t138 * t147, 0, 0; 0, -t156, t137, -(-(t142 * t153 + t146 * t155) * t150 + (t140 * t148 - t144 * t156) * t152) * t143 + t137 * t147, 0, 0; 0, t149, t139, -(t149 * t144 * t152 + (t148 * t152 * t153 - t150 * t151) * t145) * t143 + t139 * t147, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:23:37
	% EndTime: 2019-10-09 23:23:37
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (55->27), mult. (161->59), div. (0->0), fcn. (227->14), ass. (0->34)
	t201 = sin(pkin(14));
	t204 = sin(pkin(6));
	t223 = t201 * t204;
	t203 = sin(pkin(7));
	t222 = t203 * t204;
	t208 = cos(pkin(6));
	t221 = t203 * t208;
	t205 = cos(pkin(14));
	t220 = t205 * t204;
	t207 = cos(pkin(7));
	t214 = cos(qJ(2));
	t219 = t207 * t214;
	t211 = sin(qJ(2));
	t218 = t208 * t211;
	t217 = t208 * t214;
	t197 = -t201 * t211 + t205 * t217;
	t216 = t197 * t207 - t203 * t220;
	t199 = -t201 * t217 - t205 * t211;
	t215 = t199 * t207 + t201 * t222;
	t213 = cos(qJ(3));
	t212 = cos(qJ(4));
	t210 = sin(qJ(3));
	t209 = sin(qJ(4));
	t206 = cos(pkin(8));
	t202 = sin(pkin(8));
	t200 = -t201 * t218 + t205 * t214;
	t198 = t201 * t214 + t205 * t218;
	t196 = t208 * t207 - t214 * t222;
	t195 = -t199 * t203 + t207 * t223;
	t194 = -t197 * t203 - t207 * t220;
	t193 = t213 * t221 + (-t210 * t211 + t213 * t219) * t204;
	t192 = -t200 * t210 + t215 * t213;
	t191 = -t198 * t210 + t216 * t213;
	t1 = [0, t223, t195, -t192 * t202 + t195 * t206, (t200 * t213 + t215 * t210) * t209 + (-t192 * t206 - t195 * t202) * t212, 0; 0, -t220, t194, -t191 * t202 + t194 * t206, (t198 * t213 + t216 * t210) * t209 + (-t191 * t206 - t194 * t202) * t212, 0; 0, t208, t196, -t193 * t202 + t196 * t206, (t210 * t221 + (t210 * t219 + t211 * t213) * t204) * t209 + (-t193 * t206 - t196 * t202) * t212, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:23:37
	% EndTime: 2019-10-09 23:23:38
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (107->33), mult. (311->71), div. (0->0), fcn. (432->16), ass. (0->45)
	t255 = sin(pkin(14));
	t258 = sin(pkin(6));
	t282 = t255 * t258;
	t257 = sin(pkin(7));
	t281 = t257 * t258;
	t262 = cos(pkin(6));
	t280 = t257 * t262;
	t259 = cos(pkin(14));
	t279 = t259 * t258;
	t261 = cos(pkin(7));
	t270 = cos(qJ(2));
	t278 = t261 * t270;
	t266 = sin(qJ(2));
	t277 = t262 * t266;
	t276 = t262 * t270;
	t252 = t255 * t270 + t259 * t277;
	t265 = sin(qJ(3));
	t269 = cos(qJ(3));
	t251 = -t255 * t266 + t259 * t276;
	t272 = t251 * t261 - t257 * t279;
	t242 = -t252 * t265 + t272 * t269;
	t248 = -t251 * t257 - t261 * t279;
	t256 = sin(pkin(8));
	t260 = cos(pkin(8));
	t275 = t242 * t260 + t248 * t256;
	t254 = -t255 * t277 + t259 * t270;
	t253 = -t255 * t276 - t259 * t266;
	t271 = t253 * t261 + t255 * t281;
	t244 = -t254 * t265 + t271 * t269;
	t249 = -t253 * t257 + t261 * t282;
	t274 = t244 * t260 + t249 * t256;
	t246 = t269 * t280 + (-t265 * t266 + t269 * t278) * t258;
	t250 = t262 * t261 - t270 * t281;
	t273 = t246 * t260 + t250 * t256;
	t268 = cos(qJ(4));
	t267 = cos(qJ(5));
	t264 = sin(qJ(4));
	t263 = sin(qJ(5));
	t247 = t265 * t280 + (t265 * t278 + t266 * t269) * t258;
	t245 = t254 * t269 + t271 * t265;
	t243 = t252 * t269 + t272 * t265;
	t241 = -t246 * t256 + t250 * t260;
	t240 = -t244 * t256 + t249 * t260;
	t239 = -t242 * t256 + t248 * t260;
	t1 = [0, t282, t249, t240, t245 * t264 - t274 * t268, (t245 * t268 + t274 * t264) * t263 - t240 * t267; 0, -t279, t248, t239, t243 * t264 - t275 * t268, (t243 * t268 + t275 * t264) * t263 - t239 * t267; 0, t262, t250, t241, t247 * t264 - t273 * t268, (t247 * t268 + t273 * t264) * t263 - t241 * t267;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end