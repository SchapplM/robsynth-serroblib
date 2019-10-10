% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRRR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:11
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRPRRR14_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR14_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobig_rot_sym_varpar: pkin has to be [14x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:10:58
	% EndTime: 2019-10-10 11:10:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:10:58
	% EndTime: 2019-10-10 11:10:58
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:10:58
	% EndTime: 2019-10-10 11:10:58
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t43 = sin(pkin(6));
	t1 = [0, sin(qJ(1)) * t43, 0, 0, 0, 0; 0, -cos(qJ(1)) * t43, 0, 0, 0, 0; 1, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:10:58
	% EndTime: 2019-10-10 11:10:58
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t79 = sin(pkin(6));
	t1 = [0, sin(qJ(1)) * t79, 0, 0, 0, 0; 0, -cos(qJ(1)) * t79, 0, 0, 0, 0; 1, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:10:59
	% EndTime: 2019-10-10 11:10:59
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (19->17), mult. (54->37), div. (0->0), fcn. (78->12), ass. (0->21)
	t150 = sin(pkin(6));
	t156 = sin(qJ(1));
	t164 = t156 * t150;
	t155 = sin(qJ(2));
	t163 = t156 * t155;
	t157 = cos(qJ(2));
	t162 = t156 * t157;
	t158 = cos(qJ(1));
	t161 = t158 * t150;
	t160 = t158 * t155;
	t159 = t158 * t157;
	t154 = cos(pkin(6));
	t153 = cos(pkin(7));
	t152 = cos(pkin(8));
	t151 = cos(pkin(14));
	t149 = sin(pkin(7));
	t148 = sin(pkin(8));
	t147 = sin(pkin(14));
	t146 = -t154 * t162 - t160;
	t145 = t154 * t159 - t163;
	t1 = [0, t164, 0, -(-(-t154 * t163 + t159) * t147 + (t146 * t153 + t149 * t164) * t151) * t148 + (-t146 * t149 + t153 * t164) * t152, 0, 0; 0, -t161, 0, -(-(t154 * t160 + t162) * t147 + (t145 * t153 - t149 * t161) * t151) * t148 + (-t145 * t149 - t153 * t161) * t152, 0, 0; 1, t154, 0, -(t154 * t149 * t151 + (t151 * t153 * t157 - t147 * t155) * t150) * t148 + (-t150 * t157 * t149 + t154 * t153) * t152, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:11:00
	% EndTime: 2019-10-10 11:11:00
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (50->27), mult. (146->57), div. (0->0), fcn. (204->14), ass. (0->35)
	t222 = sin(pkin(7));
	t227 = cos(pkin(6));
	t243 = t222 * t227;
	t226 = cos(pkin(7));
	t232 = cos(qJ(2));
	t242 = t226 * t232;
	t223 = sin(pkin(6));
	t230 = sin(qJ(1));
	t241 = t230 * t223;
	t229 = sin(qJ(2));
	t240 = t230 * t229;
	t239 = t230 * t232;
	t233 = cos(qJ(1));
	t238 = t233 * t223;
	t237 = t233 * t229;
	t236 = t233 * t232;
	t216 = t227 * t236 - t240;
	t235 = t216 * t226 - t222 * t238;
	t218 = -t227 * t239 - t237;
	t234 = t218 * t226 + t222 * t241;
	t231 = cos(qJ(4));
	t228 = sin(qJ(4));
	t225 = cos(pkin(8));
	t224 = cos(pkin(14));
	t221 = sin(pkin(8));
	t220 = sin(pkin(14));
	t219 = -t227 * t240 + t236;
	t217 = t227 * t237 + t239;
	t215 = -t223 * t232 * t222 + t227 * t226;
	t214 = -t218 * t222 + t226 * t241;
	t213 = -t216 * t222 - t226 * t238;
	t212 = t224 * t243 + (-t220 * t229 + t224 * t242) * t223;
	t211 = -t219 * t220 + t234 * t224;
	t210 = -t217 * t220 + t235 * t224;
	t1 = [0, t241, 0, -t211 * t221 + t214 * t225, (t219 * t224 + t234 * t220) * t228 + (-t211 * t225 - t214 * t221) * t231, 0; 0, -t238, 0, -t210 * t221 + t213 * t225, (t217 * t224 + t235 * t220) * t228 + (-t210 * t225 - t213 * t221) * t231, 0; 1, t227, 0, -t212 * t221 + t215 * t225, (t223 * t229 * t224 + (t223 * t242 + t243) * t220) * t228 + (-t212 * t225 - t215 * t221) * t231, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:11:01
	% EndTime: 2019-10-10 11:11:01
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (102->33), mult. (296->69), div. (0->0), fcn. (409->16), ass. (0->46)
	t274 = sin(pkin(7));
	t279 = cos(pkin(6));
	t300 = t274 * t279;
	t278 = cos(pkin(7));
	t286 = cos(qJ(2));
	t299 = t278 * t286;
	t275 = sin(pkin(6));
	t283 = sin(qJ(1));
	t298 = t283 * t275;
	t282 = sin(qJ(2));
	t297 = t283 * t282;
	t296 = t283 * t286;
	t287 = cos(qJ(1));
	t295 = t287 * t275;
	t294 = t287 * t282;
	t293 = t287 * t286;
	t269 = t279 * t294 + t296;
	t272 = sin(pkin(14));
	t276 = cos(pkin(14));
	t268 = t279 * t293 - t297;
	t289 = t268 * t278 - t274 * t295;
	t259 = -t269 * t272 + t289 * t276;
	t265 = -t268 * t274 - t278 * t295;
	t273 = sin(pkin(8));
	t277 = cos(pkin(8));
	t292 = t259 * t277 + t265 * t273;
	t271 = -t279 * t297 + t293;
	t270 = -t279 * t296 - t294;
	t288 = t270 * t278 + t274 * t298;
	t261 = -t271 * t272 + t288 * t276;
	t266 = -t270 * t274 + t278 * t298;
	t291 = t261 * t277 + t266 * t273;
	t263 = t276 * t300 + (-t272 * t282 + t276 * t299) * t275;
	t267 = -t275 * t286 * t274 + t279 * t278;
	t290 = t263 * t277 + t267 * t273;
	t285 = cos(qJ(4));
	t284 = cos(qJ(5));
	t281 = sin(qJ(4));
	t280 = sin(qJ(5));
	t264 = t275 * t282 * t276 + (t275 * t299 + t300) * t272;
	t262 = t271 * t276 + t288 * t272;
	t260 = t269 * t276 + t289 * t272;
	t258 = -t263 * t273 + t267 * t277;
	t257 = -t261 * t273 + t266 * t277;
	t256 = -t259 * t273 + t265 * t277;
	t1 = [0, t298, 0, t257, t262 * t281 - t291 * t285, (t262 * t285 + t291 * t281) * t280 - t257 * t284; 0, -t295, 0, t256, t260 * t281 - t292 * t285, (t260 * t285 + t292 * t281) * t280 - t256 * t284; 1, t279, 0, t258, t264 * t281 - t290 * t285, (t264 * t285 + t290 * t281) * t280 - t258 * t284;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end