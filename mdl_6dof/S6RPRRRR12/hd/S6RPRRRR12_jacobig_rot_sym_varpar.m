% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRRR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:16
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RPRRRR12_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR12_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRR12_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_jacobig_rot_sym_varpar: pkin has to be [14x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:15:46
	% EndTime: 2019-10-10 09:15:46
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:15:46
	% EndTime: 2019-10-10 09:15:46
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:15:46
	% EndTime: 2019-10-10 09:15:46
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:15:46
	% EndTime: 2019-10-10 09:15:46
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (5->5), mult. (15->13), div. (0->0), fcn. (23->8), ass. (0->11)
	t78 = sin(pkin(6));
	t80 = cos(pkin(7));
	t85 = t78 * t80;
	t79 = cos(pkin(14));
	t81 = cos(pkin(6));
	t84 = t79 * t81;
	t83 = cos(qJ(1));
	t82 = sin(qJ(1));
	t77 = sin(pkin(7));
	t76 = sin(pkin(14));
	t1 = [0, 0, -(-t83 * t76 - t82 * t84) * t77 + t82 * t85, 0, 0, 0; 0, 0, -(-t82 * t76 + t83 * t84) * t77 - t83 * t85, 0, 0, 0; 1, 0, -t78 * t79 * t77 + t81 * t80, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:15:47
	% EndTime: 2019-10-10 09:15:47
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (23->16), mult. (67->37), div. (0->0), fcn. (96->12), ass. (0->24)
	t142 = sin(pkin(6));
	t148 = sin(qJ(1));
	t156 = t142 * t148;
	t150 = cos(qJ(1));
	t155 = t142 * t150;
	t139 = sin(pkin(14));
	t154 = t148 * t139;
	t143 = cos(pkin(14));
	t153 = t148 * t143;
	t152 = t150 * t139;
	t151 = t150 * t143;
	t149 = cos(qJ(3));
	t147 = sin(qJ(3));
	t146 = cos(pkin(6));
	t145 = cos(pkin(7));
	t144 = cos(pkin(8));
	t141 = sin(pkin(7));
	t140 = sin(pkin(8));
	t138 = -t146 * t153 - t152;
	t137 = t146 * t151 - t154;
	t136 = -t142 * t143 * t141 + t146 * t145;
	t135 = -t138 * t141 + t145 * t156;
	t134 = -t137 * t141 - t145 * t155;
	t1 = [0, 0, t135, -(-(-t146 * t154 + t151) * t147 + (t138 * t145 + t141 * t156) * t149) * t140 + t135 * t144, 0, 0; 0, 0, t134, -(-(t146 * t152 + t153) * t147 + (t137 * t145 - t141 * t155) * t149) * t140 + t134 * t144, 0, 0; 1, 0, t136, -(t146 * t141 * t149 + (t143 * t145 * t149 - t139 * t147) * t142) * t140 + t136 * t144, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:15:47
	% EndTime: 2019-10-10 09:15:48
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (54->26), mult. (159->57), div. (0->0), fcn. (222->14), ass. (0->35)
	t197 = sin(pkin(7));
	t202 = cos(pkin(6));
	t218 = t197 * t202;
	t198 = sin(pkin(6));
	t205 = sin(qJ(1));
	t217 = t198 * t205;
	t208 = cos(qJ(1));
	t216 = t198 * t208;
	t199 = cos(pkin(14));
	t201 = cos(pkin(7));
	t215 = t199 * t201;
	t195 = sin(pkin(14));
	t214 = t205 * t195;
	t213 = t205 * t199;
	t212 = t208 * t195;
	t211 = t208 * t199;
	t191 = t202 * t211 - t214;
	t210 = t191 * t201 - t197 * t216;
	t193 = -t202 * t213 - t212;
	t209 = t193 * t201 + t197 * t217;
	t207 = cos(qJ(3));
	t206 = cos(qJ(4));
	t204 = sin(qJ(3));
	t203 = sin(qJ(4));
	t200 = cos(pkin(8));
	t196 = sin(pkin(8));
	t194 = -t202 * t214 + t211;
	t192 = t202 * t212 + t213;
	t190 = -t198 * t199 * t197 + t202 * t201;
	t189 = -t193 * t197 + t201 * t217;
	t188 = -t191 * t197 - t201 * t216;
	t187 = t207 * t218 + (-t195 * t204 + t207 * t215) * t198;
	t186 = -t194 * t204 + t209 * t207;
	t185 = -t192 * t204 + t210 * t207;
	t1 = [0, 0, t189, -t186 * t196 + t189 * t200, (t194 * t207 + t209 * t204) * t203 + (-t186 * t200 - t189 * t196) * t206, 0; 0, 0, t188, -t185 * t196 + t188 * t200, (t192 * t207 + t210 * t204) * t203 + (-t185 * t200 - t188 * t196) * t206, 0; 1, 0, t190, -t187 * t196 + t190 * t200, (t204 * t218 + (t195 * t207 + t204 * t215) * t198) * t203 + (-t187 * t200 - t190 * t196) * t206, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:15:49
	% EndTime: 2019-10-10 09:15:49
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (106->32), mult. (309->69), div. (0->0), fcn. (427->16), ass. (0->46)
	t259 = sin(pkin(7));
	t264 = cos(pkin(6));
	t285 = t259 * t264;
	t260 = sin(pkin(6));
	t268 = sin(qJ(1));
	t284 = t260 * t268;
	t272 = cos(qJ(1));
	t283 = t260 * t272;
	t261 = cos(pkin(14));
	t263 = cos(pkin(7));
	t282 = t261 * t263;
	t257 = sin(pkin(14));
	t281 = t268 * t257;
	t280 = t268 * t261;
	t279 = t272 * t257;
	t278 = t272 * t261;
	t254 = t264 * t279 + t280;
	t267 = sin(qJ(3));
	t271 = cos(qJ(3));
	t253 = t264 * t278 - t281;
	t274 = t253 * t263 - t259 * t283;
	t244 = -t254 * t267 + t274 * t271;
	t250 = -t253 * t259 - t263 * t283;
	t258 = sin(pkin(8));
	t262 = cos(pkin(8));
	t277 = t244 * t262 + t250 * t258;
	t256 = -t264 * t281 + t278;
	t255 = -t264 * t280 - t279;
	t273 = t255 * t263 + t259 * t284;
	t246 = -t256 * t267 + t273 * t271;
	t251 = -t255 * t259 + t263 * t284;
	t276 = t246 * t262 + t251 * t258;
	t248 = t271 * t285 + (-t257 * t267 + t271 * t282) * t260;
	t252 = -t260 * t261 * t259 + t264 * t263;
	t275 = t248 * t262 + t252 * t258;
	t270 = cos(qJ(4));
	t269 = cos(qJ(5));
	t266 = sin(qJ(4));
	t265 = sin(qJ(5));
	t249 = t267 * t285 + (t257 * t271 + t267 * t282) * t260;
	t247 = t256 * t271 + t273 * t267;
	t245 = t254 * t271 + t274 * t267;
	t243 = -t248 * t258 + t252 * t262;
	t242 = -t246 * t258 + t251 * t262;
	t241 = -t244 * t258 + t250 * t262;
	t1 = [0, 0, t251, t242, t247 * t266 - t276 * t270, (t247 * t270 + t276 * t266) * t265 - t242 * t269; 0, 0, t250, t241, t245 * t266 - t277 * t270, (t245 * t270 + t277 * t266) * t265 - t241 * t269; 1, 0, t252, t243, t249 * t266 - t275 * t270, (t249 * t270 + t275 * t266) * t265 - t243 * t269;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end