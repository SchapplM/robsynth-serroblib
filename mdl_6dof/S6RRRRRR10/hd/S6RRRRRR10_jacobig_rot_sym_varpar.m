% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRRRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 13:35
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRRRR10_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRR10_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_jacobig_rot_sym_varpar: pkin has to be [14x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:34:53
	% EndTime: 2019-10-10 13:34:53
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:34:53
	% EndTime: 2019-10-10 13:34:53
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:34:53
	% EndTime: 2019-10-10 13:34:53
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t43 = sin(pkin(6));
	t1 = [0, sin(qJ(1)) * t43, 0, 0, 0, 0; 0, -cos(qJ(1)) * t43, 0, 0, 0, 0; 1, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:34:53
	% EndTime: 2019-10-10 13:34:53
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (6->6), mult. (17->14), div. (0->0), fcn. (28->8), ass. (0->12)
	t94 = cos(pkin(6));
	t97 = cos(qJ(2));
	t101 = t94 * t97;
	t92 = sin(pkin(6));
	t96 = sin(qJ(1));
	t100 = t96 * t92;
	t98 = cos(qJ(1));
	t99 = t98 * t92;
	t95 = sin(qJ(2));
	t93 = cos(pkin(7));
	t91 = sin(pkin(7));
	t1 = [0, t100, -(-t96 * t101 - t98 * t95) * t91 + t93 * t100, 0, 0, 0; 0, -t99, -(t98 * t101 - t96 * t95) * t91 - t93 * t99, 0, 0, 0; 1, t94, -t92 * t97 * t91 + t94 * t93, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:34:54
	% EndTime: 2019-10-10 13:34:54
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (24->17), mult. (69->37), div. (0->0), fcn. (101->12), ass. (0->24)
	t175 = sin(pkin(6));
	t181 = sin(qJ(1));
	t190 = t181 * t175;
	t180 = sin(qJ(2));
	t189 = t181 * t180;
	t183 = cos(qJ(2));
	t188 = t181 * t183;
	t184 = cos(qJ(1));
	t187 = t184 * t175;
	t186 = t184 * t180;
	t185 = t184 * t183;
	t182 = cos(qJ(3));
	t179 = sin(qJ(3));
	t178 = cos(pkin(6));
	t177 = cos(pkin(7));
	t176 = cos(pkin(8));
	t174 = sin(pkin(7));
	t173 = sin(pkin(8));
	t172 = -t178 * t188 - t186;
	t171 = t178 * t185 - t189;
	t170 = -t175 * t183 * t174 + t178 * t177;
	t169 = -t172 * t174 + t177 * t190;
	t168 = -t171 * t174 - t177 * t187;
	t1 = [0, t190, t169, -(-(-t178 * t189 + t185) * t179 + (t172 * t177 + t174 * t190) * t182) * t173 + t169 * t176, 0, 0; 0, -t187, t168, -(-(t178 * t186 + t188) * t179 + (t171 * t177 - t174 * t187) * t182) * t173 + t168 * t176, 0, 0; 1, t178, t170, -(t178 * t174 * t182 + (t177 * t182 * t183 - t179 * t180) * t175) * t173 + t170 * t176, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:34:55
	% EndTime: 2019-10-10 13:34:55
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (55->27), mult. (161->57), div. (0->0), fcn. (227->14), ass. (0->35)
	t241 = sin(pkin(7));
	t245 = cos(pkin(6));
	t263 = t241 * t245;
	t244 = cos(pkin(7));
	t252 = cos(qJ(2));
	t262 = t244 * t252;
	t242 = sin(pkin(6));
	t249 = sin(qJ(1));
	t261 = t249 * t242;
	t248 = sin(qJ(2));
	t260 = t249 * t248;
	t259 = t249 * t252;
	t253 = cos(qJ(1));
	t258 = t253 * t242;
	t257 = t253 * t248;
	t256 = t253 * t252;
	t236 = t245 * t256 - t260;
	t255 = t236 * t244 - t241 * t258;
	t238 = -t245 * t259 - t257;
	t254 = t238 * t244 + t241 * t261;
	t251 = cos(qJ(3));
	t250 = cos(qJ(4));
	t247 = sin(qJ(3));
	t246 = sin(qJ(4));
	t243 = cos(pkin(8));
	t240 = sin(pkin(8));
	t239 = -t245 * t260 + t256;
	t237 = t245 * t257 + t259;
	t235 = -t242 * t252 * t241 + t245 * t244;
	t234 = -t238 * t241 + t244 * t261;
	t233 = -t236 * t241 - t244 * t258;
	t232 = t251 * t263 + (-t247 * t248 + t251 * t262) * t242;
	t231 = -t239 * t247 + t254 * t251;
	t230 = -t237 * t247 + t255 * t251;
	t1 = [0, t261, t234, -t231 * t240 + t234 * t243, (t239 * t251 + t254 * t247) * t246 + (-t231 * t243 - t234 * t240) * t250, 0; 0, -t258, t233, -t230 * t240 + t233 * t243, (t237 * t251 + t255 * t247) * t246 + (-t230 * t243 - t233 * t240) * t250, 0; 1, t245, t235, -t232 * t240 + t235 * t243, (t247 * t263 + (t247 * t262 + t248 * t251) * t242) * t246 + (-t232 * t243 - t235 * t240) * t250, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:34:57
	% EndTime: 2019-10-10 13:34:57
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (107->33), mult. (311->69), div. (0->0), fcn. (432->16), ass. (0->46)
	t308 = sin(pkin(7));
	t312 = cos(pkin(6));
	t335 = t308 * t312;
	t311 = cos(pkin(7));
	t321 = cos(qJ(2));
	t334 = t311 * t321;
	t309 = sin(pkin(6));
	t317 = sin(qJ(1));
	t333 = t317 * t309;
	t316 = sin(qJ(2));
	t332 = t317 * t316;
	t331 = t317 * t321;
	t322 = cos(qJ(1));
	t330 = t322 * t309;
	t329 = t322 * t316;
	t328 = t322 * t321;
	t304 = t312 * t329 + t331;
	t315 = sin(qJ(3));
	t320 = cos(qJ(3));
	t303 = t312 * t328 - t332;
	t324 = t303 * t311 - t308 * t330;
	t294 = -t304 * t315 + t320 * t324;
	t300 = -t303 * t308 - t311 * t330;
	t307 = sin(pkin(8));
	t310 = cos(pkin(8));
	t327 = t294 * t310 + t300 * t307;
	t306 = -t312 * t332 + t328;
	t305 = -t312 * t331 - t329;
	t323 = t305 * t311 + t308 * t333;
	t296 = -t306 * t315 + t320 * t323;
	t301 = -t305 * t308 + t311 * t333;
	t326 = t296 * t310 + t301 * t307;
	t298 = t320 * t335 + (-t315 * t316 + t320 * t334) * t309;
	t302 = -t308 * t309 * t321 + t311 * t312;
	t325 = t298 * t310 + t302 * t307;
	t319 = cos(qJ(4));
	t318 = cos(qJ(5));
	t314 = sin(qJ(4));
	t313 = sin(qJ(5));
	t299 = t315 * t335 + (t315 * t334 + t316 * t320) * t309;
	t297 = t306 * t320 + t315 * t323;
	t295 = t304 * t320 + t315 * t324;
	t293 = -t298 * t307 + t302 * t310;
	t292 = -t296 * t307 + t301 * t310;
	t291 = -t294 * t307 + t300 * t310;
	t1 = [0, t333, t301, t292, t297 * t314 - t319 * t326, (t297 * t319 + t314 * t326) * t313 - t292 * t318; 0, -t330, t300, t291, t295 * t314 - t319 * t327, (t295 * t319 + t314 * t327) * t313 - t291 * t318; 1, t312, t302, t293, t299 * t314 - t319 * t325, (t299 * t319 + t314 * t325) * t313 - t293 * t318;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end