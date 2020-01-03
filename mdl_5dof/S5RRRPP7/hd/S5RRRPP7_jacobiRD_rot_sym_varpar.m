% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RRRPP7_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP7_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP7_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPP7_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP7_jacobiRD_rot_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:06:38
	% EndTime: 2019-12-31 21:06:39
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:06:39
	% EndTime: 2019-12-31 21:06:39
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0; -t31, 0, 0, 0, 0; 0, 0, 0, 0, 0; t31, 0, 0, 0, 0; -t30, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:06:38
	% EndTime: 2019-12-31 21:06:39
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->9), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->14)
	t34 = sin(qJ(1));
	t41 = qJD(1) * t34;
	t36 = cos(qJ(1));
	t40 = qJD(1) * t36;
	t33 = sin(qJ(2));
	t39 = qJD(2) * t33;
	t35 = cos(qJ(2));
	t38 = qJD(2) * t35;
	t37 = qJD(2) * t36;
	t32 = t34 * t39 - t35 * t40;
	t31 = t33 * t40 + t34 * t38;
	t30 = t33 * t37 + t35 * t41;
	t29 = t33 * t41 - t35 * t37;
	t1 = [t32, t29, 0, 0, 0; -t30, -t31, 0, 0, 0; 0, -t39, 0, 0, 0; t31, t30, 0, 0, 0; t29, t32, 0, 0, 0; 0, -t38, 0, 0, 0; -t41, 0, 0, 0, 0; t40, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:06:40
	% EndTime: 2019-12-31 21:06:40
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (48->26), mult. (173->50), div. (0->0), fcn. (173->6), ass. (0->33)
	t232 = cos(qJ(3));
	t234 = cos(qJ(1));
	t256 = t232 * t234;
	t231 = sin(qJ(1));
	t255 = qJD(1) * t231;
	t233 = cos(qJ(2));
	t254 = qJD(1) * t233;
	t253 = qJD(1) * t234;
	t230 = sin(qJ(2));
	t252 = qJD(2) * t230;
	t251 = qJD(2) * t233;
	t250 = qJD(2) * t234;
	t229 = sin(qJ(3));
	t249 = qJD(3) * t229;
	t248 = qJD(3) * t230;
	t247 = qJD(3) * t233;
	t246 = t232 * t252;
	t245 = t232 * t248;
	t244 = t231 * t252;
	t243 = t231 * t251;
	t242 = t230 * t250;
	t241 = t233 * t250;
	t240 = -qJD(1) + t247;
	t239 = -qJD(3) + t254;
	t238 = t240 * t229;
	t237 = t230 * t253 + t243;
	t236 = -t230 * t255 + t241;
	t235 = t239 * t231 + t242;
	t228 = -t239 * t256 + (t238 + t246) * t231;
	t227 = t240 * t232 * t231 + (t239 * t234 - t244) * t229;
	t226 = t235 * t232 + t234 * t238;
	t225 = t235 * t229 - t240 * t256;
	t1 = [t228, -t232 * t241 + (t232 * t255 + t234 * t249) * t230, t225, 0, 0; -t226, -t232 * t243 + (t231 * t249 - t232 * t253) * t230, -t227, 0, 0; 0, -t229 * t247 - t246, -t229 * t251 - t245, 0, 0; t227, t236 * t229 + t234 * t245, t226, 0, 0; t225, t237 * t229 + t231 * t245, t228, 0, 0; 0, t229 * t252 - t232 * t247, t229 * t248 - t232 * t251, 0, 0; -t237, -t231 * t254 - t242, 0, 0, 0; t236, t233 * t253 - t244, 0, 0, 0; 0, t251, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:06:40
	% EndTime: 2019-12-31 21:06:40
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (48->26), mult. (173->50), div. (0->0), fcn. (173->6), ass. (0->33)
	t296 = cos(qJ(1));
	t295 = cos(qJ(2));
	t315 = qJD(1) * t295;
	t300 = -qJD(3) + t315;
	t318 = t300 * t296;
	t293 = sin(qJ(1));
	t292 = sin(qJ(2));
	t311 = qJD(2) * t296;
	t303 = t292 * t311;
	t317 = t300 * t293 + t303;
	t316 = qJD(1) * t293;
	t314 = qJD(1) * t296;
	t313 = qJD(2) * t292;
	t312 = qJD(2) * t295;
	t291 = sin(qJ(3));
	t310 = qJD(3) * t291;
	t309 = qJD(3) * t292;
	t308 = qJD(3) * t295;
	t294 = cos(qJ(3));
	t307 = t294 * t313;
	t306 = t294 * t309;
	t305 = t293 * t313;
	t304 = t293 * t312;
	t302 = t295 * t311;
	t301 = qJD(1) - t308;
	t299 = t301 * t296;
	t298 = -t292 * t314 - t304;
	t297 = t292 * t316 - t302;
	t290 = t294 * t318 + (t301 * t291 - t307) * t293;
	t289 = t301 * t294 * t293 + (t305 - t318) * t291;
	t288 = t291 * t299 - t317 * t294;
	t287 = t317 * t291 + t294 * t299;
	t1 = [-t290, -t294 * t302 + (t294 * t316 + t296 * t310) * t292, t287, 0, 0; t288, -t294 * t304 + (t293 * t310 - t294 * t314) * t292, t289, 0, 0; 0, -t291 * t308 - t307, -t291 * t312 - t306, 0, 0; t298, -t293 * t315 - t303, 0, 0, 0; -t297, t295 * t314 - t305, 0, 0, 0; 0, t312, 0, 0, 0; t289, t297 * t291 - t296 * t306, t288, 0, 0; -t287, t298 * t291 - t293 * t306, t290, 0, 0; 0, -t291 * t313 + t294 * t308, -t291 * t309 + t294 * t312, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:06:40
	% EndTime: 2019-12-31 21:06:40
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (49->26), mult. (173->50), div. (0->0), fcn. (173->6), ass. (0->33)
	t267 = cos(qJ(1));
	t266 = cos(qJ(2));
	t286 = qJD(1) * t266;
	t271 = -qJD(3) + t286;
	t289 = t271 * t267;
	t264 = sin(qJ(1));
	t263 = sin(qJ(2));
	t282 = qJD(2) * t267;
	t274 = t263 * t282;
	t288 = t271 * t264 + t274;
	t287 = qJD(1) * t264;
	t285 = qJD(1) * t267;
	t284 = qJD(2) * t263;
	t283 = qJD(2) * t266;
	t262 = sin(qJ(3));
	t281 = qJD(3) * t262;
	t280 = qJD(3) * t263;
	t279 = qJD(3) * t266;
	t265 = cos(qJ(3));
	t278 = t265 * t284;
	t277 = t265 * t280;
	t276 = t264 * t284;
	t275 = t264 * t283;
	t273 = t266 * t282;
	t272 = qJD(1) - t279;
	t270 = t272 * t267;
	t269 = t263 * t285 + t275;
	t268 = t263 * t287 - t273;
	t261 = t265 * t289 + (t272 * t262 - t278) * t264;
	t260 = t272 * t265 * t264 + (t276 - t289) * t262;
	t259 = t262 * t270 - t288 * t265;
	t258 = t288 * t262 + t265 * t270;
	t1 = [-t261, -t265 * t273 + (t265 * t287 + t267 * t281) * t263, t258, 0, 0; t259, -t265 * t275 + (t264 * t281 - t265 * t285) * t263, t260, 0, 0; 0, -t262 * t279 - t278, -t262 * t283 - t277, 0, 0; t260, t268 * t262 - t267 * t277, t259, 0, 0; -t258, -t269 * t262 - t264 * t277, t261, 0, 0; 0, -t262 * t284 + t265 * t279, -t262 * t280 + t265 * t283, 0, 0; t269, t264 * t286 + t274, 0, 0, 0; t268, -t266 * t285 + t276, 0, 0, 0; 0, -t283, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end