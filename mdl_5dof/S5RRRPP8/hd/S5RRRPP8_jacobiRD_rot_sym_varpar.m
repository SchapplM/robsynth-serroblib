% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRRPP8
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
% Datum: 2019-12-29 19:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RRRPP8_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP8_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP8_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPP8_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP8_jacobiRD_rot_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:53:57
	% EndTime: 2019-12-29 19:53:57
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:53:57
	% EndTime: 2019-12-29 19:53:57
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0; -t31, 0, 0, 0, 0; 0, 0, 0, 0, 0; t31, 0, 0, 0, 0; -t30, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:53:58
	% EndTime: 2019-12-29 19:53:58
	% DurationCPUTime: 0.05s
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
	% StartTime: 2019-12-29 19:53:59
	% EndTime: 2019-12-29 19:53:59
	% DurationCPUTime: 0.17s
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
	% StartTime: 2019-12-29 19:54:00
	% EndTime: 2019-12-29 19:54:00
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (48->26), mult. (173->50), div. (0->0), fcn. (173->6), ass. (0->33)
	t290 = sin(qJ(1));
	t292 = cos(qJ(2));
	t313 = qJD(1) * t292;
	t298 = -qJD(3) + t313;
	t289 = sin(qJ(2));
	t293 = cos(qJ(1));
	t309 = qJD(2) * t293;
	t301 = t289 * t309;
	t315 = t298 * t290 + t301;
	t314 = qJD(1) * t290;
	t312 = qJD(1) * t293;
	t311 = qJD(2) * t289;
	t310 = qJD(2) * t292;
	t288 = sin(qJ(3));
	t308 = qJD(3) * t288;
	t307 = qJD(3) * t289;
	t306 = qJD(3) * t292;
	t291 = cos(qJ(3));
	t305 = t291 * t311;
	t304 = t291 * t307;
	t303 = t290 * t311;
	t302 = t290 * t310;
	t300 = t292 * t309;
	t299 = -qJD(1) + t306;
	t297 = t299 * t293;
	t296 = t298 * t293;
	t295 = -t289 * t312 - t302;
	t294 = t289 * t314 - t300;
	t287 = t291 * t296 + (-t299 * t288 - t305) * t290;
	t286 = t299 * t291 * t290 + (t296 - t303) * t288;
	t285 = t288 * t297 + t315 * t291;
	t284 = -t315 * t288 + t291 * t297;
	t1 = [t295, -t290 * t313 - t301, 0, 0, 0; -t294, t292 * t312 - t303, 0, 0, 0; 0, t310, 0, 0, 0; t287, t291 * t300 + (-t291 * t314 - t293 * t308) * t289, t284, 0, 0; t285, t291 * t302 + (-t290 * t308 + t291 * t312) * t289, t286, 0, 0; 0, t288 * t306 + t305, t288 * t310 + t304, 0, 0; -t286, t294 * t288 - t293 * t304, -t285, 0, 0; t284, t295 * t288 - t290 * t304, t287, 0, 0; 0, -t288 * t311 + t291 * t306, -t288 * t307 + t291 * t310, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:54:00
	% EndTime: 2019-12-29 19:54:01
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (48->26), mult. (173->50), div. (0->0), fcn. (173->6), ass. (0->33)
	t292 = cos(qJ(1));
	t291 = cos(qJ(2));
	t311 = qJD(1) * t291;
	t296 = -qJD(3) + t311;
	t314 = t296 * t292;
	t289 = sin(qJ(1));
	t288 = sin(qJ(2));
	t307 = qJD(2) * t292;
	t299 = t288 * t307;
	t313 = t296 * t289 + t299;
	t312 = qJD(1) * t289;
	t310 = qJD(1) * t292;
	t309 = qJD(2) * t288;
	t308 = qJD(2) * t291;
	t287 = sin(qJ(3));
	t306 = qJD(3) * t287;
	t305 = qJD(3) * t288;
	t304 = qJD(3) * t291;
	t290 = cos(qJ(3));
	t303 = t290 * t309;
	t302 = t290 * t305;
	t301 = t289 * t309;
	t300 = t289 * t308;
	t298 = t291 * t307;
	t297 = qJD(1) - t304;
	t295 = t297 * t292;
	t294 = -t288 * t310 - t300;
	t293 = t288 * t312 - t298;
	t286 = t290 * t314 + (t297 * t287 - t303) * t289;
	t285 = t297 * t290 * t289 + (t301 - t314) * t287;
	t284 = t287 * t295 - t313 * t290;
	t283 = t313 * t287 + t290 * t295;
	t1 = [t294, -t289 * t311 - t299, 0, 0, 0; -t293, t291 * t310 - t301, 0, 0, 0; 0, t308, 0, 0, 0; t285, t293 * t287 - t292 * t302, t284, 0, 0; -t283, t294 * t287 - t289 * t302, t286, 0, 0; 0, -t287 * t309 + t290 * t304, -t287 * t305 + t290 * t308, 0, 0; -t286, -t290 * t298 + (t290 * t312 + t292 * t306) * t288, t283, 0, 0; t284, -t290 * t300 + (t289 * t306 - t290 * t310) * t288, t285, 0, 0; 0, -t287 * t304 - t303, -t287 * t308 - t302, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end