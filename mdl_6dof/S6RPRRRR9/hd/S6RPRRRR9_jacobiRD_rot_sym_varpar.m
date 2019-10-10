% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:09
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPRRRR9_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR9_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR9_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRR9_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR9_jacobiRD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:09:39
	% EndTime: 2019-10-10 09:09:39
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:09:39
	% EndTime: 2019-10-10 09:09:39
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0, 0; -t31, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t31, 0, 0, 0, 0, 0; -t30, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:09:39
	% EndTime: 2019-10-10 09:09:39
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (1->1), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t13 = qJD(1) * sin(qJ(1));
	t11 = qJD(1) * cos(qJ(1));
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t11, 0, 0, 0, 0, 0; t13, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t13, 0, 0, 0, 0, 0; t11, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:09:39
	% EndTime: 2019-10-10 09:09:39
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (11->9), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->14)
	t35 = sin(qJ(1));
	t42 = qJD(1) * t35;
	t37 = cos(qJ(1));
	t41 = qJD(1) * t37;
	t34 = sin(qJ(3));
	t40 = qJD(3) * t34;
	t36 = cos(qJ(3));
	t39 = qJD(3) * t36;
	t38 = qJD(3) * t37;
	t33 = -t35 * t40 + t36 * t41;
	t32 = t34 * t41 + t35 * t39;
	t31 = t34 * t38 + t36 * t42;
	t30 = -t34 * t42 + t36 * t38;
	t1 = [t30, 0, t33, 0, 0, 0; t32, 0, t31, 0, 0, 0; 0, 0, -t39, 0, 0, 0; -t31, 0, -t32, 0, 0, 0; t33, 0, t30, 0, 0, 0; 0, 0, t40, 0, 0, 0; -t41, 0, 0, 0, 0, 0; -t42, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:09:40
	% EndTime: 2019-10-10 09:09:40
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (49->24), mult. (173->45), div. (0->0), fcn. (173->6), ass. (0->32)
	t239 = cos(qJ(1));
	t235 = sin(qJ(3));
	t244 = qJD(1) * t235 + qJD(4);
	t260 = t244 * t239;
	t236 = sin(qJ(1));
	t238 = cos(qJ(3));
	t254 = qJD(3) * t239;
	t246 = t238 * t254;
	t259 = t244 * t236 - t246;
	t258 = qJD(1) * t236;
	t257 = qJD(1) * t239;
	t256 = qJD(3) * t235;
	t255 = qJD(3) * t238;
	t253 = qJD(4) * t235;
	t252 = qJD(4) * t238;
	t251 = t238 * t257;
	t237 = cos(qJ(4));
	t250 = t237 * t255;
	t249 = t236 * t255;
	t234 = sin(qJ(4));
	t248 = t234 * t252;
	t247 = t237 * t252;
	t245 = -qJD(1) - t253;
	t243 = t245 * t239;
	t242 = t236 * t256 - t251;
	t241 = t235 * t254 + t238 * t258;
	t240 = t237 * t256 + t248;
	t233 = t237 * t260 + (t245 * t234 + t250) * t236;
	t232 = t245 * t237 * t236 + (-t249 - t260) * t234;
	t231 = t234 * t243 - t259 * t237;
	t230 = t259 * t234 + t237 * t243;
	t1 = [t231, 0, -t240 * t236 + t237 * t251, t232, 0, 0; t233, 0, t241 * t237 + t239 * t248, -t230, 0, 0; 0, 0, t234 * t253 - t250, t234 * t256 - t247, 0, 0; t230, 0, t242 * t234 - t236 * t247, -t233, 0, 0; t232, 0, -t241 * t234 + t239 * t247, t231, 0, 0; 0, 0, t234 * t255 + t237 * t253, t240, 0, 0; t241, 0, t235 * t257 + t249, 0, 0, 0; t242, 0, t235 * t258 - t246, 0, 0, 0; 0, 0, -t256, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:09:40
	% EndTime: 2019-10-10 09:09:41
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (185->30), mult. (233->47), div. (0->0), fcn. (233->6), ass. (0->34)
	t297 = sin(qJ(1));
	t294 = qJD(4) + qJD(5);
	t296 = sin(qJ(3));
	t302 = qJD(1) * t296 + t294;
	t298 = cos(qJ(3));
	t299 = cos(qJ(1));
	t310 = qJD(3) * t299;
	t304 = t298 * t310;
	t318 = t302 * t297 - t304;
	t311 = qJD(3) * t298;
	t305 = t297 * t311;
	t317 = t302 * t299 + t305;
	t316 = t294 * t296;
	t315 = t294 * t298;
	t314 = qJD(1) * t297;
	t313 = qJD(1) * t299;
	t312 = qJD(3) * t296;
	t295 = qJ(4) + qJ(5);
	t292 = sin(t295);
	t309 = t292 * t316;
	t308 = t292 * t315;
	t293 = cos(t295);
	t307 = t293 * t315;
	t306 = t299 * t294 * t293;
	t303 = -qJD(1) - t316;
	t301 = -t297 * t312 + t298 * t313;
	t300 = t296 * t310 + t298 * t314;
	t286 = t293 * t312 + t308;
	t285 = t292 * t312 - t307;
	t284 = -t292 * t314 + t317 * t293 - t297 * t309;
	t283 = t303 * t297 * t293 - t317 * t292;
	t282 = t303 * t299 * t292 - t318 * t293;
	t281 = t318 * t292 - t293 * t313 - t296 * t306;
	t1 = [t282, 0, t301 * t293 - t297 * t308, t283, t283, 0; t284, 0, t300 * t293 + t299 * t308, -t281, -t281, 0; 0, 0, -t293 * t311 + t309, t285, t285, 0; t281, 0, -t301 * t292 - t297 * t307, -t284, -t284, 0; t283, 0, -t300 * t292 + t298 * t306, t282, t282, 0; 0, 0, t292 * t311 + t293 * t316, t286, t286, 0; t300, 0, t296 * t313 + t305, 0, 0, 0; -t301, 0, t296 * t314 - t304, 0, 0, 0; 0, 0, -t312, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:09:41
	% EndTime: 2019-10-10 09:09:41
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (389->32), mult. (293->47), div. (0->0), fcn. (293->6), ass. (0->34)
	t316 = sin(qJ(1));
	t313 = qJD(4) + qJD(5) + qJD(6);
	t315 = sin(qJ(3));
	t321 = qJD(1) * t315 + t313;
	t317 = cos(qJ(3));
	t318 = cos(qJ(1));
	t329 = qJD(3) * t318;
	t323 = t317 * t329;
	t337 = t321 * t316 - t323;
	t330 = qJD(3) * t317;
	t324 = t316 * t330;
	t336 = t321 * t318 + t324;
	t335 = t313 * t315;
	t334 = t313 * t317;
	t333 = qJD(1) * t316;
	t332 = qJD(1) * t318;
	t331 = qJD(3) * t315;
	t314 = qJ(4) + qJ(5) + qJ(6);
	t311 = sin(t314);
	t328 = t311 * t335;
	t327 = t311 * t334;
	t312 = cos(t314);
	t326 = t312 * t334;
	t325 = t318 * t313 * t312;
	t322 = -qJD(1) - t335;
	t320 = -t316 * t331 + t317 * t332;
	t319 = t315 * t329 + t317 * t333;
	t305 = t312 * t331 + t327;
	t304 = t311 * t331 - t326;
	t303 = -t311 * t333 + t336 * t312 - t316 * t328;
	t302 = t322 * t316 * t312 - t336 * t311;
	t301 = t322 * t318 * t311 - t337 * t312;
	t300 = t337 * t311 - t312 * t332 - t315 * t325;
	t1 = [t301, 0, t320 * t312 - t316 * t327, t302, t302, t302; t303, 0, t319 * t312 + t318 * t327, -t300, -t300, -t300; 0, 0, -t312 * t330 + t328, t304, t304, t304; t300, 0, -t320 * t311 - t316 * t326, -t303, -t303, -t303; t302, 0, -t319 * t311 + t317 * t325, t301, t301, t301; 0, 0, t311 * t330 + t312 * t335, t305, t305, t305; t319, 0, t315 * t332 + t324, 0, 0, 0; -t320, 0, t315 * t333 - t323, 0, 0, 0; 0, 0, -t331, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end