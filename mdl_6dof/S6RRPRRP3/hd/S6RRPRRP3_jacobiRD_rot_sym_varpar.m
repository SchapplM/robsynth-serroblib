% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPRRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:33
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRPRRP3_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP3_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP3_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRP3_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP3_jacobiRD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:33:36
	% EndTime: 2019-10-10 10:33:36
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:33:36
	% EndTime: 2019-10-10 10:33:36
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
	% StartTime: 2019-10-10 10:33:36
	% EndTime: 2019-10-10 10:33:36
	% DurationCPUTime: 0.04s
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
	t1 = [t32, t29, 0, 0, 0, 0; -t30, -t31, 0, 0, 0, 0; 0, -t39, 0, 0, 0, 0; t31, t30, 0, 0, 0, 0; t29, t32, 0, 0, 0, 0; 0, -t38, 0, 0, 0, 0; -t41, 0, 0, 0, 0, 0; t40, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:33:36
	% EndTime: 2019-10-10 10:33:36
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->10), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t50 = sin(qJ(1));
	t55 = qJD(1) * t50;
	t51 = cos(qJ(1));
	t54 = qJD(1) * t51;
	t53 = qJD(2) * t50;
	t52 = qJD(2) * t51;
	t49 = qJ(2) + pkin(10);
	t48 = cos(t49);
	t47 = sin(t49);
	t46 = t47 * t53 - t48 * t54;
	t45 = t47 * t54 + t48 * t53;
	t44 = t47 * t52 + t48 * t55;
	t43 = t47 * t55 - t48 * t52;
	t1 = [t46, t43, 0, 0, 0, 0; -t44, -t45, 0, 0, 0, 0; 0, -qJD(2) * t47, 0, 0, 0, 0; t45, t44, 0, 0, 0, 0; t43, t46, 0, 0, 0, 0; 0, -qJD(2) * t48, 0, 0, 0, 0; -t55, 0, 0, 0, 0, 0; t54, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:33:38
	% EndTime: 2019-10-10 10:33:38
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (101->28), mult. (173->57), div. (0->0), fcn. (173->6), ass. (0->31)
	t256 = cos(qJ(4));
	t257 = cos(qJ(1));
	t276 = t256 * t257;
	t255 = sin(qJ(1));
	t275 = qJD(1) * t255;
	t274 = qJD(1) * t257;
	t254 = sin(qJ(4));
	t273 = qJD(2) * t254;
	t272 = qJD(2) * t255;
	t271 = qJD(2) * t256;
	t270 = qJD(2) * t257;
	t269 = qJD(4) * t254;
	t268 = qJD(4) * t256;
	t267 = qJD(4) * t257;
	t253 = qJ(2) + pkin(10);
	t251 = sin(t253);
	t266 = t251 * t271;
	t265 = t251 * t272;
	t252 = cos(t253);
	t264 = t252 * t272;
	t263 = t251 * t270;
	t262 = t252 * t270;
	t261 = qJD(4) * t252 - qJD(1);
	t260 = qJD(1) * t252 - qJD(4);
	t259 = t261 * t254;
	t258 = t260 * t255 + t263;
	t250 = -t260 * t276 + (t259 + t266) * t255;
	t249 = t261 * t256 * t255 + (t260 * t257 - t265) * t254;
	t248 = t258 * t256 + t257 * t259;
	t247 = t258 * t254 - t261 * t276;
	t1 = [t250, -t256 * t262 + (t254 * t267 + t256 * t275) * t251, 0, t247, 0, 0; -t248, -t256 * t264 + (t255 * t269 - t256 * t274) * t251, 0, -t249, 0, 0; 0, -t252 * t269 - t266, 0, -t251 * t268 - t252 * t273, 0, 0; t249, t254 * t262 + (-t254 * t275 + t256 * t267) * t251, 0, t248, 0, 0; t247, t254 * t264 + (t254 * t274 + t255 * t268) * t251, 0, t250, 0, 0; 0, t251 * t273 - t252 * t268, 0, t251 * t269 - t252 * t271, 0, 0; -t251 * t274 - t264, -t252 * t275 - t263, 0, 0, 0, 0; -t251 * t275 + t262, t252 * t274 - t265, 0, 0, 0, 0; 0, qJD(2) * t252, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:33:38
	% EndTime: 2019-10-10 10:33:38
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (253->31), mult. (233->57), div. (0->0), fcn. (233->6), ass. (0->36)
	t318 = qJ(4) + qJ(5);
	t314 = sin(t318);
	t316 = qJD(4) + qJD(5);
	t340 = t314 * t316;
	t315 = cos(t318);
	t339 = t315 * t316;
	t319 = sin(qJ(1));
	t338 = t316 * t319;
	t320 = cos(qJ(1));
	t337 = t316 * t320;
	t336 = qJD(1) * t319;
	t335 = qJD(1) * t320;
	t334 = qJD(2) * t314;
	t333 = qJD(2) * t315;
	t332 = qJD(2) * t319;
	t331 = qJD(2) * t320;
	t317 = qJ(2) + pkin(10);
	t312 = sin(t317);
	t330 = t312 * t332;
	t313 = cos(t317);
	t329 = t313 * t332;
	t328 = t312 * t331;
	t327 = t313 * t331;
	t326 = t313 * t316 - qJD(1);
	t325 = qJD(1) * t313 - t316;
	t324 = t314 * t326;
	t323 = t314 * t337 + t315 * t336;
	t322 = t314 * t335 + t315 * t338;
	t321 = t325 * t319 + t328;
	t308 = t312 * t340 - t313 * t333;
	t307 = -t312 * t339 - t313 * t334;
	t306 = t319 * t324 + (-t325 * t320 + t330) * t315;
	t305 = t322 * t313 - t314 * t330 - t323;
	t304 = t321 * t315 + t320 * t324;
	t303 = -t326 * t320 * t315 + t321 * t314;
	t1 = [t306, t323 * t312 - t315 * t327, 0, t303, t303, 0; -t304, -t315 * t329 + (t314 * t338 - t315 * t335) * t312, 0, -t305, -t305, 0; 0, -t312 * t333 - t313 * t340, 0, t307, t307, 0; t305, t314 * t327 + (-t314 * t336 + t315 * t337) * t312, 0, t304, t304, 0; t303, t322 * t312 + t314 * t329, 0, t306, t306, 0; 0, t312 * t334 - t313 * t339, 0, t308, t308, 0; -t312 * t335 - t329, -t313 * t336 - t328, 0, 0, 0, 0; -t312 * t336 + t327, t313 * t335 - t330, 0, 0, 0, 0; 0, qJD(2) * t313, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:33:38
	% EndTime: 2019-10-10 10:33:38
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (253->31), mult. (233->57), div. (0->0), fcn. (233->6), ass. (0->36)
	t327 = qJ(4) + qJ(5);
	t323 = sin(t327);
	t325 = qJD(4) + qJD(5);
	t349 = t323 * t325;
	t324 = cos(t327);
	t348 = t324 * t325;
	t328 = sin(qJ(1));
	t347 = t325 * t328;
	t329 = cos(qJ(1));
	t346 = t325 * t329;
	t345 = qJD(1) * t328;
	t344 = qJD(1) * t329;
	t343 = qJD(2) * t323;
	t342 = qJD(2) * t324;
	t341 = qJD(2) * t328;
	t340 = qJD(2) * t329;
	t326 = qJ(2) + pkin(10);
	t321 = sin(t326);
	t339 = t321 * t341;
	t322 = cos(t326);
	t338 = t322 * t341;
	t337 = t321 * t340;
	t336 = t322 * t340;
	t335 = t322 * t325 - qJD(1);
	t334 = qJD(1) * t322 - t325;
	t333 = t323 * t335;
	t332 = t323 * t346 + t324 * t345;
	t331 = t323 * t344 + t324 * t347;
	t330 = t334 * t328 + t337;
	t317 = t321 * t349 - t322 * t342;
	t316 = -t321 * t348 - t322 * t343;
	t315 = t328 * t333 + (-t334 * t329 + t339) * t324;
	t314 = t331 * t322 - t323 * t339 - t332;
	t313 = t330 * t324 + t329 * t333;
	t312 = -t335 * t329 * t324 + t330 * t323;
	t1 = [t315, t332 * t321 - t324 * t336, 0, t312, t312, 0; -t313, -t324 * t338 + (t323 * t347 - t324 * t344) * t321, 0, -t314, -t314, 0; 0, -t321 * t342 - t322 * t349, 0, t316, t316, 0; t314, t323 * t336 + (-t323 * t345 + t324 * t346) * t321, 0, t313, t313, 0; t312, t331 * t321 + t323 * t338, 0, t315, t315, 0; 0, t321 * t343 - t322 * t348, 0, t317, t317, 0; -t321 * t344 - t338, -t322 * t345 - t337, 0, 0, 0, 0; -t321 * t345 + t336, t322 * t344 - t339, 0, 0, 0, 0; 0, qJD(2) * t322, 0, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end