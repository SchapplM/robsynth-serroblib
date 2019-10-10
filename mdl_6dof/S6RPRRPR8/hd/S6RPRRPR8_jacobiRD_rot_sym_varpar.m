% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:35
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPRRPR8_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR8_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR8_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR8_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR8_jacobiRD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:35:37
	% EndTime: 2019-10-10 01:35:37
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:35:37
	% EndTime: 2019-10-10 01:35:37
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0, 0; -t31, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t31, 0, 0, 0, 0, 0; -t30, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:35:37
	% EndTime: 2019-10-10 01:35:38
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (1->1), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t13 = qJD(1) * sin(qJ(1));
	t11 = qJD(1) * cos(qJ(1));
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t11, 0, 0, 0, 0, 0; t13, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t13, 0, 0, 0, 0, 0; t11, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:35:38
	% EndTime: 2019-10-10 01:35:38
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
	% StartTime: 2019-10-10 01:35:39
	% EndTime: 2019-10-10 01:35:39
	% DurationCPUTime: 0.13s
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
	% StartTime: 2019-10-10 01:35:39
	% EndTime: 2019-10-10 01:35:39
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (109->25), mult. (173->44), div. (0->0), fcn. (173->6), ass. (0->31)
	t265 = sin(qJ(1));
	t264 = sin(qJ(3));
	t272 = qJD(1) * t264 + qJD(4);
	t266 = cos(qJ(3));
	t267 = cos(qJ(1));
	t280 = qJD(3) * t267;
	t274 = t266 * t280;
	t286 = t272 * t265 - t274;
	t281 = qJD(3) * t266;
	t277 = t265 * t281;
	t285 = t272 * t267 + t277;
	t284 = qJD(1) * t265;
	t283 = qJD(1) * t267;
	t282 = qJD(3) * t264;
	t279 = qJD(4) * t264;
	t278 = qJD(4) * t266;
	t263 = qJ(4) + pkin(10);
	t261 = sin(t263);
	t276 = t261 * t278;
	t262 = cos(t263);
	t275 = t262 * t278;
	t273 = -qJD(1) - t279;
	t271 = t273 * t265;
	t270 = t273 * t267;
	t269 = -t265 * t282 + t266 * t283;
	t268 = t264 * t280 + t266 * t284;
	t260 = t261 * t271 + t285 * t262;
	t259 = -t285 * t261 + t262 * t271;
	t258 = t261 * t270 - t286 * t262;
	t257 = t286 * t261 + t262 * t270;
	t1 = [t258, 0, t269 * t262 - t265 * t276, t259, 0, 0; t260, 0, t268 * t262 + t267 * t276, -t257, 0, 0; 0, 0, t261 * t279 - t262 * t281, t261 * t282 - t275, 0, 0; t257, 0, -t269 * t261 - t265 * t275, -t260, 0, 0; t259, 0, -t268 * t261 + t267 * t275, t258, 0, 0; 0, 0, t261 * t281 + t262 * t279, t262 * t282 + t276, 0, 0; t268, 0, t264 * t283 + t277, 0, 0, 0; -t269, 0, t264 * t284 - t274, 0, 0, 0; 0, 0, -t282, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:35:39
	% EndTime: 2019-10-10 01:35:39
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (269->30), mult. (233->47), div. (0->0), fcn. (233->6), ass. (0->34)
	t316 = sin(qJ(1));
	t314 = qJD(4) + qJD(6);
	t315 = sin(qJ(3));
	t321 = qJD(1) * t315 + t314;
	t317 = cos(qJ(3));
	t318 = cos(qJ(1));
	t329 = qJD(3) * t318;
	t323 = t317 * t329;
	t337 = t316 * t321 - t323;
	t330 = qJD(3) * t317;
	t324 = t316 * t330;
	t336 = t318 * t321 + t324;
	t335 = t314 * t315;
	t334 = t314 * t317;
	t333 = qJD(1) * t316;
	t332 = qJD(1) * t318;
	t331 = qJD(3) * t315;
	t313 = qJ(4) + pkin(10) + qJ(6);
	t311 = sin(t313);
	t328 = t311 * t335;
	t327 = t311 * t334;
	t312 = cos(t313);
	t326 = t312 * t334;
	t325 = t318 * t314 * t312;
	t322 = -qJD(1) - t335;
	t320 = -t316 * t331 + t317 * t332;
	t319 = t315 * t329 + t317 * t333;
	t305 = t312 * t331 + t327;
	t304 = t311 * t331 - t326;
	t303 = -t311 * t333 + t336 * t312 - t316 * t328;
	t302 = t322 * t316 * t312 - t336 * t311;
	t301 = t322 * t318 * t311 - t337 * t312;
	t300 = t337 * t311 - t312 * t332 - t315 * t325;
	t1 = [t301, 0, t312 * t320 - t316 * t327, t302, 0, t302; t303, 0, t312 * t319 + t318 * t327, -t300, 0, -t300; 0, 0, -t312 * t330 + t328, t304, 0, t304; t300, 0, -t311 * t320 - t316 * t326, -t303, 0, -t303; t302, 0, -t311 * t319 + t317 * t325, t301, 0, t301; 0, 0, t311 * t330 + t312 * t335, t305, 0, t305; t319, 0, t315 * t332 + t324, 0, 0, 0; -t320, 0, t315 * t333 - t323, 0, 0, 0; 0, 0, -t331, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end