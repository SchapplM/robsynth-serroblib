% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRPRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:36
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPRPRP5_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP5_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP5_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRP5_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP5_jacobiRD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:35:58
	% EndTime: 2019-10-10 00:35:58
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:35:58
	% EndTime: 2019-10-10 00:35:58
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
	% StartTime: 2019-10-10 00:35:58
	% EndTime: 2019-10-10 00:35:58
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (3->3), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->5)
	t22 = qJD(1) * sin(qJ(1));
	t21 = qJD(1) * cos(qJ(1));
	t18 = cos(pkin(9));
	t17 = sin(pkin(9));
	t1 = [-t18 * t21, 0, 0, 0, 0, 0; -t18 * t22, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t17 * t21, 0, 0, 0, 0, 0; t17 * t22, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t22, 0, 0, 0, 0, 0; t21, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:35:58
	% EndTime: 2019-10-10 00:35:58
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->10), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t43 = sin(qJ(1));
	t48 = qJD(1) * t43;
	t44 = cos(qJ(1));
	t47 = qJD(1) * t44;
	t46 = qJD(3) * t43;
	t45 = qJD(3) * t44;
	t42 = pkin(9) + qJ(3);
	t41 = cos(t42);
	t40 = sin(t42);
	t39 = t40 * t46 - t41 * t47;
	t38 = t40 * t47 + t41 * t46;
	t37 = t40 * t45 + t41 * t48;
	t36 = t40 * t48 - t41 * t45;
	t1 = [t39, 0, t36, 0, 0, 0; -t37, 0, -t38, 0, 0, 0; 0, 0, -qJD(3) * t40, 0, 0, 0; t38, 0, t37, 0, 0, 0; t36, 0, t39, 0, 0, 0; 0, 0, -qJD(3) * t41, 0, 0, 0; -t48, 0, 0, 0, 0, 0; t47, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:35:59
	% EndTime: 2019-10-10 00:35:59
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (44->18), mult. (77->36), div. (0->0), fcn. (77->6), ass. (0->21)
	t194 = sin(pkin(10));
	t196 = sin(qJ(1));
	t210 = t194 * t196;
	t197 = cos(qJ(1));
	t209 = t194 * t197;
	t195 = cos(pkin(10));
	t208 = t195 * t196;
	t207 = t195 * t197;
	t206 = qJD(1) * t196;
	t205 = qJD(1) * t197;
	t193 = pkin(9) + qJ(3);
	t191 = sin(t193);
	t204 = qJD(3) * t191;
	t203 = qJD(3) * t196;
	t202 = qJD(3) * t197;
	t201 = t191 * t203;
	t200 = t191 * t202;
	t192 = cos(t193);
	t199 = t191 * t205 + t192 * t203;
	t198 = t191 * t206 - t192 * t202;
	t1 = [t195 * t201 + (-t192 * t207 - t210) * qJD(1), 0, t198 * t195, 0, 0, 0; -t195 * t200 + (-t192 * t208 + t209) * qJD(1), 0, -t199 * t195, 0, 0, 0; 0, 0, -t195 * t204, 0, 0, 0; -t194 * t201 + (t192 * t209 - t208) * qJD(1), 0, -t198 * t194, 0, 0, 0; t194 * t200 + (t192 * t210 + t207) * qJD(1), 0, t199 * t194, 0, 0, 0; 0, 0, t194 * t204, 0, 0, 0; -t199, 0, -t192 * t206 - t200, 0, 0, 0; -t198, 0, t192 * t205 - t201, 0, 0, 0; 0, 0, qJD(3) * t192, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:35:59
	% EndTime: 2019-10-10 00:35:59
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (161->27), mult. (173->50), div. (0->0), fcn. (173->6), ass. (0->35)
	t266 = cos(qJ(1));
	t264 = pkin(9) + qJ(3);
	t262 = cos(t264);
	t278 = qJD(5) * t262;
	t272 = -qJD(1) + t278;
	t288 = t266 * t272;
	t271 = qJD(1) * t262 - qJD(5);
	t260 = sin(t264);
	t265 = sin(qJ(1));
	t282 = qJD(3) * t265;
	t276 = t260 * t282;
	t287 = t271 * t266 - t276;
	t286 = qJD(1) * t265;
	t285 = qJD(1) * t266;
	t284 = qJD(3) * t260;
	t283 = qJD(3) * t262;
	t281 = qJD(3) * t266;
	t263 = pkin(10) + qJ(5);
	t259 = sin(t263);
	t280 = qJD(5) * t259;
	t279 = qJD(5) * t260;
	t261 = cos(t263);
	t277 = t261 * t279;
	t275 = t262 * t282;
	t274 = t260 * t281;
	t273 = t262 * t281;
	t270 = t272 * t265;
	t269 = t260 * t285 + t275;
	t268 = -t260 * t286 + t273;
	t267 = t271 * t265 + t274;
	t258 = t259 * t270 - t287 * t261;
	t257 = t287 * t259 + t261 * t270;
	t256 = t259 * t288 + t267 * t261;
	t255 = t267 * t259 - t261 * t288;
	t1 = [t258, 0, -t261 * t273 + (t261 * t286 + t266 * t280) * t260, 0, t255, 0; -t256, 0, -t261 * t275 + (-t261 * t285 + t265 * t280) * t260, 0, -t257, 0; 0, 0, -t259 * t278 - t261 * t284, 0, -t259 * t283 - t277, 0; t257, 0, t268 * t259 + t266 * t277, 0, t256, 0; t255, 0, t269 * t259 + t265 * t277, 0, t258, 0; 0, 0, t259 * t284 - t261 * t278, 0, t259 * t279 - t261 * t283, 0; -t269, 0, -t262 * t286 - t274, 0, 0, 0; t268, 0, t262 * t285 - t276, 0, 0, 0; 0, 0, t283, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:36:00
	% EndTime: 2019-10-10 00:36:00
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (161->27), mult. (173->50), div. (0->0), fcn. (173->6), ass. (0->35)
	t327 = sin(qJ(1));
	t326 = pkin(9) + qJ(3);
	t324 = cos(t326);
	t333 = qJD(1) * t324 - qJD(5);
	t322 = sin(t326);
	t328 = cos(qJ(1));
	t343 = qJD(3) * t328;
	t336 = t322 * t343;
	t350 = t333 * t327 + t336;
	t344 = qJD(3) * t327;
	t338 = t322 * t344;
	t349 = t333 * t328 - t338;
	t348 = qJD(1) * t327;
	t347 = qJD(1) * t328;
	t346 = qJD(3) * t322;
	t345 = qJD(3) * t324;
	t325 = pkin(10) + qJ(5);
	t321 = sin(t325);
	t342 = qJD(5) * t321;
	t341 = qJD(5) * t322;
	t340 = qJD(5) * t324;
	t323 = cos(t325);
	t339 = t323 * t341;
	t337 = t324 * t344;
	t335 = t324 * t343;
	t334 = qJD(1) - t340;
	t332 = t334 * t327;
	t331 = t334 * t328;
	t330 = -t322 * t347 - t337;
	t329 = t322 * t348 - t335;
	t320 = t321 * t332 + t349 * t323;
	t319 = -t349 * t321 + t323 * t332;
	t318 = t321 * t331 - t350 * t323;
	t317 = t350 * t321 + t323 * t331;
	t1 = [-t320, 0, -t323 * t335 + (t323 * t348 + t328 * t342) * t322, 0, t317, 0; t318, 0, -t323 * t337 + (-t323 * t347 + t327 * t342) * t322, 0, t319, 0; 0, 0, -t321 * t340 - t323 * t346, 0, -t321 * t345 - t339, 0; t330, 0, -t324 * t348 - t336, 0, 0, 0; -t329, 0, t324 * t347 - t338, 0, 0, 0; 0, 0, t345, 0, 0, 0; t319, 0, t329 * t321 - t328 * t339, 0, t318, 0; -t317, 0, t330 * t321 - t327 * t339, 0, t320, 0; 0, 0, -t321 * t346 + t323 * t340, 0, -t321 * t341 + t323 * t345, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end