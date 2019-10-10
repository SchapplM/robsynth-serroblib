% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRPRR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:58
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPRPRR8_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR8_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR8_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR8_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR8_jacobiRD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:58:22
	% EndTime: 2019-10-10 00:58:22
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:58:22
	% EndTime: 2019-10-10 00:58:22
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
	% StartTime: 2019-10-10 00:58:22
	% EndTime: 2019-10-10 00:58:22
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
	% StartTime: 2019-10-10 00:58:22
	% EndTime: 2019-10-10 00:58:22
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
	% StartTime: 2019-10-10 00:58:22
	% EndTime: 2019-10-10 00:58:22
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->10), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t47 = sin(qJ(1));
	t52 = qJD(1) * t47;
	t48 = cos(qJ(1));
	t51 = qJD(1) * t48;
	t50 = qJD(3) * t47;
	t49 = qJD(3) * t48;
	t46 = qJ(3) + pkin(10);
	t45 = cos(t46);
	t44 = sin(t46);
	t43 = -t44 * t50 + t45 * t51;
	t42 = t44 * t51 + t45 * t50;
	t41 = t44 * t49 + t45 * t52;
	t40 = -t44 * t52 + t45 * t49;
	t1 = [t40, 0, t43, 0, 0, 0; t42, 0, t41, 0, 0, 0; 0, 0, -qJD(3) * t45, 0, 0, 0; -t41, 0, -t42, 0, 0, 0; t43, 0, t40, 0, 0, 0; 0, 0, qJD(3) * t44, 0, 0, 0; -t51, 0, 0, 0, 0, 0; -t52, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:58:23
	% EndTime: 2019-10-10 00:58:23
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (102->29), mult. (173->57), div. (0->0), fcn. (173->6), ass. (0->31)
	t259 = cos(qJ(1));
	t255 = qJ(3) + pkin(10);
	t253 = sin(t255);
	t261 = qJD(1) * t253 + qJD(5);
	t278 = t261 * t259;
	t257 = sin(qJ(1));
	t254 = cos(t255);
	t271 = qJD(3) * t259;
	t263 = t254 * t271;
	t277 = t261 * t257 - t263;
	t276 = qJD(1) * t257;
	t275 = qJD(1) * t259;
	t256 = sin(qJ(5));
	t274 = qJD(3) * t256;
	t273 = qJD(3) * t257;
	t258 = cos(qJ(5));
	t272 = qJD(3) * t258;
	t270 = qJD(5) * t256;
	t269 = qJD(5) * t258;
	t268 = qJD(5) * t259;
	t267 = t254 * t272;
	t266 = t253 * t273;
	t265 = t254 * t273;
	t264 = t253 * t271;
	t262 = -qJD(5) * t253 - qJD(1);
	t260 = t262 * t259;
	t252 = t258 * t278 + (t262 * t256 + t267) * t257;
	t251 = t262 * t258 * t257 + (-t265 - t278) * t256;
	t250 = t256 * t260 - t277 * t258;
	t249 = t277 * t256 + t258 * t260;
	t1 = [t250, 0, -t258 * t266 + (-t257 * t270 + t258 * t275) * t254, 0, t251, 0; t252, 0, t258 * t264 + (t256 * t268 + t258 * t276) * t254, 0, -t249, 0; 0, 0, t253 * t270 - t267, 0, t253 * t274 - t254 * t269, 0; t249, 0, t256 * t266 + (-t256 * t275 - t257 * t269) * t254, 0, -t252, 0; t251, 0, -t256 * t264 + (-t256 * t276 + t258 * t268) * t254, 0, t250, 0; 0, 0, t253 * t269 + t254 * t274, 0, t253 * t272 + t254 * t270, 0; t254 * t276 + t264, 0, t253 * t275 + t265, 0, 0, 0; -t254 * t275 + t266, 0, t253 * t276 - t263, 0, 0, 0; 0, 0, -qJD(3) * t253, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:58:23
	% EndTime: 2019-10-10 00:58:23
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (254->34), mult. (233->57), div. (0->0), fcn. (233->6), ass. (0->38)
	t325 = sin(qJ(1));
	t323 = qJ(3) + pkin(10);
	t318 = sin(t323);
	t322 = qJD(5) + qJD(6);
	t327 = qJD(1) * t318 + t322;
	t319 = cos(t323);
	t326 = cos(qJ(1));
	t337 = qJD(3) * t326;
	t329 = t319 * t337;
	t348 = t327 * t325 - t329;
	t338 = qJD(3) * t325;
	t331 = t319 * t338;
	t347 = t327 * t326 + t331;
	t324 = qJ(5) + qJ(6);
	t320 = sin(t324);
	t346 = t320 * t322;
	t321 = cos(t324);
	t345 = t321 * t322;
	t344 = t322 * t325;
	t343 = t322 * t326;
	t342 = qJD(1) * t325;
	t341 = qJD(1) * t326;
	t340 = qJD(3) * t320;
	t339 = qJD(3) * t321;
	t336 = t320 * t344;
	t335 = t321 * t343;
	t334 = t320 * t342;
	t333 = t321 * t341;
	t332 = t318 * t338;
	t330 = t318 * t337;
	t328 = -t318 * t322 - qJD(1);
	t312 = t318 * t339 + t319 * t346;
	t311 = t318 * t340 - t319 * t345;
	t310 = -t318 * t336 + t347 * t321 - t334;
	t309 = t328 * t325 * t321 - t347 * t320;
	t308 = t328 * t326 * t320 - t348 * t321;
	t307 = -t318 * t335 + t348 * t320 - t333;
	t1 = [t308, 0, -t321 * t332 + (t333 - t336) * t319, 0, t309, t309; t310, 0, t321 * t330 + (t320 * t343 + t321 * t342) * t319, 0, -t307, -t307; 0, 0, t318 * t346 - t319 * t339, 0, t311, t311; t307, 0, t320 * t332 + (-t320 * t341 - t321 * t344) * t319, 0, -t310, -t310; t309, 0, -t320 * t330 + (-t334 + t335) * t319, 0, t308, t308; 0, 0, t318 * t345 + t319 * t340, 0, t312, t312; t319 * t342 + t330, 0, t318 * t341 + t331, 0, 0, 0; -t319 * t341 + t332, 0, t318 * t342 - t329, 0, 0, 0; 0, 0, -qJD(3) * t318, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end