% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPPRRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:49
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPPRRP2_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP2_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP2_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRP2_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP2_jacobiRD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:49:30
	% EndTime: 2019-10-09 23:49:30
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:49:30
	% EndTime: 2019-10-09 23:49:30
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0, 0; -t31, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t31, 0, 0, 0, 0, 0; -t30, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:49:30
	% EndTime: 2019-10-09 23:49:30
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (7->4), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->4)
	t35 = qJ(1) + pkin(9);
	t37 = qJD(1) * sin(t35);
	t36 = qJD(1) * cos(t35);
	t1 = [-t36, 0, 0, 0, 0, 0; -t37, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t37, 0, 0, 0, 0, 0; -t36, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:49:30
	% EndTime: 2019-10-09 23:49:30
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (9->4), mult. (10->8), div. (0->0), fcn. (10->4), ass. (0->6)
	t26 = qJD(1) * sin(pkin(10));
	t25 = qJD(1) * cos(pkin(10));
	t22 = qJ(1) + pkin(9);
	t21 = cos(t22);
	t20 = sin(t22);
	t1 = [-t21 * t25, 0, 0, 0, 0, 0; -t20 * t25, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t21 * t26, 0, 0, 0, 0, 0; t20 * t26, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -qJD(1) * t20, 0, 0, 0, 0, 0; qJD(1) * t21, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:49:30
	% EndTime: 2019-10-09 23:49:30
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (47->11), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->16)
	t50 = qJ(1) + pkin(9);
	t46 = sin(t50);
	t55 = qJD(1) * t46;
	t48 = cos(t50);
	t54 = qJD(1) * t48;
	t49 = pkin(10) + qJ(4);
	t45 = sin(t49);
	t53 = qJD(4) * t45;
	t47 = cos(t49);
	t52 = qJD(4) * t47;
	t51 = qJD(4) * t48;
	t44 = t46 * t53 - t47 * t54;
	t43 = t45 * t54 + t46 * t52;
	t42 = t45 * t51 + t47 * t55;
	t41 = t45 * t55 - t47 * t51;
	t1 = [t44, 0, 0, t41, 0, 0; -t42, 0, 0, -t43, 0, 0; 0, 0, 0, -t53, 0, 0; t43, 0, 0, t42, 0, 0; t41, 0, 0, t44, 0, 0; 0, 0, 0, -t52, 0, 0; -t55, 0, 0, 0, 0, 0; t54, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:49:31
	% EndTime: 2019-10-09 23:49:32
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (161->28), mult. (173->56), div. (0->0), fcn. (173->6), ass. (0->34)
	t256 = cos(qJ(5));
	t253 = pkin(10) + qJ(4);
	t251 = cos(t253);
	t274 = qJD(1) * t251;
	t259 = -qJD(5) + t274;
	t277 = t256 * t259;
	t260 = qJD(5) * t251 - qJD(1);
	t249 = sin(t253);
	t255 = sin(qJ(5));
	t268 = qJD(4) * t255;
	t264 = t249 * t268;
	t276 = t260 * t256 - t264;
	t254 = qJ(1) + pkin(9);
	t250 = sin(t254);
	t275 = t249 * t250;
	t252 = cos(t254);
	t273 = qJD(1) * t252;
	t272 = qJD(1) * t255;
	t271 = qJD(1) * t256;
	t270 = qJD(4) * t251;
	t269 = qJD(4) * t252;
	t267 = qJD(4) * t256;
	t266 = qJD(5) * t255;
	t265 = qJD(5) * t256;
	t263 = t249 * t267;
	t262 = t250 * t270;
	t261 = t251 * t269;
	t258 = t259 * t255;
	t257 = t260 * t255 + t263;
	t248 = t257 * t250 - t252 * t277;
	t247 = t276 * t250 + t252 * t258;
	t246 = t250 * t277 + t257 * t252;
	t245 = t250 * t258 - t276 * t252;
	t1 = [t248, 0, 0, -t256 * t261 + (t250 * t271 + t252 * t266) * t249, t245, 0; -t246, 0, 0, -t256 * t262 + (t250 * t266 - t252 * t271) * t249, -t247, 0; 0, 0, 0, -t251 * t266 - t263, -t249 * t265 - t251 * t268, 0; t247, 0, 0, t255 * t261 + (-t250 * t272 + t252 * t265) * t249, t246, 0; t245, 0, 0, t255 * t262 + (t250 * t265 + t252 * t272) * t249, t248, 0; 0, 0, 0, -t251 * t265 + t264, t249 * t266 - t251 * t267, 0; -t249 * t273 - t262, 0, 0, -t249 * t269 - t250 * t274, 0, 0; -qJD(1) * t275 + t261, 0, 0, -qJD(4) * t275 + t251 * t273, 0, 0; 0, 0, 0, t270, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:49:32
	% EndTime: 2019-10-09 23:49:32
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (161->28), mult. (173->56), div. (0->0), fcn. (173->6), ass. (0->34)
	t317 = sin(qJ(5));
	t315 = pkin(10) + qJ(4);
	t313 = cos(t315);
	t336 = qJD(1) * t313;
	t321 = -qJD(5) + t336;
	t339 = t317 * t321;
	t318 = cos(qJ(5));
	t338 = t318 * t321;
	t311 = sin(t315);
	t316 = qJ(1) + pkin(9);
	t312 = sin(t316);
	t337 = t311 * t312;
	t314 = cos(t316);
	t335 = qJD(1) * t314;
	t334 = qJD(1) * t317;
	t333 = qJD(1) * t318;
	t332 = qJD(4) * t313;
	t331 = qJD(4) * t314;
	t330 = qJD(4) * t317;
	t329 = qJD(4) * t318;
	t328 = qJD(5) * t317;
	t327 = qJD(5) * t318;
	t326 = t311 * t330;
	t325 = t311 * t329;
	t324 = t312 * t332;
	t323 = t313 * t331;
	t322 = -qJD(5) * t313 + qJD(1);
	t320 = t318 * t322 + t326;
	t319 = t317 * t322 - t325;
	t310 = t312 * t319 + t314 * t338;
	t309 = t312 * t320 - t314 * t339;
	t308 = -t312 * t338 + t314 * t319;
	t307 = t312 * t339 + t314 * t320;
	t1 = [-t310, 0, 0, -t318 * t323 + (t312 * t333 + t314 * t328) * t311, t307, 0; t308, 0, 0, -t318 * t324 + (t312 * t328 - t314 * t333) * t311, t309, 0; 0, 0, 0, -t313 * t328 - t325, -t311 * t327 - t313 * t330, 0; -t311 * t335 - t324, 0, 0, -t311 * t331 - t312 * t336, 0, 0; -qJD(1) * t337 + t323, 0, 0, -qJD(4) * t337 + t313 * t335, 0, 0; 0, 0, 0, t332, 0, 0; t309, 0, 0, -t317 * t323 + (t312 * t334 - t314 * t327) * t311, t308, 0; -t307, 0, 0, -t317 * t324 + (-t312 * t327 - t314 * t334) * t311, t310, 0; 0, 0, 0, t313 * t327 - t326, -t311 * t328 + t313 * t329, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end