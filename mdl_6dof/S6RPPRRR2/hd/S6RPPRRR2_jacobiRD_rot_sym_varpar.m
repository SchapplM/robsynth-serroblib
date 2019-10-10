% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPPRRR2
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:03
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPPRRR2_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR2_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR2_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRR2_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR2_jacobiRD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:03:06
	% EndTime: 2019-10-10 00:03:06
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:03:06
	% EndTime: 2019-10-10 00:03:06
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
	% StartTime: 2019-10-10 00:03:06
	% EndTime: 2019-10-10 00:03:06
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (7->4), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->4)
	t35 = qJ(1) + pkin(10);
	t37 = qJD(1) * sin(t35);
	t36 = qJD(1) * cos(t35);
	t1 = [-t36, 0, 0, 0, 0, 0; -t37, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t37, 0, 0, 0, 0, 0; -t36, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:03:06
	% EndTime: 2019-10-10 00:03:06
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (9->4), mult. (10->8), div. (0->0), fcn. (10->4), ass. (0->6)
	t26 = qJD(1) * sin(pkin(11));
	t25 = qJD(1) * cos(pkin(11));
	t22 = qJ(1) + pkin(10);
	t21 = cos(t22);
	t20 = sin(t22);
	t1 = [-t21 * t25, 0, 0, 0, 0, 0; -t20 * t25, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t21 * t26, 0, 0, 0, 0, 0; t20 * t26, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -qJD(1) * t20, 0, 0, 0, 0, 0; qJD(1) * t21, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:03:06
	% EndTime: 2019-10-10 00:03:06
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (47->11), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->16)
	t50 = qJ(1) + pkin(10);
	t46 = sin(t50);
	t55 = qJD(1) * t46;
	t48 = cos(t50);
	t54 = qJD(1) * t48;
	t49 = pkin(11) + qJ(4);
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
	% StartTime: 2019-10-10 00:03:08
	% EndTime: 2019-10-10 00:03:08
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (161->28), mult. (173->56), div. (0->0), fcn. (173->6), ass. (0->34)
	t256 = cos(qJ(5));
	t253 = pkin(11) + qJ(4);
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
	t254 = qJ(1) + pkin(10);
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
	% StartTime: 2019-10-10 00:03:08
	% EndTime: 2019-10-10 00:03:08
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (333->32), mult. (233->56), div. (0->0), fcn. (233->6), ass. (0->39)
	t317 = pkin(11) + qJ(4);
	t311 = sin(t317);
	t319 = qJ(1) + pkin(10);
	t312 = sin(t319);
	t342 = t311 * t312;
	t320 = qJ(5) + qJ(6);
	t315 = sin(t320);
	t318 = qJD(5) + qJD(6);
	t341 = t315 * t318;
	t316 = cos(t320);
	t340 = t316 * t318;
	t313 = cos(t317);
	t339 = qJD(1) * t313;
	t314 = cos(t319);
	t338 = qJD(1) * t314;
	t337 = qJD(1) * t315;
	t336 = qJD(1) * t316;
	t335 = qJD(4) * t313;
	t334 = qJD(4) * t314;
	t333 = qJD(4) * t315;
	t332 = qJD(4) * t316;
	t331 = t311 * t333;
	t330 = t311 * t332;
	t329 = qJD(4) * t342;
	t328 = t312 * t335;
	t327 = t313 * t334;
	t326 = t313 * t318 - qJD(1);
	t325 = -t318 + t339;
	t324 = t312 * t325;
	t323 = t312 * t340 + t314 * t337;
	t322 = t312 * t336 + t314 * t341;
	t321 = t326 * t315 + t330;
	t307 = t311 * t341 - t313 * t332;
	t306 = -t311 * t340 - t313 * t333;
	t305 = -t325 * t316 * t314 + t321 * t312;
	t304 = t323 * t313 - t315 * t329 - t322;
	t303 = t321 * t314 + t316 * t324;
	t302 = t315 * t324 + (-t326 * t316 + t331) * t314;
	t1 = [t305, 0, 0, t322 * t311 - t316 * t327, t302, t302; -t303, 0, 0, -t316 * t328 + (t312 * t341 - t314 * t336) * t311, -t304, -t304; 0, 0, 0, -t313 * t341 - t330, t306, t306; t304, 0, 0, t315 * t327 + (-t312 * t337 + t314 * t340) * t311, t303, t303; t302, 0, 0, t323 * t311 + t315 * t328, t305, t305; 0, 0, 0, -t313 * t340 + t331, t307, t307; -t311 * t338 - t328, 0, 0, -t311 * t334 - t312 * t339, 0, 0; -qJD(1) * t342 + t327, 0, 0, t313 * t338 - t329, 0, 0; 0, 0, 0, t335, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end