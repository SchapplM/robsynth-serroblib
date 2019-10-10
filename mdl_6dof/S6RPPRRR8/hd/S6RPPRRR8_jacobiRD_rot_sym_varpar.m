% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPPRRR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:13
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPPRRR8_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR8_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR8_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRR8_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR8_jacobiRD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:13:28
	% EndTime: 2019-10-10 00:13:28
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:13:28
	% EndTime: 2019-10-10 00:13:28
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
	% StartTime: 2019-10-10 00:13:28
	% EndTime: 2019-10-10 00:13:28
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
	% StartTime: 2019-10-10 00:13:28
	% EndTime: 2019-10-10 00:13:28
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (4->4), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->5)
	t24 = qJD(1) * sin(qJ(1));
	t23 = qJD(1) * cos(qJ(1));
	t20 = cos(pkin(10));
	t19 = sin(pkin(10));
	t1 = [-t19 * t24, 0, 0, 0, 0, 0; t19 * t23, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t20 * t24, 0, 0, 0, 0, 0; t20 * t23, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t23, 0, 0, 0, 0, 0; -t24, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:13:28
	% EndTime: 2019-10-10 00:13:28
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->10), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t44 = sin(qJ(1));
	t49 = qJD(1) * t44;
	t45 = cos(qJ(1));
	t48 = qJD(1) * t45;
	t47 = qJD(4) * t44;
	t46 = qJD(4) * t45;
	t43 = pkin(10) + qJ(4);
	t42 = cos(t43);
	t41 = sin(t43);
	t40 = -t41 * t47 + t42 * t48;
	t39 = t41 * t48 + t42 * t47;
	t38 = t41 * t46 + t42 * t49;
	t37 = -t41 * t49 + t42 * t46;
	t1 = [t37, 0, 0, t40, 0, 0; t39, 0, 0, t38, 0, 0; 0, 0, 0, -qJD(4) * t42, 0, 0; -t38, 0, 0, -t39, 0, 0; t40, 0, 0, t37, 0, 0; 0, 0, 0, qJD(4) * t41, 0, 0; -t48, 0, 0, 0, 0, 0; -t49, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:13:30
	% EndTime: 2019-10-10 00:13:30
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (102->29), mult. (173->57), div. (0->0), fcn. (173->6), ass. (0->31)
	t249 = cos(qJ(1));
	t245 = pkin(10) + qJ(4);
	t243 = sin(t245);
	t251 = qJD(1) * t243 + qJD(5);
	t268 = t251 * t249;
	t247 = sin(qJ(1));
	t244 = cos(t245);
	t261 = qJD(4) * t249;
	t253 = t244 * t261;
	t267 = t251 * t247 - t253;
	t266 = qJD(1) * t247;
	t265 = qJD(1) * t249;
	t246 = sin(qJ(5));
	t264 = qJD(4) * t246;
	t263 = qJD(4) * t247;
	t248 = cos(qJ(5));
	t262 = qJD(4) * t248;
	t260 = qJD(5) * t246;
	t259 = qJD(5) * t248;
	t258 = qJD(5) * t249;
	t257 = t244 * t262;
	t256 = t243 * t263;
	t255 = t244 * t263;
	t254 = t243 * t261;
	t252 = -qJD(5) * t243 - qJD(1);
	t250 = t252 * t249;
	t242 = t248 * t268 + (t252 * t246 + t257) * t247;
	t241 = t252 * t248 * t247 + (-t255 - t268) * t246;
	t240 = t246 * t250 - t267 * t248;
	t239 = t267 * t246 + t248 * t250;
	t1 = [t240, 0, 0, -t248 * t256 + (-t247 * t260 + t248 * t265) * t244, t241, 0; t242, 0, 0, t248 * t254 + (t246 * t258 + t248 * t266) * t244, -t239, 0; 0, 0, 0, t243 * t260 - t257, t243 * t264 - t244 * t259, 0; t239, 0, 0, t246 * t256 + (-t246 * t265 - t247 * t259) * t244, -t242, 0; t241, 0, 0, -t246 * t254 + (-t246 * t266 + t248 * t258) * t244, t240, 0; 0, 0, 0, t243 * t259 + t244 * t264, t243 * t262 + t244 * t260, 0; t244 * t266 + t254, 0, 0, t243 * t265 + t255, 0, 0; -t244 * t265 + t256, 0, 0, t243 * t266 - t253, 0, 0; 0, 0, 0, -qJD(4) * t243, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:13:30
	% EndTime: 2019-10-10 00:13:30
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (254->34), mult. (233->57), div. (0->0), fcn. (233->6), ass. (0->38)
	t314 = sin(qJ(1));
	t311 = pkin(10) + qJ(4);
	t307 = sin(t311);
	t312 = qJD(5) + qJD(6);
	t316 = qJD(1) * t307 + t312;
	t308 = cos(t311);
	t315 = cos(qJ(1));
	t326 = qJD(4) * t315;
	t318 = t308 * t326;
	t337 = t316 * t314 - t318;
	t327 = qJD(4) * t314;
	t320 = t308 * t327;
	t336 = t316 * t315 + t320;
	t313 = qJ(5) + qJ(6);
	t309 = sin(t313);
	t335 = t309 * t312;
	t310 = cos(t313);
	t334 = t310 * t312;
	t333 = t312 * t314;
	t332 = t312 * t315;
	t331 = qJD(1) * t314;
	t330 = qJD(1) * t315;
	t329 = qJD(4) * t309;
	t328 = qJD(4) * t310;
	t325 = t309 * t333;
	t324 = t310 * t332;
	t323 = t309 * t331;
	t322 = t310 * t330;
	t321 = t307 * t327;
	t319 = t307 * t326;
	t317 = -t307 * t312 - qJD(1);
	t301 = t307 * t328 + t308 * t335;
	t300 = t307 * t329 - t308 * t334;
	t299 = -t307 * t325 + t336 * t310 - t323;
	t298 = t317 * t314 * t310 - t336 * t309;
	t297 = t317 * t315 * t309 - t337 * t310;
	t296 = -t307 * t324 + t337 * t309 - t322;
	t1 = [t297, 0, 0, -t310 * t321 + (t322 - t325) * t308, t298, t298; t299, 0, 0, t310 * t319 + (t309 * t332 + t310 * t331) * t308, -t296, -t296; 0, 0, 0, t307 * t335 - t308 * t328, t300, t300; t296, 0, 0, t309 * t321 + (-t309 * t330 - t310 * t333) * t308, -t299, -t299; t298, 0, 0, -t309 * t319 + (-t323 + t324) * t308, t297, t297; 0, 0, 0, t307 * t334 + t308 * t329, t301, t301; t308 * t331 + t319, 0, 0, t307 * t330 + t320, 0, 0; -t308 * t330 + t321, 0, 0, t307 * t331 - t318, 0, 0; 0, 0, 0, -qJD(4) * t307, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end