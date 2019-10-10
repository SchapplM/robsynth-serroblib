% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRRPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:25
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPRRPR2_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR2_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR2_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR2_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR2_jacobiRD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:25:14
	% EndTime: 2019-10-10 01:25:14
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:25:14
	% EndTime: 2019-10-10 01:25:14
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
	% StartTime: 2019-10-10 01:25:14
	% EndTime: 2019-10-10 01:25:14
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
	% StartTime: 2019-10-10 01:25:14
	% EndTime: 2019-10-10 01:25:14
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->10), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t41 = sin(qJ(3));
	t46 = qJD(1) * t41;
	t42 = cos(qJ(3));
	t45 = qJD(1) * t42;
	t44 = qJD(3) * t41;
	t43 = qJD(3) * t42;
	t40 = qJ(1) + pkin(10);
	t39 = cos(t40);
	t38 = sin(t40);
	t37 = t38 * t44 - t39 * t45;
	t36 = t38 * t43 + t39 * t46;
	t35 = t38 * t45 + t39 * t44;
	t34 = t38 * t46 - t39 * t43;
	t1 = [t37, 0, t34, 0, 0, 0; -t35, 0, -t36, 0, 0, 0; 0, 0, -t44, 0, 0, 0; t36, 0, t35, 0, 0, 0; t34, 0, t37, 0, 0, 0; 0, 0, -t43, 0, 0, 0; -qJD(1) * t38, 0, 0, 0, 0, 0; qJD(1) * t39, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:25:15
	% EndTime: 2019-10-10 01:25:15
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (108->24), mult. (173->44), div. (0->0), fcn. (173->6), ass. (0->30)
	t245 = cos(qJ(4));
	t246 = cos(qJ(3));
	t261 = qJD(1) * t246;
	t251 = -qJD(4) + t261;
	t264 = t245 * t251;
	t257 = qJD(4) * t246;
	t252 = -qJD(1) + t257;
	t243 = sin(qJ(4));
	t244 = sin(qJ(3));
	t260 = qJD(3) * t244;
	t254 = t243 * t260;
	t263 = t252 * t245 - t254;
	t262 = qJD(1) * t244;
	t259 = qJD(3) * t246;
	t258 = qJD(4) * t244;
	t256 = t243 * t262;
	t255 = t245 * t262;
	t253 = t245 * t260;
	t250 = t251 * t243;
	t249 = t243 * t258 - t245 * t259;
	t248 = t243 * t259 + t245 * t258;
	t247 = t252 * t243 + t253;
	t242 = qJ(1) + pkin(10);
	t241 = cos(t242);
	t240 = sin(t242);
	t239 = t247 * t240 - t241 * t264;
	t238 = t263 * t240 + t241 * t250;
	t237 = t240 * t264 + t247 * t241;
	t236 = t240 * t250 - t263 * t241;
	t1 = [t239, 0, t240 * t255 + t249 * t241, t236, 0, 0; -t237, 0, t249 * t240 - t241 * t255, -t238, 0, 0; 0, 0, -t243 * t257 - t253, -t248, 0, 0; t238, 0, -t240 * t256 + t248 * t241, t237, 0, 0; t236, 0, t248 * t240 + t241 * t256, t239, 0, 0; 0, 0, -t245 * t257 + t254, t249, 0, 0; -t240 * t259 - t241 * t262, 0, -t240 * t261 - t241 * t260, 0, 0, 0; -t240 * t262 + t241 * t259, 0, -t240 * t260 + t241 * t261, 0, 0, 0; 0, 0, t259, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:25:15
	% EndTime: 2019-10-10 01:25:15
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (168->27), mult. (173->43), div. (0->0), fcn. (173->6), ass. (0->33)
	t265 = qJ(4) + pkin(11);
	t263 = cos(t265);
	t266 = qJ(1) + pkin(10);
	t264 = cos(t266);
	t288 = t263 * t264;
	t267 = sin(qJ(3));
	t287 = qJD(1) * t267;
	t268 = cos(qJ(3));
	t286 = qJD(1) * t268;
	t285 = qJD(3) * t267;
	t284 = qJD(3) * t268;
	t283 = qJD(4) * t267;
	t282 = qJD(4) * t268;
	t281 = t264 * t287;
	t280 = t263 * t285;
	t262 = sin(t266);
	t279 = t262 * t285;
	t278 = t264 * t285;
	t261 = sin(t265);
	t277 = t261 * t283;
	t276 = t263 * t283;
	t275 = -qJD(1) + t282;
	t274 = -qJD(4) + t286;
	t273 = t275 * t261;
	t272 = t262 * t284 + t281;
	t271 = t262 * t287 - t264 * t284;
	t270 = -t263 * t284 + t277;
	t269 = t274 * t262 + t278;
	t260 = -t274 * t288 + (t273 + t280) * t262;
	t259 = t275 * t263 * t262 + (t274 * t264 - t279) * t261;
	t258 = t269 * t263 + t264 * t273;
	t257 = t269 * t261 - t275 * t288;
	t1 = [t260, 0, t271 * t263 + t264 * t277, t257, 0, 0; -t258, 0, t270 * t262 - t263 * t281, -t259, 0, 0; 0, 0, -t261 * t282 - t280, -t261 * t284 - t276, 0, 0; t259, 0, -t271 * t261 + t264 * t276, t258, 0, 0; t257, 0, t272 * t261 + t262 * t276, t260, 0, 0; 0, 0, t261 * t285 - t263 * t282, t270, 0, 0; -t272, 0, -t262 * t286 - t278, 0, 0, 0; -t271, 0, t264 * t286 - t279, 0, 0, 0; 0, 0, t284, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:25:15
	% EndTime: 2019-10-10 01:25:16
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (348->31), mult. (233->48), div. (0->0), fcn. (233->6), ass. (0->34)
	t304 = qJ(4) + pkin(11) + qJ(6);
	t300 = sin(t304);
	t306 = qJ(1) + pkin(10);
	t303 = cos(t306);
	t326 = t300 * t303;
	t301 = cos(t304);
	t302 = sin(t306);
	t325 = t301 * t302;
	t305 = qJD(4) + qJD(6);
	t307 = sin(qJ(3));
	t324 = t305 * t307;
	t308 = cos(qJ(3));
	t323 = t305 * t308;
	t322 = qJD(1) * t307;
	t321 = qJD(1) * t308;
	t320 = qJD(3) * t307;
	t319 = qJD(3) * t308;
	t318 = t300 * t324;
	t317 = t301 * t324;
	t316 = t302 * t320;
	t315 = t303 * t320;
	t314 = -qJD(1) + t323;
	t313 = -t305 + t321;
	t312 = t300 * t314;
	t311 = t302 * t319 + t303 * t322;
	t310 = t302 * t322 - t303 * t319;
	t309 = t313 * t302 + t315;
	t296 = -t301 * t319 + t318;
	t295 = -t300 * t319 - t317;
	t294 = t302 * t312 + (-t313 * t303 + t316) * t301;
	t293 = -t300 * t316 - t305 * t326 - qJD(1) * t325 + (qJD(1) * t326 + t305 * t325) * t308;
	t292 = t309 * t301 + t303 * t312;
	t291 = -t314 * t303 * t301 + t309 * t300;
	t1 = [t294, 0, t310 * t301 + t303 * t318, t291, 0, t291; -t292, 0, -t311 * t301 + t302 * t318, -t293, 0, -t293; 0, 0, -t300 * t323 - t301 * t320, t295, 0, t295; t293, 0, -t310 * t300 + t303 * t317, t292, 0, t292; t291, 0, t311 * t300 + t302 * t317, t294, 0, t294; 0, 0, t300 * t320 - t301 * t323, t296, 0, t296; -t311, 0, -t302 * t321 - t315, 0, 0, 0; -t310, 0, t303 * t321 - t316, 0, 0, 0; 0, 0, t319, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end