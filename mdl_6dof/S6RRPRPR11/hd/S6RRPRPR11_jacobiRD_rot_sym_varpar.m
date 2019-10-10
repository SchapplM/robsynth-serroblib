% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPRPR11
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:22
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRPRPR11_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR11_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR11_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR11_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR11_jacobiRD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:22:42
	% EndTime: 2019-10-10 10:22:43
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:22:43
	% EndTime: 2019-10-10 10:22:43
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
	% StartTime: 2019-10-10 10:22:43
	% EndTime: 2019-10-10 10:22:43
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
	% StartTime: 2019-10-10 10:22:44
	% EndTime: 2019-10-10 10:22:44
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->7), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->14)
	t161 = sin(qJ(1));
	t168 = qJD(1) * t161;
	t163 = cos(qJ(1));
	t167 = qJD(1) * t163;
	t160 = sin(qJ(2));
	t166 = qJD(2) * t160;
	t162 = cos(qJ(2));
	t165 = qJD(2) * t162;
	t164 = qJD(2) * t163;
	t159 = -t161 * t166 + t162 * t167;
	t158 = t160 * t167 + t161 * t165;
	t157 = t160 * t164 + t162 * t168;
	t156 = -t160 * t168 + t162 * t164;
	t1 = [-t168, 0, 0, 0, 0, 0; t167, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t159, t156, 0, 0, 0, 0; t157, t158, 0, 0, 0, 0; 0, t166, 0, 0, 0, 0; -t158, -t157, 0, 0, 0, 0; t156, t159, 0, 0, 0, 0; 0, t165, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:22:44
	% EndTime: 2019-10-10 10:22:44
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (49->25), mult. (173->44), div. (0->0), fcn. (173->6), ass. (0->33)
	t237 = cos(qJ(4));
	t235 = sin(qJ(2));
	t254 = qJD(4) * t235;
	t246 = qJD(1) + t254;
	t261 = t237 * t246;
	t234 = sin(qJ(4));
	t260 = t246 * t234;
	t236 = sin(qJ(1));
	t259 = qJD(1) * t236;
	t239 = cos(qJ(1));
	t258 = qJD(1) * t239;
	t257 = qJD(2) * t235;
	t238 = cos(qJ(2));
	t256 = qJD(2) * t238;
	t255 = qJD(2) * t239;
	t253 = qJD(4) * t238;
	t252 = t238 * t258;
	t251 = t237 * t256;
	t250 = t236 * t256;
	t249 = t234 * t253;
	t248 = t237 * t253;
	t247 = t238 * t255;
	t245 = -qJD(1) * t235 - qJD(4);
	t244 = t245 * t239;
	t243 = -t236 * t257 + t252;
	t242 = -t235 * t255 - t238 * t259;
	t241 = t237 * t257 + t249;
	t240 = t245 * t236 + t247;
	t233 = t240 * t234 + t239 * t261;
	t232 = t240 * t237 - t239 * t260;
	t231 = -t236 * t261 + (t244 - t250) * t234;
	t230 = t237 * t244 + (-t251 + t260) * t236;
	t1 = [t231, t242 * t234 + t239 * t248, 0, t232, 0, 0; t233, t243 * t234 + t236 * t248, 0, -t230, 0, 0; 0, t234 * t256 + t237 * t254, 0, t241, 0, 0; t230, t242 * t237 - t239 * t249, 0, -t233, 0, 0; t232, -t241 * t236 + t237 * t252, 0, t231, 0, 0; 0, -t234 * t254 + t251, 0, -t234 * t257 + t248, 0, 0; -t243, t235 * t259 - t247, 0, 0, 0, 0; t242, -t235 * t258 - t250, 0, 0, 0, 0; 0, -t257, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:22:44
	% EndTime: 2019-10-10 10:22:44
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (109->25), mult. (173->44), div. (0->0), fcn. (173->6), ass. (0->31)
	t261 = sin(qJ(1));
	t260 = sin(qJ(2));
	t275 = qJD(4) * t260;
	t269 = qJD(1) + t275;
	t282 = t261 * t269;
	t263 = cos(qJ(1));
	t281 = t263 * t269;
	t280 = qJD(1) * t261;
	t279 = qJD(1) * t263;
	t278 = qJD(2) * t260;
	t262 = cos(qJ(2));
	t277 = qJD(2) * t262;
	t276 = qJD(2) * t263;
	t274 = qJD(4) * t262;
	t273 = t261 * t277;
	t259 = qJ(4) + pkin(10);
	t257 = sin(t259);
	t272 = t257 * t274;
	t258 = cos(t259);
	t271 = t258 * t274;
	t270 = t262 * t276;
	t268 = -qJD(1) * t260 - qJD(4);
	t267 = -t261 * t278 + t262 * t279;
	t266 = -t260 * t276 - t262 * t280;
	t265 = t268 * t263 - t273;
	t264 = t268 * t261 + t270;
	t256 = t264 * t257 + t258 * t281;
	t255 = -t257 * t281 + t264 * t258;
	t254 = t265 * t257 - t258 * t282;
	t253 = t257 * t282 + t265 * t258;
	t1 = [t254, t266 * t257 + t263 * t271, 0, t255, 0, 0; t256, t267 * t257 + t261 * t271, 0, -t253, 0, 0; 0, t257 * t277 + t258 * t275, 0, t258 * t278 + t272, 0, 0; t253, t266 * t258 - t263 * t272, 0, -t256, 0, 0; t255, t267 * t258 - t261 * t272, 0, t254, 0, 0; 0, -t257 * t275 + t258 * t277, 0, -t257 * t278 + t271, 0, 0; -t267, t260 * t280 - t270, 0, 0, 0, 0; t266, -t260 * t279 - t273, 0, 0, 0, 0; 0, -t278, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:22:44
	% EndTime: 2019-10-10 10:22:44
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (269->31), mult. (233->47), div. (0->0), fcn. (233->6), ass. (0->36)
	t311 = qJ(4) + pkin(10) + qJ(6);
	t310 = cos(t311);
	t312 = qJD(4) + qJD(6);
	t313 = sin(qJ(2));
	t335 = t312 * t313;
	t320 = qJD(1) + t335;
	t338 = t320 * t310;
	t314 = sin(qJ(1));
	t315 = cos(qJ(2));
	t330 = qJD(2) * t315;
	t322 = t314 * t330;
	t316 = cos(qJ(1));
	t332 = qJD(1) * t316;
	t337 = -t313 * t332 - t322;
	t309 = sin(t311);
	t336 = t309 * t314;
	t334 = t312 * t315;
	t333 = qJD(1) * t314;
	t331 = qJD(2) * t313;
	t329 = qJD(2) * t316;
	t328 = t309 * t334;
	t327 = t310 * t334;
	t326 = t316 * t312 * t310;
	t325 = t313 * t333;
	t323 = t309 * t330;
	t321 = t315 * t329;
	t319 = -qJD(1) * t313 - t312;
	t318 = -t314 * t331 + t315 * t332;
	t317 = -t313 * t329 - t315 * t333;
	t303 = t310 * t331 + t328;
	t302 = -t309 * t331 + t327;
	t301 = -t309 * t325 - t312 * t336 + (t323 + t338) * t316;
	t300 = -t320 * t316 * t309 + (t319 * t314 + t321) * t310;
	t299 = -t314 * t338 + (t319 * t316 - t322) * t309;
	t298 = t337 * t310 + t320 * t336 - t326;
	t1 = [t299, t317 * t309 + t315 * t326, 0, t300, 0, t300; t301, t318 * t309 + t314 * t327, 0, -t298, 0, -t298; 0, t310 * t335 + t323, 0, t303, 0, t303; t298, t317 * t310 - t316 * t328, 0, -t301, 0, -t301; t300, t318 * t310 - t314 * t328, 0, t299, 0, t299; 0, -t309 * t335 + t310 * t330, 0, t302, 0, t302; -t318, -t321 + t325, 0, 0, 0, 0; t317, t337, 0, 0, 0, 0; 0, -t331, 0, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end