% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 19:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RRRPP6_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP6_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP6_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPP6_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP6_jacobiRD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:48:17
	% EndTime: 2019-12-29 19:48:17
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:48:22
	% EndTime: 2019-12-29 19:48:22
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0; -t31, 0, 0, 0, 0; 0, 0, 0, 0, 0; t31, 0, 0, 0, 0; -t30, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:48:17
	% EndTime: 2019-12-29 19:48:17
	% DurationCPUTime: 0.05s
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
	t1 = [t32, t29, 0, 0, 0; -t30, -t31, 0, 0, 0; 0, -t39, 0, 0, 0; t31, t30, 0, 0, 0; t29, t32, 0, 0, 0; 0, -t38, 0, 0, 0; -t41, 0, 0, 0, 0; t40, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:48:19
	% EndTime: 2019-12-29 19:48:19
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (48->26), mult. (173->50), div. (0->0), fcn. (173->6), ass. (0->33)
	t232 = cos(qJ(3));
	t234 = cos(qJ(1));
	t256 = t232 * t234;
	t231 = sin(qJ(1));
	t255 = qJD(1) * t231;
	t233 = cos(qJ(2));
	t254 = qJD(1) * t233;
	t253 = qJD(1) * t234;
	t230 = sin(qJ(2));
	t252 = qJD(2) * t230;
	t251 = qJD(2) * t233;
	t250 = qJD(2) * t234;
	t229 = sin(qJ(3));
	t249 = qJD(3) * t229;
	t248 = qJD(3) * t230;
	t247 = qJD(3) * t233;
	t246 = t232 * t252;
	t245 = t232 * t248;
	t244 = t231 * t252;
	t243 = t231 * t251;
	t242 = t230 * t250;
	t241 = t233 * t250;
	t240 = -qJD(1) + t247;
	t239 = -qJD(3) + t254;
	t238 = t240 * t229;
	t237 = t230 * t253 + t243;
	t236 = -t230 * t255 + t241;
	t235 = t239 * t231 + t242;
	t228 = -t239 * t256 + (t238 + t246) * t231;
	t227 = t240 * t232 * t231 + (t239 * t234 - t244) * t229;
	t226 = t235 * t232 + t234 * t238;
	t225 = t235 * t229 - t240 * t256;
	t1 = [t228, -t232 * t241 + (t232 * t255 + t234 * t249) * t230, t225, 0, 0; -t226, -t232 * t243 + (t231 * t249 - t232 * t253) * t230, -t227, 0, 0; 0, -t229 * t247 - t246, -t229 * t251 - t245, 0, 0; t227, t236 * t229 + t234 * t245, t226, 0, 0; t225, t237 * t229 + t231 * t245, t228, 0, 0; 0, t229 * t252 - t232 * t247, t229 * t248 - t232 * t251, 0, 0; -t237, -t231 * t254 - t242, 0, 0, 0; t236, t233 * t253 - t244, 0, 0, 0; 0, t251, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:48:19
	% EndTime: 2019-12-29 19:48:19
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (108->25), mult. (173->44), div. (0->0), fcn. (173->6), ass. (0->31)
	t260 = cos(qJ(1));
	t259 = cos(qJ(2));
	t271 = qJD(3) * t259;
	t266 = -qJD(1) + t271;
	t279 = t260 * t266;
	t277 = qJD(1) * t259;
	t265 = -qJD(3) + t277;
	t258 = sin(qJ(1));
	t257 = sin(qJ(2));
	t275 = qJD(2) * t257;
	t268 = t258 * t275;
	t278 = t265 * t260 - t268;
	t276 = qJD(1) * t260;
	t274 = qJD(2) * t259;
	t273 = qJD(2) * t260;
	t272 = qJD(3) * t257;
	t256 = qJ(3) + pkin(8);
	t254 = sin(t256);
	t270 = t254 * t272;
	t255 = cos(t256);
	t269 = t255 * t272;
	t267 = t257 * t273;
	t264 = t266 * t258;
	t263 = t257 * t276 + t258 * t274;
	t262 = qJD(1) * t258 * t257 - t259 * t273;
	t261 = t265 * t258 + t267;
	t253 = t254 * t264 - t278 * t255;
	t252 = t278 * t254 + t255 * t264;
	t251 = t254 * t279 + t261 * t255;
	t250 = t261 * t254 - t255 * t279;
	t1 = [t253, t262 * t255 + t260 * t270, t250, 0, 0; -t251, -t263 * t255 + t258 * t270, -t252, 0, 0; 0, -t254 * t271 - t255 * t275, -t254 * t274 - t269, 0, 0; t252, -t262 * t254 + t260 * t269, t251, 0, 0; t250, t263 * t254 + t258 * t269, t253, 0, 0; 0, t254 * t275 - t255 * t271, -t255 * t274 + t270, 0, 0; -t263, -t258 * t277 - t267, 0, 0, 0; -t262, t259 * t276 - t268, 0, 0, 0; 0, t274, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:48:20
	% EndTime: 2019-12-29 19:48:20
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (108->24), mult. (173->44), div. (0->0), fcn. (173->6), ass. (0->31)
	t320 = sin(qJ(1));
	t321 = cos(qJ(2));
	t339 = qJD(1) * t321;
	t327 = -qJD(3) + t339;
	t319 = sin(qJ(2));
	t322 = cos(qJ(1));
	t335 = qJD(2) * t322;
	t329 = t319 * t335;
	t341 = t327 * t320 + t329;
	t337 = qJD(2) * t319;
	t330 = t320 * t337;
	t340 = t327 * t322 - t330;
	t338 = qJD(1) * t322;
	t336 = qJD(2) * t321;
	t334 = qJD(3) * t319;
	t333 = qJD(3) * t321;
	t318 = qJ(3) + pkin(8);
	t316 = sin(t318);
	t332 = t316 * t334;
	t317 = cos(t318);
	t331 = t317 * t334;
	t328 = qJD(1) - t333;
	t326 = t328 * t320;
	t325 = t328 * t322;
	t324 = -t319 * t338 - t320 * t336;
	t323 = qJD(1) * t320 * t319 - t321 * t335;
	t315 = t316 * t326 + t340 * t317;
	t314 = -t340 * t316 + t317 * t326;
	t313 = t316 * t325 - t341 * t317;
	t312 = t341 * t316 + t317 * t325;
	t1 = [-t315, t323 * t317 + t322 * t332, t312, 0, 0; t313, t324 * t317 + t320 * t332, t314, 0, 0; 0, -t316 * t333 - t317 * t337, -t316 * t336 - t331, 0, 0; t324, -t320 * t339 - t329, 0, 0, 0; -t323, t321 * t338 - t330, 0, 0, 0; 0, t336, 0, 0, 0; t314, t323 * t316 - t322 * t331, t313, 0, 0; -t312, t324 * t316 - t320 * t331, t315, 0, 0; 0, -t316 * t337 + t317 * t333, t317 * t336 - t332, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end