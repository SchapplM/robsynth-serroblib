% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRPRP8
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:41
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPRPRP8_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP8_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP8_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRP8_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP8_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:41:05
	% EndTime: 2019-10-10 00:41:05
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:41:05
	% EndTime: 2019-10-10 00:41:05
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
	% StartTime: 2019-10-10 00:41:05
	% EndTime: 2019-10-10 00:41:05
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
	% StartTime: 2019-10-10 00:41:05
	% EndTime: 2019-10-10 00:41:05
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
	% StartTime: 2019-10-10 00:41:05
	% EndTime: 2019-10-10 00:41:05
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->10), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t47 = sin(qJ(1));
	t52 = qJD(1) * t47;
	t48 = cos(qJ(1));
	t51 = qJD(1) * t48;
	t50 = qJD(3) * t47;
	t49 = qJD(3) * t48;
	t46 = qJ(3) + pkin(9);
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
	% StartTime: 2019-10-10 00:41:06
	% EndTime: 2019-10-10 00:41:07
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (102->29), mult. (173->57), div. (0->0), fcn. (173->6), ass. (0->31)
	t259 = cos(qJ(1));
	t255 = qJ(3) + pkin(9);
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
	t267 = t253 * t271;
	t266 = t254 * t272;
	t265 = t253 * t273;
	t264 = t254 * t273;
	t262 = -qJD(5) * t253 - qJD(1);
	t260 = t262 * t259;
	t252 = t258 * t278 + (t262 * t256 + t266) * t257;
	t251 = t262 * t258 * t257 + (-t264 - t278) * t256;
	t250 = t256 * t260 - t277 * t258;
	t249 = t277 * t256 + t258 * t260;
	t1 = [t250, 0, -t258 * t265 + (-t257 * t270 + t258 * t275) * t254, 0, t251, 0; t252, 0, t258 * t267 + (t256 * t268 + t258 * t276) * t254, 0, -t249, 0; 0, 0, t253 * t270 - t266, 0, t253 * t274 - t254 * t269, 0; t249, 0, t256 * t265 + (-t256 * t275 - t257 * t269) * t254, 0, -t252, 0; t251, 0, -t256 * t267 + (-t256 * t276 + t258 * t268) * t254, 0, t250, 0; 0, 0, t253 * t269 + t254 * t274, 0, t253 * t272 + t254 * t270, 0; t254 * t276 + t267, 0, t253 * t275 + t264, 0, 0, 0; -t254 * t275 + t265, 0, t253 * t276 - t263, 0, 0, 0; 0, 0, -qJD(3) * t253, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:41:07
	% EndTime: 2019-10-10 00:41:07
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (102->29), mult. (173->57), div. (0->0), fcn. (173->6), ass. (0->31)
	t302 = sin(qJ(1));
	t300 = qJ(3) + pkin(9);
	t298 = sin(t300);
	t307 = qJD(1) * t298 + qJD(5);
	t299 = cos(t300);
	t304 = cos(qJ(1));
	t317 = qJD(3) * t304;
	t309 = t299 * t317;
	t323 = t307 * t302 - t309;
	t322 = qJD(1) * t302;
	t321 = qJD(1) * t304;
	t301 = sin(qJ(5));
	t320 = qJD(3) * t301;
	t319 = qJD(3) * t302;
	t303 = cos(qJ(5));
	t318 = qJD(3) * t303;
	t316 = qJD(5) * t301;
	t315 = qJD(5) * t303;
	t314 = qJD(5) * t304;
	t313 = t299 * t318;
	t312 = t298 * t319;
	t311 = t299 * t319;
	t310 = t298 * t317;
	t308 = qJD(5) * t298 + qJD(1);
	t306 = t308 * t304;
	t305 = t307 * t304;
	t297 = t303 * t305 + (-t308 * t301 + t313) * t302;
	t296 = t308 * t303 * t302 + (t305 + t311) * t301;
	t295 = t301 * t306 + t323 * t303;
	t294 = -t323 * t301 + t303 * t306;
	t1 = [-t295, 0, -t303 * t312 + (-t302 * t316 + t303 * t321) * t299, 0, -t296, 0; t297, 0, t303 * t310 + (t301 * t314 + t303 * t322) * t299, 0, t294, 0; 0, 0, t298 * t316 - t313, 0, t298 * t320 - t299 * t315, 0; t299 * t322 + t310, 0, t298 * t321 + t311, 0, 0, 0; -t299 * t321 + t312, 0, t298 * t322 - t309, 0, 0, 0; 0, 0, -qJD(3) * t298, 0, 0, 0; t294, 0, -t301 * t312 + (t301 * t321 + t302 * t315) * t299, 0, t297, 0; t296, 0, t301 * t310 + (t301 * t322 - t303 * t314) * t299, 0, t295, 0; 0, 0, -t298 * t315 - t299 * t320, 0, -t298 * t318 - t299 * t316, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end