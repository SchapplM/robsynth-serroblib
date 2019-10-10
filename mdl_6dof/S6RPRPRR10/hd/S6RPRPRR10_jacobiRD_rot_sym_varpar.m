% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRPRR10
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
% Datum: 2019-10-10 01:02
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPRPRR10_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR10_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR10_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR10_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR10_jacobiRD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:02:11
	% EndTime: 2019-10-10 01:02:11
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:02:11
	% EndTime: 2019-10-10 01:02:11
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0, 0; -t31, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t31, 0, 0, 0, 0, 0; -t30, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:02:11
	% EndTime: 2019-10-10 01:02:11
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (1->1), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t13 = qJD(1) * sin(qJ(1));
	t11 = qJD(1) * cos(qJ(1));
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t11, 0, 0, 0, 0, 0; t13, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t13, 0, 0, 0, 0, 0; t11, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:02:11
	% EndTime: 2019-10-10 01:02:11
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
	% StartTime: 2019-10-10 01:02:12
	% EndTime: 2019-10-10 01:02:12
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (18->17), mult. (77->37), div. (0->0), fcn. (77->6), ass. (0->18)
	t184 = sin(qJ(3));
	t185 = sin(qJ(1));
	t198 = t184 * t185;
	t187 = cos(qJ(1));
	t197 = t184 * t187;
	t196 = qJD(1) * t185;
	t195 = qJD(1) * t187;
	t194 = qJD(3) * t184;
	t186 = cos(qJ(3));
	t193 = qJD(3) * t186;
	t192 = qJD(3) * t187;
	t191 = t185 * t193;
	t190 = t186 * t192;
	t189 = -t185 * t194 + t186 * t195;
	t188 = t184 * t192 + t186 * t196;
	t183 = cos(pkin(10));
	t182 = sin(pkin(10));
	t1 = [t183 * t190 + (-t182 * t187 - t183 * t198) * qJD(1), 0, t189 * t183, 0, 0, 0; t183 * t191 + (-t182 * t185 + t183 * t197) * qJD(1), 0, t188 * t183, 0, 0, 0; 0, 0, -t183 * t193, 0, 0, 0; -t182 * t190 + (t182 * t198 - t183 * t187) * qJD(1), 0, -t189 * t182, 0, 0, 0; -t182 * t191 + (-t182 * t197 - t183 * t185) * qJD(1), 0, -t188 * t182, 0, 0, 0; 0, 0, t182 * t193, 0, 0, 0; t188, 0, t184 * t195 + t191, 0, 0, 0; -t189, 0, t184 * t196 - t190, 0, 0, 0; 0, 0, -t194, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:02:12
	% EndTime: 2019-10-10 01:02:12
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (109->25), mult. (173->44), div. (0->0), fcn. (173->6), ass. (0->31)
	t254 = sin(qJ(1));
	t253 = sin(qJ(3));
	t261 = qJD(1) * t253 + qJD(5);
	t255 = cos(qJ(3));
	t256 = cos(qJ(1));
	t269 = qJD(3) * t256;
	t263 = t255 * t269;
	t275 = t261 * t254 - t263;
	t270 = qJD(3) * t255;
	t266 = t254 * t270;
	t274 = t261 * t256 + t266;
	t273 = qJD(1) * t254;
	t272 = qJD(1) * t256;
	t271 = qJD(3) * t253;
	t268 = qJD(5) * t253;
	t267 = qJD(5) * t255;
	t252 = pkin(10) + qJ(5);
	t250 = sin(t252);
	t265 = t250 * t267;
	t251 = cos(t252);
	t264 = t251 * t267;
	t262 = -qJD(1) - t268;
	t260 = t262 * t254;
	t259 = t262 * t256;
	t258 = -t254 * t271 + t255 * t272;
	t257 = t253 * t269 + t255 * t273;
	t249 = t250 * t260 + t274 * t251;
	t248 = -t274 * t250 + t251 * t260;
	t247 = t250 * t259 - t275 * t251;
	t246 = t275 * t250 + t251 * t259;
	t1 = [t247, 0, t258 * t251 - t254 * t265, 0, t248, 0; t249, 0, t257 * t251 + t256 * t265, 0, -t246, 0; 0, 0, t250 * t268 - t251 * t270, 0, t250 * t271 - t264, 0; t246, 0, -t258 * t250 - t254 * t264, 0, -t249, 0; t248, 0, -t257 * t250 + t256 * t264, 0, t247, 0; 0, 0, t250 * t270 + t251 * t268, 0, t251 * t271 + t265, 0; t257, 0, t253 * t272 + t266, 0, 0, 0; -t258, 0, t253 * t273 - t263, 0, 0, 0; 0, 0, -t271, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:02:12
	% EndTime: 2019-10-10 01:02:13
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (269->30), mult. (233->47), div. (0->0), fcn. (233->6), ass. (0->34)
	t313 = sin(qJ(1));
	t311 = qJD(5) + qJD(6);
	t312 = sin(qJ(3));
	t318 = qJD(1) * t312 + t311;
	t314 = cos(qJ(3));
	t315 = cos(qJ(1));
	t326 = qJD(3) * t315;
	t320 = t314 * t326;
	t334 = t318 * t313 - t320;
	t327 = qJD(3) * t314;
	t321 = t313 * t327;
	t333 = t318 * t315 + t321;
	t332 = t311 * t312;
	t331 = t311 * t314;
	t330 = qJD(1) * t313;
	t329 = qJD(1) * t315;
	t328 = qJD(3) * t312;
	t310 = pkin(10) + qJ(5) + qJ(6);
	t308 = sin(t310);
	t325 = t308 * t332;
	t324 = t308 * t331;
	t309 = cos(t310);
	t323 = t309 * t331;
	t322 = t315 * t311 * t309;
	t319 = -qJD(1) - t332;
	t317 = -t313 * t328 + t314 * t329;
	t316 = t312 * t326 + t314 * t330;
	t302 = t309 * t328 + t324;
	t301 = t308 * t328 - t323;
	t300 = -t308 * t330 + t333 * t309 - t313 * t325;
	t299 = t319 * t313 * t309 - t333 * t308;
	t298 = t319 * t315 * t308 - t334 * t309;
	t297 = t334 * t308 - t309 * t329 - t312 * t322;
	t1 = [t298, 0, t317 * t309 - t313 * t324, 0, t299, t299; t300, 0, t316 * t309 + t315 * t324, 0, -t297, -t297; 0, 0, -t309 * t327 + t325, 0, t301, t301; t297, 0, -t317 * t308 - t313 * t323, 0, -t300, -t300; t299, 0, -t316 * t308 + t314 * t322, 0, t298, t298; 0, 0, t308 * t327 + t309 * t332, 0, t302, t302; t316, 0, t312 * t329 + t321, 0, 0, 0; -t317, 0, t312 * t330 - t320, 0, 0, 0; 0, 0, -t328, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end