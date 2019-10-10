% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPPRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3,theta5]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:34
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPPRPR1_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR1_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR1_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRPR1_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRPR1_jacobiRD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:34:17
	% EndTime: 2019-10-09 23:34:17
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:34:17
	% EndTime: 2019-10-09 23:34:17
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
	% StartTime: 2019-10-09 23:34:17
	% EndTime: 2019-10-09 23:34:17
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
	% StartTime: 2019-10-09 23:34:17
	% EndTime: 2019-10-09 23:34:17
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
	% StartTime: 2019-10-09 23:34:17
	% EndTime: 2019-10-09 23:34:18
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
	% StartTime: 2019-10-09 23:34:19
	% EndTime: 2019-10-09 23:34:19
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (72->19), mult. (77->35), div. (0->0), fcn. (77->6), ass. (0->22)
	t206 = qJ(1) + pkin(9);
	t202 = sin(t206);
	t207 = sin(pkin(11));
	t221 = t202 * t207;
	t208 = cos(pkin(11));
	t220 = t202 * t208;
	t204 = cos(t206);
	t219 = t204 * t207;
	t218 = t204 * t208;
	t217 = qJD(1) * t202;
	t216 = qJD(1) * t204;
	t205 = pkin(10) + qJ(4);
	t201 = sin(t205);
	t215 = qJD(4) * t201;
	t203 = cos(t205);
	t214 = qJD(4) * t203;
	t213 = qJD(4) * t204;
	t212 = t202 * t215;
	t211 = t201 * t213;
	t210 = t201 * t216 + t202 * t214;
	t209 = t201 * t217 - t203 * t213;
	t1 = [t208 * t212 + (-t203 * t218 - t221) * qJD(1), 0, 0, t209 * t208, 0, 0; -t208 * t211 + (-t203 * t220 + t219) * qJD(1), 0, 0, -t210 * t208, 0, 0; 0, 0, 0, -t208 * t215, 0, 0; -t207 * t212 + (t203 * t219 - t220) * qJD(1), 0, 0, -t209 * t207, 0, 0; t207 * t211 + (t203 * t221 + t218) * qJD(1), 0, 0, t210 * t207, 0, 0; 0, 0, 0, t207 * t215, 0, 0; -t210, 0, 0, -t203 * t217 - t211, 0, 0; -t209, 0, 0, t203 * t216 - t212, 0, 0; 0, 0, 0, t214, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:34:19
	% EndTime: 2019-10-09 23:34:19
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (221->29), mult. (173->50), div. (0->0), fcn. (173->6), ass. (0->36)
	t270 = pkin(11) + qJ(6);
	t267 = cos(t270);
	t272 = qJ(1) + pkin(9);
	t269 = cos(t272);
	t294 = t267 * t269;
	t266 = sin(t272);
	t293 = qJD(1) * t266;
	t271 = pkin(10) + qJ(4);
	t268 = cos(t271);
	t292 = qJD(1) * t268;
	t291 = qJD(1) * t269;
	t265 = sin(t271);
	t290 = qJD(4) * t265;
	t289 = qJD(4) * t268;
	t288 = qJD(4) * t269;
	t264 = sin(t270);
	t287 = qJD(6) * t264;
	t286 = qJD(6) * t265;
	t285 = qJD(6) * t268;
	t284 = t267 * t290;
	t283 = t267 * t286;
	t282 = t266 * t290;
	t281 = t266 * t289;
	t280 = t265 * t288;
	t279 = t268 * t288;
	t278 = -qJD(1) + t285;
	t277 = -qJD(6) + t292;
	t276 = t278 * t264;
	t275 = t265 * t291 + t281;
	t274 = -t265 * t293 + t279;
	t273 = t277 * t266 + t280;
	t263 = -t277 * t294 + (t276 + t284) * t266;
	t262 = t278 * t267 * t266 + (t277 * t269 - t282) * t264;
	t261 = t273 * t267 + t269 * t276;
	t260 = t273 * t264 - t278 * t294;
	t1 = [t263, 0, 0, -t267 * t279 + (t267 * t293 + t269 * t287) * t265, 0, t260; -t261, 0, 0, -t267 * t281 + (t266 * t287 - t267 * t291) * t265, 0, -t262; 0, 0, 0, -t264 * t285 - t284, 0, -t264 * t289 - t283; t262, 0, 0, t274 * t264 + t269 * t283, 0, t261; t260, 0, 0, t275 * t264 + t266 * t283, 0, t263; 0, 0, 0, t264 * t290 - t267 * t285, 0, t264 * t286 - t267 * t289; -t275, 0, 0, -t266 * t292 - t280, 0, 0; t274, 0, 0, t268 * t291 - t282, 0, 0; 0, 0, 0, t289, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end