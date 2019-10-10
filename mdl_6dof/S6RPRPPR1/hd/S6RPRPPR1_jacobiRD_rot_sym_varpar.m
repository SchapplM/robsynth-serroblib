% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRPPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4,theta5]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:15
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPRPPR1_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR1_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR1_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPPR1_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPPR1_jacobiRD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:15:13
	% EndTime: 2019-10-10 00:15:13
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:15:13
	% EndTime: 2019-10-10 00:15:13
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
	% StartTime: 2019-10-10 00:15:13
	% EndTime: 2019-10-10 00:15:13
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
	% StartTime: 2019-10-10 00:15:13
	% EndTime: 2019-10-10 00:15:13
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->10), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t41 = sin(qJ(3));
	t46 = qJD(1) * t41;
	t42 = cos(qJ(3));
	t45 = qJD(1) * t42;
	t44 = qJD(3) * t41;
	t43 = qJD(3) * t42;
	t40 = qJ(1) + pkin(9);
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
	% StartTime: 2019-10-10 00:15:14
	% EndTime: 2019-10-10 00:15:14
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (47->11), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->16)
	t57 = qJ(1) + pkin(9);
	t53 = sin(t57);
	t62 = qJD(1) * t53;
	t55 = cos(t57);
	t61 = qJD(1) * t55;
	t56 = qJ(3) + pkin(10);
	t52 = sin(t56);
	t60 = qJD(3) * t52;
	t54 = cos(t56);
	t59 = qJD(3) * t54;
	t58 = qJD(3) * t55;
	t51 = t53 * t60 - t54 * t61;
	t50 = t52 * t61 + t53 * t59;
	t49 = t52 * t58 + t54 * t62;
	t48 = t52 * t62 - t54 * t58;
	t1 = [t51, 0, t48, 0, 0, 0; -t49, 0, -t50, 0, 0, 0; 0, 0, -t60, 0, 0, 0; t50, 0, t49, 0, 0, 0; t48, 0, t51, 0, 0, 0; 0, 0, -t59, 0, 0, 0; -t62, 0, 0, 0, 0, 0; t61, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:15:15
	% EndTime: 2019-10-10 00:15:15
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (72->19), mult. (77->35), div. (0->0), fcn. (77->6), ass. (0->22)
	t215 = qJ(1) + pkin(9);
	t211 = sin(t215);
	t216 = sin(pkin(11));
	t230 = t211 * t216;
	t217 = cos(pkin(11));
	t229 = t211 * t217;
	t213 = cos(t215);
	t228 = t213 * t216;
	t227 = t213 * t217;
	t226 = qJD(1) * t211;
	t225 = qJD(1) * t213;
	t214 = qJ(3) + pkin(10);
	t210 = sin(t214);
	t224 = qJD(3) * t210;
	t212 = cos(t214);
	t223 = qJD(3) * t212;
	t222 = qJD(3) * t213;
	t221 = t211 * t224;
	t220 = t210 * t222;
	t219 = t210 * t225 + t211 * t223;
	t218 = t210 * t226 - t212 * t222;
	t1 = [t217 * t221 + (-t212 * t227 - t230) * qJD(1), 0, t218 * t217, 0, 0, 0; -t217 * t220 + (-t212 * t229 + t228) * qJD(1), 0, -t219 * t217, 0, 0, 0; 0, 0, -t217 * t224, 0, 0, 0; -t216 * t221 + (t212 * t228 - t229) * qJD(1), 0, -t218 * t216, 0, 0, 0; t216 * t220 + (t212 * t230 + t227) * qJD(1), 0, t219 * t216, 0, 0, 0; 0, 0, t216 * t224, 0, 0, 0; -t219, 0, -t212 * t226 - t220, 0, 0, 0; -t218, 0, t212 * t225 - t221, 0, 0, 0; 0, 0, t223, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:15:15
	% EndTime: 2019-10-10 00:15:15
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (221->29), mult. (173->50), div. (0->0), fcn. (173->6), ass. (0->36)
	t277 = pkin(11) + qJ(6);
	t274 = cos(t277);
	t279 = qJ(1) + pkin(9);
	t276 = cos(t279);
	t301 = t274 * t276;
	t273 = sin(t279);
	t300 = qJD(1) * t273;
	t278 = qJ(3) + pkin(10);
	t275 = cos(t278);
	t299 = qJD(1) * t275;
	t298 = qJD(1) * t276;
	t272 = sin(t278);
	t297 = qJD(3) * t272;
	t296 = qJD(3) * t275;
	t295 = qJD(3) * t276;
	t271 = sin(t277);
	t294 = qJD(6) * t271;
	t293 = qJD(6) * t272;
	t292 = qJD(6) * t275;
	t291 = t274 * t297;
	t290 = t274 * t293;
	t289 = t273 * t297;
	t288 = t273 * t296;
	t287 = t272 * t295;
	t286 = t275 * t295;
	t285 = -qJD(1) + t292;
	t284 = -qJD(6) + t299;
	t283 = t285 * t271;
	t282 = t272 * t298 + t288;
	t281 = -t272 * t300 + t286;
	t280 = t284 * t273 + t287;
	t270 = -t284 * t301 + (t283 + t291) * t273;
	t269 = t285 * t274 * t273 + (t284 * t276 - t289) * t271;
	t268 = t280 * t274 + t276 * t283;
	t267 = t280 * t271 - t285 * t301;
	t1 = [t270, 0, -t274 * t286 + (t274 * t300 + t276 * t294) * t272, 0, 0, t267; -t268, 0, -t274 * t288 + (t273 * t294 - t274 * t298) * t272, 0, 0, -t269; 0, 0, -t271 * t292 - t291, 0, 0, -t271 * t296 - t290; t269, 0, t281 * t271 + t276 * t290, 0, 0, t268; t267, 0, t282 * t271 + t273 * t290, 0, 0, t270; 0, 0, t271 * t297 - t274 * t292, 0, 0, t271 * t293 - t274 * t296; -t282, 0, -t273 * t299 - t287, 0, 0, 0; t281, 0, t275 * t298 - t289, 0, 0, 0; 0, 0, t296, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end