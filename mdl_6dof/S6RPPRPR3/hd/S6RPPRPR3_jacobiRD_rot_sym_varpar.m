% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPPRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta5]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:37
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPPRPR3_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR3_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR3_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRPR3_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR3_jacobiRD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:37:42
	% EndTime: 2019-10-09 23:37:43
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:37:43
	% EndTime: 2019-10-09 23:37:43
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
	% StartTime: 2019-10-09 23:37:43
	% EndTime: 2019-10-09 23:37:43
	% DurationCPUTime: 0.05s
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
	% StartTime: 2019-10-09 23:37:43
	% EndTime: 2019-10-09 23:37:43
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (5->2), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->4)
	t16 = qJ(1) + pkin(9);
	t17 = qJD(1) * sin(t16);
	t14 = qJD(1) * cos(t16);
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t14, 0, 0, 0, 0, 0; t17, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t17, 0, 0, 0, 0, 0; t14, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:37:43
	% EndTime: 2019-10-09 23:37:43
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->10), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t42 = sin(qJ(4));
	t47 = qJD(1) * t42;
	t43 = cos(qJ(4));
	t46 = qJD(1) * t43;
	t45 = qJD(4) * t42;
	t44 = qJD(4) * t43;
	t41 = qJ(1) + pkin(9);
	t40 = cos(t41);
	t39 = sin(t41);
	t38 = -t39 * t45 + t40 * t46;
	t37 = t39 * t44 + t40 * t47;
	t36 = t39 * t46 + t40 * t45;
	t35 = -t39 * t47 + t40 * t44;
	t1 = [t35, 0, 0, t38, 0, 0; t37, 0, 0, t36, 0, 0; 0, 0, 0, -t44, 0, 0; -t36, 0, 0, -t37, 0, 0; t38, 0, 0, t35, 0, 0; 0, 0, 0, t45, 0, 0; -qJD(1) * t40, 0, 0, 0, 0, 0; -qJD(1) * t39, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:37:43
	% EndTime: 2019-10-09 23:37:43
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (47->11), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->16)
	t54 = qJ(1) + pkin(9);
	t50 = sin(t54);
	t59 = qJD(1) * t50;
	t52 = cos(t54);
	t58 = qJD(1) * t52;
	t53 = qJ(4) + pkin(10);
	t49 = sin(t53);
	t57 = qJD(4) * t49;
	t51 = cos(t53);
	t56 = qJD(4) * t51;
	t55 = qJD(4) * t52;
	t48 = -t50 * t57 + t51 * t58;
	t47 = t49 * t58 + t50 * t56;
	t46 = t49 * t55 + t51 * t59;
	t45 = -t49 * t59 + t51 * t55;
	t1 = [t45, 0, 0, t48, 0, 0; t47, 0, 0, t46, 0, 0; 0, 0, 0, -t56, 0, 0; -t46, 0, 0, -t47, 0, 0; t48, 0, 0, t45, 0, 0; 0, 0, 0, t57, 0, 0; -t58, 0, 0, 0, 0, 0; -t59, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:37:44
	% EndTime: 2019-10-09 23:37:44
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (162->27), mult. (173->51), div. (0->0), fcn. (173->6), ass. (0->33)
	t267 = sin(qJ(6));
	t265 = qJ(4) + pkin(10);
	t261 = sin(t265);
	t273 = qJD(1) * t261 + qJD(6);
	t288 = t267 * t273;
	t268 = cos(qJ(6));
	t287 = t268 * t273;
	t266 = qJ(1) + pkin(9);
	t262 = sin(t266);
	t286 = qJD(1) * t262;
	t264 = cos(t266);
	t285 = qJD(1) * t264;
	t284 = qJD(4) * t262;
	t283 = qJD(4) * t264;
	t282 = qJD(4) * t267;
	t281 = qJD(4) * t268;
	t280 = qJD(6) * t267;
	t279 = qJD(6) * t268;
	t263 = cos(t265);
	t278 = t263 * t285;
	t277 = t263 * t282;
	t276 = t263 * t281;
	t275 = t261 * t283;
	t274 = -qJD(6) * t261 - qJD(1);
	t272 = t261 * t281 + t263 * t280;
	t271 = t261 * t282 - t263 * t279;
	t270 = t274 * t268 - t277;
	t269 = t274 * t267 + t276;
	t260 = t269 * t262 + t264 * t287;
	t259 = t270 * t262 - t264 * t288;
	t258 = -t262 * t287 + t269 * t264;
	t257 = t262 * t288 + t270 * t264;
	t1 = [t258, 0, 0, -t272 * t262 + t268 * t278, 0, t259; t260, 0, 0, t268 * t275 + (t264 * t280 + t268 * t286) * t263, 0, -t257; 0, 0, 0, t261 * t280 - t276, 0, t271; t257, 0, 0, t271 * t262 - t267 * t278, 0, -t260; t259, 0, 0, -t267 * t275 + (t264 * t279 - t267 * t286) * t263, 0, t258; 0, 0, 0, t261 * t279 + t277, 0, t272; t263 * t286 + t275, 0, 0, t261 * t285 + t263 * t284, 0, 0; t261 * t284 - t278, 0, 0, t261 * t286 - t263 * t283, 0, 0; 0, 0, 0, -qJD(4) * t261, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end