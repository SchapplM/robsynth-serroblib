% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRPPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:17
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPRPPR2_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR2_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR2_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPPR2_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR2_jacobiRD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:17:00
	% EndTime: 2019-10-10 00:17:00
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:17:00
	% EndTime: 2019-10-10 00:17:00
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0, 0; -t31, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t31, 0, 0, 0, 0, 0; -t30, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:17:00
	% EndTime: 2019-10-10 00:17:00
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
	% StartTime: 2019-10-10 00:17:00
	% EndTime: 2019-10-10 00:17:00
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
	% StartTime: 2019-10-10 00:17:00
	% EndTime: 2019-10-10 00:17:00
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
	% StartTime: 2019-10-10 00:17:01
	% EndTime: 2019-10-10 00:17:01
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (45->9), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->16)
	t190 = qJ(1) + pkin(9);
	t186 = sin(t190);
	t195 = qJD(1) * t186;
	t188 = cos(t190);
	t194 = qJD(1) * t188;
	t189 = qJ(3) + pkin(10);
	t185 = sin(t189);
	t193 = qJD(3) * t185;
	t187 = cos(t189);
	t192 = qJD(3) * t187;
	t191 = qJD(3) * t188;
	t184 = -t186 * t193 + t187 * t194;
	t183 = t185 * t194 + t186 * t192;
	t182 = t185 * t191 + t187 * t195;
	t181 = -t185 * t195 + t187 * t191;
	t1 = [-t195, 0, 0, 0, 0, 0; t194, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t184, 0, t181, 0, 0, 0; t182, 0, t183, 0, 0, 0; 0, 0, t193, 0, 0, 0; -t183, 0, -t182, 0, 0, 0; t181, 0, t184, 0, 0, 0; 0, 0, t192, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:17:01
	% EndTime: 2019-10-10 00:17:02
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (162->27), mult. (173->51), div. (0->0), fcn. (173->6), ass. (0->33)
	t269 = sin(qJ(6));
	t267 = qJ(3) + pkin(10);
	t263 = sin(t267);
	t276 = qJD(6) * t263 + qJD(1);
	t265 = cos(t267);
	t270 = cos(qJ(6));
	t283 = qJD(3) * t270;
	t278 = t265 * t283;
	t290 = t276 * t269 - t278;
	t284 = qJD(3) * t269;
	t279 = t265 * t284;
	t289 = t276 * t270 + t279;
	t268 = qJ(1) + pkin(9);
	t264 = sin(t268);
	t288 = qJD(1) * t264;
	t266 = cos(t268);
	t287 = qJD(1) * t266;
	t286 = qJD(3) * t264;
	t285 = qJD(3) * t266;
	t282 = qJD(6) * t269;
	t281 = qJD(6) * t270;
	t280 = t265 * t287;
	t277 = t263 * t285;
	t275 = -qJD(1) * t263 - qJD(6);
	t274 = t275 * t269;
	t273 = t275 * t270;
	t272 = t263 * t283 + t265 * t282;
	t271 = -t263 * t284 + t265 * t281;
	t262 = t264 * t274 + t289 * t266;
	t261 = t264 * t273 - t290 * t266;
	t260 = -t289 * t264 + t266 * t274;
	t259 = t290 * t264 + t266 * t273;
	t1 = [t260, 0, -t269 * t277 + (t266 * t281 - t269 * t288) * t265, 0, 0, t261; t262, 0, t271 * t264 + t269 * t280, 0, 0, -t259; 0, 0, t263 * t281 + t279, 0, 0, t272; t259, 0, -t270 * t277 + (-t266 * t282 - t270 * t288) * t265, 0, 0, -t262; t261, 0, -t272 * t264 + t270 * t280, 0, 0, t260; 0, 0, -t263 * t282 + t278, 0, 0, t271; t263 * t286 - t280, 0, t263 * t288 - t265 * t285, 0, 0, 0; -t265 * t288 - t277, 0, -t263 * t287 - t265 * t286, 0, 0, 0; 0, 0, -qJD(3) * t263, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end