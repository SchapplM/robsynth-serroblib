% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPPPRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2,theta4]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:27
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPPPRR2_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR2_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR2_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPPRR2_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR2_jacobiRD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:27:33
	% EndTime: 2019-10-09 23:27:33
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:27:33
	% EndTime: 2019-10-09 23:27:33
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
	% StartTime: 2019-10-09 23:27:33
	% EndTime: 2019-10-09 23:27:33
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
	% StartTime: 2019-10-09 23:27:33
	% EndTime: 2019-10-09 23:27:33
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
	% StartTime: 2019-10-09 23:27:33
	% EndTime: 2019-10-09 23:27:33
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (10->5), mult. (10->8), div. (0->0), fcn. (10->4), ass. (0->6)
	t28 = qJD(1) * sin(pkin(10));
	t27 = qJD(1) * cos(pkin(10));
	t24 = qJ(1) + pkin(9);
	t23 = cos(t24);
	t22 = sin(t24);
	t1 = [-t22 * t28, 0, 0, 0, 0, 0; t23 * t28, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t22 * t27, 0, 0, 0, 0, 0; t23 * t27, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -qJD(1) * t23, 0, 0, 0, 0, 0; -qJD(1) * t22, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:27:34
	% EndTime: 2019-10-09 23:27:34
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (47->11), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->16)
	t51 = qJ(1) + pkin(9);
	t47 = sin(t51);
	t56 = qJD(1) * t47;
	t49 = cos(t51);
	t55 = qJD(1) * t49;
	t50 = pkin(10) + qJ(5);
	t46 = sin(t50);
	t54 = qJD(5) * t46;
	t48 = cos(t50);
	t53 = qJD(5) * t48;
	t52 = qJD(5) * t49;
	t45 = -t47 * t54 + t48 * t55;
	t44 = t46 * t55 + t47 * t53;
	t43 = t46 * t52 + t48 * t56;
	t42 = -t46 * t56 + t48 * t52;
	t1 = [t42, 0, 0, 0, t45, 0; t44, 0, 0, 0, t43, 0; 0, 0, 0, 0, -t53, 0; -t43, 0, 0, 0, -t44, 0; t45, 0, 0, 0, t42, 0; 0, 0, 0, 0, t54, 0; -t55, 0, 0, 0, 0, 0; -t56, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:27:35
	% EndTime: 2019-10-09 23:27:35
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (162->27), mult. (173->51), div. (0->0), fcn. (173->6), ass. (0->33)
	t255 = sin(qJ(6));
	t253 = pkin(10) + qJ(5);
	t249 = sin(t253);
	t261 = qJD(1) * t249 + qJD(6);
	t276 = t255 * t261;
	t256 = cos(qJ(6));
	t275 = t256 * t261;
	t254 = qJ(1) + pkin(9);
	t250 = sin(t254);
	t274 = qJD(1) * t250;
	t252 = cos(t254);
	t273 = qJD(1) * t252;
	t272 = qJD(5) * t250;
	t271 = qJD(5) * t252;
	t270 = qJD(5) * t255;
	t269 = qJD(5) * t256;
	t268 = qJD(6) * t255;
	t267 = qJD(6) * t256;
	t251 = cos(t253);
	t266 = t251 * t273;
	t265 = t251 * t270;
	t264 = t251 * t269;
	t263 = t249 * t271;
	t262 = -qJD(6) * t249 - qJD(1);
	t260 = t249 * t269 + t251 * t268;
	t259 = t249 * t270 - t251 * t267;
	t258 = t262 * t256 - t265;
	t257 = t262 * t255 + t264;
	t248 = t257 * t250 + t252 * t275;
	t247 = t258 * t250 - t252 * t276;
	t246 = -t250 * t275 + t257 * t252;
	t245 = t250 * t276 + t258 * t252;
	t1 = [t246, 0, 0, 0, -t260 * t250 + t256 * t266, t247; t248, 0, 0, 0, t256 * t263 + (t252 * t268 + t256 * t274) * t251, -t245; 0, 0, 0, 0, t249 * t268 - t264, t259; t245, 0, 0, 0, t259 * t250 - t255 * t266, -t248; t247, 0, 0, 0, -t255 * t263 + (t252 * t267 - t255 * t274) * t251, t246; 0, 0, 0, 0, t249 * t267 + t265, t260; t251 * t274 + t263, 0, 0, 0, t249 * t273 + t251 * t272, 0; t249 * t272 - t266, 0, 0, 0, t249 * t274 - t251 * t271, 0; 0, 0, 0, 0, -qJD(5) * t249, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end