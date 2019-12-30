% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPPRR6
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 16:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RPPRR6_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR6_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR6_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPRR6_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR6_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:13:23
	% EndTime: 2019-12-29 16:13:23
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:13:23
	% EndTime: 2019-12-29 16:13:23
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
	% StartTime: 2019-12-29 16:13:23
	% EndTime: 2019-12-29 16:13:23
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (7->4), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->4)
	t35 = qJ(1) + pkin(8);
	t37 = qJD(1) * sin(t35);
	t36 = qJD(1) * cos(t35);
	t1 = [-t36, 0, 0, 0, 0; -t37, 0, 0, 0, 0; 0, 0, 0, 0, 0; t37, 0, 0, 0, 0; -t36, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:13:23
	% EndTime: 2019-12-29 16:13:23
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->4), mult. (10->8), div. (0->0), fcn. (10->4), ass. (0->6)
	t26 = qJD(1) * sin(pkin(9));
	t25 = qJD(1) * cos(pkin(9));
	t22 = qJ(1) + pkin(8);
	t21 = cos(t22);
	t20 = sin(t22);
	t1 = [-t21 * t25, 0, 0, 0, 0; -t20 * t25, 0, 0, 0, 0; 0, 0, 0, 0, 0; t21 * t26, 0, 0, 0, 0; t20 * t26, 0, 0, 0, 0; 0, 0, 0, 0, 0; -qJD(1) * t20, 0, 0, 0, 0; qJD(1) * t21, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:13:24
	% EndTime: 2019-12-29 16:13:24
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (47->11), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->16)
	t50 = qJ(1) + pkin(8);
	t46 = sin(t50);
	t55 = qJD(1) * t46;
	t48 = cos(t50);
	t54 = qJD(1) * t48;
	t49 = pkin(9) + qJ(4);
	t45 = sin(t49);
	t53 = qJD(4) * t45;
	t47 = cos(t49);
	t52 = qJD(4) * t47;
	t51 = qJD(4) * t48;
	t44 = t46 * t53 - t47 * t54;
	t43 = t45 * t54 + t46 * t52;
	t42 = t45 * t51 + t47 * t55;
	t41 = t45 * t55 - t47 * t51;
	t1 = [t44, 0, 0, t41, 0; -t42, 0, 0, -t43, 0; 0, 0, 0, -t53, 0; t43, 0, 0, t42, 0; t41, 0, 0, t44, 0; 0, 0, 0, -t52, 0; -t55, 0, 0, 0, 0; t54, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:13:26
	% EndTime: 2019-12-29 16:13:26
	% DurationCPUTime: 0.30s
	% Computational Cost: add. (161->28), mult. (173->56), div. (0->0), fcn. (173->6), ass. (0->34)
	t256 = cos(qJ(5));
	t253 = pkin(9) + qJ(4);
	t251 = cos(t253);
	t274 = qJD(1) * t251;
	t259 = -qJD(5) + t274;
	t277 = t256 * t259;
	t260 = qJD(5) * t251 - qJD(1);
	t249 = sin(t253);
	t255 = sin(qJ(5));
	t268 = qJD(4) * t255;
	t264 = t249 * t268;
	t276 = t260 * t256 - t264;
	t254 = qJ(1) + pkin(8);
	t250 = sin(t254);
	t275 = t249 * t250;
	t252 = cos(t254);
	t273 = qJD(1) * t252;
	t272 = qJD(1) * t255;
	t271 = qJD(1) * t256;
	t270 = qJD(4) * t251;
	t269 = qJD(4) * t252;
	t267 = qJD(4) * t256;
	t266 = qJD(5) * t255;
	t265 = qJD(5) * t256;
	t263 = t249 * t267;
	t262 = t250 * t270;
	t261 = t251 * t269;
	t258 = t259 * t255;
	t257 = t260 * t255 + t263;
	t248 = t257 * t250 - t252 * t277;
	t247 = t276 * t250 + t252 * t258;
	t246 = t250 * t277 + t257 * t252;
	t245 = t250 * t258 - t276 * t252;
	t1 = [t248, 0, 0, -t256 * t261 + (t250 * t271 + t252 * t266) * t249, t245; -t246, 0, 0, -t256 * t262 + (t250 * t266 - t252 * t271) * t249, -t247; 0, 0, 0, -t251 * t266 - t263, -t249 * t265 - t251 * t268; t247, 0, 0, t255 * t261 + (-t250 * t272 + t252 * t265) * t249, t246; t245, 0, 0, t255 * t262 + (t250 * t265 + t252 * t272) * t249, t248; 0, 0, 0, -t251 * t265 + t264, t249 * t266 - t251 * t267; -t249 * t273 - t262, 0, 0, -t249 * t269 - t250 * t274, 0; -qJD(1) * t275 + t261, 0, 0, -qJD(4) * t275 + t251 * t273, 0; 0, 0, 0, t270, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end