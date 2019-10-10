% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRPRR7
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
% Datum: 2019-10-10 00:56
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPRPRR7_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR7_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR7_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR7_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR7_jacobiRD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:56:39
	% EndTime: 2019-10-10 00:56:39
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:56:39
	% EndTime: 2019-10-10 00:56:39
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0, 0; -t31, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t31, 0, 0, 0, 0, 0; -t30, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:56:39
	% EndTime: 2019-10-10 00:56:39
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
	% StartTime: 2019-10-10 00:56:39
	% EndTime: 2019-10-10 00:56:39
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
	% StartTime: 2019-10-10 00:56:39
	% EndTime: 2019-10-10 00:56:39
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->10), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t47 = sin(qJ(1));
	t52 = qJD(1) * t47;
	t48 = cos(qJ(1));
	t51 = qJD(1) * t48;
	t50 = qJD(3) * t47;
	t49 = qJD(3) * t48;
	t46 = qJ(3) + pkin(10);
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
	% StartTime: 2019-10-10 00:56:39
	% EndTime: 2019-10-10 00:56:39
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (88->13), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t78 = qJ(3) + pkin(10) + qJ(5);
	t77 = cos(t78);
	t79 = qJD(3) + qJD(5);
	t86 = t79 * t77;
	t80 = sin(qJ(1));
	t85 = t79 * t80;
	t81 = cos(qJ(1));
	t84 = t79 * t81;
	t83 = qJD(1) * t80;
	t82 = qJD(1) * t81;
	t76 = sin(t78);
	t75 = t79 * t76;
	t74 = -t76 * t85 + t77 * t82;
	t73 = t76 * t82 + t77 * t85;
	t72 = t76 * t84 + t77 * t83;
	t71 = -t76 * t83 + t77 * t84;
	t1 = [t71, 0, t74, 0, t74, 0; t73, 0, t72, 0, t72, 0; 0, 0, -t86, 0, -t86, 0; -t72, 0, -t73, 0, -t73, 0; t74, 0, t71, 0, t71, 0; 0, 0, t75, 0, t75, 0; -t82, 0, 0, 0, 0, 0; -t83, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:56:40
	% EndTime: 2019-10-10 00:56:41
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (240->31), mult. (226->57), div. (0->0), fcn. (226->6), ass. (0->41)
	t334 = cos(qJ(1));
	t329 = qJ(3) + pkin(10) + qJ(5);
	t327 = sin(t329);
	t336 = qJD(1) * t327 + qJD(6);
	t354 = t336 * t334;
	t332 = sin(qJ(1));
	t328 = cos(t329);
	t330 = qJD(3) + qJD(5);
	t348 = t330 * t334;
	t338 = t328 * t348;
	t353 = t336 * t332 - t338;
	t352 = t330 * t327;
	t331 = sin(qJ(6));
	t351 = t330 * t331;
	t350 = t330 * t332;
	t333 = cos(qJ(6));
	t349 = t330 * t333;
	t347 = qJD(1) * t332;
	t346 = qJD(1) * t334;
	t345 = qJD(6) * t331;
	t344 = qJD(6) * t333;
	t343 = qJD(6) * t334;
	t342 = t328 * t349;
	t341 = t327 * t350;
	t340 = t328 * t350;
	t339 = t327 * t348;
	t337 = -qJD(6) * t327 - qJD(1);
	t335 = t337 * t334;
	t326 = t327 * t346 + t340;
	t325 = t327 * t347 - t338;
	t324 = t327 * t345 - t342;
	t323 = t327 * t344 + t328 * t351;
	t322 = -t333 * t341 + (-t332 * t345 + t333 * t346) * t328;
	t321 = t331 * t341 + (-t331 * t346 - t332 * t344) * t328;
	t320 = t333 * t339 + (t331 * t343 + t333 * t347) * t328;
	t319 = -t331 * t339 + (-t331 * t347 + t333 * t343) * t328;
	t318 = t333 * t354 + (t337 * t331 + t342) * t332;
	t317 = t337 * t333 * t332 + (-t340 - t354) * t331;
	t316 = t331 * t335 - t353 * t333;
	t315 = t353 * t331 + t333 * t335;
	t1 = [t316, 0, t322, 0, t322, t317; t318, 0, t320, 0, t320, -t315; 0, 0, t324, 0, t324, t327 * t351 - t328 * t344; t315, 0, t321, 0, t321, -t318; t317, 0, t319, 0, t319, t316; 0, 0, t323, 0, t323, t327 * t349 + t328 * t345; t328 * t347 + t339, 0, t326, 0, t326, 0; -t328 * t346 + t341, 0, t325, 0, t325, 0; 0, 0, -t352, 0, -t352, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end