% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPPRRR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:11
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPPRRR7_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR7_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR7_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRR7_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR7_jacobiRD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:11:46
	% EndTime: 2019-10-10 00:11:46
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:11:46
	% EndTime: 2019-10-10 00:11:46
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
	% StartTime: 2019-10-10 00:11:46
	% EndTime: 2019-10-10 00:11:46
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
	% StartTime: 2019-10-10 00:11:46
	% EndTime: 2019-10-10 00:11:46
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (4->4), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->5)
	t24 = qJD(1) * sin(qJ(1));
	t23 = qJD(1) * cos(qJ(1));
	t20 = cos(pkin(10));
	t19 = sin(pkin(10));
	t1 = [-t19 * t24, 0, 0, 0, 0, 0; t19 * t23, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t20 * t24, 0, 0, 0, 0, 0; t20 * t23, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t23, 0, 0, 0, 0, 0; -t24, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:11:46
	% EndTime: 2019-10-10 00:11:46
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->10), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t44 = sin(qJ(1));
	t49 = qJD(1) * t44;
	t45 = cos(qJ(1));
	t48 = qJD(1) * t45;
	t47 = qJD(4) * t44;
	t46 = qJD(4) * t45;
	t43 = pkin(10) + qJ(4);
	t42 = cos(t43);
	t41 = sin(t43);
	t40 = -t41 * t47 + t42 * t48;
	t39 = t41 * t48 + t42 * t47;
	t38 = t41 * t46 + t42 * t49;
	t37 = -t41 * t49 + t42 * t46;
	t1 = [t37, 0, 0, t40, 0, 0; t39, 0, 0, t38, 0, 0; 0, 0, 0, -qJD(4) * t42, 0, 0; -t38, 0, 0, -t39, 0, 0; t40, 0, 0, t37, 0, 0; 0, 0, 0, qJD(4) * t41, 0, 0; -t48, 0, 0, 0, 0, 0; -t49, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:11:46
	% EndTime: 2019-10-10 00:11:46
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (88->13), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t75 = pkin(10) + qJ(4) + qJ(5);
	t74 = cos(t75);
	t76 = qJD(4) + qJD(5);
	t83 = t76 * t74;
	t77 = sin(qJ(1));
	t82 = t76 * t77;
	t78 = cos(qJ(1));
	t81 = t76 * t78;
	t80 = qJD(1) * t77;
	t79 = qJD(1) * t78;
	t73 = sin(t75);
	t72 = t76 * t73;
	t71 = -t73 * t82 + t74 * t79;
	t70 = t73 * t79 + t74 * t82;
	t69 = t73 * t81 + t74 * t80;
	t68 = -t73 * t80 + t74 * t81;
	t1 = [t68, 0, 0, t71, t71, 0; t70, 0, 0, t69, t69, 0; 0, 0, 0, -t83, -t83, 0; -t69, 0, 0, -t70, -t70, 0; t71, 0, 0, t68, t68, 0; 0, 0, 0, t72, t72, 0; -t79, 0, 0, 0, 0, 0; -t80, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:11:48
	% EndTime: 2019-10-10 00:11:48
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (240->31), mult. (226->57), div. (0->0), fcn. (226->6), ass. (0->41)
	t332 = cos(qJ(1));
	t327 = pkin(10) + qJ(4) + qJ(5);
	t325 = sin(t327);
	t334 = qJD(1) * t325 + qJD(6);
	t352 = t334 * t332;
	t330 = sin(qJ(1));
	t326 = cos(t327);
	t328 = qJD(4) + qJD(5);
	t346 = t328 * t332;
	t336 = t326 * t346;
	t351 = t334 * t330 - t336;
	t350 = t328 * t325;
	t329 = sin(qJ(6));
	t349 = t328 * t329;
	t348 = t328 * t330;
	t331 = cos(qJ(6));
	t347 = t328 * t331;
	t345 = qJD(1) * t330;
	t344 = qJD(1) * t332;
	t343 = qJD(6) * t329;
	t342 = qJD(6) * t331;
	t341 = qJD(6) * t332;
	t340 = t326 * t347;
	t339 = t325 * t348;
	t338 = t326 * t348;
	t337 = t325 * t346;
	t335 = -qJD(6) * t325 - qJD(1);
	t333 = t335 * t332;
	t324 = t325 * t344 + t338;
	t323 = t325 * t345 - t336;
	t322 = t325 * t343 - t340;
	t321 = t325 * t342 + t326 * t349;
	t320 = -t331 * t339 + (-t330 * t343 + t331 * t344) * t326;
	t319 = t329 * t339 + (-t329 * t344 - t330 * t342) * t326;
	t318 = t331 * t337 + (t329 * t341 + t331 * t345) * t326;
	t317 = -t329 * t337 + (-t329 * t345 + t331 * t341) * t326;
	t316 = t331 * t352 + (t335 * t329 + t340) * t330;
	t315 = t335 * t331 * t330 + (-t338 - t352) * t329;
	t314 = t329 * t333 - t351 * t331;
	t313 = t351 * t329 + t331 * t333;
	t1 = [t314, 0, 0, t320, t320, t315; t316, 0, 0, t318, t318, -t313; 0, 0, 0, t322, t322, t325 * t349 - t326 * t342; t313, 0, 0, t319, t319, -t316; t315, 0, 0, t317, t317, t314; 0, 0, 0, t321, t321, t325 * t347 + t326 * t343; t326 * t345 + t337, 0, 0, t324, t324, 0; -t326 * t344 + t339, 0, 0, t323, t323, 0; 0, 0, 0, -t350, -t350, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end