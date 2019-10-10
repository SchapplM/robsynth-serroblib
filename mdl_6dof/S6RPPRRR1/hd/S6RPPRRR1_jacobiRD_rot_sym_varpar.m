% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPPRRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:01
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPPRRR1_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR1_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR1_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRR1_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR1_jacobiRD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:01:23
	% EndTime: 2019-10-10 00:01:23
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:01:23
	% EndTime: 2019-10-10 00:01:23
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
	% StartTime: 2019-10-10 00:01:23
	% EndTime: 2019-10-10 00:01:23
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (7->4), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->4)
	t35 = qJ(1) + pkin(10);
	t37 = qJD(1) * sin(t35);
	t36 = qJD(1) * cos(t35);
	t1 = [-t36, 0, 0, 0, 0, 0; -t37, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t37, 0, 0, 0, 0, 0; -t36, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:01:23
	% EndTime: 2019-10-10 00:01:23
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (9->4), mult. (10->8), div. (0->0), fcn. (10->4), ass. (0->6)
	t26 = qJD(1) * sin(pkin(11));
	t25 = qJD(1) * cos(pkin(11));
	t22 = qJ(1) + pkin(10);
	t21 = cos(t22);
	t20 = sin(t22);
	t1 = [-t21 * t25, 0, 0, 0, 0, 0; -t20 * t25, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t21 * t26, 0, 0, 0, 0, 0; t20 * t26, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -qJD(1) * t20, 0, 0, 0, 0, 0; qJD(1) * t21, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:01:23
	% EndTime: 2019-10-10 00:01:23
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (47->11), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->16)
	t50 = qJ(1) + pkin(10);
	t46 = sin(t50);
	t55 = qJD(1) * t46;
	t48 = cos(t50);
	t54 = qJD(1) * t48;
	t49 = pkin(11) + qJ(4);
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
	% StartTime: 2019-10-10 00:01:24
	% EndTime: 2019-10-10 00:01:24
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (115->15), mult. (54->12), div. (0->0), fcn. (54->4), ass. (0->16)
	t79 = pkin(11) + qJ(4) + qJ(5);
	t75 = sin(t79);
	t80 = qJD(4) + qJD(5);
	t85 = t80 * t75;
	t76 = cos(t79);
	t84 = t80 * t76;
	t81 = qJ(1) + pkin(10);
	t77 = sin(t81);
	t83 = qJD(1) * t77;
	t78 = cos(t81);
	t82 = qJD(1) * t78;
	t74 = -t76 * t82 + t77 * t85;
	t73 = t75 * t82 + t77 * t84;
	t72 = t76 * t83 + t78 * t85;
	t71 = t75 * t83 - t78 * t84;
	t1 = [t74, 0, 0, t71, t71, 0; -t72, 0, 0, -t73, -t73, 0; 0, 0, 0, -t85, -t85, 0; t73, 0, 0, t72, t72, 0; t71, 0, 0, t74, t74, 0; 0, 0, 0, -t84, -t84, 0; -t83, 0, 0, 0, 0, 0; t82, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:01:25
	% EndTime: 2019-10-10 00:01:25
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (314->29), mult. (226->56), div. (0->0), fcn. (226->6), ass. (0->42)
	t347 = cos(qJ(6));
	t343 = pkin(11) + qJ(4) + qJ(5);
	t340 = cos(t343);
	t350 = qJD(1) * t340 - qJD(6);
	t366 = t347 * t350;
	t351 = qJD(6) * t340 - qJD(1);
	t339 = sin(t343);
	t344 = qJD(4) + qJD(5);
	t346 = sin(qJ(6));
	t363 = t344 * t346;
	t355 = t339 * t363;
	t365 = t351 * t347 - t355;
	t364 = t339 * t344;
	t338 = t344 * t340;
	t362 = t344 * t347;
	t345 = qJ(1) + pkin(10);
	t341 = sin(t345);
	t361 = qJD(1) * t341;
	t342 = cos(t345);
	t360 = qJD(1) * t342;
	t359 = qJD(1) * t346;
	t358 = qJD(1) * t347;
	t357 = qJD(6) * t346;
	t356 = qJD(6) * t347;
	t354 = t339 * t362;
	t353 = t340 * t363;
	t352 = t340 * t362;
	t349 = t350 * t346;
	t348 = t351 * t346 + t354;
	t337 = -t340 * t357 - t354;
	t336 = -t340 * t356 + t355;
	t335 = t340 * t360 - t341 * t364;
	t334 = -t340 * t361 - t342 * t364;
	t333 = -t341 * t352 + (t341 * t357 - t342 * t358) * t339;
	t332 = t341 * t353 + (t341 * t356 + t342 * t359) * t339;
	t331 = -t342 * t352 + (t341 * t358 + t342 * t357) * t339;
	t330 = t342 * t353 + (-t341 * t359 + t342 * t356) * t339;
	t329 = t348 * t341 - t342 * t366;
	t328 = t365 * t341 + t342 * t349;
	t327 = t341 * t366 + t348 * t342;
	t326 = t341 * t349 - t365 * t342;
	t1 = [t329, 0, 0, t331, t331, t326; -t327, 0, 0, t333, t333, -t328; 0, 0, 0, t337, t337, -t339 * t356 - t353; t328, 0, 0, t330, t330, t327; t326, 0, 0, t332, t332, t329; 0, 0, 0, t336, t336, t339 * t357 - t352; -t341 * t338 - t339 * t360, 0, 0, t334, t334, 0; t342 * t338 - t339 * t361, 0, 0, t335, t335, 0; 0, 0, 0, t338, t338, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end