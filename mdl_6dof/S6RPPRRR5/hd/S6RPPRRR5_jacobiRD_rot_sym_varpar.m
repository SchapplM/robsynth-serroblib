% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPPRRR5
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:08
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPPRRR5_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR5_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR5_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRR5_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR5_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:08:19
	% EndTime: 2019-10-10 00:08:19
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:08:19
	% EndTime: 2019-10-10 00:08:19
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
	% StartTime: 2019-10-10 00:08:19
	% EndTime: 2019-10-10 00:08:19
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
	% StartTime: 2019-10-10 00:08:19
	% EndTime: 2019-10-10 00:08:19
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t16 = qJD(1) * sin(qJ(1));
	t15 = qJD(1) * cos(qJ(1));
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t16, 0, 0, 0, 0, 0; t15, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t15, 0, 0, 0, 0, 0; -t16, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:08:19
	% EndTime: 2019-10-10 00:08:19
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (10->8), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->14)
	t35 = sin(qJ(1));
	t42 = qJD(1) * t35;
	t37 = cos(qJ(1));
	t41 = qJD(1) * t37;
	t34 = sin(qJ(4));
	t40 = qJD(4) * t34;
	t36 = cos(qJ(4));
	t39 = qJD(4) * t36;
	t38 = qJD(4) * t37;
	t33 = -t35 * t40 + t36 * t41;
	t32 = -t34 * t41 - t35 * t39;
	t31 = -t34 * t38 - t36 * t42;
	t30 = t34 * t42 - t36 * t38;
	t1 = [t32, 0, 0, t31, 0, 0; -t30, 0, 0, t33, 0, 0; 0, 0, 0, -t39, 0, 0; -t33, 0, 0, t30, 0, 0; t31, 0, 0, t32, 0, 0; 0, 0, 0, t40, 0, 0; t42, 0, 0, 0, 0, 0; -t41, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:08:19
	% EndTime: 2019-10-10 00:08:19
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (59->11), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t74 = qJ(4) + qJ(5);
	t72 = cos(t74);
	t73 = qJD(4) + qJD(5);
	t81 = t73 * t72;
	t75 = sin(qJ(1));
	t80 = t73 * t75;
	t76 = cos(qJ(1));
	t79 = t73 * t76;
	t78 = qJD(1) * t75;
	t77 = qJD(1) * t76;
	t71 = sin(t74);
	t70 = t73 * t71;
	t69 = -t71 * t80 + t72 * t77;
	t68 = -t71 * t77 - t72 * t80;
	t67 = -t71 * t79 - t72 * t78;
	t66 = t71 * t78 - t72 * t79;
	t1 = [t68, 0, 0, t67, t67, 0; -t66, 0, 0, t69, t69, 0; 0, 0, 0, -t81, -t81, 0; -t69, 0, 0, t66, t66, 0; t67, 0, 0, t68, t68, 0; 0, 0, 0, t70, t70, 0; t78, 0, 0, 0, 0, 0; -t77, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:08:20
	% EndTime: 2019-10-10 00:08:20
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (166->31), mult. (226->57), div. (0->0), fcn. (226->6), ass. (0->41)
	t328 = qJ(4) + qJ(5);
	t325 = sin(t328);
	t327 = qJD(4) + qJD(5);
	t352 = t327 * t325;
	t329 = sin(qJ(6));
	t351 = t327 * t329;
	t330 = sin(qJ(1));
	t350 = t327 * t330;
	t331 = cos(qJ(6));
	t349 = t327 * t331;
	t332 = cos(qJ(1));
	t348 = t327 * t332;
	t347 = t331 * t332;
	t346 = qJD(1) * t330;
	t345 = qJD(1) * t332;
	t344 = qJD(6) * t329;
	t343 = qJD(6) * t331;
	t342 = qJD(6) * t332;
	t341 = t325 * t348;
	t326 = cos(t328);
	t340 = t326 * t349;
	t339 = t325 * t350;
	t338 = t326 * t350;
	t337 = t326 * t348;
	t336 = qJD(6) * t325 + qJD(1);
	t335 = qJD(1) * t325 + qJD(6);
	t334 = t336 * t329;
	t333 = t335 * t330 - t337;
	t324 = t325 * t345 + t338;
	t323 = -t325 * t346 + t337;
	t322 = t325 * t344 - t340;
	t321 = t325 * t343 + t326 * t351;
	t320 = -t331 * t339 + (-t330 * t344 + t331 * t345) * t326;
	t319 = t329 * t339 + (-t329 * t345 - t330 * t343) * t326;
	t318 = -t331 * t341 + (-t329 * t342 - t331 * t346) * t326;
	t317 = t329 * t341 + (t329 * t346 - t331 * t342) * t326;
	t316 = -t335 * t347 + (t334 - t340) * t330;
	t315 = t336 * t331 * t330 + (t335 * t332 + t338) * t329;
	t314 = t333 * t331 + t332 * t334;
	t313 = t333 * t329 - t336 * t347;
	t1 = [t316, 0, 0, t318, t318, t313; -t314, 0, 0, t320, t320, -t315; 0, 0, 0, t322, t322, t325 * t351 - t326 * t343; t315, 0, 0, t317, t317, t314; t313, 0, 0, t319, t319, t316; 0, 0, 0, t321, t321, t325 * t349 + t326 * t344; t326 * t345 - t339, 0, 0, t323, t323, 0; t326 * t346 + t341, 0, 0, t324, t324, 0; 0, 0, 0, -t352, -t352, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end