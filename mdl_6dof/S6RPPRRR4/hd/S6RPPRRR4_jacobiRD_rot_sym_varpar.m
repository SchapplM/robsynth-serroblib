% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPPRRR4
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
% Datum: 2019-10-10 00:06
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPPRRR4_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR4_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR4_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRR4_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR4_jacobiRD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:06:31
	% EndTime: 2019-10-10 00:06:31
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:06:31
	% EndTime: 2019-10-10 00:06:31
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
	% StartTime: 2019-10-10 00:06:31
	% EndTime: 2019-10-10 00:06:32
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t14 = qJD(1) * sin(qJ(1));
	t13 = qJD(1) * cos(qJ(1));
	t1 = [-t13, 0, 0, 0, 0, 0; -t14, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t14, 0, 0, 0, 0, 0; t13, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:06:32
	% EndTime: 2019-10-10 00:06:32
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (4->3), mult. (16->6), div. (0->0), fcn. (16->4), ass. (0->7)
	t64 = cos(qJ(1));
	t63 = sin(qJ(1));
	t62 = cos(pkin(10));
	t61 = sin(pkin(10));
	t60 = (t61 * t64 - t62 * t63) * qJD(1);
	t59 = (t61 * t63 + t62 * t64) * qJD(1);
	t1 = [-t59, 0, 0, 0, 0, 0; t60, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t60, 0, 0, 0, 0, 0; t59, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:06:32
	% EndTime: 2019-10-10 00:06:32
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (26->10), mult. (82->16), div. (0->0), fcn. (90->6), ass. (0->17)
	t67 = sin(qJ(4));
	t74 = qJD(4) * t67;
	t69 = cos(qJ(4));
	t73 = qJD(4) * t69;
	t65 = sin(pkin(10));
	t66 = cos(pkin(10));
	t68 = sin(qJ(1));
	t70 = cos(qJ(1));
	t62 = t70 * t65 - t68 * t66;
	t63 = t68 * t65 + t70 * t66;
	t60 = t63 * qJD(1);
	t72 = -t60 * t67 + t62 * t73;
	t71 = -t60 * t69 - t62 * t74;
	t61 = t62 * qJD(1);
	t59 = t61 * t69 - t63 * t74;
	t58 = -t61 * t67 - t63 * t73;
	t1 = [t71, 0, 0, t58, 0, 0; t59, 0, 0, t72, 0, 0; 0, 0, 0, t74, 0, 0; -t72, 0, 0, -t59, 0, 0; t58, 0, 0, t71, 0, 0; 0, 0, 0, t73, 0, 0; -t61, 0, 0, 0, 0, 0; -t60, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:06:33
	% EndTime: 2019-10-10 00:06:33
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (109->25), mult. (317->48), div. (0->0), fcn. (353->8), ass. (0->31)
	t299 = sin(pkin(10));
	t300 = cos(pkin(10));
	t301 = sin(qJ(1));
	t302 = cos(qJ(1));
	t273 = t302 * t299 - t301 * t300;
	t278 = sin(qJ(5));
	t280 = cos(qJ(5));
	t272 = -t301 * t299 - t302 * t300;
	t270 = t272 * qJD(1);
	t281 = cos(qJ(4));
	t279 = sin(qJ(4));
	t296 = qJD(4) * t279;
	t288 = -t270 * t281 + t273 * t296;
	t284 = -qJD(5) * t272 + t288;
	t271 = t273 * qJD(1);
	t293 = qJD(5) * t281;
	t289 = t273 * t293 + t271;
	t303 = t284 * t278 - t289 * t280;
	t298 = t278 * t279;
	t297 = t279 * t280;
	t295 = qJD(4) * t281;
	t294 = qJD(5) * t279;
	t290 = t272 * t293 + t270;
	t287 = t271 * t281 + t272 * t296;
	t286 = -t278 * t294 + t280 * t295;
	t285 = t278 * t295 + t280 * t294;
	t283 = qJD(5) * t273 + t287;
	t282 = -t289 * t278 - t284 * t280;
	t269 = t290 * t278 + t283 * t280;
	t268 = -t283 * t278 + t290 * t280;
	t1 = [t282, 0, 0, -t271 * t297 + t286 * t272, t268, 0; t269, 0, 0, t270 * t297 + t286 * t273, -t303, 0; 0, 0, 0, t278 * t293 + t280 * t296, t285, 0; t303, 0, 0, t271 * t298 - t285 * t272, -t269, 0; t268, 0, 0, -t270 * t298 - t285 * t273, t282, 0; 0, 0, 0, -t278 * t296 + t280 * t293, t286, 0; t270 * t279 + t273 * t295, 0, 0, t287, 0, 0; t271 * t279 - t272 * t295, 0, 0, t288, 0, 0; 0, 0, 0, -t295, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:06:33
	% EndTime: 2019-10-10 00:06:33
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (265->30), mult. (425->49), div. (0->0), fcn. (473->8), ass. (0->33)
	t365 = sin(pkin(10));
	t366 = cos(pkin(10));
	t367 = sin(qJ(1));
	t368 = cos(qJ(1));
	t338 = t368 * t365 - t367 * t366;
	t345 = qJD(5) + qJD(6);
	t336 = t338 * qJD(1);
	t337 = -t367 * t365 - t368 * t366;
	t348 = cos(qJ(4));
	t347 = sin(qJ(4));
	t358 = qJD(4) * t347;
	t350 = t336 * t348 + t337 * t358;
	t370 = -t338 * t345 - t350;
	t346 = qJ(5) + qJ(6);
	t343 = sin(t346);
	t344 = cos(t346);
	t335 = t337 * qJD(1);
	t351 = -t335 * t348 + t338 * t358;
	t349 = -t337 * t345 + t351;
	t359 = t345 * t348;
	t352 = t338 * t359 + t336;
	t369 = t349 * t343 - t352 * t344;
	t362 = t343 * t347;
	t361 = t344 * t347;
	t360 = t345 * t347;
	t357 = qJD(4) * t348;
	t356 = t343 * t359;
	t334 = -t343 * t360 + t344 * t357;
	t333 = t343 * t357 + t344 * t360;
	t328 = -t352 * t343 - t349 * t344;
	t330 = t335 * t343 + t337 * t356 - t370 * t344;
	t329 = (t337 * t359 + t335) * t344 + t370 * t343;
	t1 = [t328, 0, 0, t334 * t337 - t336 * t361, t329, t329; t330, 0, 0, t334 * t338 + t335 * t361, -t369, -t369; 0, 0, 0, t344 * t358 + t356, t333, t333; t369, 0, 0, -t333 * t337 + t336 * t362, -t330, -t330; t329, 0, 0, -t333 * t338 - t335 * t362, t328, t328; 0, 0, 0, -t343 * t358 + t344 * t359, t334, t334; t335 * t347 + t338 * t357, 0, 0, t350, 0, 0; t336 * t347 - t337 * t357, 0, 0, t351, 0, 0; 0, 0, 0, -t357, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end