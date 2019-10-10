% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPPRRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:52
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPPRRP4_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP4_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP4_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRP4_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP4_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:52:56
	% EndTime: 2019-10-09 23:52:56
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:52:56
	% EndTime: 2019-10-09 23:52:56
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
	% StartTime: 2019-10-09 23:52:56
	% EndTime: 2019-10-09 23:52:56
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
	% StartTime: 2019-10-09 23:52:56
	% EndTime: 2019-10-09 23:52:56
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (4->3), mult. (16->6), div. (0->0), fcn. (16->4), ass. (0->7)
	t64 = cos(qJ(1));
	t63 = sin(qJ(1));
	t62 = cos(pkin(9));
	t61 = sin(pkin(9));
	t60 = (t61 * t64 - t62 * t63) * qJD(1);
	t59 = (t61 * t63 + t62 * t64) * qJD(1);
	t1 = [-t59, 0, 0, 0, 0, 0; t60, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t60, 0, 0, 0, 0, 0; t59, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:52:56
	% EndTime: 2019-10-09 23:52:56
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (26->10), mult. (82->16), div. (0->0), fcn. (90->6), ass. (0->17)
	t67 = sin(qJ(4));
	t74 = qJD(4) * t67;
	t69 = cos(qJ(4));
	t73 = qJD(4) * t69;
	t65 = sin(pkin(9));
	t66 = cos(pkin(9));
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
	% StartTime: 2019-10-09 23:52:58
	% EndTime: 2019-10-09 23:52:58
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (109->25), mult. (317->48), div. (0->0), fcn. (353->8), ass. (0->31)
	t299 = sin(pkin(9));
	t300 = cos(pkin(9));
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
	% StartTime: 2019-10-09 23:52:59
	% EndTime: 2019-10-09 23:52:59
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (109->26), mult. (317->48), div. (0->0), fcn. (353->8), ass. (0->31)
	t367 = sin(pkin(9));
	t368 = cos(pkin(9));
	t369 = sin(qJ(1));
	t370 = cos(qJ(1));
	t341 = t370 * t367 - t369 * t368;
	t346 = sin(qJ(5));
	t348 = cos(qJ(5));
	t340 = -t369 * t367 - t370 * t368;
	t338 = t340 * qJD(1);
	t349 = cos(qJ(4));
	t347 = sin(qJ(4));
	t364 = qJD(4) * t347;
	t356 = -t338 * t349 + t341 * t364;
	t352 = -qJD(5) * t340 + t356;
	t339 = t341 * qJD(1);
	t361 = qJD(5) * t349;
	t357 = t341 * t361 + t339;
	t371 = t357 * t346 + t352 * t348;
	t366 = t346 * t347;
	t365 = t347 * t348;
	t363 = qJD(4) * t349;
	t362 = qJD(5) * t347;
	t358 = t340 * t361 + t338;
	t355 = t339 * t349 + t340 * t364;
	t354 = -t346 * t362 + t348 * t363;
	t353 = t346 * t363 + t348 * t362;
	t351 = qJD(5) * t341 + t355;
	t350 = -t352 * t346 + t357 * t348;
	t337 = t358 * t346 + t351 * t348;
	t336 = t351 * t346 - t358 * t348;
	t1 = [-t371, 0, 0, -t339 * t365 + t354 * t340, -t336, 0; t337, 0, 0, t338 * t365 + t354 * t341, t350, 0; 0, 0, 0, t346 * t361 + t348 * t364, t353, 0; t338 * t347 + t341 * t363, 0, 0, t355, 0, 0; t339 * t347 - t340 * t363, 0, 0, t356, 0, 0; 0, 0, 0, -t363, 0, 0; t350, 0, 0, -t339 * t366 + t353 * t340, t337, 0; t336, 0, 0, t338 * t366 + t353 * t341, t371, 0; 0, 0, 0, t346 * t364 - t348 * t361, -t354, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end