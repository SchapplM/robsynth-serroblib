% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRPPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4,theta5]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:17
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRRPPR1_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR1_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR1_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPPR1_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR1_jacobiRD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:17:04
	% EndTime: 2019-10-10 11:17:04
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:17:04
	% EndTime: 2019-10-10 11:17:04
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
	% StartTime: 2019-10-10 11:17:04
	% EndTime: 2019-10-10 11:17:04
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (11->9), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->14)
	t34 = sin(qJ(1));
	t41 = qJD(1) * t34;
	t36 = cos(qJ(1));
	t40 = qJD(1) * t36;
	t33 = sin(qJ(2));
	t39 = qJD(2) * t33;
	t35 = cos(qJ(2));
	t38 = qJD(2) * t35;
	t37 = qJD(2) * t36;
	t32 = t34 * t39 - t35 * t40;
	t31 = t33 * t40 + t34 * t38;
	t30 = t33 * t37 + t35 * t41;
	t29 = t33 * t41 - t35 * t37;
	t1 = [t32, t29, 0, 0, 0, 0; -t30, -t31, 0, 0, 0, 0; 0, -t39, 0, 0, 0, 0; t31, t30, 0, 0, 0, 0; t29, t32, 0, 0, 0, 0; 0, -t38, 0, 0, 0, 0; -t41, 0, 0, 0, 0, 0; t40, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:17:04
	% EndTime: 2019-10-10 11:17:04
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (61->14), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t71 = qJ(2) + qJ(3);
	t68 = sin(t71);
	t70 = qJD(2) + qJD(3);
	t79 = t70 * t68;
	t69 = cos(t71);
	t78 = t70 * t69;
	t72 = sin(qJ(1));
	t77 = t70 * t72;
	t73 = cos(qJ(1));
	t76 = t70 * t73;
	t75 = qJD(1) * t72;
	t74 = qJD(1) * t73;
	t67 = t68 * t77 - t69 * t74;
	t66 = t68 * t74 + t69 * t77;
	t65 = t68 * t76 + t69 * t75;
	t64 = t68 * t75 - t69 * t76;
	t1 = [t67, t64, t64, 0, 0, 0; -t65, -t66, -t66, 0, 0, 0; 0, -t79, -t79, 0, 0, 0; t66, t65, t65, 0, 0, 0; t64, t67, t67, 0, 0, 0; 0, -t78, -t78, 0, 0, 0; -t75, 0, 0, 0, 0, 0; t74, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:17:04
	% EndTime: 2019-10-10 11:17:04
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (89->14), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t85 = qJ(2) + qJ(3) + pkin(10);
	t83 = sin(t85);
	t86 = qJD(2) + qJD(3);
	t94 = t86 * t83;
	t84 = cos(t85);
	t93 = t86 * t84;
	t87 = sin(qJ(1));
	t92 = t86 * t87;
	t88 = cos(qJ(1));
	t91 = t86 * t88;
	t90 = qJD(1) * t87;
	t89 = qJD(1) * t88;
	t82 = t83 * t92 - t84 * t89;
	t81 = t83 * t89 + t84 * t92;
	t80 = t83 * t91 + t84 * t90;
	t79 = t83 * t90 - t84 * t91;
	t1 = [t82, t79, t79, 0, 0, 0; -t80, -t81, -t81, 0, 0, 0; 0, -t94, -t94, 0, 0, 0; t81, t80, t80, 0, 0, 0; t79, t82, t82, 0, 0, 0; 0, -t93, -t93, 0, 0, 0; -t90, 0, 0, 0, 0, 0; t89, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:17:06
	% EndTime: 2019-10-10 11:17:06
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (132->22), mult. (114->36), div. (0->0), fcn. (114->6), ass. (0->31)
	t292 = qJ(2) + qJ(3) + pkin(10);
	t290 = sin(t292);
	t293 = qJD(2) + qJD(3);
	t311 = t290 * t293;
	t296 = sin(qJ(1));
	t310 = t293 * t296;
	t297 = cos(qJ(1));
	t309 = t293 * t297;
	t294 = sin(pkin(11));
	t308 = t294 * t296;
	t307 = t294 * t297;
	t295 = cos(pkin(11));
	t306 = t295 * t296;
	t305 = t295 * t297;
	t304 = qJD(1) * t296;
	t303 = qJD(1) * t297;
	t302 = t295 * t311;
	t301 = t290 * t310;
	t300 = t290 * t309;
	t291 = cos(t292);
	t299 = t290 * t303 + t291 * t310;
	t298 = t290 * t304 - t291 * t309;
	t289 = t293 * t291;
	t288 = t294 * t311;
	t287 = t291 * t303 - t301;
	t286 = -t291 * t304 - t300;
	t285 = t299 * t295;
	t284 = t299 * t294;
	t283 = t298 * t295;
	t282 = t298 * t294;
	t1 = [t295 * t301 + (-t291 * t305 - t308) * qJD(1), t283, t283, 0, 0, 0; -t295 * t300 + (-t291 * t306 + t307) * qJD(1), -t285, -t285, 0, 0, 0; 0, -t302, -t302, 0, 0, 0; -t294 * t301 + (t291 * t307 - t306) * qJD(1), -t282, -t282, 0, 0, 0; t294 * t300 + (t291 * t308 + t305) * qJD(1), t284, t284, 0, 0, 0; 0, t288, t288, 0, 0, 0; -t299, t286, t286, 0, 0, 0; -t298, t287, t287, 0, 0, 0; 0, t289, t289, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:17:06
	% EndTime: 2019-10-10 11:17:06
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (314->29), mult. (226->56), div. (0->0), fcn. (226->6), ass. (0->42)
	t364 = cos(qJ(1));
	t360 = qJ(2) + qJ(3) + pkin(10);
	t357 = cos(t360);
	t368 = qJD(6) * t357 - qJD(1);
	t383 = t364 * t368;
	t367 = qJD(1) * t357 - qJD(6);
	t356 = sin(t360);
	t362 = qJD(2) + qJD(3);
	t363 = sin(qJ(1));
	t380 = t362 * t363;
	t372 = t356 * t380;
	t382 = t367 * t364 - t372;
	t381 = t356 * t362;
	t355 = t362 * t357;
	t379 = t362 * t364;
	t378 = qJD(1) * t363;
	t377 = qJD(1) * t364;
	t361 = pkin(11) + qJ(6);
	t358 = sin(t361);
	t376 = qJD(6) * t358;
	t359 = cos(t361);
	t375 = qJD(6) * t359;
	t374 = qJD(6) * t363;
	t373 = qJD(6) * t364;
	t371 = t357 * t380;
	t370 = t356 * t379;
	t369 = t357 * t379;
	t366 = t368 * t363;
	t365 = t367 * t363 + t370;
	t354 = t357 * t377 - t372;
	t353 = -t357 * t378 - t370;
	t352 = -t357 * t376 - t359 * t381;
	t351 = -t357 * t375 + t358 * t381;
	t350 = -t359 * t371 + (t358 * t374 - t359 * t377) * t356;
	t349 = t358 * t371 + (t358 * t377 + t359 * t374) * t356;
	t348 = -t359 * t369 + (t358 * t373 + t359 * t378) * t356;
	t347 = t358 * t369 + (-t358 * t378 + t359 * t373) * t356;
	t346 = t358 * t366 - t382 * t359;
	t345 = t382 * t358 + t359 * t366;
	t344 = t358 * t383 + t365 * t359;
	t343 = t365 * t358 - t359 * t383;
	t1 = [t346, t348, t348, 0, 0, t343; -t344, t350, t350, 0, 0, -t345; 0, t352, t352, 0, 0, -t358 * t355 - t356 * t375; t345, t347, t347, 0, 0, t344; t343, t349, t349, 0, 0, t346; 0, t351, t351, 0, 0, -t359 * t355 + t356 * t376; -t356 * t377 - t371, t353, t353, 0, 0, 0; -t356 * t378 + t369, t354, t354, 0, 0, 0; 0, t355, t355, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end