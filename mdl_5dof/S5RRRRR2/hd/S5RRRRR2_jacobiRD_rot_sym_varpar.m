% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRRRR2
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RRRRR2_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR2_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR2_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR2_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5RRRRR2_jacobiRD_rot_sym_varpar: pkin has to be [2x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:54:12
	% EndTime: 2019-12-05 18:54:12
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:54:12
	% EndTime: 2019-12-05 18:54:12
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0; -t31, 0, 0, 0, 0; 0, 0, 0, 0, 0; t31, 0, 0, 0, 0; -t30, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:54:13
	% EndTime: 2019-12-05 18:54:13
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (22->8), mult. (8->2), div. (0->0), fcn. (8->2), ass. (0->5)
	t47 = qJD(1) + qJD(2);
	t48 = qJ(1) + qJ(2);
	t49 = t47 * cos(t48);
	t44 = t47 * sin(t48);
	t1 = [-t49, -t49, 0, 0, 0; -t44, -t44, 0, 0, 0; 0, 0, 0, 0, 0; t44, t44, 0, 0, 0; -t49, -t49, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:54:13
	% EndTime: 2019-12-05 18:54:13
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (60->13), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t62 = qJ(1) + qJ(2);
	t59 = sin(t62);
	t61 = qJD(1) + qJD(2);
	t69 = t61 * t59;
	t63 = sin(qJ(3));
	t68 = t61 * t63;
	t64 = cos(qJ(3));
	t67 = t61 * t64;
	t66 = qJD(3) * t63;
	t65 = qJD(3) * t64;
	t60 = cos(t62);
	t58 = t61 * t60;
	t57 = t59 * t66 - t60 * t67;
	t56 = t59 * t65 + t60 * t68;
	t55 = t59 * t67 + t60 * t66;
	t54 = t59 * t68 - t60 * t65;
	t1 = [t57, t57, t54, 0, 0; -t55, -t55, -t56, 0, 0; 0, 0, -t66, 0, 0; t56, t56, t55, 0, 0; t54, t54, t57, 0, 0; 0, 0, -t65, 0, 0; -t69, -t69, 0, 0, 0; t58, t58, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:54:13
	% EndTime: 2019-12-05 18:54:13
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (134->18), mult. (72->12), div. (0->0), fcn. (72->4), ass. (0->17)
	t95 = qJ(3) + qJ(4);
	t89 = sin(t95);
	t93 = qJD(3) + qJD(4);
	t99 = t93 * t89;
	t91 = cos(t95);
	t98 = t93 * t91;
	t96 = qJ(1) + qJ(2);
	t90 = sin(t96);
	t94 = qJD(1) + qJD(2);
	t97 = t94 * t90;
	t92 = cos(t96);
	t88 = t94 * t92;
	t87 = -t91 * t88 + t90 * t99;
	t86 = t89 * t88 + t90 * t98;
	t85 = t91 * t97 + t92 * t99;
	t84 = t89 * t97 - t92 * t98;
	t1 = [t87, t87, t84, t84, 0; -t85, -t85, -t86, -t86, 0; 0, 0, -t99, -t99, 0; t86, t86, t85, t85, 0; t84, t84, t87, t87, 0; 0, 0, -t98, -t98, 0; -t97, -t97, 0, 0, 0; t88, t88, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:54:14
	% EndTime: 2019-12-05 18:54:14
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (344->32), mult. (286->56), div. (0->0), fcn. (286->6), ass. (0->45)
	t349 = cos(qJ(5));
	t346 = qJ(3) + qJ(4);
	t342 = cos(t346);
	t345 = qJD(1) + qJD(2);
	t353 = qJD(5) * t342 - t345;
	t340 = sin(t346);
	t344 = qJD(3) + qJD(4);
	t348 = sin(qJ(5));
	t364 = t344 * t348;
	t358 = t340 * t364;
	t368 = t353 * t349 - t358;
	t367 = t340 * t344;
	t366 = t342 * t345;
	t347 = qJ(1) + qJ(2);
	t343 = cos(t347);
	t365 = t343 * t345;
	t339 = t344 * t342;
	t363 = t344 * t349;
	t362 = t345 * t348;
	t361 = t345 * t349;
	t360 = qJD(5) * t348;
	t359 = qJD(5) * t349;
	t357 = t340 * t363;
	t356 = t342 * t364;
	t355 = t342 * t363;
	t354 = -qJD(5) + t366;
	t352 = t354 * t348;
	t341 = sin(t347);
	t351 = t341 * t361 + t343 * t360;
	t350 = t341 * t359 + t343 * t362;
	t336 = -t342 * t360 - t357;
	t335 = -t342 * t359 + t358;
	t334 = -t341 * t367 + t342 * t365;
	t333 = -t341 * t339 - t340 * t365;
	t332 = -t341 * t366 - t343 * t367;
	t331 = -t345 * t341 * t340 + t343 * t339;
	t330 = -t341 * t355 + (t341 * t360 - t343 * t361) * t340;
	t329 = t350 * t340 + t341 * t356;
	t328 = t351 * t340 - t343 * t355;
	t327 = t343 * t356 + (-t341 * t362 + t343 * t359) * t340;
	t326 = -t354 * t349 * t343 + (t353 * t348 + t357) * t341;
	t325 = t368 * t341 + t343 * t352;
	t324 = t351 * t342 + t343 * t357 - t350;
	t323 = t341 * t352 - t368 * t343;
	t1 = [t326, t326, t328, t328, t323; -t324, -t324, t330, t330, -t325; 0, 0, t336, t336, -t340 * t359 - t356; t325, t325, t327, t327, t324; t323, t323, t329, t329, t326; 0, 0, t335, t335, t340 * t360 - t355; t333, t333, t332, t332, 0; t331, t331, t334, t334, 0; 0, 0, t339, t339, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end