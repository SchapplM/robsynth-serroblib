% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRRPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:33
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPRRPR7_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR7_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR7_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR7_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR7_jacobiRD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:33:56
	% EndTime: 2019-10-10 01:33:56
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:33:56
	% EndTime: 2019-10-10 01:33:56
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
	% StartTime: 2019-10-10 01:33:56
	% EndTime: 2019-10-10 01:33:56
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
	% StartTime: 2019-10-10 01:33:56
	% EndTime: 2019-10-10 01:33:56
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
	% StartTime: 2019-10-10 01:33:56
	% EndTime: 2019-10-10 01:33:56
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (60->13), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t71 = qJ(3) + qJ(4);
	t69 = cos(t71);
	t70 = qJD(3) + qJD(4);
	t78 = t70 * t69;
	t72 = sin(qJ(1));
	t77 = t70 * t72;
	t73 = cos(qJ(1));
	t76 = t70 * t73;
	t75 = qJD(1) * t72;
	t74 = qJD(1) * t73;
	t68 = sin(t71);
	t67 = t70 * t68;
	t66 = -t68 * t77 + t69 * t74;
	t65 = t68 * t74 + t69 * t77;
	t64 = t68 * t76 + t69 * t75;
	t63 = -t68 * t75 + t69 * t76;
	t1 = [t63, 0, t66, t66, 0, 0; t65, 0, t64, t64, 0, 0; 0, 0, -t78, -t78, 0, 0; -t64, 0, -t65, -t65, 0, 0; t66, 0, t63, t63, 0, 0; 0, 0, t67, t67, 0, 0; -t74, 0, 0, 0, 0, 0; -t75, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:33:56
	% EndTime: 2019-10-10 01:33:56
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (88->13), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t90 = qJ(3) + qJ(4) + pkin(10);
	t89 = cos(t90);
	t91 = qJD(3) + qJD(4);
	t98 = t91 * t89;
	t92 = sin(qJ(1));
	t97 = t91 * t92;
	t93 = cos(qJ(1));
	t96 = t91 * t93;
	t95 = qJD(1) * t92;
	t94 = qJD(1) * t93;
	t88 = sin(t90);
	t87 = t91 * t88;
	t86 = -t88 * t97 + t89 * t94;
	t85 = t88 * t94 + t89 * t97;
	t84 = t88 * t96 + t89 * t95;
	t83 = -t88 * t95 + t89 * t96;
	t1 = [t83, 0, t86, t86, 0, 0; t85, 0, t84, t84, 0, 0; 0, 0, -t98, -t98, 0, 0; -t84, 0, -t85, -t85, 0, 0; t86, 0, t83, t83, 0, 0; 0, 0, t87, t87, 0, 0; -t94, 0, 0, 0, 0, 0; -t95, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:33:57
	% EndTime: 2019-10-10 01:33:58
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (240->31), mult. (226->57), div. (0->0), fcn. (226->6), ass. (0->41)
	t345 = cos(qJ(1));
	t340 = qJ(3) + qJ(4) + pkin(10);
	t338 = sin(t340);
	t347 = qJD(1) * t338 + qJD(6);
	t365 = t347 * t345;
	t343 = sin(qJ(1));
	t339 = cos(t340);
	t341 = qJD(3) + qJD(4);
	t359 = t341 * t345;
	t349 = t339 * t359;
	t364 = t347 * t343 - t349;
	t363 = t341 * t338;
	t342 = sin(qJ(6));
	t362 = t341 * t342;
	t361 = t341 * t343;
	t344 = cos(qJ(6));
	t360 = t341 * t344;
	t358 = qJD(1) * t343;
	t357 = qJD(1) * t345;
	t356 = qJD(6) * t342;
	t355 = qJD(6) * t344;
	t354 = qJD(6) * t345;
	t353 = t339 * t360;
	t352 = t338 * t361;
	t351 = t339 * t361;
	t350 = t338 * t359;
	t348 = -qJD(6) * t338 - qJD(1);
	t346 = t348 * t345;
	t337 = t338 * t357 + t351;
	t336 = t338 * t358 - t349;
	t335 = t338 * t356 - t353;
	t334 = t338 * t355 + t339 * t362;
	t333 = -t344 * t352 + (-t343 * t356 + t344 * t357) * t339;
	t332 = t342 * t352 + (-t342 * t357 - t343 * t355) * t339;
	t331 = t344 * t350 + (t342 * t354 + t344 * t358) * t339;
	t330 = -t342 * t350 + (-t342 * t358 + t344 * t354) * t339;
	t329 = t344 * t365 + (t348 * t342 + t353) * t343;
	t328 = t348 * t344 * t343 + (-t351 - t365) * t342;
	t327 = t342 * t346 - t364 * t344;
	t326 = t364 * t342 + t344 * t346;
	t1 = [t327, 0, t333, t333, 0, t328; t329, 0, t331, t331, 0, -t326; 0, 0, t335, t335, 0, t338 * t362 - t339 * t355; t326, 0, t332, t332, 0, -t329; t328, 0, t330, t330, 0, t327; 0, 0, t334, t334, 0, t338 * t360 + t339 * t356; t339 * t358 + t350, 0, t337, t337, 0, 0; -t339 * t357 + t352, 0, t336, t336, 0, 0; 0, 0, -t363, -t363, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end