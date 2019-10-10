% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:35
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRRRPR1_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR1_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR1_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPR1_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR1_jacobiRD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:35:07
	% EndTime: 2019-10-10 12:35:07
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:35:07
	% EndTime: 2019-10-10 12:35:07
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
	% StartTime: 2019-10-10 12:35:07
	% EndTime: 2019-10-10 12:35:07
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
	% StartTime: 2019-10-10 12:35:07
	% EndTime: 2019-10-10 12:35:07
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
	% StartTime: 2019-10-10 12:35:07
	% EndTime: 2019-10-10 12:35:07
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (143->17), mult. (72->14), div. (0->0), fcn. (72->4), ass. (0->17)
	t90 = qJ(2) + qJ(3) + qJ(4);
	t87 = sin(t90);
	t89 = qJD(2) + qJD(3) + qJD(4);
	t98 = t89 * t87;
	t88 = cos(t90);
	t97 = t89 * t88;
	t91 = sin(qJ(1));
	t96 = t89 * t91;
	t92 = cos(qJ(1));
	t95 = t89 * t92;
	t94 = qJD(1) * t91;
	t93 = qJD(1) * t92;
	t86 = t87 * t96 - t88 * t93;
	t85 = t87 * t93 + t88 * t96;
	t84 = t87 * t95 + t88 * t94;
	t83 = t87 * t94 - t88 * t95;
	t1 = [t86, t83, t83, t83, 0, 0; -t84, -t85, -t85, -t85, 0, 0; 0, -t98, -t98, -t98, 0, 0; t85, t84, t84, t84, 0, 0; t83, t86, t86, t86, 0, 0; 0, -t97, -t97, -t97, 0, 0; -t94, 0, 0, 0, 0, 0; t93, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:35:07
	% EndTime: 2019-10-10 12:35:07
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (181->17), mult. (72->14), div. (0->0), fcn. (72->4), ass. (0->17)
	t102 = qJ(2) + qJ(3) + qJ(4) + pkin(11);
	t100 = sin(t102);
	t103 = qJD(2) + qJD(3) + qJD(4);
	t111 = t103 * t100;
	t101 = cos(t102);
	t110 = t103 * t101;
	t104 = sin(qJ(1));
	t109 = t103 * t104;
	t105 = cos(qJ(1));
	t108 = t103 * t105;
	t107 = qJD(1) * t104;
	t106 = qJD(1) * t105;
	t99 = t100 * t109 - t101 * t106;
	t98 = t100 * t106 + t101 * t109;
	t97 = t100 * t108 + t101 * t107;
	t96 = t100 * t107 - t101 * t108;
	t1 = [t99, t96, t96, t96, 0, 0; -t97, -t98, -t98, -t98, 0, 0; 0, -t111, -t111, -t111, 0, 0; t98, t97, t97, t97, 0, 0; t96, t99, t99, t99, 0, 0; 0, -t110, -t110, -t110, 0, 0; -t107, 0, 0, 0, 0, 0; t106, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:35:09
	% EndTime: 2019-10-10 12:35:09
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (435->29), mult. (279->57), div. (0->0), fcn. (279->6), ass. (0->41)
	t359 = qJD(2) + qJD(3) + qJD(4);
	t360 = sin(qJ(6));
	t382 = t359 * t360;
	t361 = sin(qJ(1));
	t381 = t359 * t361;
	t362 = cos(qJ(6));
	t380 = t359 * t362;
	t363 = cos(qJ(1));
	t379 = t359 * t363;
	t378 = t362 * t363;
	t377 = qJD(1) * t361;
	t376 = qJD(1) * t363;
	t375 = qJD(6) * t360;
	t374 = qJD(6) * t362;
	t373 = qJD(6) * t363;
	t358 = qJ(2) + qJ(3) + qJ(4) + pkin(11);
	t356 = sin(t358);
	t372 = t356 * t380;
	t371 = t356 * t381;
	t357 = cos(t358);
	t370 = t357 * t381;
	t369 = t356 * t379;
	t368 = t357 * t379;
	t367 = qJD(6) * t357 - qJD(1);
	t366 = qJD(1) * t357 - qJD(6);
	t365 = t367 * t360;
	t364 = t366 * t361 + t369;
	t355 = t359 * t357;
	t354 = t357 * t376 - t371;
	t353 = -t357 * t377 - t369;
	t352 = -t357 * t375 - t372;
	t351 = t356 * t382 - t357 * t374;
	t350 = -t362 * t370 + (t361 * t375 - t362 * t376) * t356;
	t349 = t360 * t370 + (t360 * t376 + t361 * t374) * t356;
	t348 = -t362 * t368 + (t360 * t373 + t362 * t377) * t356;
	t347 = t360 * t368 + (-t360 * t377 + t362 * t373) * t356;
	t346 = -t366 * t378 + (t365 + t372) * t361;
	t345 = t367 * t362 * t361 + (t366 * t363 - t371) * t360;
	t344 = t364 * t362 + t363 * t365;
	t343 = t364 * t360 - t367 * t378;
	t1 = [t346, t348, t348, t348, 0, t343; -t344, t350, t350, t350, 0, -t345; 0, t352, t352, t352, 0, -t356 * t374 - t357 * t382; t345, t347, t347, t347, 0, t344; t343, t349, t349, t349, 0, t346; 0, t351, t351, t351, 0, t356 * t375 - t357 * t380; -t356 * t376 - t370, t353, t353, t353, 0, 0; -t356 * t377 + t368, t354, t354, t354, 0, 0; 0, t355, t355, t355, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end