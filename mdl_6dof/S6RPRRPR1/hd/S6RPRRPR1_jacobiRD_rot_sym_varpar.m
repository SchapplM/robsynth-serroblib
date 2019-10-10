% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:23
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPRRPR1_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR1_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR1_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR1_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR1_jacobiRD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:23:31
	% EndTime: 2019-10-10 01:23:31
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:23:31
	% EndTime: 2019-10-10 01:23:31
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
	% StartTime: 2019-10-10 01:23:31
	% EndTime: 2019-10-10 01:23:31
	% DurationCPUTime: 0.06s
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
	% StartTime: 2019-10-10 01:23:31
	% EndTime: 2019-10-10 01:23:31
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->10), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t41 = sin(qJ(3));
	t46 = qJD(1) * t41;
	t42 = cos(qJ(3));
	t45 = qJD(1) * t42;
	t44 = qJD(3) * t41;
	t43 = qJD(3) * t42;
	t40 = qJ(1) + pkin(10);
	t39 = cos(t40);
	t38 = sin(t40);
	t37 = t38 * t44 - t39 * t45;
	t36 = t38 * t43 + t39 * t46;
	t35 = t38 * t45 + t39 * t44;
	t34 = t38 * t46 - t39 * t43;
	t1 = [t37, 0, t34, 0, 0, 0; -t35, 0, -t36, 0, 0, 0; 0, 0, -t44, 0, 0, 0; t36, 0, t35, 0, 0, 0; t34, 0, t37, 0, 0, 0; 0, 0, -t43, 0, 0, 0; -qJD(1) * t38, 0, 0, 0, 0, 0; qJD(1) * t39, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:23:31
	% EndTime: 2019-10-10 01:23:32
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (87->15), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->16)
	t77 = qJ(3) + qJ(4);
	t73 = sin(t77);
	t75 = qJD(3) + qJD(4);
	t81 = t75 * t73;
	t74 = cos(t77);
	t80 = t75 * t74;
	t79 = qJD(1) * t73;
	t78 = qJD(1) * t74;
	t76 = qJ(1) + pkin(10);
	t72 = cos(t76);
	t71 = sin(t76);
	t70 = t71 * t81 - t72 * t78;
	t69 = t71 * t80 + t72 * t79;
	t68 = t71 * t78 + t72 * t81;
	t67 = t71 * t79 - t72 * t80;
	t1 = [t70, 0, t67, t67, 0, 0; -t68, 0, -t69, -t69, 0, 0; 0, 0, -t81, -t81, 0, 0; t69, 0, t68, t68, 0, 0; t67, 0, t70, t70, 0, 0; 0, 0, -t80, -t80, 0, 0; -qJD(1) * t71, 0, 0, 0, 0, 0; qJD(1) * t72, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:23:32
	% EndTime: 2019-10-10 01:23:32
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (115->15), mult. (54->12), div. (0->0), fcn. (54->4), ass. (0->16)
	t92 = qJ(3) + qJ(4) + pkin(11);
	t88 = sin(t92);
	t93 = qJD(3) + qJD(4);
	t98 = t93 * t88;
	t89 = cos(t92);
	t97 = t93 * t89;
	t94 = qJ(1) + pkin(10);
	t90 = sin(t94);
	t96 = qJD(1) * t90;
	t91 = cos(t94);
	t95 = qJD(1) * t91;
	t87 = -t89 * t95 + t90 * t98;
	t86 = t88 * t95 + t90 * t97;
	t85 = t89 * t96 + t91 * t98;
	t84 = t88 * t96 - t91 * t97;
	t1 = [t87, 0, t84, t84, 0, 0; -t85, 0, -t86, -t86, 0, 0; 0, 0, -t98, -t98, 0, 0; t86, 0, t85, t85, 0, 0; t84, 0, t87, t87, 0, 0; 0, 0, -t97, -t97, 0, 0; -t96, 0, 0, 0, 0, 0; t95, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:23:33
	% EndTime: 2019-10-10 01:23:33
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (314->29), mult. (226->56), div. (0->0), fcn. (226->6), ass. (0->42)
	t358 = cos(qJ(6));
	t354 = qJ(3) + qJ(4) + pkin(11);
	t351 = cos(t354);
	t361 = qJD(1) * t351 - qJD(6);
	t377 = t358 * t361;
	t362 = qJD(6) * t351 - qJD(1);
	t350 = sin(t354);
	t355 = qJD(3) + qJD(4);
	t357 = sin(qJ(6));
	t374 = t355 * t357;
	t366 = t350 * t374;
	t376 = t362 * t358 - t366;
	t375 = t350 * t355;
	t349 = t355 * t351;
	t373 = t355 * t358;
	t356 = qJ(1) + pkin(10);
	t352 = sin(t356);
	t372 = qJD(1) * t352;
	t353 = cos(t356);
	t371 = qJD(1) * t353;
	t370 = qJD(1) * t357;
	t369 = qJD(1) * t358;
	t368 = qJD(6) * t357;
	t367 = qJD(6) * t358;
	t365 = t350 * t373;
	t364 = t351 * t374;
	t363 = t351 * t373;
	t360 = t361 * t357;
	t359 = t362 * t357 + t365;
	t348 = -t351 * t368 - t365;
	t347 = -t351 * t367 + t366;
	t346 = t351 * t371 - t352 * t375;
	t345 = -t351 * t372 - t353 * t375;
	t344 = -t352 * t363 + (t352 * t368 - t353 * t369) * t350;
	t343 = t352 * t364 + (t352 * t367 + t353 * t370) * t350;
	t342 = -t353 * t363 + (t352 * t369 + t353 * t368) * t350;
	t341 = t353 * t364 + (-t352 * t370 + t353 * t367) * t350;
	t340 = t359 * t352 - t353 * t377;
	t339 = t376 * t352 + t353 * t360;
	t338 = t352 * t377 + t359 * t353;
	t337 = t352 * t360 - t376 * t353;
	t1 = [t340, 0, t342, t342, 0, t337; -t338, 0, t344, t344, 0, -t339; 0, 0, t348, t348, 0, -t350 * t367 - t364; t339, 0, t341, t341, 0, t338; t337, 0, t343, t343, 0, t340; 0, 0, t347, t347, 0, t350 * t368 - t363; -t352 * t349 - t350 * t371, 0, t345, t345, 0, 0; t353 * t349 - t350 * t372, 0, t346, t346, 0, 0; 0, 0, t349, t349, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end