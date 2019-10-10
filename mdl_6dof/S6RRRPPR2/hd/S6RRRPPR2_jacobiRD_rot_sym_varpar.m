% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRPPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:18
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRRPPR2_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR2_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR2_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPPR2_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR2_jacobiRD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:18:48
	% EndTime: 2019-10-10 11:18:48
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:18:48
	% EndTime: 2019-10-10 11:18:48
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
	% StartTime: 2019-10-10 11:18:48
	% EndTime: 2019-10-10 11:18:48
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
	% StartTime: 2019-10-10 11:18:48
	% EndTime: 2019-10-10 11:18:48
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
	% StartTime: 2019-10-10 11:18:48
	% EndTime: 2019-10-10 11:18:48
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
	% StartTime: 2019-10-10 11:18:49
	% EndTime: 2019-10-10 11:18:49
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (85->10), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t251 = qJD(2) + qJD(3);
	t252 = sin(qJ(1));
	t257 = t251 * t252;
	t253 = cos(qJ(1));
	t256 = t251 * t253;
	t255 = qJD(1) * t252;
	t254 = qJD(1) * t253;
	t250 = qJ(2) + qJ(3) + pkin(10);
	t249 = cos(t250);
	t248 = sin(t250);
	t247 = t251 * t249;
	t246 = t251 * t248;
	t245 = -t248 * t257 + t249 * t254;
	t244 = t248 * t254 + t249 * t257;
	t243 = t248 * t256 + t249 * t255;
	t242 = -t248 * t255 + t249 * t256;
	t1 = [-t255, 0, 0, 0, 0, 0; t254, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t245, t242, t242, 0, 0, 0; t243, t244, t244, 0, 0, 0; 0, t246, t246, 0, 0, 0; -t244, -t243, -t243, 0, 0, 0; t242, t245, t245, 0, 0, 0; 0, t247, t247, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:18:50
	% EndTime: 2019-10-10 11:18:50
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (240->31), mult. (226->56), div. (0->0), fcn. (226->6), ass. (0->42)
	t349 = cos(qJ(6));
	t345 = qJ(2) + qJ(3) + pkin(10);
	t343 = sin(t345);
	t354 = qJD(6) * t343 + qJD(1);
	t371 = t349 * t354;
	t347 = sin(qJ(6));
	t370 = t354 * t347;
	t346 = qJD(2) + qJD(3);
	t369 = t346 * t343;
	t368 = t346 * t347;
	t348 = sin(qJ(1));
	t367 = t346 * t348;
	t366 = t346 * t349;
	t350 = cos(qJ(1));
	t365 = t346 * t350;
	t364 = qJD(1) * t348;
	t363 = qJD(1) * t350;
	t362 = qJD(6) * t347;
	t361 = qJD(6) * t349;
	t360 = qJD(6) * t350;
	t344 = cos(t345);
	t359 = t344 * t366;
	t358 = t343 * t367;
	t357 = t344 * t367;
	t356 = t343 * t365;
	t355 = t344 * t365;
	t353 = -qJD(1) * t343 - qJD(6);
	t352 = t353 * t350;
	t351 = t353 * t348 + t355;
	t342 = -t343 * t363 - t357;
	t341 = t343 * t364 - t355;
	t340 = -t343 * t362 + t359;
	t339 = t343 * t361 + t344 * t368;
	t338 = -t349 * t358 + (-t348 * t362 + t349 * t363) * t344;
	t337 = -t347 * t358 + (t347 * t363 + t348 * t361) * t344;
	t336 = -t349 * t356 + (-t347 * t360 - t349 * t364) * t344;
	t335 = -t347 * t356 + (-t347 * t364 + t349 * t360) * t344;
	t334 = t351 * t347 + t350 * t371;
	t333 = t351 * t349 - t350 * t370;
	t332 = -t348 * t371 + (t352 - t357) * t347;
	t331 = t349 * t352 + (-t359 + t370) * t348;
	t1 = [t332, t335, t335, 0, 0, t333; t334, t337, t337, 0, 0, -t331; 0, t339, t339, 0, 0, t343 * t366 + t344 * t362; t331, t336, t336, 0, 0, -t334; t333, t338, t338, 0, 0, t332; 0, t340, t340, 0, 0, -t343 * t368 + t344 * t361; -t344 * t363 + t358, t341, t341, 0, 0, 0; -t344 * t364 - t356, t342, t342, 0, 0, 0; 0, -t369, -t369, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end