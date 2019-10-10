% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRPPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:20
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRRPPR3_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR3_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR3_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPPR3_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPPR3_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:20:30
	% EndTime: 2019-10-10 11:20:30
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:20:30
	% EndTime: 2019-10-10 11:20:30
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
	% StartTime: 2019-10-10 11:20:30
	% EndTime: 2019-10-10 11:20:30
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
	% StartTime: 2019-10-10 11:20:30
	% EndTime: 2019-10-10 11:20:30
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
	% StartTime: 2019-10-10 11:20:31
	% EndTime: 2019-10-10 11:20:31
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (59->11), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t226 = qJ(2) + qJ(3);
	t223 = sin(t226);
	t225 = qJD(2) + qJD(3);
	t233 = t225 * t223;
	t227 = sin(qJ(1));
	t232 = t225 * t227;
	t228 = cos(qJ(1));
	t231 = t225 * t228;
	t230 = qJD(1) * t227;
	t229 = qJD(1) * t228;
	t224 = cos(t226);
	t222 = t225 * t224;
	t221 = -t223 * t232 + t224 * t229;
	t220 = -t223 * t229 - t224 * t232;
	t219 = -t223 * t231 - t224 * t230;
	t218 = t223 * t230 - t224 * t231;
	t1 = [-t221, t218, t218, 0, 0, 0; t219, t220, t220, 0, 0, 0; 0, -t233, -t233, 0, 0, 0; -t230, 0, 0, 0, 0, 0; t229, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t220, t219, t219, 0, 0, 0; -t218, t221, t221, 0, 0, 0; 0, t222, t222, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:20:30
	% EndTime: 2019-10-10 11:20:31
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (57->10), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t120 = qJD(2) + qJD(3);
	t122 = sin(qJ(1));
	t127 = t120 * t122;
	t123 = cos(qJ(1));
	t126 = t120 * t123;
	t125 = qJD(1) * t122;
	t124 = qJD(1) * t123;
	t121 = qJ(2) + qJ(3);
	t119 = cos(t121);
	t118 = sin(t121);
	t117 = t120 * t119;
	t116 = t120 * t118;
	t115 = -t118 * t127 + t119 * t124;
	t114 = t118 * t124 + t119 * t127;
	t113 = t118 * t126 + t119 * t125;
	t112 = -t118 * t125 + t119 * t126;
	t1 = [-t114, -t113, -t113, 0, 0, 0; t112, t115, t115, 0, 0, 0; 0, t117, t117, 0, 0, 0; t115, t112, t112, 0, 0, 0; t113, t114, t114, 0, 0, 0; 0, t116, t116, 0, 0, 0; t125, 0, 0, 0, 0, 0; -t124, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:20:32
	% EndTime: 2019-10-10 11:20:32
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (166->31), mult. (226->57), div. (0->0), fcn. (226->6), ass. (0->41)
	t336 = qJ(2) + qJ(3);
	t333 = sin(t336);
	t335 = qJD(2) + qJD(3);
	t360 = t335 * t333;
	t337 = sin(qJ(6));
	t359 = t335 * t337;
	t338 = sin(qJ(1));
	t358 = t335 * t338;
	t339 = cos(qJ(6));
	t357 = t335 * t339;
	t340 = cos(qJ(1));
	t356 = t335 * t340;
	t355 = t339 * t340;
	t354 = qJD(1) * t338;
	t353 = qJD(1) * t340;
	t352 = qJD(6) * t337;
	t351 = qJD(6) * t339;
	t350 = qJD(6) * t340;
	t334 = cos(t336);
	t349 = t334 * t357;
	t348 = t333 * t358;
	t347 = t334 * t358;
	t346 = t333 * t356;
	t345 = t334 * t356;
	t344 = qJD(6) * t333 + qJD(1);
	t343 = qJD(1) * t333 + qJD(6);
	t342 = t344 * t337;
	t341 = t343 * t338 - t345;
	t332 = -t333 * t353 - t347;
	t331 = t333 * t354 - t345;
	t330 = -t333 * t352 + t349;
	t329 = -t333 * t351 - t334 * t359;
	t328 = -t339 * t348 + (-t338 * t352 + t339 * t353) * t334;
	t327 = t337 * t348 + (-t337 * t353 - t338 * t351) * t334;
	t326 = -t339 * t346 + (-t337 * t350 - t339 * t354) * t334;
	t325 = t337 * t346 + (t337 * t354 - t339 * t350) * t334;
	t324 = -t343 * t355 + (t342 - t349) * t338;
	t323 = t344 * t339 * t338 + (t343 * t340 + t347) * t337;
	t322 = t341 * t339 + t340 * t342;
	t321 = t341 * t337 - t344 * t355;
	t1 = [t324, t326, t326, 0, 0, t321; -t322, t328, t328, 0, 0, -t323; 0, t330, t330, 0, 0, -t333 * t359 + t334 * t351; t323, t325, t325, 0, 0, t322; t321, t327, t327, 0, 0, t324; 0, t329, t329, 0, 0, -t333 * t357 - t334 * t352; -t334 * t353 + t348, t331, t331, 0, 0, 0; -t334 * t354 - t346, t332, t332, 0, 0, 0; 0, -t360, -t360, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end