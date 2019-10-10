% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PRRRR1
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
% Datum: 2019-10-09 20:54
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5PRRRR1_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR1_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRR1_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_jacobiRD_rot_sym_varpar: pkin has to be [2x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:54:42
	% EndTime: 2019-10-09 20:54:42
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:54:42
	% EndTime: 2019-10-09 20:54:42
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:54:42
	% EndTime: 2019-10-09 20:54:42
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t10 = qJD(2) * sin(qJ(2));
	t9 = qJD(2) * cos(qJ(2));
	t1 = [0, -t9, 0, 0, 0; 0, 0, 0, 0, 0; 0, -t10, 0, 0, 0; 0, t10, 0, 0, 0; 0, 0, 0, 0, 0; 0, -t9, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:54:42
	% EndTime: 2019-10-09 20:54:42
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->7), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->14)
	t81 = sin(qJ(2));
	t88 = qJD(2) * t81;
	t83 = cos(qJ(2));
	t87 = qJD(2) * t83;
	t80 = sin(qJ(3));
	t86 = qJD(3) * t80;
	t82 = cos(qJ(3));
	t85 = qJD(3) * t82;
	t84 = qJD(3) * t83;
	t79 = t81 * t86 - t82 * t87;
	t78 = t80 * t87 + t81 * t85;
	t77 = t80 * t84 + t82 * t88;
	t76 = t80 * t88 - t82 * t84;
	t1 = [0, t79, t76, 0, 0; 0, 0, t86, 0, 0; 0, -t77, -t78, 0, 0; 0, t78, t77, 0, 0; 0, 0, t85, 0, 0; 0, t76, t79, 0, 0; 0, -t88, 0, 0, 0; 0, 0, 0, 0, 0; 0, t87, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:54:42
	% EndTime: 2019-10-09 20:54:42
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (57->10), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t122 = qJD(3) + qJD(4);
	t124 = sin(qJ(2));
	t129 = t122 * t124;
	t125 = cos(qJ(2));
	t128 = t122 * t125;
	t127 = qJD(2) * t124;
	t126 = qJD(2) * t125;
	t123 = qJ(3) + qJ(4);
	t121 = cos(t123);
	t120 = sin(t123);
	t119 = t122 * t121;
	t118 = t122 * t120;
	t117 = t120 * t129 - t121 * t126;
	t116 = t120 * t126 + t121 * t129;
	t115 = t120 * t128 + t121 * t127;
	t114 = t120 * t127 - t121 * t128;
	t1 = [0, t117, t114, t114, 0; 0, 0, t118, t118, 0; 0, -t115, -t116, -t116, 0; 0, t116, t115, t115, 0; 0, 0, t119, t119, 0; 0, t114, t117, t117, 0; 0, -t127, 0, 0, 0; 0, 0, 0, 0, 0; 0, t126, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:54:43
	% EndTime: 2019-10-09 20:54:43
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (166->31), mult. (226->57), div. (0->0), fcn. (226->6), ass. (0->41)
	t311 = qJ(3) + qJ(4);
	t309 = cos(t311);
	t310 = qJD(3) + qJD(4);
	t335 = t310 * t309;
	t312 = sin(qJ(5));
	t334 = t310 * t312;
	t313 = sin(qJ(2));
	t333 = t310 * t313;
	t314 = cos(qJ(5));
	t332 = t310 * t314;
	t315 = cos(qJ(2));
	t331 = t310 * t315;
	t330 = t314 * t315;
	t329 = qJD(2) * t313;
	t328 = qJD(2) * t315;
	t327 = qJD(5) * t312;
	t326 = qJD(5) * t314;
	t325 = qJD(5) * t315;
	t308 = sin(t311);
	t324 = t308 * t332;
	t323 = t308 * t333;
	t322 = t309 * t333;
	t321 = t308 * t331;
	t320 = t309 * t331;
	t319 = qJD(5) * t309 - qJD(2);
	t318 = qJD(2) * t309 - qJD(5);
	t317 = t319 * t312;
	t316 = t318 * t313 + t321;
	t307 = t309 * t328 - t323;
	t306 = -t309 * t329 - t321;
	t305 = t309 * t327 + t324;
	t304 = -t308 * t334 + t309 * t326;
	t303 = -t314 * t322 + (t313 * t327 - t314 * t328) * t308;
	t302 = t312 * t322 + (t312 * t328 + t313 * t326) * t308;
	t301 = -t314 * t320 + (t312 * t325 + t314 * t329) * t308;
	t300 = t312 * t320 + (-t312 * t329 + t314 * t325) * t308;
	t299 = -t318 * t330 + (t317 + t324) * t313;
	t298 = t319 * t314 * t313 + (t318 * t315 - t323) * t312;
	t297 = t316 * t314 + t315 * t317;
	t296 = t316 * t312 - t319 * t330;
	t1 = [0, t299, t301, t301, t296; 0, 0, t305, t305, t308 * t326 + t309 * t334; 0, -t297, t303, t303, -t298; 0, t298, t300, t300, t297; 0, 0, t304, t304, -t308 * t327 + t309 * t332; 0, t296, t302, t302, t299; 0, -t308 * t328 - t322, t306, t306, 0; 0, 0, -t335, -t335, 0; 0, -t308 * t329 + t320, t307, t307, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end