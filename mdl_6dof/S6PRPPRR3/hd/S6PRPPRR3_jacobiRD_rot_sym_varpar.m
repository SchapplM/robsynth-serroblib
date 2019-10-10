% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:28
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6PRPPRR3_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR3_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPPRR3_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_jacobiRD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:28:08
	% EndTime: 2019-10-09 21:28:08
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:28:08
	% EndTime: 2019-10-09 21:28:08
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:28:08
	% EndTime: 2019-10-09 21:28:08
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (6->6), mult. (24->17), div. (0->0), fcn. (24->6), ass. (0->9)
	t61 = cos(pkin(6));
	t62 = sin(qJ(2));
	t66 = t61 * t62;
	t63 = cos(qJ(2));
	t65 = t61 * t63;
	t64 = qJD(2) * sin(pkin(6));
	t60 = cos(pkin(10));
	t58 = sin(pkin(10));
	t1 = [0, (t58 * t66 - t60 * t63) * qJD(2), 0, 0, 0, 0; 0, (-t58 * t63 - t60 * t66) * qJD(2), 0, 0, 0, 0; 0, -t62 * t64, 0, 0, 0, 0; 0, (t58 * t65 + t60 * t62) * qJD(2), 0, 0, 0, 0; 0, (t58 * t62 - t60 * t65) * qJD(2), 0, 0, 0, 0; 0, -t63 * t64, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:28:09
	% EndTime: 2019-10-09 21:28:09
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (5->5), mult. (24->17), div. (0->0), fcn. (24->6), ass. (0->9)
	t150 = cos(pkin(6));
	t151 = sin(qJ(2));
	t155 = t150 * t151;
	t152 = cos(qJ(2));
	t154 = t150 * t152;
	t153 = qJD(2) * sin(pkin(6));
	t149 = cos(pkin(10));
	t147 = sin(pkin(10));
	t1 = [0, (t147 * t155 - t149 * t152) * qJD(2), 0, 0, 0, 0; 0, (-t147 * t152 - t149 * t155) * qJD(2), 0, 0, 0, 0; 0, -t151 * t153, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, (-t147 * t154 - t149 * t151) * qJD(2), 0, 0, 0, 0; 0, (-t147 * t151 + t149 * t154) * qJD(2), 0, 0, 0, 0; 0, t152 * t153, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:28:08
	% EndTime: 2019-10-09 21:28:08
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (14->10), mult. (60->29), div. (0->0), fcn. (60->8), ass. (0->15)
	t122 = cos(pkin(6));
	t123 = sin(qJ(2));
	t127 = t122 * t123;
	t124 = cos(qJ(2));
	t126 = t122 * t124;
	t125 = qJD(2) * sin(pkin(6));
	t121 = cos(pkin(10));
	t120 = cos(pkin(11));
	t118 = sin(pkin(10));
	t117 = sin(pkin(11));
	t116 = (t118 * t127 - t121 * t124) * qJD(2);
	t115 = (-t118 * t126 - t121 * t123) * qJD(2);
	t114 = (-t118 * t124 - t121 * t127) * qJD(2);
	t113 = (-t118 * t123 + t121 * t126) * qJD(2);
	t1 = [0, t115 * t117 + t116 * t120, 0, 0, 0, 0; 0, t113 * t117 + t114 * t120, 0, 0, 0, 0; 0, (t117 * t124 - t120 * t123) * t125, 0, 0, 0, 0; 0, t115 * t120 - t116 * t117, 0, 0, 0, 0; 0, t113 * t120 - t114 * t117, 0, 0, 0, 0; 0, (t117 * t123 + t120 * t124) * t125, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:28:10
	% EndTime: 2019-10-09 21:28:10
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (81->36), mult. (282->83), div. (0->0), fcn. (310->10), ass. (0->42)
	t301 = sin(pkin(11));
	t308 = sin(qJ(2));
	t322 = t301 * t308;
	t303 = sin(pkin(6));
	t307 = sin(qJ(5));
	t321 = t303 * t307;
	t309 = cos(qJ(5));
	t320 = t303 * t309;
	t306 = cos(pkin(6));
	t319 = t306 * t308;
	t310 = cos(qJ(2));
	t318 = t306 * t310;
	t317 = qJD(2) * t303;
	t316 = qJD(2) * t310;
	t315 = qJD(5) * t307;
	t314 = qJD(5) * t309;
	t302 = sin(pkin(10));
	t313 = t302 * t319;
	t305 = cos(pkin(10));
	t311 = -t302 * t308 + t305 * t318;
	t290 = t311 * qJD(2);
	t295 = t302 * t310 + t305 * t319;
	t291 = t295 * qJD(2);
	t304 = cos(pkin(11));
	t278 = t290 * t304 + t291 * t301;
	t296 = t302 * t318 + t305 * t308;
	t292 = t296 * qJD(2);
	t293 = -qJD(2) * t313 + t305 * t316;
	t280 = -t292 * t304 + t293 * t301;
	t312 = t301 * t310 - t304 * t308;
	t297 = t305 * t310 - t313;
	t289 = t312 * t303;
	t288 = (t304 * t310 + t322) * t303;
	t287 = t312 * t317;
	t286 = -t303 * t304 * t316 - t317 * t322;
	t285 = t296 * t301 + t297 * t304;
	t284 = -t296 * t304 + t297 * t301;
	t283 = t295 * t304 - t301 * t311;
	t282 = t295 * t301 + t304 * t311;
	t281 = -t292 * t301 - t293 * t304;
	t279 = t290 * t301 - t291 * t304;
	t1 = [0, t281 * t309 - t284 * t315, 0, 0, -t280 * t307 + (-t285 * t309 + t302 * t321) * qJD(5), 0; 0, t279 * t309 - t282 * t315, 0, 0, -t278 * t307 + (-t283 * t309 - t305 * t321) * qJD(5), 0; 0, t287 * t309 - t288 * t315, 0, 0, t286 * t307 + (t289 * t309 + t306 * t307) * qJD(5), 0; 0, -t281 * t307 - t284 * t314, 0, 0, -t280 * t309 + (t285 * t307 + t302 * t320) * qJD(5), 0; 0, -t279 * t307 - t282 * t314, 0, 0, -t278 * t309 + (t283 * t307 - t305 * t320) * qJD(5), 0; 0, -t287 * t307 - t288 * t314, 0, 0, t286 * t309 + (-t289 * t307 + t306 * t309) * qJD(5), 0; 0, -t280, 0, 0, 0, 0; 0, -t278, 0, 0, 0, 0; 0, t286, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:28:11
	% EndTime: 2019-10-09 21:28:12
	% DurationCPUTime: 0.42s
	% Computational Cost: add. (303->70), mult. (966->146), div. (0->0), fcn. (1104->12), ass. (0->65)
	t478 = sin(pkin(11));
	t481 = cos(pkin(11));
	t486 = sin(qJ(2));
	t489 = cos(qJ(2));
	t517 = -t478 * t489 + t481 * t486;
	t516 = t478 * t486;
	t480 = sin(pkin(6));
	t485 = sin(qJ(5));
	t514 = t480 * t485;
	t488 = cos(qJ(5));
	t513 = t480 * t488;
	t483 = cos(pkin(6));
	t511 = t483 * t486;
	t510 = t483 * t489;
	t509 = qJD(2) * t480;
	t508 = qJD(2) * t489;
	t507 = qJD(5) * t485;
	t506 = qJD(5) * t488;
	t484 = sin(qJ(6));
	t505 = qJD(6) * t484;
	t487 = cos(qJ(6));
	t504 = qJD(6) * t487;
	t503 = qJD(6) * t488;
	t479 = sin(pkin(10));
	t502 = t479 * t511;
	t482 = cos(pkin(10));
	t493 = -t479 * t486 + t482 * t510;
	t465 = t493 * qJD(2);
	t470 = t479 * t489 + t482 * t511;
	t466 = t470 * qJD(2);
	t501 = t465 * t478 - t466 * t481;
	t471 = t479 * t510 + t482 * t486;
	t467 = t471 * qJD(2);
	t468 = -qJD(2) * t502 + t482 * t508;
	t500 = -t467 * t478 - t468 * t481;
	t499 = t470 * t478 + t481 * t493;
	t472 = t482 * t489 - t502;
	t498 = -t471 * t481 + t472 * t478;
	t440 = t465 * t481 + t466 * t478;
	t497 = -t499 * t503 - t440;
	t444 = -t467 * t481 + t468 * t478;
	t496 = -t498 * t503 - t444;
	t459 = -t480 * t481 * t508 - t509 * t516;
	t463 = (t481 * t489 + t516) * t480;
	t495 = -t463 * t503 + t459;
	t464 = t517 * t480;
	t456 = t464 * t488 - t483 * t485;
	t455 = -t464 * t485 - t483 * t488;
	t450 = t470 * t481 - t478 * t493;
	t454 = t471 * t478 + t472 * t481;
	t435 = -t450 * t485 + t482 * t513;
	t436 = t450 * t488 + t482 * t514;
	t437 = -t454 * t485 - t479 * t513;
	t494 = -t454 * t488 + t479 * t514;
	t492 = qJD(6) * t450 - t488 * t501 + t499 * t507;
	t491 = qJD(6) * t454 - t488 * t500 + t498 * t507;
	t460 = t517 * t509;
	t490 = qJD(6) * t464 + t460 * t488 + t463 * t507;
	t434 = t455 * qJD(5) - t459 * t488;
	t433 = -t456 * qJD(5) + t459 * t485;
	t432 = t437 * qJD(5) + t444 * t488;
	t431 = t494 * qJD(5) - t444 * t485;
	t430 = t435 * qJD(5) + t440 * t488;
	t429 = -t436 * qJD(5) - t440 * t485;
	t1 = [0, t496 * t484 - t491 * t487, 0, 0, t431 * t487 - t437 * t505, -t432 * t484 + t500 * t487 + (-t484 * t498 + t487 * t494) * qJD(6); 0, t497 * t484 - t492 * t487, 0, 0, t429 * t487 - t435 * t505, -t430 * t484 + t501 * t487 + (-t436 * t487 - t484 * t499) * qJD(6); 0, t495 * t484 - t490 * t487, 0, 0, t433 * t487 - t455 * t505, -t434 * t484 - t460 * t487 + (-t456 * t487 - t463 * t484) * qJD(6); 0, t491 * t484 + t496 * t487, 0, 0, -t431 * t484 - t437 * t504, -t432 * t487 - t500 * t484 + (-t484 * t494 - t487 * t498) * qJD(6); 0, t492 * t484 + t497 * t487, 0, 0, -t429 * t484 - t435 * t504, -t430 * t487 - t501 * t484 + (t436 * t484 - t487 * t499) * qJD(6); 0, t490 * t484 + t495 * t487, 0, 0, -t433 * t484 - t455 * t504, -t434 * t487 + t460 * t484 + (t456 * t484 - t463 * t487) * qJD(6); 0, t485 * t500 + t498 * t506, 0, 0, t432, 0; 0, t485 * t501 + t499 * t506, 0, 0, t430, 0; 0, -t460 * t485 + t463 * t506, 0, 0, t434, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end