% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRPRRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:01
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6PRPRRR5_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR5_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR5_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRR5_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR5_jacobiRD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:01:16
	% EndTime: 2019-10-09 22:01:16
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:01:16
	% EndTime: 2019-10-09 22:01:16
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:01:16
	% EndTime: 2019-10-09 22:01:16
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (6->6), mult. (24->17), div. (0->0), fcn. (24->6), ass. (0->9)
	t61 = cos(pkin(6));
	t62 = sin(qJ(2));
	t66 = t61 * t62;
	t63 = cos(qJ(2));
	t65 = t61 * t63;
	t64 = qJD(2) * sin(pkin(6));
	t60 = cos(pkin(11));
	t58 = sin(pkin(11));
	t1 = [0, (t58 * t66 - t60 * t63) * qJD(2), 0, 0, 0, 0; 0, (-t58 * t63 - t60 * t66) * qJD(2), 0, 0, 0, 0; 0, -t62 * t64, 0, 0, 0, 0; 0, (t58 * t65 + t60 * t62) * qJD(2), 0, 0, 0, 0; 0, (t58 * t62 - t60 * t65) * qJD(2), 0, 0, 0, 0; 0, -t63 * t64, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:01:17
	% EndTime: 2019-10-09 22:01:17
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (4->4), mult. (24->17), div. (0->0), fcn. (24->6), ass. (0->9)
	t139 = cos(pkin(6));
	t140 = sin(qJ(2));
	t144 = t139 * t140;
	t141 = cos(qJ(2));
	t143 = t139 * t141;
	t142 = qJD(2) * sin(pkin(6));
	t138 = cos(pkin(11));
	t136 = sin(pkin(11));
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, (-t136 * t144 + t138 * t141) * qJD(2), 0, 0, 0, 0; 0, (t136 * t141 + t138 * t144) * qJD(2), 0, 0, 0, 0; 0, t140 * t142, 0, 0, 0, 0; 0, (-t136 * t143 - t138 * t140) * qJD(2), 0, 0, 0, 0; 0, (-t136 * t140 + t138 * t143) * qJD(2), 0, 0, 0, 0; 0, t141 * t142, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:01:17
	% EndTime: 2019-10-09 22:01:17
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (37->25), mult. (140->59), div. (0->0), fcn. (148->8), ass. (0->26)
	t209 = sin(pkin(6));
	t212 = sin(qJ(4));
	t225 = t209 * t212;
	t214 = cos(qJ(4));
	t224 = t209 * t214;
	t215 = cos(qJ(2));
	t223 = t209 * t215;
	t211 = cos(pkin(6));
	t213 = sin(qJ(2));
	t222 = t211 * t213;
	t221 = t211 * t215;
	t220 = qJD(2) * t215;
	t219 = qJD(4) * t212;
	t218 = qJD(4) * t214;
	t217 = t209 * qJD(2) * t213;
	t208 = sin(pkin(11));
	t210 = cos(pkin(11));
	t216 = -t208 * t213 + t210 * t221;
	t205 = t208 * t215 + t210 * t222;
	t206 = t208 * t221 + t210 * t213;
	t207 = -t208 * t222 + t210 * t215;
	t203 = t207 * qJD(2);
	t202 = t206 * qJD(2);
	t201 = t205 * qJD(2);
	t200 = t216 * qJD(2);
	t1 = [0, -t202 * t212 + t207 * t218, 0, t203 * t214 + (-t206 * t212 - t208 * t224) * qJD(4), 0, 0; 0, t200 * t212 + t205 * t218, 0, t201 * t214 + (t210 * t224 + t212 * t216) * qJD(4), 0, 0; 0, (t212 * t220 + t213 * t218) * t209, 0, t214 * t217 + (-t211 * t214 + t212 * t223) * qJD(4), 0, 0; 0, -t202 * t214 - t207 * t219, 0, -t203 * t212 + (-t206 * t214 + t208 * t225) * qJD(4), 0, 0; 0, t200 * t214 - t205 * t219, 0, -t201 * t212 + (-t210 * t225 + t214 * t216) * qJD(4), 0, 0; 0, (-t213 * t219 + t214 * t220) * t209, 0, -t212 * t217 + (t211 * t212 + t214 * t223) * qJD(4), 0, 0; 0, -t203, 0, 0, 0, 0; 0, -t201, 0, 0, 0, 0; 0, -t217, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:01:17
	% EndTime: 2019-10-09 22:01:17
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (135->24), mult. (212->50), div. (0->0), fcn. (224->8), ass. (0->37)
	t271 = qJ(4) + qJ(5);
	t268 = sin(t271);
	t270 = qJD(4) + qJD(5);
	t290 = t268 * t270;
	t269 = cos(t271);
	t289 = t269 * t270;
	t273 = sin(pkin(6));
	t288 = t270 * t273;
	t276 = sin(qJ(2));
	t287 = t270 * t276;
	t275 = cos(pkin(6));
	t286 = t275 * t276;
	t277 = cos(qJ(2));
	t285 = t275 * t277;
	t284 = qJD(2) * t277;
	t283 = t277 * t288;
	t282 = t273 * qJD(2) * t276;
	t272 = sin(pkin(11));
	t274 = cos(pkin(11));
	t265 = t272 * t277 + t274 * t286;
	t261 = t265 * qJD(2);
	t281 = t274 * t288 + t261;
	t267 = -t272 * t286 + t274 * t277;
	t263 = t267 * qJD(2);
	t280 = t272 * t288 - t263;
	t279 = -t272 * t276 + t274 * t285;
	t266 = t272 * t285 + t274 * t276;
	t278 = -t270 * t275 + t282;
	t262 = t266 * qJD(2);
	t260 = t279 * qJD(2);
	t259 = t268 * t283 + t269 * t278;
	t258 = -t268 * t278 + t269 * t283;
	t257 = t269 * t281 + t279 * t290;
	t256 = -t268 * t281 + t279 * t289;
	t255 = -t266 * t290 - t269 * t280;
	t254 = -t266 * t289 + t268 * t280;
	t1 = [0, -t262 * t268 + t267 * t289, 0, t255, t255, 0; 0, t260 * t268 + t265 * t289, 0, t257, t257, 0; 0, (t268 * t284 + t269 * t287) * t273, 0, t259, t259, 0; 0, -t262 * t269 - t267 * t290, 0, t254, t254, 0; 0, t260 * t269 - t265 * t290, 0, t256, t256, 0; 0, (-t268 * t287 + t269 * t284) * t273, 0, t258, t258, 0; 0, -t263, 0, 0, 0, 0; 0, -t261, 0, 0, 0, 0; 0, -t282, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:01:19
	% EndTime: 2019-10-09 22:01:20
	% DurationCPUTime: 0.36s
	% Computational Cost: add. (388->65), mult. (672->132), div. (0->0), fcn. (736->10), ass. (0->66)
	t437 = qJ(4) + qJ(5);
	t434 = sin(t437);
	t445 = cos(qJ(2));
	t474 = (qJD(2) * t434 + qJD(6)) * t445;
	t436 = qJD(4) + qJD(5);
	t473 = t434 * t436;
	t435 = cos(t437);
	t472 = t435 * t436;
	t443 = sin(qJ(2));
	t471 = t436 * t443;
	t438 = sin(pkin(11));
	t439 = sin(pkin(6));
	t470 = t438 * t439;
	t440 = cos(pkin(11));
	t469 = t439 * t440;
	t468 = t439 * t443;
	t467 = t439 * t445;
	t441 = cos(pkin(6));
	t466 = t441 * t443;
	t465 = t441 * t445;
	t444 = cos(qJ(6));
	t464 = t443 * t444;
	t463 = qJD(2) * t445;
	t462 = qJD(6) * t434;
	t442 = sin(qJ(6));
	t461 = qJD(6) * t442;
	t460 = qJD(6) * t444;
	t459 = t438 * t466;
	t458 = t434 * t467;
	t457 = t435 * t467;
	t456 = t439 * t463;
	t455 = -qJD(2) - t462;
	t429 = t438 * t445 + t440 * t466;
	t425 = t429 * qJD(2);
	t453 = t436 * t469 + t425;
	t427 = -qJD(2) * t459 + t440 * t463;
	t452 = -t436 * t470 + t427;
	t451 = -t429 * t462 - t425;
	t431 = t440 * t445 - t459;
	t450 = -t431 * t462 - t427;
	t449 = -t438 * t443 + t440 * t465;
	t430 = t438 * t465 + t440 * t443;
	t448 = qJD(2) * t468 - t436 * t441;
	t424 = t449 * qJD(2);
	t447 = qJD(6) * t449 + t424 * t434 + t429 * t472;
	t426 = t430 * qJD(2);
	t446 = -qJD(6) * t430 - t426 * t434 + t431 * t472;
	t423 = t441 * t435 - t458;
	t422 = -t441 * t434 - t457;
	t421 = -t434 * t449 - t435 * t469;
	t420 = t434 * t469 - t435 * t449;
	t419 = t430 * t434 + t435 * t470;
	t418 = t430 * t435 - t434 * t470;
	t417 = t448 * t435 + t436 * t458;
	t416 = t448 * t434 - t436 * t457;
	t415 = t453 * t435 + t449 * t473;
	t414 = t453 * t434 - t449 * t472;
	t413 = -t430 * t473 + t452 * t435;
	t412 = t430 * t472 + t452 * t434;
	t411 = t417 * t444 - t422 * t461;
	t410 = -t417 * t442 - t422 * t460;
	t409 = t415 * t444 - t420 * t461;
	t408 = -t415 * t442 - t420 * t460;
	t407 = t413 * t444 - t418 * t461;
	t406 = -t413 * t442 - t418 * t460;
	t1 = [0, t450 * t442 + t446 * t444, 0, t407, t407, -t412 * t442 - t426 * t444 + (-t419 * t444 - t431 * t442) * qJD(6); 0, t451 * t442 + t447 * t444, 0, t409, t409, -t414 * t442 + t424 * t444 + (-t421 * t444 - t429 * t442) * qJD(6); 0, (t444 * t474 + (t455 * t442 + t444 * t472) * t443) * t439, 0, t411, t411, t444 * t456 - t416 * t442 + (-t423 * t444 - t442 * t468) * qJD(6); 0, -t446 * t442 + t450 * t444, 0, t406, t406, -t412 * t444 + t426 * t442 + (t419 * t442 - t431 * t444) * qJD(6); 0, -t447 * t442 + t451 * t444, 0, t408, t408, -t414 * t444 - t424 * t442 + (t421 * t442 - t429 * t444) * qJD(6); 0, (t455 * t464 + (-t435 * t471 - t474) * t442) * t439, 0, t410, t410, -t442 * t456 - t416 * t444 + (t423 * t442 - t439 * t464) * qJD(6); 0, t426 * t435 + t431 * t473, 0, t412, t412, 0; 0, -t424 * t435 + t429 * t473, 0, t414, t414, 0; 0, (t434 * t471 - t435 * t463) * t439, 0, t416, t416, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end