% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRPRRR3
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:57
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6PRPRRR3_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR3_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR3_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRR3_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR3_jacobiRD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:57:34
	% EndTime: 2019-10-09 21:57:35
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:57:34
	% EndTime: 2019-10-09 21:57:34
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:57:35
	% EndTime: 2019-10-09 21:57:35
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
	% StartTime: 2019-10-09 21:57:35
	% EndTime: 2019-10-09 21:57:35
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (9->7), mult. (42->23), div. (0->0), fcn. (42->8), ass. (0->14)
	t176 = cos(pkin(6));
	t177 = sin(qJ(2));
	t182 = t176 * t177;
	t178 = cos(qJ(2));
	t181 = t176 * t178;
	t180 = qJD(2) * sin(pkin(6));
	t179 = t177 * t180;
	t175 = cos(pkin(11));
	t174 = cos(pkin(12));
	t172 = sin(pkin(11));
	t171 = sin(pkin(12));
	t170 = (t172 * t182 - t175 * t178) * qJD(2);
	t169 = (-t172 * t178 - t175 * t182) * qJD(2);
	t1 = [0, t170 * t174, 0, 0, 0, 0; 0, t169 * t174, 0, 0, 0, 0; 0, -t174 * t179, 0, 0, 0, 0; 0, -t170 * t171, 0, 0, 0, 0; 0, -t169 * t171, 0, 0, 0, 0; 0, t171 * t179, 0, 0, 0, 0; 0, (-t172 * t181 - t175 * t177) * qJD(2), 0, 0, 0, 0; 0, (-t172 * t177 + t175 * t181) * qJD(2), 0, 0, 0, 0; 0, t178 * t180, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:57:36
	% EndTime: 2019-10-09 21:57:36
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (66->23), mult. (140->60), div. (0->0), fcn. (148->8), ass. (0->28)
	t224 = sin(pkin(11));
	t225 = sin(pkin(6));
	t240 = t224 * t225;
	t226 = cos(pkin(11));
	t239 = t225 * t226;
	t228 = sin(qJ(2));
	t238 = t225 * t228;
	t227 = cos(pkin(6));
	t237 = t227 * t228;
	t229 = cos(qJ(2));
	t236 = t227 * t229;
	t235 = qJD(2) * t228;
	t223 = pkin(12) + qJ(4);
	t221 = sin(t223);
	t234 = qJD(4) * t221;
	t222 = cos(t223);
	t233 = qJD(4) * t222;
	t232 = qJD(4) * t229;
	t231 = t225 * qJD(2) * t229;
	t217 = -t224 * t228 + t226 * t236;
	t218 = t224 * t229 + t226 * t237;
	t219 = -t224 * t236 - t226 * t228;
	t230 = t224 * t237 - t226 * t229;
	t216 = t230 * qJD(2);
	t215 = t219 * qJD(2);
	t214 = t218 * qJD(2);
	t213 = t217 * qJD(2);
	t1 = [0, t216 * t222 - t219 * t234, 0, -t215 * t221 + (-t221 * t240 + t222 * t230) * qJD(4), 0, 0; 0, -t214 * t222 - t217 * t234, 0, -t213 * t221 + (-t218 * t222 + t221 * t239) * qJD(4), 0, 0; 0, (-t221 * t232 - t222 * t235) * t225, 0, -t221 * t231 + (-t221 * t227 - t222 * t238) * qJD(4), 0, 0; 0, -t216 * t221 - t219 * t233, 0, -t215 * t222 + (-t221 * t230 - t222 * t240) * qJD(4), 0, 0; 0, t214 * t221 - t217 * t233, 0, -t213 * t222 + (t218 * t221 + t222 * t239) * qJD(4), 0, 0; 0, (t221 * t235 - t222 * t232) * t225, 0, -t222 * t231 + (t221 * t238 - t222 * t227) * qJD(4), 0, 0; 0, t215, 0, 0, 0, 0; 0, t213, 0, 0, 0, 0; 0, t231, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:57:36
	% EndTime: 2019-10-09 21:57:36
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (182->21), mult. (212->50), div. (0->0), fcn. (224->8), ass. (0->37)
	t275 = pkin(12) + qJ(4) + qJ(5);
	t273 = sin(t275);
	t276 = qJD(4) + qJD(5);
	t295 = t273 * t276;
	t274 = cos(t275);
	t294 = t274 * t276;
	t278 = sin(pkin(6));
	t293 = t276 * t278;
	t282 = cos(qJ(2));
	t292 = t276 * t282;
	t280 = cos(pkin(6));
	t281 = sin(qJ(2));
	t291 = t280 * t281;
	t290 = t280 * t282;
	t289 = qJD(2) * t281;
	t288 = t281 * t293;
	t287 = t278 * qJD(2) * t282;
	t277 = sin(pkin(11));
	t279 = cos(pkin(11));
	t269 = -t277 * t281 + t279 * t290;
	t265 = t269 * qJD(2);
	t286 = t279 * t293 - t265;
	t271 = -t277 * t290 - t279 * t281;
	t267 = t271 * qJD(2);
	t285 = -t277 * t293 - t267;
	t270 = t277 * t282 + t279 * t291;
	t284 = t277 * t291 - t279 * t282;
	t283 = -t276 * t280 - t287;
	t268 = t284 * qJD(2);
	t266 = t270 * qJD(2);
	t264 = t273 * t288 + t283 * t274;
	t263 = t283 * t273 - t274 * t288;
	t262 = t285 * t274 - t284 * t295;
	t261 = t285 * t273 + t284 * t294;
	t260 = t270 * t295 + t286 * t274;
	t259 = -t270 * t294 + t286 * t273;
	t1 = [0, t268 * t274 - t271 * t295, 0, t261, t261, 0; 0, -t266 * t274 - t269 * t295, 0, t259, t259, 0; 0, (-t273 * t292 - t274 * t289) * t278, 0, t263, t263, 0; 0, -t268 * t273 - t271 * t294, 0, t262, t262, 0; 0, t266 * t273 - t269 * t294, 0, t260, t260, 0; 0, (t273 * t289 - t274 * t292) * t278, 0, t264, t264, 0; 0, t267, 0, 0, 0, 0; 0, t265, 0, 0, 0, 0; 0, t287, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:57:37
	% EndTime: 2019-10-09 21:57:38
	% DurationCPUTime: 0.36s
	% Computational Cost: add. (520->64), mult. (672->131), div. (0->0), fcn. (736->10), ass. (0->66)
	t445 = pkin(12) + qJ(4) + qJ(5);
	t443 = sin(t445);
	t444 = cos(t445);
	t452 = sin(qJ(2));
	t446 = qJD(4) + qJD(5);
	t454 = cos(qJ(2));
	t481 = t446 * t454;
	t484 = (qJD(2) * t444 - qJD(6)) * t452 + t443 * t481;
	t483 = t443 * t446;
	t482 = t444 * t446;
	t447 = sin(pkin(11));
	t448 = sin(pkin(6));
	t480 = t447 * t448;
	t449 = cos(pkin(11));
	t479 = t448 * t449;
	t478 = t448 * t452;
	t477 = t448 * t454;
	t450 = cos(pkin(6));
	t476 = t450 * t452;
	t475 = t450 * t454;
	t474 = qJD(2) * t452;
	t473 = qJD(2) * t454;
	t472 = qJD(6) * t444;
	t451 = sin(qJ(6));
	t471 = qJD(6) * t451;
	t453 = cos(qJ(6));
	t470 = qJD(6) * t453;
	t468 = t447 * t476;
	t467 = t443 * t478;
	t466 = t444 * t478;
	t465 = t448 * t474;
	t458 = -t447 * t452 + t449 * t475;
	t433 = t458 * qJD(2);
	t463 = t446 * t479 - t433;
	t439 = t447 * t475 + t449 * t452;
	t435 = t439 * qJD(2);
	t462 = t446 * t480 - t435;
	t461 = -t458 * t472 + t433;
	t460 = t439 * t472 - t435;
	t459 = (qJD(2) - t472) * t454;
	t438 = t447 * t454 + t449 * t476;
	t457 = t446 * t450 + t448 * t473;
	t434 = t438 * qJD(2);
	t456 = qJD(6) * t438 - t434 * t444 - t458 * t483;
	t436 = -qJD(2) * t468 + t449 * t473;
	t440 = t449 * t454 - t468;
	t455 = qJD(6) * t440 - t436 * t444 + t439 * t483;
	t432 = t450 * t443 + t466;
	t431 = t450 * t444 - t467;
	t430 = t440 * t444 + t443 * t480;
	t429 = -t440 * t443 + t444 * t480;
	t428 = t438 * t444 - t443 * t479;
	t427 = -t438 * t443 - t444 * t479;
	t426 = t457 * t444 - t446 * t467;
	t425 = -t457 * t443 - t446 * t466;
	t424 = -t440 * t483 + t462 * t444;
	t423 = -t440 * t482 - t462 * t443;
	t422 = -t438 * t483 - t463 * t444;
	t421 = -t438 * t482 + t463 * t443;
	t420 = t425 * t453 - t431 * t471;
	t419 = -t425 * t451 - t431 * t470;
	t418 = t423 * t453 - t429 * t471;
	t417 = -t423 * t451 - t429 * t470;
	t416 = t421 * t453 - t427 * t471;
	t415 = -t421 * t451 - t427 * t470;
	t1 = [0, t460 * t451 + t455 * t453, 0, t418, t418, -t424 * t451 + t436 * t453 + (-t430 * t453 - t439 * t451) * qJD(6); 0, t461 * t451 + t456 * t453, 0, t416, t416, -t422 * t451 + t434 * t453 + (-t428 * t453 + t451 * t458) * qJD(6); 0, (t451 * t459 - t484 * t453) * t448, 0, t420, t420, t453 * t465 - t426 * t451 + (-t432 * t453 + t451 * t477) * qJD(6); 0, -t455 * t451 + t460 * t453, 0, t417, t417, -t424 * t453 - t436 * t451 + (t430 * t451 - t439 * t453) * qJD(6); 0, -t456 * t451 + t461 * t453, 0, t415, t415, -t422 * t453 - t434 * t451 + (t428 * t451 + t453 * t458) * qJD(6); 0, (t484 * t451 + t453 * t459) * t448, 0, t419, t419, -t451 * t465 - t426 * t453 + (t432 * t451 + t453 * t477) * qJD(6); 0, -t436 * t443 - t439 * t482, 0, t424, t424, 0; 0, -t434 * t443 + t458 * t482, 0, t422, t422, 0; 0, (-t443 * t474 + t444 * t481) * t448, 0, t426, t426, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end