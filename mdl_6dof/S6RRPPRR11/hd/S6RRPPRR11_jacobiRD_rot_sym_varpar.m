% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta4]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:53
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRPPRR11_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR11_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR11_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR11_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR11_jacobiRD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:53:50
	% EndTime: 2019-10-10 09:53:50
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:53:50
	% EndTime: 2019-10-10 09:53:51
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
	% StartTime: 2019-10-10 09:53:51
	% EndTime: 2019-10-10 09:53:51
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (27->13), mult. (88->22), div. (0->0), fcn. (88->6), ass. (0->21)
	t158 = sin(qJ(2));
	t159 = sin(qJ(1));
	t171 = t158 * t159;
	t161 = cos(qJ(1));
	t170 = t158 * t161;
	t160 = cos(qJ(2));
	t169 = t159 * t160;
	t168 = t160 * t161;
	t156 = sin(pkin(6));
	t167 = qJD(1) * t156;
	t166 = qJD(2) * t156;
	t157 = cos(pkin(6));
	t165 = -t157 * t168 + t171;
	t164 = t157 * t169 + t170;
	t163 = t157 * t170 + t169;
	t162 = t157 * t171 - t168;
	t155 = t162 * qJD(1) + t165 * qJD(2);
	t154 = t164 * qJD(1) + t163 * qJD(2);
	t153 = t163 * qJD(1) + t164 * qJD(2);
	t152 = t165 * qJD(1) + t162 * qJD(2);
	t1 = [t155, t152, 0, 0, 0, 0; -t153, -t154, 0, 0, 0, 0; 0, -t158 * t166, 0, 0, 0, 0; t154, t153, 0, 0, 0, 0; t152, t155, 0, 0, 0, 0; 0, -t160 * t166, 0, 0, 0, 0; -t159 * t167, 0, 0, 0, 0, 0; t161 * t167, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:53:51
	% EndTime: 2019-10-10 09:53:52
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (25->11), mult. (88->22), div. (0->0), fcn. (88->6), ass. (0->21)
	t201 = sin(qJ(2));
	t202 = sin(qJ(1));
	t214 = t201 * t202;
	t204 = cos(qJ(1));
	t213 = t201 * t204;
	t203 = cos(qJ(2));
	t212 = t202 * t203;
	t211 = t203 * t204;
	t199 = sin(pkin(6));
	t210 = qJD(1) * t199;
	t209 = qJD(2) * t199;
	t200 = cos(pkin(6));
	t208 = t200 * t211 - t214;
	t207 = t200 * t212 + t213;
	t206 = t200 * t213 + t212;
	t205 = -t200 * t214 + t211;
	t198 = t205 * qJD(1) + t208 * qJD(2);
	t197 = t207 * qJD(1) + t206 * qJD(2);
	t196 = t206 * qJD(1) + t207 * qJD(2);
	t195 = t208 * qJD(1) + t205 * qJD(2);
	t1 = [-t202 * t210, 0, 0, 0, 0, 0; t204 * t210, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t198, t195, 0, 0, 0, 0; t196, t197, 0, 0, 0, 0; 0, t201 * t209, 0, 0, 0, 0; -t197, -t196, 0, 0, 0, 0; t195, t198, 0, 0, 0, 0; 0, t203 * t209, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:53:52
	% EndTime: 2019-10-10 09:53:52
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (41->16), mult. (148->39), div. (0->0), fcn. (148->8), ass. (0->25)
	t234 = sin(qJ(2));
	t235 = sin(qJ(1));
	t249 = t234 * t235;
	t237 = cos(qJ(1));
	t248 = t234 * t237;
	t236 = cos(qJ(2));
	t247 = t235 * t236;
	t246 = t236 * t237;
	t231 = sin(pkin(6));
	t245 = qJD(1) * t231;
	t244 = qJD(2) * t236;
	t233 = cos(pkin(6));
	t243 = t233 * t246;
	t242 = t235 * t245;
	t241 = t237 * t245;
	t240 = t231 * t244;
	t239 = -t233 * t247 - t248;
	t238 = -t233 * t248 - t247;
	t232 = cos(pkin(11));
	t230 = sin(pkin(11));
	t227 = (t243 - t249) * qJD(2) + (-t233 * t249 + t246) * qJD(1);
	t226 = t239 * qJD(1) + t238 * qJD(2);
	t225 = t238 * qJD(1) + t239 * qJD(2);
	t224 = -qJD(1) * t243 - t237 * t244 + (qJD(2) * t233 + qJD(1)) * t249;
	t1 = [t226 * t230 - t232 * t242, t225 * t230, 0, 0, 0, 0; -t224 * t230 + t232 * t241, t227 * t230, 0, 0, 0, 0; 0, t230 * t240, 0, 0, 0, 0; t226 * t232 + t230 * t242, t225 * t232, 0, 0, 0, 0; -t224 * t232 - t230 * t241, t227 * t232, 0, 0, 0, 0; 0, t232 * t240, 0, 0, 0, 0; -t227, t224, 0, 0, 0, 0; t225, t226, 0, 0, 0, 0; 0, -t231 * qJD(2) * t234, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:53:52
	% EndTime: 2019-10-10 09:53:52
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (145->38), mult. (310->71), div. (0->0), fcn. (322->8), ass. (0->38)
	t302 = cos(pkin(6));
	t303 = sin(qJ(2));
	t306 = cos(qJ(1));
	t321 = t306 * t303;
	t304 = sin(qJ(1));
	t305 = cos(qJ(2));
	t322 = t304 * t305;
	t292 = t302 * t321 + t322;
	t293 = t302 * t322 + t321;
	t289 = t293 * qJD(1) + t292 * qJD(2);
	t300 = pkin(11) + qJ(5);
	t298 = sin(t300);
	t299 = cos(t300);
	t320 = t306 * t305;
	t314 = t302 * t320;
	t323 = t304 * t303;
	t309 = t314 - t323;
	t301 = sin(pkin(6));
	t319 = qJD(1) * t301;
	t313 = t304 * t319;
	t324 = t301 * t306;
	t327 = (t298 * t309 + t299 * t324) * qJD(5) + t289 * t299 - t298 * t313;
	t326 = t301 * t304;
	t325 = t301 * t305;
	t318 = qJD(2) * t305;
	t317 = qJD(5) * t298;
	t316 = qJD(5) * t299;
	t315 = qJD(5) * t303;
	t312 = t306 * t319;
	t311 = t301 * qJD(2) * t303;
	t294 = -t302 * t323 + t320;
	t307 = -t299 * t313 - t289 * t298 + (-t298 * t324 + t299 * t309) * qJD(5);
	t290 = t294 * qJD(1) + t309 * qJD(2);
	t288 = -t292 * qJD(1) - t293 * qJD(2);
	t287 = -qJD(1) * t314 - t306 * t318 + (qJD(2) * t302 + qJD(1)) * t323;
	t286 = t299 * t312 - t287 * t298 + (t293 * t299 - t298 * t326) * qJD(5);
	t285 = -t298 * t312 - t287 * t299 + (-t293 * t298 - t299 * t326) * qJD(5);
	t1 = [t307, t288 * t298 + t294 * t316, 0, 0, t285, 0; t286, t290 * t298 + t292 * t316, 0, 0, t327, 0; 0, (t298 * t318 + t299 * t315) * t301, 0, 0, t299 * t311 + (t298 * t325 - t299 * t302) * qJD(5), 0; -t327, t288 * t299 - t294 * t317, 0, 0, -t286, 0; t285, t290 * t299 - t292 * t317, 0, 0, t307, 0; 0, (-t298 * t315 + t299 * t318) * t301, 0, 0, -t298 * t311 + (t298 * t302 + t299 * t325) * qJD(5), 0; -t290, t287, 0, 0, 0, 0; t288, -t289, 0, 0, 0, 0; 0, -t311, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:53:55
	% EndTime: 2019-10-10 09:53:55
	% DurationCPUTime: 0.56s
	% Computational Cost: add. (424->73), mult. (898->141), div. (0->0), fcn. (964->10), ass. (0->65)
	t486 = cos(pkin(6));
	t488 = sin(qJ(2));
	t492 = cos(qJ(1));
	t518 = t492 * t488;
	t489 = sin(qJ(1));
	t491 = cos(qJ(2));
	t519 = t489 * t491;
	t475 = t486 * t518 + t519;
	t476 = t486 * t519 + t518;
	t463 = t476 * qJD(1) + t475 * qJD(2);
	t517 = t492 * t491;
	t520 = t489 * t488;
	t474 = -t486 * t517 + t520;
	t484 = pkin(11) + qJ(5);
	t482 = sin(t484);
	t483 = cos(t484);
	t485 = sin(pkin(6));
	t522 = t485 * t492;
	t469 = t474 * t483 + t482 * t522;
	t516 = qJD(1) * t485;
	t505 = t489 * t516;
	t456 = t469 * qJD(5) + t463 * t482 + t483 * t505;
	t507 = t486 * t520;
	t515 = qJD(2) * t488;
	t464 = -qJD(1) * t507 - t489 * t515 + (qJD(2) * t486 + qJD(1)) * t517;
	t506 = t483 * t522;
	t470 = -t474 * t482 + t506;
	t487 = sin(qJ(6));
	t490 = cos(qJ(6));
	t533 = -t456 * t490 + (-t470 * t487 - t475 * t490) * qJD(6) - t464 * t487;
	t532 = (t470 * t490 - t475 * t487) * qJD(6) - t456 * t487 + t464 * t490;
	t529 = (qJD(2) * t482 + qJD(6)) * t491;
	t524 = t485 * t489;
	t523 = t485 * t491;
	t521 = t488 * t490;
	t514 = qJD(2) * t491;
	t513 = qJD(5) * t482;
	t512 = qJD(5) * t483;
	t511 = qJD(5) * t488;
	t510 = qJD(6) * t482;
	t509 = qJD(6) * t487;
	t508 = qJD(6) * t490;
	t504 = t492 * t516;
	t503 = t485 * t515;
	t502 = t485 * t514;
	t501 = -qJD(2) - t510;
	t496 = t507 - t517;
	t461 = t474 * qJD(1) + t496 * qJD(2);
	t499 = t496 * t510 + t461;
	t498 = -t475 * t510 - t463;
	t468 = t476 * t482 + t483 * t524;
	t467 = t476 * t483 - t482 * t524;
	t472 = -t486 * t482 - t483 * t523;
	t497 = t482 * t523 - t486 * t483;
	t462 = -t475 * qJD(1) - t476 * qJD(2);
	t494 = -qJD(6) * t476 + t462 * t482 - t496 * t512;
	t493 = -qJD(6) * t474 + t464 * t482 + t475 * t512;
	t455 = t463 * t483 + qJD(5) * t506 + (-qJD(5) * t474 - t505) * t482;
	t466 = t497 * qJD(5) + t483 * t503;
	t465 = t472 * qJD(5) + t482 * t503;
	t459 = t467 * qJD(5) - t461 * t482 + t483 * t504;
	t458 = t468 * qJD(5) + t461 * t483 + t482 * t504;
	t454 = t459 * t490 + t462 * t487 + (-t468 * t487 - t490 * t496) * qJD(6);
	t453 = -t459 * t487 + t462 * t490 + (-t468 * t490 + t487 * t496) * qJD(6);
	t1 = [t533, t499 * t487 + t494 * t490, 0, 0, -t458 * t490 - t467 * t509, t453; t454, t498 * t487 + t493 * t490, 0, 0, t455 * t490 - t469 * t509, t532; 0, (t490 * t529 + (t501 * t487 + t490 * t512) * t488) * t485, 0, 0, t466 * t490 - t472 * t509, t490 * t502 - t465 * t487 + (-t485 * t487 * t488 + t490 * t497) * qJD(6); -t532, -t494 * t487 + t499 * t490, 0, 0, t458 * t487 - t467 * t508, -t454; t453, -t493 * t487 + t498 * t490, 0, 0, -t455 * t487 - t469 * t508, t533; 0, (t501 * t521 + (-t483 * t511 - t529) * t487) * t485, 0, 0, -t466 * t487 - t472 * t508, -t487 * t502 - t465 * t490 + (-t485 * t521 - t487 * t497) * qJD(6); t455, -t462 * t483 - t496 * t513, 0, 0, t459, 0; t458, -t464 * t483 + t475 * t513, 0, 0, t456, 0; 0, (t482 * t511 - t483 * t514) * t485, 0, 0, t465, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end