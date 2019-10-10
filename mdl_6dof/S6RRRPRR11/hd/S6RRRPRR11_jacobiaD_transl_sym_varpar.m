% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRPRR11
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:10
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRR11_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR11_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR11_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR11_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR11_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR11_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:10:18
	% EndTime: 2019-10-10 12:10:18
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:10:18
	% EndTime: 2019-10-10 12:10:18
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:10:19
	% EndTime: 2019-10-10 12:10:19
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (35->18), mult. (110->35), div. (0->0), fcn. (94->6), ass. (0->20)
	t136 = sin(pkin(6));
	t151 = t136 * (pkin(8) + r_i_i_C(3));
	t138 = sin(qJ(2));
	t139 = sin(qJ(1));
	t149 = t138 * t139;
	t141 = cos(qJ(1));
	t148 = t138 * t141;
	t140 = cos(qJ(2));
	t147 = t139 * t140;
	t146 = t140 * t141;
	t137 = cos(pkin(6));
	t145 = -t137 * t146 + t149;
	t144 = t137 * t147 + t148;
	t143 = t137 * t148 + t147;
	t142 = t137 * t149 - t146;
	t135 = t142 * qJD(1) + t145 * qJD(2);
	t134 = t144 * qJD(1) + t143 * qJD(2);
	t133 = t143 * qJD(1) + t144 * qJD(2);
	t132 = t145 * qJD(1) + t142 * qJD(2);
	t1 = [t135 * r_i_i_C(1) + t134 * r_i_i_C(2) + (-pkin(1) * t141 - t139 * t151) * qJD(1), t132 * r_i_i_C(1) + t133 * r_i_i_C(2), 0, 0, 0, 0; -t133 * r_i_i_C(1) + t132 * r_i_i_C(2) + (-pkin(1) * t139 + t141 * t151) * qJD(1), -t134 * r_i_i_C(1) + t135 * r_i_i_C(2), 0, 0, 0, 0; 0, (-r_i_i_C(1) * t138 - r_i_i_C(2) * t140) * t136 * qJD(2), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:10:20
	% EndTime: 2019-10-10 12:10:20
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (145->51), mult. (441->91), div. (0->0), fcn. (412->8), ass. (0->38)
	t242 = cos(pkin(6));
	t244 = sin(qJ(2));
	t245 = sin(qJ(1));
	t264 = t245 * t244;
	t258 = t242 * t264;
	t247 = cos(qJ(2));
	t248 = cos(qJ(1));
	t261 = t248 * t247;
	t233 = -qJD(1) * t258 - qJD(2) * t264 + (qJD(2) * t242 + qJD(1)) * t261;
	t241 = sin(pkin(6));
	t265 = t241 * t248;
	t256 = qJD(3) * t265;
	t269 = t233 - t256;
	t268 = -r_i_i_C(3) - pkin(9);
	t267 = t241 * t245;
	t246 = cos(qJ(3));
	t266 = t241 * t246;
	t263 = t245 * t247;
	t262 = t248 * t244;
	t260 = qJD(1) * t241;
	t235 = t242 * t262 + t263;
	t259 = qJD(3) * t235;
	t257 = t248 * t260;
	t243 = sin(qJ(3));
	t255 = -r_i_i_C(1) * t243 - r_i_i_C(2) * t246;
	t254 = t246 * r_i_i_C(1) - t243 * r_i_i_C(2) + pkin(2);
	t253 = t242 * t261 - t264;
	t252 = t242 * t263 + t262;
	t251 = t258 - t261;
	t250 = qJD(3) * t255;
	t249 = t245 * t260 - t259;
	t238 = t246 * t256;
	t232 = t252 * qJD(1) + t235 * qJD(2);
	t231 = t235 * qJD(1) + t252 * qJD(2);
	t230 = -t253 * qJD(1) + t251 * qJD(2);
	t229 = t243 * t257 - t231 * t246 + (t243 * t251 + t245 * t266) * qJD(3);
	t228 = t246 * t257 + t231 * t243 + (-t243 * t267 + t246 * t251) * qJD(3);
	t1 = [(-t233 * t246 + t243 * t259 + t238) * r_i_i_C(1) + (t269 * t243 + t246 * t259) * r_i_i_C(2) - t233 * pkin(2) + t268 * t232 + (-t248 * pkin(1) + (-pkin(8) + t255) * t267) * qJD(1), t254 * t230 + t268 * t231 - t252 * t250, t228 * r_i_i_C(1) - t229 * r_i_i_C(2), 0, 0, 0; -t231 * pkin(2) + t229 * r_i_i_C(1) + t228 * r_i_i_C(2) + t268 * t230 + (-pkin(1) * t245 + pkin(8) * t265) * qJD(1), -t254 * t232 - t268 * t233 + t253 * t250, t238 * r_i_i_C(2) + (t249 * r_i_i_C(1) - t233 * r_i_i_C(2)) * t246 + (-t269 * r_i_i_C(1) - t249 * r_i_i_C(2)) * t243, 0, 0, 0; 0, (t247 * t250 + (-t254 * t244 - t268 * t247) * qJD(2)) * t241, t255 * t247 * t241 * qJD(2) + ((-t242 * t243 - t244 * t266) * r_i_i_C(1) + (t241 * t243 * t244 - t242 * t246) * r_i_i_C(2)) * qJD(3), 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:10:20
	% EndTime: 2019-10-10 12:10:21
	% DurationCPUTime: 0.29s
	% Computational Cost: add. (278->57), mult. (828->96), div. (0->0), fcn. (802->8), ass. (0->43)
	t297 = cos(pkin(6));
	t299 = sin(qJ(2));
	t300 = sin(qJ(1));
	t324 = t300 * t299;
	t318 = t297 * t324;
	t302 = cos(qJ(2));
	t303 = cos(qJ(1));
	t321 = t303 * t302;
	t283 = -qJD(1) * t318 - qJD(2) * t324 + (qJD(2) * t297 + qJD(1)) * t321;
	t296 = sin(pkin(6));
	t325 = t296 * t303;
	t334 = -qJD(3) * t325 + t283;
	t322 = t303 * t299;
	t323 = t300 * t302;
	t287 = t297 * t322 + t323;
	t320 = qJD(1) * t296;
	t333 = -qJD(3) * t287 + t300 * t320;
	t298 = sin(qJ(3));
	t301 = cos(qJ(3));
	t332 = t333 * t298 + t334 * t301;
	t328 = r_i_i_C(3) + qJ(4);
	t329 = -r_i_i_C(1) - pkin(3);
	t331 = t328 * t298 - t329 * t301 + pkin(2);
	t330 = pkin(9) + r_i_i_C(2);
	t327 = t296 * t300;
	t326 = t296 * t301;
	t316 = t303 * t320;
	t315 = qJD(2) * t296 * t302;
	t307 = t318 - t321;
	t312 = t298 * t307 + t300 * t326;
	t311 = t298 * t327 - t301 * t307;
	t310 = t297 * t298 + t299 * t326;
	t309 = t297 * t321 - t324;
	t308 = t297 * t323 + t322;
	t276 = t334 * t298 - t333 * t301;
	t304 = qJD(4) * t298 + (t329 * t298 + t328 * t301) * qJD(3);
	t284 = t310 * qJD(3) + t298 * t315;
	t282 = t308 * qJD(1) + t287 * qJD(2);
	t281 = t287 * qJD(1) + t308 * qJD(2);
	t280 = -t309 * qJD(1) + t307 * qJD(2);
	t275 = t312 * qJD(3) - t281 * t301 + t298 * t316;
	t274 = t311 * qJD(3) - t281 * t298 - t301 * t316;
	t1 = [-(t287 * t298 + t301 * t325) * qJD(4) - t283 * pkin(2) - t330 * t282 + t329 * t332 - t328 * t276 + (-t303 * pkin(1) - pkin(8) * t327) * qJD(1), t331 * t280 - t330 * t281 - t304 * t308, t311 * qJD(4) + t329 * t274 + t328 * t275, t274, 0, 0; -t312 * qJD(4) - t281 * pkin(2) - t330 * t280 - t329 * t275 + t328 * t274 + (-t300 * pkin(1) + pkin(8) * t325) * qJD(1), -t282 * t331 + t330 * t283 + t304 * t309, -(-t287 * t301 + t298 * t325) * qJD(4) + t328 * t332 + t329 * t276, t276, 0, 0; 0, (t304 * t302 + (-t299 * t331 + t330 * t302) * qJD(2)) * t296, t310 * qJD(4) + t328 * (t301 * t315 + (-t296 * t298 * t299 + t297 * t301) * qJD(3)) + t329 * t284, t284, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:10:20
	% EndTime: 2019-10-10 12:10:21
	% DurationCPUTime: 0.68s
	% Computational Cost: add. (623->93), mult. (1848->153), div. (0->0), fcn. (1896->10), ass. (0->59)
	t337 = cos(pkin(6));
	t341 = sin(qJ(1));
	t340 = sin(qJ(2));
	t374 = t341 * t340;
	t366 = t337 * t374;
	t369 = qJD(2) * t340;
	t344 = cos(qJ(2));
	t345 = cos(qJ(1));
	t371 = t345 * t344;
	t312 = -qJD(1) * t366 - t341 * t369 + (qJD(2) * t337 + qJD(1)) * t371;
	t336 = sin(pkin(6));
	t375 = t336 * t345;
	t387 = -qJD(3) * t375 + t312;
	t372 = t345 * t340;
	t373 = t341 * t344;
	t324 = t337 * t372 + t373;
	t370 = qJD(1) * t336;
	t386 = -qJD(3) * t324 + t341 * t370;
	t339 = sin(qJ(3));
	t343 = cos(qJ(3));
	t315 = t324 * t339 + t343 * t375;
	t316 = t324 * t343 - t339 * t375;
	t338 = sin(qJ(5));
	t342 = cos(qJ(5));
	t382 = t315 * t338 + t316 * t342;
	t383 = t315 * t342 - t316 * t338;
	t385 = (t382 * r_i_i_C(1) + t383 * r_i_i_C(2)) * qJD(5);
	t321 = t336 * t340 * t339 - t337 * t343;
	t376 = t336 * t343;
	t322 = t337 * t339 + t340 * t376;
	t384 = ((t321 * t338 + t322 * t342) * r_i_i_C(1) + (t321 * t342 - t322 * t338) * r_i_i_C(2)) * qJD(5);
	t306 = t386 * t339 + t387 * t343;
	t378 = pkin(3) + pkin(4);
	t353 = t342 * r_i_i_C(1) - t338 * r_i_i_C(2) + t378;
	t354 = t338 * r_i_i_C(1) + t342 * r_i_i_C(2) + qJ(4);
	t379 = t354 * t339 + t353 * t343 + pkin(2);
	t377 = t336 * t341;
	t367 = r_i_i_C(3) + pkin(10) - pkin(9);
	t364 = t345 * t370;
	t363 = qJD(2) * t336 * t344;
	t349 = t366 - t371;
	t320 = t339 * t377 - t343 * t349;
	t352 = t339 * t349 + t341 * t376;
	t358 = -t320 * t338 - t342 * t352;
	t357 = t320 * t342 - t338 * t352;
	t351 = t337 * t371 - t374;
	t350 = t337 * t373 + t372;
	t305 = t387 * t339 - t386 * t343;
	t346 = t339 * qJD(4) + ((-t338 * t343 + t339 * t342) * r_i_i_C(1) + (-t338 * t339 - t342 * t343) * r_i_i_C(2)) * qJD(5) + (-t353 * t339 + t354 * t343) * qJD(3);
	t314 = -t321 * qJD(3) + t343 * t363;
	t313 = t322 * qJD(3) + t339 * t363;
	t311 = t350 * qJD(1) + t324 * qJD(2);
	t310 = t324 * qJD(1) + t350 * qJD(2);
	t309 = -t351 * qJD(1) + t349 * qJD(2);
	t304 = t352 * qJD(3) - t310 * t343 + t339 * t364;
	t303 = t320 * qJD(3) - t310 * t339 - t343 * t364;
	t302 = t358 * qJD(5) + t303 * t338 + t304 * t342;
	t301 = -t357 * qJD(5) + t303 * t342 - t304 * t338;
	t1 = [-t312 * pkin(2) - t315 * qJD(4) - t354 * t305 - t353 * t306 + (-t383 * r_i_i_C(1) + t382 * r_i_i_C(2)) * qJD(5) + (-t345 * pkin(1) - pkin(8) * t377) * qJD(1) + t367 * t311, t379 * t309 + t367 * t310 - t346 * t350, t320 * qJD(4) + t354 * t304 - t353 * t303 + (t357 * r_i_i_C(1) + t358 * r_i_i_C(2)) * qJD(5), t303, t301 * r_i_i_C(1) - t302 * r_i_i_C(2), 0; -t310 * pkin(2) + t302 * r_i_i_C(1) + t301 * r_i_i_C(2) + t303 * qJ(4) - t352 * qJD(4) + t378 * t304 + (-pkin(1) * t341 + pkin(8) * t375) * qJD(1) + t367 * t309, -t311 * t379 - t367 * t312 + t346 * t351, t316 * qJD(4) - t353 * t305 + t354 * t306 + t385, t305, (t305 * t342 - t306 * t338) * r_i_i_C(1) + (-t305 * t338 - t306 * t342) * r_i_i_C(2) - t385, 0; 0, (-t379 * t369 + (-t367 * qJD(2) + t346) * t344) * t336, t322 * qJD(4) - t353 * t313 + t354 * t314 + t384, t313, (t313 * t342 - t314 * t338) * r_i_i_C(1) + (-t313 * t338 - t314 * t342) * r_i_i_C(2) - t384, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:10:25
	% EndTime: 2019-10-10 12:10:27
	% DurationCPUTime: 1.69s
	% Computational Cost: add. (1501->145), mult. (4424->242), div. (0->0), fcn. (4722->12), ass. (0->84)
	t625 = sin(qJ(2));
	t626 = sin(qJ(1));
	t630 = cos(qJ(2));
	t676 = cos(pkin(6));
	t677 = cos(qJ(1));
	t654 = t676 * t677;
	t607 = t625 * t654 + t626 * t630;
	t624 = sin(qJ(3));
	t629 = cos(qJ(3));
	t621 = sin(pkin(6));
	t665 = t621 * t677;
	t597 = t607 * t624 + t629 * t665;
	t598 = t607 * t629 - t624 * t665;
	t623 = sin(qJ(5));
	t628 = cos(qJ(5));
	t581 = t597 * t623 + t598 * t628;
	t606 = t626 * t625 - t630 * t654;
	t622 = sin(qJ(6));
	t627 = cos(qJ(6));
	t705 = t581 * t622 + t606 * t627;
	t704 = t581 * t627 - t606 * t622;
	t672 = t621 * t626;
	t698 = qJD(1) * t672 - qJD(3) * t607;
	t659 = t626 * t676;
	t655 = t625 * t659;
	t661 = t677 * qJD(1);
	t669 = qJD(2) * t625;
	t594 = -qJD(1) * t655 - t626 * t669 + (qJD(2) * t654 + t661) * t630;
	t699 = -qJD(3) * t665 + t594;
	t573 = t699 * t624 - t698 * t629;
	t574 = t698 * t624 + t699 * t629;
	t703 = t581 * qJD(5) - t573 * t628 + t574 * t623;
	t646 = t597 * t628 - t598 * t623;
	t559 = t646 * qJD(5) + t573 * t623 + t574 * t628;
	t643 = t623 * t629 - t624 * t628;
	t686 = qJD(5) - qJD(3);
	t632 = t686 * t643;
	t634 = t677 * t625 + t630 * t659;
	t592 = t607 * qJD(1) + t634 * qJD(2);
	t635 = -t677 * t630 + t655;
	t602 = t624 * t672 - t629 * t635;
	t657 = t621 * t661;
	t571 = t602 * qJD(3) - t592 * t624 - t629 * t657;
	t671 = t621 * t629;
	t639 = t624 * t635 + t626 * t671;
	t572 = t639 * qJD(3) - t592 * t629 + t624 * t657;
	t585 = t602 * t628 - t623 * t639;
	t554 = t585 * qJD(5) - t571 * t628 + t572 * t623;
	t645 = -t602 * t623 - t628 * t639;
	t555 = t645 * qJD(5) + t571 * t623 + t572 * t628;
	t652 = -r_i_i_C(1) * t622 - r_i_i_C(2) * t627;
	t638 = qJD(6) * t652;
	t653 = -r_i_i_C(1) * t627 + r_i_i_C(2) * t622;
	t641 = pkin(5) - t653;
	t678 = -r_i_i_C(3) - pkin(11);
	t697 = -t641 * t554 - t678 * t555 + t645 * t638;
	t680 = pkin(3) + pkin(4);
	t695 = (-qJ(4) * t629 + t680 * t624) * qJD(3) - t624 * qJD(4);
	t605 = t676 * t624 + t625 * t671;
	t670 = t621 * t630;
	t662 = qJD(2) * t670;
	t595 = t605 * qJD(3) + t624 * t662;
	t604 = t621 * t625 * t624 - t676 * t629;
	t596 = -t604 * qJD(3) + t629 * t662;
	t644 = t604 * t628 - t605 * t623;
	t569 = t644 * qJD(5) + t595 * t623 + t596 * t628;
	t590 = t604 * t623 + t605 * t628;
	t694 = -t641 * (t590 * qJD(5) - t595 * t628 + t596 * t623) - t678 * t569 + t644 * t638;
	t693 = -t678 * t559 + t646 * t638 - t641 * t703;
	t642 = t623 * t624 + t628 * t629;
	t692 = t642 * t634;
	t684 = t624 * qJ(4) + t680 * t629 + pkin(2);
	t679 = -pkin(9) + pkin(10);
	t667 = qJD(6) * t642 * t670;
	t663 = t621 * t669;
	t640 = -t652 + t679;
	t631 = t686 * t642;
	t593 = t634 * qJD(1) + t607 * qJD(2);
	t591 = t606 * qJD(1) + t635 * qJD(2);
	t586 = t642 * t606;
	t578 = (-t632 * t630 - t642 * t669) * t621;
	t551 = t555 * t627 + t591 * t622 + (-t585 * t622 - t627 * t634) * qJD(6);
	t550 = -t555 * t622 + t591 * t627 + (-t585 * t627 + t622 * t634) * qJD(6);
	t1 = [-t594 * pkin(2) - t573 * qJ(4) - t597 * qJD(4) - t641 * t559 + t640 * t593 - t680 * t574 + t678 * t703 + (t705 * r_i_i_C(1) + t704 * r_i_i_C(2)) * qJD(6) + (-t677 * pkin(1) - pkin(8) * t672) * qJD(1), t641 * (t642 * t591 + t632 * t634) + t640 * t592 - t678 * (t643 * t591 - t686 * t692) + ((t622 * t692 + t627 * t635) * r_i_i_C(1) + (-t622 * t635 + t627 * t692) * r_i_i_C(2)) * qJD(6) + t695 * t634 + t684 * t591, t572 * qJ(4) + t602 * qJD(4) - t680 * t571 - t697, t571, t697, r_i_i_C(1) * t550 - r_i_i_C(2) * t551; -t592 * pkin(2) + t555 * pkin(5) + t551 * r_i_i_C(1) + t550 * r_i_i_C(2) + t571 * qJ(4) - t639 * qJD(4) + t679 * t591 + t680 * t572 - t678 * t554 + (-pkin(1) * t626 + pkin(8) * t665) * qJD(1), t641 * (-t642 * t593 + t606 * t632) - t640 * t594 - t678 * (-t643 * t593 - t631 * t606) + ((t586 * t622 - t607 * t627) * r_i_i_C(1) + (t586 * t627 + t607 * t622) * r_i_i_C(2)) * qJD(6) + t695 * t606 - t684 * t593, t574 * qJ(4) + t598 * qJD(4) - t680 * t573 - t693, t573, t693, (-t559 * t622 - t593 * t627) * r_i_i_C(1) + (-t559 * t627 + t593 * t622) * r_i_i_C(2) + (-t704 * r_i_i_C(1) + t705 * r_i_i_C(2)) * qJD(6); 0, (t578 * t627 - t622 * t667) * r_i_i_C(1) + (-t578 * t622 - t627 * t667) * r_i_i_C(2) + t578 * pkin(5) + (t678 * (-t631 * t630 + t643 * t669) + t653 * t625 * qJD(6) - t695 * t630 + (-t625 * t684 - t640 * t630) * qJD(2)) * t621, t596 * qJ(4) + t605 * qJD(4) - t680 * t595 - t694, t595, t694, (-t569 * t622 - t627 * t663) * r_i_i_C(1) + (-t569 * t627 + t622 * t663) * r_i_i_C(2) + ((-t590 * t627 - t622 * t670) * r_i_i_C(1) + (t590 * t622 - t627 * t670) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end