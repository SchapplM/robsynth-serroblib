% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPRRR1
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RPRRR1_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR1_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPRRR1_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRR1_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_jacobiaD_transl_sym_varpar: pkin has to be [1x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:10:18
	% EndTime: 2019-12-05 18:10:18
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:10:18
	% EndTime: 2019-12-05 18:10:18
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:10:18
	% EndTime: 2019-12-05 18:10:18
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (6->5), mult. (16->10), div. (0->0), fcn. (10->2), ass. (0->4)
	t8 = r_i_i_C(3) + qJ(2);
	t7 = cos(qJ(1));
	t6 = sin(qJ(1));
	t1 = [t7 * qJD(2) + (-r_i_i_C(1) * t7 - t8 * t6) * qJD(1), qJD(1) * t7, 0, 0, 0; t6 * qJD(2) + (-r_i_i_C(1) * t6 + t8 * t7) * qJD(1), qJD(1) * t6, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:10:18
	% EndTime: 2019-12-05 18:10:18
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (19->16), mult. (64->31), div. (0->0), fcn. (42->4), ass. (0->13)
	t27 = r_i_i_C(3) + qJ(2);
	t17 = sin(qJ(1));
	t26 = qJD(1) * t17;
	t19 = cos(qJ(1));
	t25 = qJD(1) * t19;
	t24 = qJD(3) * t17;
	t23 = qJD(3) * t19;
	t16 = sin(qJ(3));
	t18 = cos(qJ(3));
	t22 = -r_i_i_C(1) * t18 + r_i_i_C(2) * t16;
	t21 = r_i_i_C(1) * t16 + r_i_i_C(2) * t18;
	t20 = t21 * qJD(3);
	t1 = [t19 * qJD(2) + t21 * t24 + (-t17 * t27 + t19 * t22) * qJD(1), t25, (t16 * t23 + t18 * t26) * r_i_i_C(2) + (t16 * t26 - t18 * t23) * r_i_i_C(1), 0, 0; t17 * qJD(2) - t19 * t20 + (t17 * t22 + t19 * t27) * qJD(1), t26, (t16 * t24 - t18 * t25) * r_i_i_C(2) + (-t16 * t25 - t18 * t24) * r_i_i_C(1), 0, 0; 0, 0, -t20, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:10:19
	% EndTime: 2019-12-05 18:10:19
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (65->35), mult. (218->71), div. (0->0), fcn. (179->6), ass. (0->27)
	t200 = cos(qJ(4));
	t202 = cos(qJ(1));
	t218 = t200 * t202;
	t199 = sin(qJ(1));
	t217 = qJD(1) * t199;
	t216 = qJD(1) * t202;
	t198 = sin(qJ(3));
	t215 = qJD(3) * t198;
	t201 = cos(qJ(3));
	t214 = qJD(3) * t201;
	t213 = qJD(3) * t202;
	t212 = qJD(4) * t198;
	t211 = qJD(4) * t201;
	t210 = -qJD(1) + t211;
	t209 = qJD(1) * t201 - qJD(4);
	t197 = sin(qJ(4));
	t208 = r_i_i_C(1) * t200 - r_i_i_C(2) * t197;
	t207 = r_i_i_C(1) * t197 + r_i_i_C(2) * t200;
	t206 = t210 * t197;
	t205 = qJD(3) * t208;
	t204 = -r_i_i_C(3) * qJD(3) + t207 * qJD(4);
	t203 = t198 * t213 + t209 * t199;
	t196 = -t209 * t218 + (t200 * t215 + t206) * t199;
	t195 = t210 * t200 * t199 + (-t199 * t215 + t209 * t202) * t197;
	t194 = t203 * t200 + t202 * t206;
	t193 = t203 * t197 - t210 * t218;
	t1 = [-qJ(2) * t217 + t196 * r_i_i_C(1) + t195 * r_i_i_C(2) + t202 * qJD(2) + (-t198 * t216 - t199 * t214) * r_i_i_C(3), t216, (-r_i_i_C(3) * t217 - t202 * t205) * t201 + (t204 * t202 + t208 * t217) * t198, t193 * r_i_i_C(1) + t194 * r_i_i_C(2), 0; qJ(2) * t216 - t194 * r_i_i_C(1) + t193 * r_i_i_C(2) + t199 * qJD(2) + (-t198 * t217 + t201 * t213) * r_i_i_C(3), t217, (r_i_i_C(3) * t216 - t199 * t205) * t201 + (t204 * t199 - t208 * t216) * t198, -t195 * r_i_i_C(1) + t196 * r_i_i_C(2), 0; 0, 0, -t207 * t211 + (r_i_i_C(3) * t201 - t208 * t198) * qJD(3), (t197 * t212 - t200 * t214) * r_i_i_C(2) + (-t197 * t214 - t200 * t212) * r_i_i_C(1), 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:10:20
	% EndTime: 2019-12-05 18:10:20
	% DurationCPUTime: 0.35s
	% Computational Cost: add. (185->67), mult. (600->128), div. (0->0), fcn. (566->8), ass. (0->57)
	t311 = cos(qJ(4));
	t358 = -qJD(5) * t311 + qJD(3);
	t308 = sin(qJ(3));
	t309 = sin(qJ(1));
	t312 = cos(qJ(3));
	t329 = qJD(1) * t312 - qJD(4);
	t313 = cos(qJ(1));
	t342 = qJD(3) * t313;
	t357 = t308 * t342 + t329 * t309;
	t343 = qJD(3) * t312;
	t344 = qJD(1) * t313;
	t307 = sin(qJ(4));
	t347 = t313 * t307;
	t348 = t309 * t312;
	t318 = -qJD(5) * (t311 * t348 - t347) + t308 * t344 + t309 * t343;
	t346 = t313 * t311;
	t301 = t309 * t307 + t312 * t346;
	t332 = qJD(3) * t308 * t309;
	t333 = t307 * t348;
	t340 = qJD(4) * t311;
	t297 = t301 * qJD(1) - qJD(4) * t333 - t311 * t332 - t313 * t340;
	t338 = qJD(5) * t308;
	t325 = t309 * t338 + t297;
	t356 = t318 * r_i_i_C(1) - t325 * r_i_i_C(2);
	t306 = sin(qJ(5));
	t310 = cos(qJ(5));
	t341 = qJD(4) * t307;
	t319 = t306 * t341 + t358 * t310;
	t320 = -t358 * t306 + t310 * t341;
	t355 = t320 * r_i_i_C(1) - t319 * r_i_i_C(2) - r_i_i_C(3) * t340;
	t350 = t308 * t311;
	t353 = r_i_i_C(3) * t307;
	t354 = (-t306 * t312 + t310 * t350) * r_i_i_C(1) - (t306 * t350 + t310 * t312) * r_i_i_C(2) + t308 * t353;
	t352 = t306 * r_i_i_C(1);
	t351 = t306 * r_i_i_C(2);
	t349 = t309 * t311;
	t345 = qJD(1) * t309;
	t339 = qJD(5) * t306;
	t337 = qJD(5) * t310;
	t330 = -qJD(4) * t312 + qJD(1);
	t328 = qJD(3) * t311 - qJD(5);
	t327 = -t310 * r_i_i_C(1) + t351;
	t323 = t330 * t313;
	t295 = t307 * t323 - t357 * t311;
	t326 = t313 * t338 + t295;
	t324 = t310 * t328;
	t317 = -qJD(5) * t301 - t308 * t345 + t312 * t342;
	t316 = -r_i_i_C(1) * t324 - qJD(3) * t353 + t328 * t351;
	t315 = -t325 * r_i_i_C(1) - t318 * r_i_i_C(2);
	t314 = t355 * t308 + t316 * t312;
	t300 = -t312 * t347 + t349;
	t298 = -t333 - t346;
	t296 = t330 * t349 + (-t329 * t313 + t332) * t307;
	t294 = t357 * t307 + t311 * t323;
	t293 = t317 * t306 + t326 * t310;
	t292 = -t326 * t306 + t317 * t310;
	t1 = [t296 * r_i_i_C(3) - qJ(2) * t345 + t313 * qJD(2) - t356 * t306 + t315 * t310, t344, t314 * t313 + t354 * t345, t295 * r_i_i_C(3) + (-t294 * t306 - t300 * t337) * r_i_i_C(2) + (t294 * t310 - t300 * t339) * r_i_i_C(1), t292 * r_i_i_C(1) - t293 * r_i_i_C(2); t293 * r_i_i_C(1) + t292 * r_i_i_C(2) - t294 * r_i_i_C(3) + qJ(2) * t344 + t309 * qJD(2), t345, t314 * t309 - t354 * t344, t297 * r_i_i_C(3) + (-t296 * t306 - t298 * t337) * r_i_i_C(2) + (t296 * t310 - t298 * t339) * r_i_i_C(1), t315 * t306 + t356 * t310; 0, 0, t316 * t308 - t355 * t312, (t327 * t308 * qJD(4) + r_i_i_C(3) * t343) * t311 + (t327 * t343 + (-qJD(4) * r_i_i_C(3) + (t310 * r_i_i_C(2) + t352) * qJD(5)) * t308) * t307, (-r_i_i_C(2) * t324 - t328 * t352) * t312 + (t319 * r_i_i_C(1) + t320 * r_i_i_C(2)) * t308;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end