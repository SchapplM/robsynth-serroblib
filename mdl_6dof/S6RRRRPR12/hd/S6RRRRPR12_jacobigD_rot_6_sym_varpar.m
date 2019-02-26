% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR12
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:37
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRRPR12_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR12_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR12_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR12_jacobigD_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:37:00
% EndTime: 2019-02-26 22:37:00
% DurationCPUTime: 0.21s
% Computational Cost: add. (112->50), mult. (332->104), div. (0->0), fcn. (355->12), ass. (0->48)
t315 = sin(pkin(7));
t320 = sin(qJ(2));
t349 = t315 * t320;
t316 = sin(pkin(6));
t321 = sin(qJ(1));
t348 = t316 * t321;
t324 = cos(qJ(1));
t347 = t316 * t324;
t317 = cos(pkin(7));
t319 = sin(qJ(3));
t346 = t317 * t319;
t345 = t319 * t320;
t323 = cos(qJ(2));
t344 = t319 * t323;
t322 = cos(qJ(3));
t343 = t320 * t322;
t342 = t321 * t320;
t341 = t321 * t323;
t340 = t322 * t323;
t339 = t324 * t320;
t338 = t324 * t323;
t337 = qJD(1) * t316;
t314 = qJ(4) + pkin(13);
t312 = sin(t314);
t336 = qJD(3) * t312;
t335 = qJD(3) * t315;
t334 = t321 * t337;
t333 = t324 * t337;
t332 = t315 * t334;
t331 = t315 * t333;
t318 = cos(pkin(6));
t308 = t318 * t338 - t342;
t330 = t308 * t317 - t315 * t347;
t310 = -t318 * t341 - t339;
t329 = t310 * t317 + t315 * t348;
t328 = t317 * t344 + t343;
t309 = t318 * t339 + t341;
t327 = t318 * t342 - t338;
t326 = t309 * t322 + t330 * t319;
t325 = t329 * t319 - t322 * t327;
t313 = cos(t314);
t307 = -t327 * qJD(1) + t308 * qJD(2);
t306 = t310 * qJD(1) - t309 * qJD(2);
t305 = -t309 * qJD(1) + t310 * qJD(2);
t304 = -t308 * qJD(1) + t327 * qJD(2);
t303 = -t306 * t315 + t317 * t334;
t302 = -t304 * t315 + t317 * t333;
t1 = [0, t333, t302, t305 * t319 + (-t304 * t317 - t331) * t322 + t325 * qJD(3), 0 (t304 * t346 + t305 * t322 + t319 * t331) * t312 - t302 * t313 + (t325 * t313 + (-t310 * t315 + t317 * t348) * t312) * qJD(4) + (t319 * t327 + t329 * t322) * t336; 0, t334, t303, t307 * t319 + (-t306 * t317 - t332) * t322 + t326 * qJD(3), 0 (t306 * t346 + t307 * t322 + t319 * t332) * t312 - t303 * t313 + (t326 * t313 + (-t308 * t315 - t317 * t347) * t312) * qJD(4) + (-t309 * t319 + t330 * t322) * t336; 0, 0, t316 * qJD(2) * t349, t318 * t319 * t335 + (t328 * qJD(3) + (t317 * t343 + t344) * qJD(2)) * t316, 0 (t322 * t312 * t335 + (t313 * t315 * t319 + t312 * t317) * qJD(4)) * t318 + ((-t315 * t323 * t312 + t328 * t313) * qJD(4) + (t317 * t340 - t345) * t336 + ((-t317 * t345 + t340) * t312 - t313 * t349) * qJD(2)) * t316;];
JgD_rot  = t1;
