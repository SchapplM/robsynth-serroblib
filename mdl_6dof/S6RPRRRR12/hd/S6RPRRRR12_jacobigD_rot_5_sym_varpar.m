% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRRR12
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:21
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RPRRRR12_jacobigD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR12_jacobigD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR12_jacobigD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_jacobigD_rot_5_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:21:11
% EndTime: 2019-02-26 21:21:11
% DurationCPUTime: 0.24s
% Computational Cost: add. (114->45), mult. (395->100), div. (0->0), fcn. (432->14), ass. (0->49)
t313 = sin(pkin(14));
t316 = sin(pkin(6));
t322 = sin(qJ(3));
t325 = cos(qJ(3));
t317 = cos(pkin(14));
t319 = cos(pkin(7));
t343 = t317 * t319;
t315 = sin(pkin(7));
t320 = cos(pkin(6));
t346 = t315 * t320;
t351 = (t313 * t325 + t322 * t343) * t316 + t322 * t346;
t326 = cos(qJ(1));
t323 = sin(qJ(1));
t341 = t323 * t313;
t312 = t317 * t326 - t320 * t341;
t340 = t323 * t317;
t311 = -t313 * t326 - t320 * t340;
t345 = t316 * t323;
t333 = t311 * t319 + t315 * t345;
t350 = t312 * t325 + t333 * t322;
t342 = t320 * t326;
t310 = t313 * t342 + t340;
t309 = t317 * t342 - t341;
t344 = t316 * t326;
t334 = -t309 * t319 + t315 * t344;
t349 = -t310 * t325 + t334 * t322;
t339 = qJD(1) * t316;
t321 = sin(qJ(4));
t338 = qJD(3) * t321;
t336 = t323 * t339;
t335 = t326 * t339;
t305 = t309 * qJD(1);
t331 = -t305 * t319 + t315 * t335;
t307 = t311 * qJD(1);
t330 = t307 * t319 + t315 * t336;
t329 = -t310 * t322 - t334 * t325;
t328 = -t312 * t322 + t333 * t325;
t327 = t325 * t346 + (-t313 * t322 + t325 * t343) * t316;
t324 = cos(qJ(4));
t318 = cos(pkin(8));
t314 = sin(pkin(8));
t308 = t312 * qJD(1);
t306 = t310 * qJD(1);
t304 = -t307 * t315 + t319 * t336;
t303 = t305 * t315 + t319 * t335;
t302 = t351 * qJD(3);
t301 = t349 * qJD(3) - t308 * t322 + t330 * t325;
t300 = -t350 * qJD(3) + t306 * t322 + t331 * t325;
t1 = [0, 0, t303, -t300 * t314 + t303 * t318 (-t306 * t325 + t331 * t322) * t321 + (-t300 * t318 - t303 * t314) * t324 + t328 * t338 + (t350 * t324 + (t328 * t318 + (-t311 * t315 + t319 * t345) * t314) * t321) * qJD(4), 0; 0, 0, t304, -t301 * t314 + t304 * t318 (t308 * t325 + t330 * t322) * t321 + (-t301 * t318 - t304 * t314) * t324 + t329 * t338 + (-t349 * t324 + (t329 * t318 + (-t309 * t315 - t319 * t344) * t314) * t321) * qJD(4), 0; 0, 0, 0, t302 * t314, t302 * t318 * t324 + (t351 * t324 + (t327 * t318 + (-t316 * t317 * t315 + t320 * t319) * t314) * t321) * qJD(4) + t327 * t338, 0;];
JgD_rot  = t1;
