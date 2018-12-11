% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPRRR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-12-10 18:38
% Revision: bb42a8b95257d9bc83910d26e849f5825122f662 (2018-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRPRRR14_jacobigD_rot_4_floatb_twist_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobigD_rot_4_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14_jacobigD_rot_4_floatb_twist_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobigD_rot_4_floatb_twist_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-10 18:38:22
% EndTime: 2018-12-10 18:38:22
% DurationCPUTime: 0.12s
% Computational Cost: add. (88->34), mult. (131->64), div. (0->0), fcn. (102->18), ass. (0->38)
t343 = qJD(2) / 0.2e1;
t321 = pkin(7) + pkin(14);
t322 = pkin(7) - pkin(14);
t328 = sin(pkin(6));
t342 = t328 * (sin(t321) / 0.2e1 + sin(t322) / 0.2e1);
t341 = qJD(1) * t328;
t323 = pkin(6) + qJ(2);
t317 = sin(t323);
t340 = qJD(2) * t317;
t324 = pkin(6) - qJ(2);
t320 = cos(t324);
t339 = qJD(2) * t320;
t332 = sin(qJ(1));
t338 = qJD(2) * t332;
t334 = cos(qJ(1));
t337 = qJD(2) * t334;
t336 = t332 * t341;
t335 = t334 * t341;
t333 = cos(qJ(2));
t331 = sin(qJ(2));
t330 = cos(pkin(7));
t329 = cos(pkin(8));
t327 = sin(pkin(7));
t326 = sin(pkin(8));
t325 = sin(pkin(14));
t319 = cos(t323);
t318 = sin(t324);
t316 = t319 * t343;
t315 = t318 * t343;
t314 = t320 / 0.2e1 + t319 / 0.2e1;
t313 = t317 / 0.2e1 - t318 / 0.2e1;
t312 = cos(t322) / 0.2e1 + cos(t321) / 0.2e1;
t310 = t316 - t339 / 0.2e1;
t309 = t316 + t339 / 0.2e1;
t308 = t315 - t340 / 0.2e1;
t307 = -t333 * t338 + t334 * t308 + (-t314 * t332 - t331 * t334) * qJD(1);
t306 = -t333 * t337 - t332 * t308 + (-t314 * t334 + t331 * t332) * qJD(1);
t1 = [0, t335, 0 -(-(-t332 * t309 - t331 * t337) * t325 + t306 * t312 + (-(-t313 * t334 - t332 * t333) * t325 + t334 * t342) * qJD(1)) * t326 + (-t306 * t327 + t330 * t335) * t329, 0, 0; 0, t336, 0 -(-(t334 * t309 - t331 * t338) * t325 + t307 * t312 + (-(-t313 * t332 + t333 * t334) * t325 + t332 * t342) * qJD(1)) * t326 + (-t307 * t327 + t330 * t336) * t329, 0, 0; 0, 0, 0 -(-(t315 + t340 / 0.2e1) * t325 + t310 * t312) * t326 - t310 * t327 * t329, 0, 0;];
JgD_rot  = t1;
