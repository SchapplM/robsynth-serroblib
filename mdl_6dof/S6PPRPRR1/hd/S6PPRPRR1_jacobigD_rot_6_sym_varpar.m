% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PPRPRR1_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRPRR1_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_jacobigD_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:39:49
% EndTime: 2019-02-26 19:39:49
% DurationCPUTime: 0.16s
% Computational Cost: add. (68->33), mult. (231->74), div. (0->0), fcn. (263->14), ass. (0->34)
t241 = sin(pkin(13));
t246 = cos(pkin(13));
t252 = sin(qJ(3));
t254 = cos(qJ(3));
t256 = t252 * t241 - t254 * t246;
t263 = t256 * qJD(3);
t243 = sin(pkin(11));
t245 = sin(pkin(6));
t262 = t243 * t245;
t250 = cos(pkin(6));
t261 = t243 * t250;
t248 = cos(pkin(11));
t260 = t245 * t248;
t249 = cos(pkin(7));
t259 = t245 * t249;
t258 = t248 * t250;
t257 = t241 * t254 + t246 * t252;
t239 = t257 * qJD(3);
t253 = cos(qJ(5));
t251 = sin(qJ(5));
t247 = cos(pkin(12));
t244 = sin(pkin(7));
t242 = sin(pkin(12));
t237 = -t242 * t261 + t248 * t247;
t236 = -t248 * t242 - t247 * t261;
t235 = t242 * t258 + t243 * t247;
t234 = -t243 * t242 + t247 * t258;
t233 = t257 * t249;
t232 = t257 * t244;
t231 = t249 * t263;
t230 = t249 * t239;
t229 = t244 * t263;
t228 = t244 * t239;
t1 = [0, 0, 0, 0, t228 * t262 + t236 * t230 - t237 * t263 (-t229 * t262 - t236 * t231 - t237 * t239) * t251 + ((t232 * t262 + t236 * t233 - t237 * t256) * t253 + (-t236 * t244 + t243 * t259) * t251) * qJD(5); 0, 0, 0, 0, -t228 * t260 + t234 * t230 - t235 * t263 (t229 * t260 - t234 * t231 - t235 * t239) * t251 + ((-t232 * t260 + t234 * t233 - t235 * t256) * t253 + (-t234 * t244 - t248 * t259) * t251) * qJD(5); 0, 0, 0, 0, t250 * t228 + (t230 * t247 - t242 * t263) * t245 (-t250 * t229 + (-t231 * t247 - t239 * t242) * t245) * t251 + ((t250 * t232 + (t233 * t247 - t242 * t256) * t245) * t253 + (-t245 * t247 * t244 + t250 * t249) * t251) * qJD(5);];
JgD_rot  = t1;
