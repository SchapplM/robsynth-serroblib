% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:20
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRPRR9_jacobigD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR9_jacobigD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR9_jacobigD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR9_jacobigD_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:20:28
% EndTime: 2019-02-26 22:20:28
% DurationCPUTime: 0.15s
% Computational Cost: add. (55->25), mult. (190->56), div. (0->0), fcn. (202->12), ass. (0->37)
t235 = cos(pkin(7));
t231 = sin(pkin(13));
t234 = cos(pkin(13));
t237 = sin(qJ(3));
t240 = cos(qJ(3));
t230 = -t240 * t231 - t237 * t234;
t243 = qJD(3) * t230;
t226 = t235 * t243;
t259 = qJD(2) * t230 + t226;
t233 = sin(pkin(6));
t239 = sin(qJ(1));
t258 = t233 * t239;
t242 = cos(qJ(1));
t257 = t233 * t242;
t238 = sin(qJ(2));
t256 = t239 * t238;
t241 = cos(qJ(2));
t255 = t239 * t241;
t254 = t242 * t238;
t253 = t242 * t241;
t252 = qJD(1) * t233;
t250 = t239 * t252;
t249 = t242 * t252;
t248 = t231 * t237 - t234 * t240;
t236 = cos(pkin(6));
t247 = t236 * t253 - t256;
t246 = t236 * t255 + t254;
t245 = t236 * t254 + t255;
t244 = t236 * t256 - t253;
t232 = sin(pkin(7));
t229 = t248 * qJD(3);
t228 = t248 * t235;
t227 = t248 * t232;
t225 = t232 * t243;
t224 = -t246 * qJD(1) - t245 * qJD(2);
t223 = -t247 * qJD(1) + t244 * qJD(2);
t1 = [0, t249, -t223 * t232 + t235 * t249, 0, t244 * t229 + t223 * t228 - t225 * t258 + (t227 * t257 + t245 * t230) * qJD(1) + t259 * t246, 0; 0, t250, -t224 * t232 + t235 * t250, 0, -t245 * t229 + t224 * t228 + t225 * t257 + (t227 * t258 + t244 * t230) * qJD(1) - t259 * t247, 0; 0, 0, t233 * qJD(2) * t238 * t232, 0, -t236 * t225 + (-t226 * t241 - t229 * t238 + (-t228 * t238 - t230 * t241) * qJD(2)) * t233, 0;];
JgD_rot  = t1;
