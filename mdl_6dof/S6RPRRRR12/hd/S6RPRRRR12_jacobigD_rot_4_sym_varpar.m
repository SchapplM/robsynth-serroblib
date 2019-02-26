% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
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

function JgD_rot = S6RPRRRR12_jacobigD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR12_jacobigD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR12_jacobigD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_jacobigD_rot_4_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:21:09
% EndTime: 2019-02-26 21:21:09
% DurationCPUTime: 0.13s
% Computational Cost: add. (31->21), mult. (115->46), div. (0->0), fcn. (119->12), ass. (0->32)
t220 = sin(pkin(7));
t221 = sin(pkin(6));
t244 = t221 * t220;
t218 = sin(pkin(14));
t227 = sin(qJ(1));
t243 = t227 * t218;
t222 = cos(pkin(14));
t242 = t227 * t222;
t229 = cos(qJ(1));
t241 = t229 * t218;
t240 = t229 * t222;
t219 = sin(pkin(8));
t239 = qJD(1) * t219;
t238 = qJD(3) * t219;
t224 = cos(pkin(7));
t228 = cos(qJ(3));
t237 = t219 * t224 * t228;
t236 = t227 * t244;
t235 = t229 * t244;
t234 = qJD(1) * t221 * t224;
t225 = cos(pkin(6));
t233 = t225 * t240 - t243;
t232 = -t225 * t242 - t241;
t231 = t225 * t241 + t242;
t230 = t225 * t243 - t240;
t226 = sin(qJ(3));
t223 = cos(pkin(8));
t217 = t232 * qJD(1);
t216 = t233 * qJD(1);
t215 = -t217 * t220 + t227 * t234;
t214 = t216 * t220 + t229 * t234;
t1 = [0, 0, t214, t216 * t237 + t214 * t223 - (t230 * t228 + (-t232 * t224 - t236) * t226) * t238 - (t231 * t226 + t228 * t235) * t239, 0, 0; 0, 0, t215, -t217 * t237 + t215 * t223 - (-t231 * t228 + (-t233 * t224 + t235) * t226) * t238 - (t230 * t226 + t228 * t236) * t239, 0, 0; 0, 0, 0 -(-t220 * t225 * t226 + (-t222 * t224 * t226 - t218 * t228) * t221) * t238, 0, 0;];
JgD_rot  = t1;
