% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPRRPR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:40
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PPRRPR1_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR1_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_jacobigD_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:40:19
% EndTime: 2019-02-26 19:40:19
% DurationCPUTime: 0.08s
% Computational Cost: add. (41->23), mult. (141->56), div. (0->0), fcn. (164->12), ass. (0->29)
t223 = sin(pkin(11));
t229 = cos(pkin(6));
t245 = t223 * t229;
t224 = sin(pkin(7));
t225 = sin(pkin(6));
t244 = t224 * t225;
t243 = t224 * t229;
t228 = cos(pkin(7));
t242 = t225 * t228;
t226 = cos(pkin(12));
t241 = t226 * t228;
t227 = cos(pkin(11));
t240 = t227 * t229;
t230 = sin(qJ(4));
t239 = qJD(3) * t230;
t222 = sin(pkin(12));
t218 = -t223 * t222 + t226 * t240;
t238 = t218 * t228 - t227 * t244;
t220 = -t227 * t222 - t226 * t245;
t237 = t220 * t228 + t223 * t244;
t219 = t222 * t240 + t223 * t226;
t231 = sin(qJ(3));
t233 = cos(qJ(3));
t236 = t219 * t233 + t238 * t231;
t221 = -t222 * t245 + t227 * t226;
t235 = t221 * t233 + t237 * t231;
t234 = t231 * t243 + (t222 * t233 + t231 * t241) * t225;
t232 = cos(qJ(4));
t1 = [0, 0, 0, t235 * qJD(3), 0 (t235 * t232 + (-t220 * t224 + t223 * t242) * t230) * qJD(4) + (-t221 * t231 + t237 * t233) * t239; 0, 0, 0, t236 * qJD(3), 0 (t236 * t232 + (-t218 * t224 - t227 * t242) * t230) * qJD(4) + (-t219 * t231 + t238 * t233) * t239; 0, 0, 0, t234 * qJD(3), 0 (t234 * t232 + (-t226 * t244 + t229 * t228) * t230) * qJD(4) + (t233 * t243 + (-t222 * t231 + t233 * t241) * t225) * t239;];
JgD_rot  = t1;
