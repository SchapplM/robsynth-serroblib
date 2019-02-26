% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRRRRR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRRRR6_jacobigD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR6_jacobigD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR6_jacobigD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRRRRR6_jacobigD_rot_4_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:21:52
% EndTime: 2019-02-26 20:21:52
% DurationCPUTime: 0.10s
% Computational Cost: add. (29->20), mult. (109->51), div. (0->0), fcn. (113->12), ass. (0->27)
t222 = sin(pkin(7));
t244 = t222 * cos(pkin(8));
t223 = sin(pkin(6));
t243 = t223 * t222;
t226 = cos(pkin(7));
t230 = cos(qJ(3));
t242 = t226 * t230;
t227 = cos(pkin(6));
t229 = sin(qJ(2));
t241 = t227 * t229;
t231 = cos(qJ(2));
t240 = t227 * t231;
t228 = sin(qJ(3));
t239 = t228 * t231;
t238 = t229 * t230;
t237 = qJD(2) * t228;
t221 = sin(pkin(8));
t236 = qJD(3) * t221;
t220 = sin(pkin(14));
t224 = cos(pkin(14));
t235 = -t220 * t229 + t224 * t240;
t234 = t220 * t231 + t224 * t241;
t233 = t220 * t240 + t224 * t229;
t232 = t220 * t241 - t224 * t231;
t219 = t232 * qJD(2);
t218 = t234 * qJD(2);
t1 = [0, 0, -t219 * t222 -(t219 * t242 + t233 * t237) * t221 - t219 * t244 - (t232 * t230 + (-t220 * t243 + t233 * t226) * t228) * t236, 0, 0; 0, 0, t218 * t222 -(-t218 * t242 - t235 * t237) * t221 + t218 * t244 - (-t234 * t230 + (t224 * t243 - t235 * t226) * t228) * t236, 0, 0; 0, 0, qJD(2) * t229 * t243, t227 * t222 * t228 * t236 + (-(-t226 * t239 - t238) * t236 + (-(-t226 * t238 - t239) * t221 + t229 * t244) * qJD(2)) * t223, 0, 0;];
JgD_rot  = t1;
