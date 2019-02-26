% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR5
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:56
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRPRRR5_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR5_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR5_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR5_jacobigD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:56:20
% EndTime: 2019-02-26 21:56:20
% DurationCPUTime: 0.09s
% Computational Cost: add. (75->22), mult. (240->54), div. (0->0), fcn. (266->10), ass. (0->27)
t223 = sin(pkin(12));
t225 = cos(pkin(12));
t228 = sin(qJ(2));
t231 = cos(qJ(2));
t234 = t228 * t223 - t231 * t225;
t241 = t234 * qJD(2);
t235 = t231 * t223 + t228 * t225;
t220 = t235 * qJD(2);
t224 = sin(pkin(6));
t229 = sin(qJ(1));
t240 = t224 * t229;
t232 = cos(qJ(1));
t239 = t224 * t232;
t238 = qJD(1) * t224;
t226 = cos(pkin(6));
t218 = t235 * t226;
t237 = t232 * t218 - t229 * t234;
t236 = -t229 * t218 - t232 * t234;
t230 = cos(qJ(4));
t227 = sin(qJ(4));
t217 = t234 * t226;
t216 = t226 * t241;
t215 = t226 * t220;
t214 = t226 * qJD(4) * t227 + (t235 * qJD(4) * t230 - t227 * t241) * t224;
t213 = (-t232 * t216 - t229 * t220) * t227 + (-t227 * t239 + t237 * t230) * qJD(4) + (t236 * t227 - t230 * t240) * qJD(1);
t212 = (t229 * t216 - t232 * t220) * t227 + (t227 * t240 + t236 * t230) * qJD(4) + (-t237 * t227 - t230 * t239) * qJD(1);
t1 = [0, t232 * t238, 0, -t229 * t215 - t232 * t241 + (-t217 * t232 - t229 * t235) * qJD(1), t212, t212; 0, t229 * t238, 0, t232 * t215 - t229 * t241 + (-t217 * t229 + t232 * t235) * qJD(1), t213, t213; 0, 0, 0, t224 * t220, t214, t214;];
JgD_rot  = t1;
