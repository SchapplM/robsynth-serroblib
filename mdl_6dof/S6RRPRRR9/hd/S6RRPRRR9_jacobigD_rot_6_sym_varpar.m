% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR9
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
% Datum: 2019-02-26 21:58
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRPRRR9_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR9_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR9_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR9_jacobigD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:58:47
% EndTime: 2019-02-26 21:58:47
% DurationCPUTime: 0.07s
% Computational Cost: add. (56->18), mult. (100->40), div. (0->0), fcn. (102->8), ass. (0->27)
t199 = sin(pkin(6));
t202 = sin(qJ(1));
t217 = t199 * t202;
t204 = cos(qJ(1));
t216 = t199 * t204;
t201 = sin(qJ(2));
t215 = t201 * t202;
t214 = t201 * t204;
t203 = cos(qJ(2));
t213 = t202 * t203;
t212 = t204 * t203;
t211 = qJD(1) * t199;
t197 = pkin(12) + qJ(4) + qJ(5);
t195 = sin(t197);
t210 = qJD(2) * t195;
t209 = qJD(2) * t199;
t200 = cos(pkin(6));
t208 = t200 * t212 - t215;
t207 = t200 * t213 + t214;
t206 = t200 * t214 + t213;
t205 = -t200 * t215 + t212;
t198 = qJD(4) + qJD(5);
t196 = cos(t197);
t194 = t201 * t209;
t193 = t207 * qJD(1) + t206 * qJD(2);
t192 = t208 * qJD(1) + t205 * qJD(2);
t1 = [0, t204 * t211, 0, t192, t192 (t195 * t217 + t205 * t196) * t198 - t207 * t210 + (-t206 * t195 - t196 * t216) * qJD(1); 0, t202 * t211, 0, t193, t193 (-t195 * t216 + t206 * t196) * t198 + t208 * t210 + (t205 * t195 - t196 * t217) * qJD(1); 0, 0, 0, t194, t194, t199 * t201 * t198 * t196 + (t198 * t200 + t203 * t209) * t195;];
JgD_rot  = t1;
