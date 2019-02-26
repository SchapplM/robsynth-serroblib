% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP14
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRPRRP14_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP14_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP14_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP14_jacobigD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:53:25
% EndTime: 2019-02-26 21:53:25
% DurationCPUTime: 0.06s
% Computational Cost: add. (22->16), mult. (78->39), div. (0->0), fcn. (80->8), ass. (0->22)
t192 = sin(pkin(6));
t196 = sin(qJ(1));
t212 = t192 * t196;
t198 = cos(qJ(2));
t211 = t192 * t198;
t199 = cos(qJ(1));
t210 = t192 * t199;
t195 = sin(qJ(2));
t209 = t196 * t195;
t208 = t196 * t198;
t207 = t198 * t199;
t206 = t199 * t195;
t205 = qJD(1) * t192;
t197 = cos(qJ(4));
t204 = qJD(2) * t197;
t193 = cos(pkin(6));
t203 = t193 * t207 - t209;
t202 = t193 * t208 + t206;
t201 = t193 * t206 + t208;
t200 = -t193 * t209 + t207;
t194 = sin(qJ(4));
t1 = [0, t199 * t205, 0, -t201 * qJD(1) - t202 * qJD(2) (t202 * t194 + t197 * t212) * qJD(4) - t200 * t204 + (t194 * t210 - t203 * t197) * qJD(1), 0; 0, t196 * t205, 0, t200 * qJD(1) + t203 * qJD(2) (-t203 * t194 - t197 * t210) * qJD(4) - t201 * t204 + (t194 * t212 - t202 * t197) * qJD(1), 0; 0, 0, 0, qJD(2) * t211, -t192 * t195 * t204 + (t193 * t197 - t194 * t211) * qJD(4), 0;];
JgD_rot  = t1;
