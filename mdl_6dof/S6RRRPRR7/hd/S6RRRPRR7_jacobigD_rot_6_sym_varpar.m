% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:19
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRPRR7_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR7_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR7_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR7_jacobigD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:19:19
% EndTime: 2019-02-26 22:19:19
% DurationCPUTime: 0.07s
% Computational Cost: add. (56->18), mult. (100->40), div. (0->0), fcn. (102->8), ass. (0->27)
t200 = sin(pkin(6));
t203 = sin(qJ(1));
t218 = t200 * t203;
t205 = cos(qJ(1));
t217 = t200 * t205;
t202 = sin(qJ(2));
t216 = t202 * t203;
t215 = t202 * t205;
t204 = cos(qJ(2));
t214 = t203 * t204;
t213 = t205 * t204;
t212 = qJD(1) * t200;
t198 = qJ(3) + pkin(12) + qJ(5);
t196 = sin(t198);
t211 = qJD(2) * t196;
t210 = qJD(2) * t200;
t201 = cos(pkin(6));
t209 = t201 * t213 - t216;
t208 = t201 * t214 + t215;
t207 = t201 * t215 + t214;
t206 = -t201 * t216 + t213;
t199 = qJD(3) + qJD(5);
t197 = cos(t198);
t195 = t202 * t210;
t194 = t208 * qJD(1) + t207 * qJD(2);
t193 = t209 * qJD(1) + t206 * qJD(2);
t1 = [0, t205 * t212, t193, 0, t193 (t196 * t218 + t206 * t197) * t199 - t208 * t211 + (-t207 * t196 - t197 * t217) * qJD(1); 0, t203 * t212, t194, 0, t194 (-t196 * t217 + t207 * t197) * t199 + t209 * t211 + (t206 * t196 - t197 * t218) * qJD(1); 0, 0, t195, 0, t195, t200 * t202 * t199 * t197 + (t199 * t201 + t204 * t210) * t196;];
JgD_rot  = t1;
