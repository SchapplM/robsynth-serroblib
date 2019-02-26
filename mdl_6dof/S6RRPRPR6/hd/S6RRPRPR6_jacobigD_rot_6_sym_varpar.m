% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR6
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:40
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRPRPR6_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR6_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR6_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR6_jacobigD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:40:37
% EndTime: 2019-02-26 21:40:37
% DurationCPUTime: 0.09s
% Computational Cost: add. (45->22), mult. (147->54), div. (0->0), fcn. (162->10), ass. (0->24)
t197 = sin(pkin(11));
t199 = cos(pkin(11));
t202 = sin(qJ(2));
t205 = cos(qJ(2));
t208 = t202 * t197 - t205 * t199;
t215 = t208 * qJD(2);
t209 = t205 * t197 + t202 * t199;
t194 = t209 * qJD(2);
t198 = sin(pkin(6));
t203 = sin(qJ(1));
t214 = t198 * t203;
t206 = cos(qJ(1));
t213 = t198 * t206;
t212 = qJD(1) * t198;
t200 = cos(pkin(6));
t192 = t209 * t200;
t211 = t206 * t192 - t203 * t208;
t210 = -t203 * t192 - t206 * t208;
t204 = cos(qJ(4));
t201 = sin(qJ(4));
t191 = t208 * t200;
t190 = t200 * t215;
t189 = t200 * t194;
t1 = [0, t206 * t212, 0, -t203 * t189 - t206 * t215 + (-t191 * t206 - t203 * t209) * qJD(1), 0 (t203 * t190 - t206 * t194) * t204 + (-t210 * t201 + t204 * t214) * qJD(4) + (t201 * t213 - t211 * t204) * qJD(1); 0, t203 * t212, 0, t206 * t189 - t203 * t215 + (-t191 * t203 + t206 * t209) * qJD(1), 0 (-t206 * t190 - t203 * t194) * t204 + (-t211 * t201 - t204 * t213) * qJD(4) + (t201 * t214 + t210 * t204) * qJD(1); 0, 0, 0, t198 * t194, 0, t200 * qJD(4) * t204 + (-t209 * qJD(4) * t201 - t204 * t215) * t198;];
JgD_rot  = t1;
