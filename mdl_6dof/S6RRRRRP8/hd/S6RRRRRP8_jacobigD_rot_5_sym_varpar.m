% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRRP8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:43
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRRRP8_jacobigD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP8_jacobigD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP8_jacobigD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP8_jacobigD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:43:49
% EndTime: 2019-02-26 22:43:49
% DurationCPUTime: 0.06s
% Computational Cost: add. (45->18), mult. (100->40), div. (0->0), fcn. (102->8), ass. (0->27)
t197 = sin(pkin(6));
t200 = sin(qJ(1));
t215 = t197 * t200;
t202 = cos(qJ(1));
t214 = t197 * t202;
t199 = sin(qJ(2));
t213 = t199 * t200;
t212 = t199 * t202;
t201 = cos(qJ(2));
t211 = t200 * t201;
t210 = t202 * t201;
t209 = qJD(1) * t197;
t196 = qJ(3) + qJ(4);
t193 = sin(t196);
t208 = qJD(2) * t193;
t207 = qJD(2) * t197;
t198 = cos(pkin(6));
t206 = t198 * t210 - t213;
t205 = t198 * t211 + t212;
t204 = t198 * t212 + t211;
t203 = -t198 * t213 + t210;
t195 = qJD(3) + qJD(4);
t194 = cos(t196);
t192 = t199 * t207;
t191 = t205 * qJD(1) + t204 * qJD(2);
t190 = t206 * qJD(1) + t203 * qJD(2);
t1 = [0, t202 * t209, t190, t190 (t193 * t215 + t203 * t194) * t195 - t205 * t208 + (-t204 * t193 - t194 * t214) * qJD(1), 0; 0, t200 * t209, t191, t191 (-t193 * t214 + t204 * t194) * t195 + t206 * t208 + (t203 * t193 - t194 * t215) * qJD(1), 0; 0, 0, t192, t192, t197 * t199 * t195 * t194 + (t195 * t198 + t201 * t207) * t193, 0;];
JgD_rot  = t1;
