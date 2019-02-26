% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRP7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:12
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRPRP7_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP7_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP7_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP7_jacobigD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:12:39
% EndTime: 2019-02-26 22:12:39
% DurationCPUTime: 0.06s
% Computational Cost: add. (33->17), mult. (78->40), div. (0->0), fcn. (80->8), ass. (0->23)
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
t199 = qJ(3) + pkin(11);
t197 = sin(t199);
t211 = qJD(2) * t197;
t210 = qJD(2) * t200;
t201 = cos(pkin(6));
t209 = t201 * t213 - t216;
t208 = t201 * t214 + t215;
t207 = t201 * t215 + t214;
t206 = -t201 * t216 + t213;
t198 = cos(t199);
t1 = [0, t205 * t212, t209 * qJD(1) + t206 * qJD(2), 0 (t197 * t218 + t206 * t198) * qJD(3) - t208 * t211 + (-t207 * t197 - t198 * t217) * qJD(1), 0; 0, t203 * t212, t208 * qJD(1) + t207 * qJD(2), 0 (-t197 * t217 + t207 * t198) * qJD(3) + t209 * t211 + (t206 * t197 - t198 * t218) * qJD(1), 0; 0, 0, t202 * t210, 0, t204 * t197 * t210 + (t198 * t200 * t202 + t197 * t201) * qJD(3), 0;];
JgD_rot  = t1;
