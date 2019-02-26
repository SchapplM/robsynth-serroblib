% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRPRRR2_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR2_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR2_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR2_jacobigD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:54:12
% EndTime: 2019-02-26 19:54:12
% DurationCPUTime: 0.06s
% Computational Cost: add. (47->14), mult. (164->37), div. (0->0), fcn. (178->10), ass. (0->22)
t201 = sin(pkin(12));
t204 = cos(pkin(12));
t208 = sin(qJ(2));
t210 = cos(qJ(2));
t212 = t208 * t201 - t210 * t204;
t215 = t212 * qJD(2);
t213 = t201 * t210 + t204 * t208;
t199 = t213 * qJD(2);
t203 = sin(pkin(6));
t207 = sin(qJ(4));
t214 = t203 * t207;
t209 = cos(qJ(4));
t206 = cos(pkin(6));
t205 = cos(pkin(11));
t202 = sin(pkin(11));
t197 = t213 * t206;
t196 = t206 * t215;
t195 = t206 * t199;
t194 = t206 * qJD(4) * t207 + (t213 * qJD(4) * t209 - t207 * t215) * t203;
t193 = (t202 * t196 - t205 * t199) * t207 + ((-t202 * t197 - t205 * t212) * t209 + t202 * t214) * qJD(4);
t192 = (-t205 * t196 - t202 * t199) * t207 + ((t205 * t197 - t202 * t212) * t209 - t205 * t214) * qJD(4);
t1 = [0, 0, 0, -t202 * t195 - t205 * t215, t193, t193; 0, 0, 0, t205 * t195 - t202 * t215, t192, t192; 0, 0, 0, t203 * t199, t194, t194;];
JgD_rot  = t1;
