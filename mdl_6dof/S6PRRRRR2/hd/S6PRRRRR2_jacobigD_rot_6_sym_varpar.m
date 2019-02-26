% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:19
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRRRR2_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR2_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR2_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR2_jacobigD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:19:15
% EndTime: 2019-02-26 20:19:15
% DurationCPUTime: 0.05s
% Computational Cost: add. (54->12), mult. (96->30), div. (0->0), fcn. (100->8), ass. (0->24)
t200 = qJ(3) + qJ(4);
t197 = sin(t200);
t202 = sin(pkin(6));
t213 = t202 * t197;
t204 = cos(pkin(6));
t205 = sin(qJ(2));
t212 = t204 * t205;
t206 = cos(qJ(2));
t211 = t204 * t206;
t210 = qJD(2) * t197;
t209 = qJD(2) * t202;
t201 = sin(pkin(12));
t203 = cos(pkin(12));
t208 = t201 * t206 + t203 * t212;
t207 = -t201 * t212 + t203 * t206;
t199 = qJD(3) + qJD(4);
t198 = cos(t200);
t196 = t205 * t209;
t195 = t207 * qJD(2);
t194 = t208 * qJD(2);
t193 = t202 * t205 * t199 * t198 + (t199 * t204 + t206 * t209) * t197;
t192 = (t207 * t198 + t201 * t213) * t199 + (-t201 * t211 - t203 * t205) * t210;
t191 = (t208 * t198 - t203 * t213) * t199 + (-t201 * t205 + t203 * t211) * t210;
t1 = [0, 0, t195, t195, t192, t192; 0, 0, t194, t194, t191, t191; 0, 0, t196, t196, t193, t193;];
JgD_rot  = t1;
