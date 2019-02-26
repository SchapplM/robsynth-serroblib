% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRP10
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
% Datum: 2019-02-26 22:14
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRPRP10_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP10_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP10_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP10_jacobigD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:14:24
% EndTime: 2019-02-26 22:14:24
% DurationCPUTime: 0.06s
% Computational Cost: add. (22->16), mult. (78->40), div. (0->0), fcn. (80->8), ass. (0->21)
t195 = sin(pkin(6));
t200 = cos(qJ(3));
t214 = t195 * t200;
t202 = cos(qJ(1));
t213 = t195 * t202;
t198 = sin(qJ(2));
t199 = sin(qJ(1));
t212 = t198 * t199;
t211 = t198 * t202;
t201 = cos(qJ(2));
t210 = t199 * t201;
t209 = t202 * t201;
t208 = qJD(1) * t195;
t197 = sin(qJ(3));
t207 = qJD(2) * t197;
t196 = cos(pkin(6));
t206 = t196 * t209 - t212;
t205 = t196 * t210 + t211;
t204 = t196 * t211 + t210;
t203 = -t196 * t212 + t209;
t1 = [0, t202 * t208, t206 * qJD(1) + t203 * qJD(2), 0 (t199 * t195 * t197 + t203 * t200) * qJD(3) - t205 * t207 + (-t204 * t197 - t200 * t213) * qJD(1), 0; 0, t199 * t208, t205 * qJD(1) + t204 * qJD(2), 0 (-t197 * t213 + t204 * t200) * qJD(3) + t206 * t207 + (t203 * t197 - t199 * t214) * qJD(1), 0; 0, 0, t195 * qJD(2) * t198, 0, t195 * t201 * t207 + (t196 * t197 + t198 * t214) * qJD(3), 0;];
JgD_rot  = t1;
