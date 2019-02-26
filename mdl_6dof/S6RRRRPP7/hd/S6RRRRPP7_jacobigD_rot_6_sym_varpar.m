% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPP7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,theta5]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:29
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRRPP7_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP7_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP7_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPP7_jacobigD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:28:55
% EndTime: 2019-02-26 22:28:55
% DurationCPUTime: 0.06s
% Computational Cost: add. (22->16), mult. (78->40), div. (0->0), fcn. (80->8), ass. (0->21)
t196 = sin(pkin(6));
t201 = cos(qJ(3));
t215 = t196 * t201;
t203 = cos(qJ(1));
t214 = t196 * t203;
t199 = sin(qJ(2));
t200 = sin(qJ(1));
t213 = t199 * t200;
t212 = t199 * t203;
t202 = cos(qJ(2));
t211 = t200 * t202;
t210 = t203 * t202;
t209 = qJD(1) * t196;
t198 = sin(qJ(3));
t208 = qJD(2) * t198;
t197 = cos(pkin(6));
t207 = t197 * t210 - t213;
t206 = t197 * t211 + t212;
t205 = t197 * t212 + t211;
t204 = -t197 * t213 + t210;
t1 = [0, t203 * t209, t207 * qJD(1) + t204 * qJD(2) (t200 * t196 * t198 + t204 * t201) * qJD(3) - t206 * t208 + (-t205 * t198 - t201 * t214) * qJD(1), 0, 0; 0, t200 * t209, t206 * qJD(1) + t205 * qJD(2) (-t198 * t214 + t205 * t201) * qJD(3) + t207 * t208 + (t204 * t198 - t200 * t215) * qJD(1), 0, 0; 0, 0, t196 * qJD(2) * t199, t196 * t202 * t208 + (t197 * t198 + t199 * t215) * qJD(3), 0, 0;];
JgD_rot  = t1;
