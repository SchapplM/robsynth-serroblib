% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRPPR9_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR9_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR9_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR9_jacobigD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:08:10
% EndTime: 2019-02-26 22:08:10
% DurationCPUTime: 0.06s
% Computational Cost: add. (22->16), mult. (78->40), div. (0->0), fcn. (80->8), ass. (0->21)
t203 = sin(pkin(6));
t208 = cos(qJ(3));
t222 = t203 * t208;
t210 = cos(qJ(1));
t221 = t203 * t210;
t206 = sin(qJ(2));
t207 = sin(qJ(1));
t220 = t206 * t207;
t219 = t206 * t210;
t209 = cos(qJ(2));
t218 = t207 * t209;
t217 = t210 * t209;
t216 = qJD(1) * t203;
t205 = sin(qJ(3));
t215 = qJD(2) * t205;
t204 = cos(pkin(6));
t214 = t204 * t217 - t220;
t213 = t204 * t218 + t219;
t212 = t204 * t219 + t218;
t211 = t204 * t220 - t217;
t1 = [0, t210 * t216, t214 * qJD(1) - t211 * qJD(2), 0, 0 (-t207 * t203 * t205 + t211 * t208) * qJD(3) + t213 * t215 + (t212 * t205 + t208 * t221) * qJD(1); 0, t207 * t216, t213 * qJD(1) + t212 * qJD(2), 0, 0 (t205 * t221 - t212 * t208) * qJD(3) - t214 * t215 + (t211 * t205 + t207 * t222) * qJD(1); 0, 0, t203 * qJD(2) * t206, 0, 0, -t203 * t209 * t215 + (-t204 * t205 - t206 * t222) * qJD(3);];
JgD_rot  = t1;
